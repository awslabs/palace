#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Bounded PBS stage runner (four-edge-physics-01..-13 / gallery run_stages.py, the
executable, its hash and the MPI wrapper taken from the plan instead of constants).

Reads a plan JSON, verifies the node identity, the absence of conflicting native
processes, the frozen executable hash and every pinned input hash, samples node
memory, then runs each stage under the executable with per-stage `timeout` caps and
a global deadline.  Never reuses an existing output directory or plan directory.
Writes <plan-dir>/status.json (stage states, parsed Palace log: H1 / ND / RT counts,
PCG iterations, per-source timings, peak memory, elapsed-time report).  Runs on the
compute node only (stdlib only, no repository imports).

A multi-node plan (decision 457: Nodes > 1, Ranks = Nodes x RanksPerNode, MPIExecArguments,
NodeGuard) runs on every node of PBS_NODEFILE: the preflight verifies the node count,
checks the conflicting processes and MemAvailable >= the admission guard on EVERY node
(pbsdsh, ssh as the fallback), a memory sampler runs on every node (this script in
--node-sampler mode) and every stage records its per-node sampled peaks; the launch adds
the plan's MPI arguments (the hostfile, RanksPerNode ranks per node).  A one-node plan
runs exactly as before.

usage: run_stages.py PLAN.json
       run_stages.py --node-sampler CSV STOPFILE   (the per-node sampler, launched by the runner)
"""
import hashlib
import ipaddress
import json
import os
import re
import signal
import socket
import subprocess
import sys
import threading
import time
from pathlib import Path

PATTERNS = {
    "H1": r"H1 \(p = (\d+)\): ([0-9,]+), ND \(p = \d+\): ([0-9,]+), RT \(p = \d+\): ([0-9,]+)",
    "Elements": r"^ elements\s+\d+\s+\d+\s+\d+\s+(\d+)",
    "It": r"It (\d+)/(\d+): Index = (\d+) \(elapsed time = ([0-9.eE+\-]+) s\)",
    "PCG": r"PCG solver converged in (\d+) iterations",
    "NotConverged": r"did NOT converge|not converged",
    "SourceTiming": r"Response source timing: index=(\d+), iterations=(\d+), solve_seconds=([0-9.eE+\-]+), total_seconds=([0-9.eE+\-]+)",
    "Pairs": r"Interface response matrix: (\d+)/(\d+) basis pairs",
    "Blocks": r"Archived response block pair (\d+)/(\d+)",
    "StreamedSources": r"Archived response source (\d+)/(\d+)",
    "ReductionSamples": r"Archived response reduction: (\d+) sources, (\d+) interfaces, quadrature samples per rank min (\d+), max (\d+), total (\d+)",
    "PeakTotal": r"Estimated peak per-node memory usage is: Min\. ([0-9.]+[KMGT]?), Max\. ([0-9.]+[KMGT]?), Avg\. ([0-9.]+[KMGT]?), Total ([0-9.]+[KMGT]?)",
    "Total": r"^Total\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s*$",
}


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def meminfo():
    values = {}
    for line in Path("/proc/meminfo").read_text().splitlines():
        key, rest = line.split(":", 1)
        values[key] = int(rest.split()[0]) * 1024
    return values


def normalize_host(entry):
    entry = entry.split("/")[0].strip()
    try:
        return socket.gethostbyaddr(str(ipaddress.ip_address(entry)))[0].split(".")[0]
    except ValueError:
        return entry.split(".")[0]


def node_sampler(csv_path, stop_path, interval=5.0):
    """The per-node memory sampler (every node of a multi-node job): node-wide used bytes
    every `interval` seconds into `csv_path` until `stop_path` exists."""
    with Path(csv_path).open("w") as stream:
        stream.write("unix,host,used_bytes\n")
        host = socket.gethostname().split(".")[0]
        while not Path(stop_path).exists():
            values = meminfo()
            stream.write(f"{time.time():.0f},{host},{values['MemTotal'] - values['MemAvailable']}\n")
            stream.flush()
            time.sleep(interval)


def unique_hosts(nodefile_text):
    """The job's hosts in PBS_NODEFILE order (one entry per rank slot; duplicates dropped)."""
    hosts = []
    for entry in nodefile_text.split():
        host = normalize_host(entry)
        if host and host not in hosts:
            hosts.append(host)
    return hosts


def shutil_which(name):
    for directory in os.environ.get("PATH", "").split(os.pathsep):
        candidate = Path(directory) / name
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)
    return None


class NodeShell:
    """Runs programs on the other nodes of the job.  pbsdsh (PBS 23: `-n <vnode index>`, the
    task's output routed to the JOB's stdout, no capture option) with the vnode index of every
    host discovered by one pbsdsh pass over all vnodes (the tasks' PBS_NODENUM + hostname
    appended to a map file); the output of a command is read back from a file on the shared
    job directory.  ssh in batch mode is the fallback when no pbsdsh exists."""

    def __init__(self, plan, base, hosts):
        pbsdsh = plan.get("PBSDsh")
        self.pbsdsh = pbsdsh if pbsdsh and Path(pbsdsh).exists() else shutil_which("pbsdsh")
        self.base = Path(base)
        self.hosts = list(hosts)
        self.index = {}
        if self.pbsdsh:
            self.index = self.discover()

    def discover(self):
        map_path = self.base / "node-map.txt"
        if map_path.exists():
            map_path.unlink()
        command = [self.pbsdsh, "--", "/bin/sh", "-c", f'echo "$PBS_NODENUM $(hostname)" >> "{map_path}"']
        result = subprocess.run(command, text=True, capture_output=True, timeout=300)
        if result.returncode != 0 or not map_path.exists():
            raise SystemExit(f"pbsdsh node discovery failed: rc {result.returncode} {result.stderr.strip()[:400]}")
        index = {}
        for line in map_path.read_text().splitlines():
            parts = line.split()
            if len(parts) == 2 and parts[0].isdigit():
                host = normalize_host(parts[1])
                index[host] = min(index.get(host, int(parts[0])), int(parts[0]))
        missing = [host for host in self.hosts if host not in index]
        if missing:
            raise SystemExit(f"pbsdsh node discovery found no vnode for {missing} (map {index})")
        return index

    def prefix(self, host):
        if self.pbsdsh:
            return [self.pbsdsh, "-n", str(self.index[host]), "--"]
        return ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=20", host]

    def run(self, host, script, output_path, timeout=300):
        """Run `script` (a /bin/sh command line) on `host`; its stdout + stderr land in
        `output_path` (through the shared file system); returns the output text."""
        output_path = Path(output_path)
        if output_path.exists():
            output_path.unlink()
        command = self.prefix(host) + ["/bin/sh", "-c", f'({script}) > "{output_path}" 2>&1']
        result = subprocess.run(command, text=True, capture_output=True, timeout=timeout)
        if result.returncode != 0 or not output_path.exists():
            raise SystemExit(f"node shell failed on {host}: rc {result.returncode} {result.stderr.strip()[:400]}")
        return output_path.read_text()

    def spawn(self, host, script, log_path):
        """Start `script` on `host` in the background (its output into `log_path`)."""
        command = self.prefix(host) + ["/bin/sh", "-c", f'({script}) > "{log_path}" 2>&1']
        return subprocess.Popen(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def remote_node_preflight(shell, host, base):
    """/proc/meminfo and the process table of another node (its admission and conflict
    checks); the node shell must land on `host` (its hostname is checked)."""
    text = shell.run(host, "hostname; echo ---MEM---; cat /proc/meminfo; echo ---PS---; ps -eo pid=,comm=",
                     Path(base) / f"node-preflight-{host}.txt")
    if "---PS---" not in text:
        raise SystemExit(f"Node preflight on {host} returned no process table: {text[:400]!r}")
    reported, _, rest = text.partition("---MEM---")
    if normalize_host(reported.strip()) != host:
        raise SystemExit(f"Node shell for {host} landed on {reported.strip()!r}")
    meminfo_text, _, processes = rest.partition("---PS---")
    values = {}
    for line in meminfo_text.splitlines():
        key, _, rest_ = line.partition(":")
        if rest_.split() and rest_.split()[0].isdigit():
            values[key.strip()] = int(rest_.split()[0]) * 1024
    if "MemTotal" not in values or "MemAvailable" not in values:
        raise SystemExit(f"Node preflight on {host} returned no MemTotal / MemAvailable: {meminfo_text[:400]!r}")
    return values, processes


def conflicting_processes(processes_text):
    return [line for line in processes_text.splitlines() if line.split()
            and (line.split()[-1].startswith(("palace", "bridge-")) or line.split()[-1] in ("mpirun", "mpiexec", "prterun", "prted"))]


def per_node_peaks(sample_paths, start_epoch, end_epoch):
    """host -> the largest sampled used bytes between the two epochs (the per-node peaks of a
    stage, read from the node samplers' CSVs)."""
    peaks = {}
    for path in sample_paths:
        if not Path(path).exists():
            continue
        for line in Path(path).read_text().splitlines()[1:]:
            parts = line.split(",")
            if len(parts) != 3:
                continue
            stamp, host, used = int(parts[0]), parts[1], int(parts[2])
            if start_epoch <= stamp <= end_epoch:
                peaks[host] = max(peaks.get(host, 0), used)
    return peaks


def parse_log(text):
    parsed = {}
    match = re.search(PATTERNS["H1"], text)
    if match:
        parsed["Order"] = int(match.group(1))
        parsed["H1"] = int(match.group(2).replace(",", ""))
        parsed["ND"] = int(match.group(3).replace(",", ""))
        parsed["RT"] = int(match.group(4).replace(",", ""))
    match = re.search(PATTERNS["Elements"], text, re.M)
    if match:
        parsed["Elements"] = int(match.group(1))
    parsed["Iterations"] = [{"Step": int(a), "Of": int(b), "Index": int(c), "ElapsedAtStart": float(d)}
                            for a, b, c, d in re.findall(PATTERNS["It"], text)]
    parsed["PCG"] = [int(x) for x in re.findall(PATTERNS["PCG"], text)]
    parsed["Nonconvergence"] = [line for line in text.splitlines() if re.search(PATTERNS["NotConverged"], line, re.I)]
    parsed["SourceTiming"] = [{"Index": int(a), "Iterations": int(b), "SolveSeconds": float(c), "TotalSeconds": float(d)}
                              for a, b, c, d in re.findall(PATTERNS["SourceTiming"], text)]
    parsed["PairsProgress"] = [[int(a), int(b)] for a, b in re.findall(PATTERNS["Pairs"], text)]
    parsed["BlockPairsProgress"] = [[int(a), int(b)] for a, b in re.findall(PATTERNS["Blocks"], text)]
    # The decision-62(4) streaming reducer: one pass over the sources, per-rank sample counts.
    parsed["StreamedSourcesProgress"] = [[int(a), int(b)] for a, b in re.findall(PATTERNS["StreamedSources"], text)]
    match = re.search(PATTERNS["ReductionSamples"], text)
    if match:
        parsed["ReductionSamples"] = {"Sources": int(match.group(1)), "Interfaces": int(match.group(2)),
                                      "PerRankMin": int(match.group(3)), "PerRankMax": int(match.group(4)),
                                      "Total": int(match.group(5))}
    match = re.search(PATTERNS["PeakTotal"], text)
    if match:
        parsed["PalacePeakMemory"] = {"Min": match.group(1), "Max": match.group(2), "Avg": match.group(3), "Total": match.group(4)}
    match = re.search(PATTERNS["Total"], text, re.M)
    if match:
        parsed["PalaceTotalSeconds"] = float(match.group(2))
    report = re.search(r"Elapsed Time Report \(s\).*?^-+\nTotal.*?$", text, re.S | re.M)
    if report:
        parsed["ElapsedTimeReport"] = report.group(0)
    memory = re.search(r"Peak Memory .*?^-+\nTotal.*?$", text, re.S | re.M)
    if memory:
        parsed["PeakMemoryReport"] = memory.group(0)
    return parsed


def main(argv):
    plan_path = Path(argv[1]).resolve()
    plan = json.loads(plan_path.read_text())
    base = plan_path.parent
    status_path = base / "status.json"
    if status_path.exists():
        raise SystemExit("status.json exists: refusing to reuse this plan directory")
    binary, binary_sha256, mpiexec = Path(plan["Binary"]), plan["BinarySHA256"], plan["MPIExec"]
    job_start = time.monotonic()
    deadline = job_start + plan["DeadlineSeconds"]
    status = {"Version": 2, "Plan": str(plan_path), "Case": plan.get("Case"), "PBSJobID": os.environ.get("PBS_JOBID"),
              "Host": socket.gethostname().split(".")[0], "StartUTC": time.strftime("%FT%TZ", time.gmtime()),
              "State": "running", "Stages": []}

    def save():
        tmp = status_path.with_suffix(".tmp")
        tmp.write_text(json.dumps(status, indent=2) + "\n")
        tmp.replace(status_path)

    # Preflight: node identity, conflicting processes, executable and pinned hashes, admission
    # (a multi-node plan: every node of PBS_NODEFILE, decision 457).
    nodes = int(plan.get("Nodes", 1))
    ordered_hosts = unique_hosts(Path(os.environ["PBS_NODEFILE"]).read_text())
    hosts = set(ordered_hosts)
    if len(hosts) != nodes or status["Host"] not in hosts:
        raise SystemExit(f"PBS node identity mismatch actual={status['Host']} allocated={sorted(hosts)} plan nodes={nodes}")
    processes = subprocess.check_output(["ps", "-eo", "pid=,comm="], text=True)
    conflicts = conflicting_processes(processes)
    if conflicts:
        raise SystemExit("Conflicting native processes: " + repr(conflicts))
    memory0 = meminfo()
    preflight = {"MemTotalBytes": memory0["MemTotal"], "MemAvailableBytes": memory0["MemAvailable"],
                 "MinimumMemAvailableBytes": plan["MinimumMemAvailableBytes"], "Instance": plan.get("Instance"),
                 "Binary": str(binary), "BinarySHA256": sha(binary), "Pinned": {}, "UTC": time.strftime("%FT%TZ", time.gmtime())}
    other_hosts = [host for host in ordered_hosts if host != status["Host"]]
    node_shell = None
    if nodes > 1:
        node_shell = NodeShell(plan, base, other_hosts)
        preflight["Nodes"] = {status["Host"]: {"MemTotalBytes": memory0["MemTotal"], "MemAvailableBytes": memory0["MemAvailable"],
                                               "Conflicts": [], "Role": "runner"}}
        preflight["NodeShell"] = {"PBSDsh": node_shell.pbsdsh, "VnodeIndex": node_shell.index}
        for host in other_hosts:
            values, remote_processes = remote_node_preflight(node_shell, host, base)
            remote_conflicts = conflicting_processes(remote_processes)
            preflight["Nodes"][host] = {"MemTotalBytes": values["MemTotal"], "MemAvailableBytes": values["MemAvailable"],
                                        "Conflicts": remote_conflicts, "Role": "node"}
            if remote_conflicts:
                raise SystemExit(f"Conflicting native processes on {host}: " + repr(remote_conflicts))
        preflight["NodeGuard"] = plan.get("NodeGuard")
        status["Nodes"] = ordered_hosts
    if preflight["BinarySHA256"] != binary_sha256:
        raise SystemExit("Executable hash mismatch: " + preflight["BinarySHA256"])
    # A stage may name its own frozen executable (an executable comparison on one archive,
    # decision 62(4)): verified here like the plan's; the stage record carries it.
    preflight["StageBinaries"] = {}
    for stage in plan["Stages"]:
        if stage.get("Binary"):
            actual = sha(stage["Binary"])
            preflight["StageBinaries"][stage["Binary"]] = {"Expected": stage["BinarySHA256"], "Actual": actual}
            if actual != stage["BinarySHA256"]:
                raise SystemExit(f"Stage executable hash mismatch: {stage['Name']} {actual}")
    for path, expected in plan["PinnedSHA256"].items():
        actual = sha(path)
        preflight["Pinned"][path] = {"Expected": expected, "Actual": actual, "OK": actual == expected}
        if actual != expected:
            raise SystemExit(f"Pinned input mismatch: {path} {actual} != {expected}")
    if memory0["MemAvailable"] < plan["MinimumMemAvailableBytes"]:
        raise SystemExit(f"Admission failed: MemAvailable={memory0['MemAvailable']} < the plan's MinimumMemAvailableBytes "
                         f"{plan['MinimumMemAvailableBytes']} ({(plan.get('Instance') or {}).get('Type')})")
    for host, node in (preflight.get("Nodes") or {}).items():
        if node["MemAvailableBytes"] < plan["MinimumMemAvailableBytes"]:
            raise SystemExit(f"Admission failed on {host}: MemAvailable={node['MemAvailableBytes']} < the plan's "
                             f"MinimumMemAvailableBytes {plan['MinimumMemAvailableBytes']} (checked per node x {nodes} nodes)")
    for stage in plan["Stages"]:
        output = Path(json.loads(Path(stage["Config"]).read_text())["Problem"]["Output"])
        if output.exists():
            raise SystemExit(f"Output reuse forbidden: {output}")
    status["Preflight"] = preflight
    save()

    # Memory sampler (node-wide used = MemTotal - MemAvailable) every 5 s.
    samples_path = base / "memory-samples.csv"
    stop_sampling = threading.Event()
    current_stage = ["preflight"]
    stage_peaks = {}

    def sampler():
        with samples_path.open("w") as stream:
            stream.write("unix,stage,used_bytes\n")
            while not stop_sampling.is_set():
                values = meminfo()
                used = values["MemTotal"] - values["MemAvailable"]
                stage_peaks[current_stage[0]] = max(stage_peaks.get(current_stage[0], 0), used)
                stream.write(f"{time.time():.0f},{current_stage[0]},{used}\n")
                stream.flush()
                stop_sampling.wait(5.0)

    threading.Thread(target=sampler, daemon=True).start()
    # The other nodes' samplers (this script in --node-sampler mode through the node shell),
    # stopped by the stop file at the end; their CSVs give every stage's per-node peaks.
    node_samplers, node_sample_paths = [], []
    stop_file = base / "node-sampler.stop"
    if nodes > 1:
        if stop_file.exists():
            stop_file.unlink()
        for host in other_hosts:
            csv_path = base / f"memory-samples-{host}.csv"
            node_sample_paths.append(csv_path)
            script = f'/usr/bin/env python3 "{Path(__file__).resolve()}" --node-sampler "{csv_path}" "{stop_file}"'
            node_samplers.append(node_shell.spawn(host, script, base / f"node-sampler-{host}.log"))

    completed = set()
    for stage in plan["Stages"]:
        name = stage["Name"]
        record = {"Name": name, "Config": stage["Config"], "Environment": stage["Environment"]}
        if any(required not in completed for required in stage.get("Requires", [])):
            record["State"] = "skipped-dependency"
            status["Stages"].append(record)
            save()
            continue
        remaining = deadline - time.monotonic() - plan["DeadlineMarginSeconds"]
        cap = min(stage["CapSeconds"], remaining)
        if cap < stage["MinimumSeconds"]:
            record.update(State="skipped-no-time", RemainingSeconds=remaining)
            status["Stages"].append(record)
            save()
            continue
        environment = dict(os.environ)
        for key in list(environment):
            if key.startswith("PALACE_RESPONSE_"):
                environment.pop(key)
        environment.update(stage["Environment"])
        exports = [part for key in stage["Environment"] for part in ("-x", key)]
        log = base / f"{name}.log"
        timefile = base / f"{name}.time"
        stage_binary = Path(stage["Binary"]) if stage.get("Binary") else binary
        node_arguments = [os.path.expandvars(argument) for argument in plan.get("MPIExecArguments", [])]
        command = ["timeout", "-k", "30", str(int(cap)), "/usr/bin/time", "-v", "-o", str(timefile),
                   str(mpiexec), *node_arguments, *exports, "-n", str(plan["Ranks"]), str(stage_binary), stage["Config"]]
        record.update(Command=command, CapSeconds=cap, Log=str(log), StartUTC=time.strftime("%FT%TZ", time.gmtime()),
                      Binary=str(stage_binary), BinarySHA256=stage.get("BinarySHA256") or binary_sha256)
        current_stage[0] = name
        started = time.monotonic()
        started_epoch = time.time()
        with log.open("w") as stream:
            process = subprocess.Popen(command, env=environment, stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
            try:
                code = process.wait()
            except KeyboardInterrupt:
                os.killpg(process.pid, signal.SIGTERM)
                code = process.wait()
        wall = time.monotonic() - started
        parsed = parse_log(log.read_text(errors="replace"))
        time_text = timefile.read_text() if timefile.exists() else ""
        match = re.search(r"Maximum resident set size \(kbytes\): (\d+)", time_text)
        record.update(ReturnCode=code, WallSeconds=wall, TimedOut=code == 124,
                      MaxSingleProcessRSSBytes=int(match.group(1)) * 1024 if match else None,
                      NodePeakUsedBytesSampled=stage_peaks.get(name), Parsed=parsed,
                      EndUTC=time.strftime("%FT%TZ", time.gmtime()))
        if nodes > 1:
            peaks = per_node_peaks(node_sample_paths, int(started_epoch), int(time.time()) + 1)
            peaks[status["Host"]] = stage_peaks.get(name, 0)
            record["NodePeakUsedBytesSampledPerNode"] = peaks
            record["NodePeakUsedBytesSampled"] = max(peaks.values()) if peaks else stage_peaks.get(name)
            record["Nodes"] = nodes
        ok = code == 0 and not parsed["Nonconvergence"]
        record["State"] = "complete" if ok else ("timed-out" if code == 124 else "failed")
        if ok:
            completed.add(name)
        status["Stages"].append(record)
        save()
        current_stage[0] = "between-stages"

    stop_sampling.set()
    time.sleep(0.2)
    if node_samplers:
        stop_file.write_text(time.strftime("%FT%TZ", time.gmtime()) + "\n")
        for process in node_samplers:
            try:
                process.wait(timeout=30)
            except subprocess.TimeoutExpired:
                process.kill()
    status.update(State="complete" if all(s["State"] == "complete" for s in status["Stages"]) else "incomplete",
                  EndUTC=time.strftime("%FT%TZ", time.gmtime()), TotalSeconds=time.monotonic() - job_start,
                  StagePeakUsedBytesSampled=stage_peaks)
    if nodes > 1:
        status["StagePeakUsedBytesSampledPerNode"] = {record["Name"]: record.get("NodePeakUsedBytesSampledPerNode")
                                                      for record in status["Stages"] if record.get("NodePeakUsedBytesSampledPerNode")}
    save()
    print(json.dumps({key: value for key, value in status.items() if key != "Stages"}, indent=2))
    for record in status["Stages"]:
        print(record["Name"], record["State"], record.get("WallSeconds"), record.get("Parsed", {}).get("H1"),
              record.get("Parsed", {}).get("PCG"))
    return 0 if status["State"] == "complete" else 1


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "--node-sampler":
        node_sampler(sys.argv[2], sys.argv[3])
        sys.exit(0)
    sys.exit(main(sys.argv))
