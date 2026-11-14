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

A form-X plan (a "Launcher" block; USER decision 1674 (2)) runs every stage INSIDE the
published container image on ONE node: the preflight pins the image file's sha256 and the
sha256 of the executable inside it (apptainer exec IMAGE sha256sum BINARY); each stage runs
`apptainer exec --bind ... IMAGE <image mpirun> <MPIExecArguments> -x ... -n R <image palace>
CONFIG` with the PBS_* variables dropped from its environment (the image's PRTE must not read
the host's node file) and the stage's PALACE_RESPONSE_* variables set and forwarded with -x.
More than one node, or a stage naming its own executable, is refused.  A plan without the
block runs exactly as before.

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


def first_slot_index(nodefile_text):
    """host -> the index of its first line in PBS_NODEFILE (the pbsdsh task-slot index of that
    host's first rank)."""
    index = {}
    for slot, entry in enumerate(entry for entry in nodefile_text.split() if entry.strip()):
        host = normalize_host(entry)
        index.setdefault(host, slot)
    return index


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
    """Runs programs on the other nodes of the job.  pbsdsh (PBS 23.06: `-n <index>` counts the
    job's TASK SLOTS in PBS_NODEFILE order - one per mpiprocs rank, so the second host of a
    192-rank-per-node job is index 192, not 1 (measured: `-n 1` lands on the first host while
    the tasks' PBS_NODENUM reports the vnode 0 / 1); the task's output is routed to the JOB's
    stdout with no capture option) with every host's index = its first line in PBS_NODEFILE,
    the landed hostname verified by every command; the output of a command is read back from
    a file on the shared job directory.  ssh in batch mode is the fallback without pbsdsh."""

    def __init__(self, plan, base, hosts, nodefile_text):
        pbsdsh = plan.get("PBSDsh")
        self.pbsdsh = pbsdsh if pbsdsh and Path(pbsdsh).exists() else shutil_which("pbsdsh")
        self.base = Path(base)
        self.hosts = list(hosts)
        self.index = first_slot_index(nodefile_text)
        missing = [host for host in self.hosts if host not in self.index]
        if missing:
            raise SystemExit(f"PBS_NODEFILE has no slot for {missing} (slots {self.index})")

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


LAUNCHER_MEMORY_STOP_FRACTION = 0.85  # the default bar: MEASURED node use above 0.85 x MemTotal stops the job
LAUNCHER_PASSED_PATTERN = re.compile(r"^(OMP_[A-Z_]+|[A-Z0-9_]+_NUM_THREADS|TMPDIR)$")


def launcher_preflight(launcher, plan, *, digest=None, run=None, environment=None):
    """Form X: the image file's sha256 = Launcher ImageSHA256; the sha256 of the executable inside
    it = the block's own pin ImageBinarySHA256 AND the plan's BinarySHA256; the plan's Binary /
    MPIExec = the image's.  Records the apptainer path and version and the environment the image
    sees under --cleanenv (image_environment).  Returns the preflight record; SystemExit otherwise."""
    digest = digest or sha
    run = run or (lambda command, env=None: subprocess.check_output(command, text=True, env=env))
    if (plan["Binary"], plan["MPIExec"]) != (launcher["ImageBinary"], launcher["ImageMPIExec"]):
        raise SystemExit(f"Form-X plan names Binary {plan['Binary']} / MPIExec {plan['MPIExec']}, not the image's "
                         f"{launcher['ImageBinary']} / {launcher['ImageMPIExec']}")
    image_sha256 = digest(launcher["Image"])
    if image_sha256 != launcher["ImageSHA256"]:
        raise SystemExit(f"Image hash mismatch: {launcher['Image']} {image_sha256} != {launcher['ImageSHA256']}")
    output = run([launcher["Apptainer"], "exec", launcher["Image"], "sha256sum", launcher["ImageBinary"]])
    inner, _, path = output.strip().partition("  ")
    if path.strip() != launcher["ImageBinary"]:
        raise SystemExit(f"No sha256 of {launcher['ImageBinary']} inside the image: {output.strip()[:200]!r}")
    if inner.strip() != launcher["ImageBinarySHA256"]:
        raise SystemExit(f"Executable hash mismatch inside the image: {inner.strip()} != the Launcher pin "
                         f"{launcher['ImageBinarySHA256']}")
    record = {"Image": launcher["Image"], "ImageSHA256": image_sha256, "Apptainer": launcher["Apptainer"],
              "BinarySHA256": inner.strip(),
              "ApptainerVersion": run([launcher["Apptainer"], "--version"]).strip()}
    if environment is not None:
        probe = launcher_launch_prefix(launcher) + ["env"]
        record["ImageEnvironment"] = dict(line.split("=", 1) for line in
                                          run(probe, env=launcher_environment(environment, launcher)).splitlines()
                                          if "=" in line)
    return record


def launcher_passed(environment, launcher):
    """Form X: the variables the image receives under --cleanenv: the Launcher's Environment, the
    stage's PALACE_RESPONSE_* variables and the thread-count / TMPDIR variables of the job."""
    passed = {key: value for key, value in environment.items()
              if key.startswith("PALACE_RESPONSE_") or LAUNCHER_PASSED_PATTERN.match(key)}
    passed.update(launcher["Environment"])
    return passed


def launcher_environment(environment, launcher):
    """Form X: apptainer's own environment (the job's without PBS_*) plus APPTAINERENV_<name> for every
    variable of launcher_passed (apptainer --cleanenv passes only those into the image)."""
    stage_environment = {key: value for key, value in environment.items()
                         if not key.startswith(("PBS_", "APPTAINERENV_", "SINGULARITYENV_"))}
    stage_environment.update(launcher["Environment"])
    stage_environment.update({f"APPTAINERENV_{key}": value for key, value in launcher_passed(environment, launcher).items()})
    return stage_environment


def launcher_launch_prefix(launcher):
    """apptainer exec --cleanenv --no-home, the read-write Binds and the ReadOnlyBinds (':ro'), IMAGE."""
    binds = [part for path in launcher["Binds"] for part in ("--bind", path)]
    binds += [part for path in launcher.get("ReadOnlyBinds", []) for part in ("--bind", f"{path}:{path}:ro")]
    return [launcher["Apptainer"], "exec", "--cleanenv", "--no-home", *binds, launcher["Image"]]


def launcher_launch(launcher, *, exports, ranks, arguments):
    """Form X: launcher_launch_prefix <image MPI launcher> <MPIExecArguments> -x <Launcher
    Environment> <exports> -n RANKS <image executable> ARGUMENTS."""
    forwarded = [part for key in launcher["Environment"] for part in ("-x", key)]
    return [*launcher_launch_prefix(launcher), launcher["ImageMPIExec"], *launcher["MPIExecArguments"],
            *forwarded, *exports, "-n", str(ranks), launcher["ImageBinary"], *arguments]


def launcher_command(launcher, *, cap, timefile, exports, ranks, config):
    """Form X: the capped stage command (timeout + /usr/bin/time -v around launcher_launch)."""
    return ["timeout", "-k", "30", str(int(cap)), "/usr/bin/time", "-v", "-o", str(timefile),
            *launcher_launch(launcher, exports=exports, ranks=ranks, arguments=[config])]


def memory_stop_bar(launcher, mem_total_bytes):
    """Form X: the memory-stop bar in bytes (Launcher MemoryStopFraction, default 0.85, x MemTotal)."""
    fraction = float(launcher.get("MemoryStopFraction", LAUNCHER_MEMORY_STOP_FRACTION))
    if not 0.0 < fraction < 1.0:
        raise SystemExit(f"Launcher MemoryStopFraction {fraction} must lie in (0, 1)")
    return int(fraction * mem_total_bytes)


def memory_stop_index(used_bytes, bar_bytes):
    """The index of the first MEASURED node-use sample above the bar (None: no stop).  The live
    sampler and the replay of recorded sample blocks (HANDBOOK s18) use this one rule."""
    for index, used in enumerate(used_bytes):
        if used > bar_bytes:
            return index
    return None


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
    launcher = plan.get("Launcher")
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
    image = None
    if launcher is not None:
        if nodes > 1:
            raise SystemExit(f"Form-X plan with {nodes} nodes: the Launcher runs on one node only")
        if any(stage.get("Binary") for stage in plan["Stages"]):
            raise SystemExit("Form-X plan with a stage naming its own executable: not supported")
        image = launcher_preflight(launcher, plan, environment={key: value for key, value in os.environ.items()
                                                                 if not key.startswith("PALACE_RESPONSE_")})
    memory0 = meminfo()
    preflight = {"MemTotalBytes": memory0["MemTotal"], "MemAvailableBytes": memory0["MemAvailable"],
                 "MinimumMemAvailableBytes": plan["MinimumMemAvailableBytes"], "Instance": plan.get("Instance"),
                 "Binary": str(binary), "BinarySHA256": sha(binary) if image is None else image["BinarySHA256"], "Pinned": {},
                 "UTC": time.strftime("%FT%TZ", time.gmtime())}
    if image is not None:
        preflight["Launcher"] = image
        preflight["MemoryStopBarBytes"] = memory_stop_bar(launcher, memory0["MemTotal"])
    other_hosts = [host for host in ordered_hosts if host != status["Host"]]
    node_shell = None
    if nodes > 1:
        node_shell = NodeShell(plan, base, other_hosts, Path(os.environ["PBS_NODEFILE"]).read_text())
        preflight["Nodes"] = {status["Host"]: {"MemTotalBytes": memory0["MemTotal"], "MemAvailableBytes": memory0["MemAvailable"],
                                               "Conflicts": [], "Role": "runner"}}
        preflight["NodeShell"] = {"PBSDsh": node_shell.pbsdsh, "SlotIndex": node_shell.index}
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
    # Form X only (E-3 (a), decision 1702): the memory STOP on the MEASURED node use of these samples.
    current_process, memory_stop, stage_samples = [None], {}, []
    memory_bar = preflight.get("MemoryStopBarBytes")

    def sampler():
        with samples_path.open("w") as stream:
            stream.write("unix,stage,used_bytes\n")
            while not stop_sampling.is_set():
                values = meminfo()
                used = values["MemTotal"] - values["MemAvailable"]
                stage_peaks[current_stage[0]] = max(stage_peaks.get(current_stage[0], 0), used)
                stream.write(f"{time.time():.0f},{current_stage[0]},{used}\n")
                stream.flush()
                if memory_bar is not None:
                    stage_samples.append((round(time.time()), current_stage[0], used))
                    if not memory_stop and memory_stop_index([used], memory_bar) is not None:
                        memory_stop.update(Stage=current_stage[0], UsedBytes=used, BarBytes=memory_bar,
                                           MemTotalBytes=values["MemTotal"], UTC=time.strftime("%FT%TZ", time.gmtime()),
                                           Samples=[list(sample) for sample in stage_samples])
                        (base / "memory-stop.json").write_text(json.dumps(memory_stop, indent=1) + "\n")
                        process = current_process[0]
                        if process is not None and process.poll() is None:
                            os.killpg(process.pid, signal.SIGTERM)
                            if stop_sampling.wait(10.0) is False and process.poll() is None:
                                os.killpg(process.pid, signal.SIGKILL)
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
        if memory_stop:
            break
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
        if launcher is not None:
            environment = launcher_environment(environment, launcher)
        exports = [part for key in stage["Environment"] for part in ("-x", key)]
        log = base / f"{name}.log"
        timefile = base / f"{name}.time"
        stage_binary = Path(stage["Binary"]) if stage.get("Binary") else binary
        node_arguments = [os.path.expandvars(argument) for argument in plan.get("MPIExecArguments", [])]
        if launcher is None:
            command = ["timeout", "-k", "30", str(int(cap)), "/usr/bin/time", "-v", "-o", str(timefile),
                       str(mpiexec), *node_arguments, *exports, "-n", str(plan["Ranks"]), str(stage_binary), stage["Config"]]
        else:
            command = launcher_command(launcher, cap=cap, timefile=timefile, exports=exports, ranks=plan["Ranks"],
                                       config=stage["Config"])
        record.update(Command=command, CapSeconds=cap, Log=str(log), StartUTC=time.strftime("%FT%TZ", time.gmtime()),
                      Binary=str(stage_binary), BinarySHA256=stage.get("BinarySHA256") or binary_sha256)
        current_stage[0] = name
        started = time.monotonic()
        started_epoch = time.time()
        with log.open("w") as stream:
            process = subprocess.Popen(command, env=environment, stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
            current_process[0] = process
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
        current_process[0] = None
        if memory_stop:
            # No partial output of a stopped job is a result: this stage and the job end here.
            ok = False
            record["State"] = "memory-stop"
            record["MemoryStop"] = {key: value for key, value in memory_stop.items() if key != "Samples"}
        if ok:
            completed.add(name)
        status["Stages"].append(record)
        save()
        current_stage[0] = "between-stages"
        if memory_stop:
            break

    stop_sampling.set()
    time.sleep(0.2)
    if node_samplers:
        stop_file.write_text(time.strftime("%FT%TZ", time.gmtime()) + "\n")
        for process in node_samplers:
            try:
                process.wait(timeout=30)
            except subprocess.TimeoutExpired:
                process.kill()
    state = "complete" if all(s["State"] == "complete" for s in status["Stages"]) else "incomplete"
    status.update(State="memory-stop" if memory_stop else state,
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
    return 0 if status["State"] == "complete" else (3 if memory_stop else 1)


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "--node-sampler":
        node_sampler(sys.argv[2], sys.argv[3])
        sys.exit(0)
    sys.exit(main(sys.argv))
