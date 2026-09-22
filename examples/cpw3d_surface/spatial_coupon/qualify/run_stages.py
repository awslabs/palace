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

usage: run_stages.py PLAN.json
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

    # Preflight: node identity, conflicting processes, executable and pinned hashes, admission.
    hosts = {normalize_host(x) for x in Path(os.environ["PBS_NODEFILE"]).read_text().split() if x.strip()}
    if len(hosts) != 1 or status["Host"] not in hosts:
        raise SystemExit(f"PBS node identity mismatch actual={status['Host']} allocated={hosts}")
    processes = subprocess.check_output(["ps", "-eo", "pid=,comm="], text=True)
    conflicts = [line for line in processes.splitlines()
                 if line.split()[-1].startswith(("palace", "bridge-")) or line.split()[-1] in ("mpirun", "mpiexec", "prterun", "prted")]
    if conflicts:
        raise SystemExit("Conflicting native processes: " + repr(conflicts))
    memory0 = meminfo()
    preflight = {"MemTotalBytes": memory0["MemTotal"], "MemAvailableBytes": memory0["MemAvailable"],
                 "Binary": str(binary), "BinarySHA256": sha(binary), "Pinned": {}, "UTC": time.strftime("%FT%TZ", time.gmtime())}
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
        raise SystemExit(f"Admission failed: MemAvailable={memory0['MemAvailable']}")
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
        command = ["timeout", "-k", "30", str(int(cap)), "/usr/bin/time", "-v", "-o", str(timefile),
                   str(mpiexec), *exports, "-n", str(plan["Ranks"]), str(stage_binary), stage["Config"]]
        record.update(Command=command, CapSeconds=cap, Log=str(log), StartUTC=time.strftime("%FT%TZ", time.gmtime()),
                      Binary=str(stage_binary), BinarySHA256=stage.get("BinarySHA256") or binary_sha256)
        current_stage[0] = name
        started = time.monotonic()
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
        ok = code == 0 and not parsed["Nonconvergence"]
        record["State"] = "complete" if ok else ("timed-out" if code == 124 else "failed")
        if ok:
            completed.add(name)
        status["Stages"].append(record)
        save()
        current_stage[0] = "between-stages"

    stop_sampling.set()
    time.sleep(0.2)
    status.update(State="complete" if all(s["State"] == "complete" for s in status["Stages"]) else "incomplete",
                  EndUTC=time.strftime("%FT%TZ", time.gmtime()), TotalSeconds=time.monotonic() - job_start,
                  StagePeakUsedBytesSampled=stage_peaks)
    save()
    print(json.dumps({key: value for key, value in status.items() if key != "Stages"}, indent=2))
    for record in status["Stages"]:
        print(record["Name"], record["State"], record.get("WallSeconds"), record.get("Parsed", {}).get("H1"),
              record.get("Parsed", {}).get("PCG"))
    return 0 if status["State"] == "complete" else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
