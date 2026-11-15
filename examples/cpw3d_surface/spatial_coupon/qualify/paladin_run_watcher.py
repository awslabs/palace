# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""The unattended run watcher of form Y-b (decision 1820 (4), review r1799 MAJOR-4; decision 1822 (c)):
the wall cap and the memory bar of ONE paladin run, whatever its driver does.  Paladin has no
walltime, so without it a stage whose driver has ended would hold its node up to the PBS ceiling
(120 h).  paladin_backend.PaladinJob starts one per run (nohup, its own session) at the submit and
for an adopted run that has none; it runs on the driver's host, never on a compute node.

Every --interval seconds: `paladin palace status RUN` (status, terminal, job id) and one read of
--host (the cluster's ssh host): `qstat -xf -F json JOB` (job state, start time, exec hosts) and,
while the job runs, each exec host's /proc/meminfo (the logic of ctr_memwatch_remote_v1.sh).
  CAP:    the job running and now - start >= --cap-seconds (start = the PBS start time, else the
          submit time; a queued job has no cap clock) -> the run is ended early;
  MEMORY: the job running and MemTotal - MemAvailable >= --bar x MemTotal on ANY exec host,
          MEASURED in that read (never projected) -> the run is ended early.
  Missing reads never end a run.  The watcher ends (exit 0) when the row is terminal or the PBS job
  has ended, and when its directory is removed.
AN EARLY END = `paladin palace cancel RUN`, then the row must read terminal within --grace seconds;
a failed cancel or a row still not terminal = `qdel JOB` through --host (decision 1822 (iii): this
does not depend on the paladin session), then the PBS job must end within --grace.  All of it is in
end.json and events.log of --dir; exit 11 (cap) / 10 (memory); 3 = another watcher holds --dir (the
driver names --dir by the run id: one watcher per run); 5 = the read bound (the cap + one day) reached.

  python3 paladin_run_watcher.py --run-id R --cli PALADIN --host HOST --pbs-bin DIR --cap-seconds S
      --submitted UTC --dir DIR [--job PBS] [--bar 0.85] [--interval 30] [--grace 120]
"""
import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paladin_backend  # noqa: E402  (the same qualify/ directory)

# One read of the cluster (bash on --host, the job id and the PBS bin directory as arguments).
SAMPLE_SCRIPT = r'''set -u
JOB="$1"; PBS_BIN="$2"
read -r STATE STIME HOSTS < <(timeout 30 "$PBS_BIN/qstat" -xf -F json "$JOB" 2>/dev/null | python3 -c '
import json, sys, time
try:
    jobs = json.load(sys.stdin).get("Jobs", {})
except Exception:
    jobs = {}
if not jobs:
    print("NOTFOUND -")
    sys.exit(0)
job = next(iter(jobs.values()))
try:
    start = str(int(time.mktime(time.strptime(job["stime"], "%a %b %d %H:%M:%S %Y"))))
except Exception:
    start = "-"
hosts = []
for part in (job.get("exec_host") or "").split("+"):
    host = part.split("/")[0].strip()
    if host and host not in hosts:
        hosts.append(host)
print(job.get("job_state", "UNKNOWN"), start, " ".join(hosts))
')
echo "STATE ${STATE:-NOTFOUND}"
echo "STIME ${STIME:--}"
[ "${STATE:-}" = "R" ] || exit 0
for h in $HOSTS; do
    out=$(timeout 20 ssh -n -o BatchMode=yes -o ConnectTimeout=8 -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
          -o LogLevel=ERROR "$h" "awk '/^MemTotal:/{t=\$2} /^MemAvailable:/{a=\$2} END{print t, a}' /proc/meminfo" 2>/dev/null)
    if [ $? -eq 0 ] && [ -n "$out" ]; then echo "MEM $h $out"; else echo "MEMERR $h"; fi
done
'''
ENDED_STATES = ("F", "X")
# Decision 1830 (3): a minimum sleep between reads and a hard bound on the reads (the cap + one day, at the interval).
WATCHER_MIN_INTERVAL_SECONDS, WATCHER_LIFE_BEYOND_CAP_SECONDS = 1.0, 86400.0


def utc():
    return paladin_backend.utc_now().strftime("%Y-%m-%dT%H:%M:%SZ")


class Watcher:
    def __init__(self, args):
        self.args = args
        self.host = None if args.host == "local" else args.host
        self.job = args.job or None
        self.start_epoch = None
        self.events = os.path.join(args.dir, "events.log")

    def event(self, line):
        with open(self.events, "a") as f:
            f.write(f"{utc()} {line}\n")

    def status(self):
        """The row's {status, terminal, job_id}, or None (a missing read)."""
        try:
            return paladin_backend.paladin_call(self.args.cli, ["palace", "status", str(self.args.run_id), "--json"],
                                                timeout=120)
        except paladin_backend.PaladinError as error:
            self.event(f"READERR palace status: {str(error)[:200]}")
            return None

    def sample(self):
        """One cluster read: {State, StartEpoch, Memory [(host, total kB, available kB)]}, or None."""
        if not self.job:
            return None
        try:
            result = paladin_backend.remote_shell(self.host, f"bash -s -- {self.job} {self.args.pbs_bin}", timeout=120,
                                                  stdin=SAMPLE_SCRIPT)
        except paladin_backend.PaladinError as error:
            self.event(f"READERR cluster: {str(error)[:200]}")
            return None
        fields = {}
        memory = []
        for line in result.stdout.splitlines():
            parts = line.split()
            if len(parts) >= 2 and parts[0] in ("STATE", "STIME"):
                fields[parts[0]] = parts[1]
            elif len(parts) == 4 and parts[0] == "MEM" and parts[2].isdigit() and parts[3].isdigit():
                memory.append((parts[1], int(parts[2]), int(parts[3])))
        if "STATE" not in fields:
            self.event(f"READERR cluster: no STATE line (rc {result.returncode})")
            return None
        start = fields.get("STIME", "-")
        return {"State": fields["STATE"], "StartEpoch": int(start) if start.isdigit() else None, "Memory": memory}

    def wait_until(self, done, seconds):
        deadline = time.time() + seconds
        while True:
            if done():
                return True
            if time.time() >= deadline:
                return False
            time.sleep(min(10.0, max(1.0, self.args.interval)))

    def end_early(self, kind, reason):
        """The early end of the run (module docstring); end.json; the exit code."""
        record = {"Kind": kind, "Reason": reason, "RunID": self.args.run_id, "Job": self.job, "UTC": utc()}
        self.event(f"{kind.upper()} {reason} -> paladin palace cancel {self.args.run_id}")
        try:
            paladin_backend.paladin_call(self.args.cli, ["palace", "cancel", str(self.args.run_id), "--json"], timeout=240)
            record["Cancel"] = {"OK": True}
        except paladin_backend.PaladinError as error:
            record["Cancel"] = {"OK": False, "Error": str(error)[:300]}
        terminal = record["Cancel"]["OK"] and self.wait_until(lambda: (self.status() or {}).get("terminal") is True,
                                                              self.args.grace)
        record["TerminalAfterCancel"] = terminal
        self.event(f"cancel {'ok' if record['Cancel']['OK'] else 'FAILED'}; row terminal within {self.args.grace:.0f} s: "
                   f"{terminal}")
        if not terminal:
            if self.job:
                command = f"{self.args.pbs_bin}/qdel {self.job}"
                try:
                    result = paladin_backend.remote_shell(self.host, command, timeout=120)
                    record["Fallback"] = {"Command": command, "ReturnCode": result.returncode,
                                          "Output": (result.stdout + result.stderr).strip()[:300]}
                except paladin_backend.PaladinError as error:
                    record["Fallback"] = {"Command": command, "Error": str(error)[:300]}
                ended = self.wait_until(lambda: (self.sample() or {}).get("State") in ENDED_STATES + ("NOTFOUND",),
                                        self.args.grace)
                record["Fallback"]["JobEnded"] = ended
                self.event(f"FALLBACK {command}: rc {record['Fallback'].get('ReturnCode')}; PBS job ended within "
                           f"{self.args.grace:.0f} s: {ended}")
            else:
                record["Fallback"] = {"Error": "no PBS job id known: no qdel possible"}
                self.event("FALLBACK impossible: no PBS job id known")
        paladin_backend.write_json_atomic(os.path.join(self.args.dir, "end.json"), record)
        code = 11 if kind == "cap" else 10
        self.event(f"END exit {code}")
        return code

    def run(self):
        args = self.args
        submitted = paladin_backend.parse_time(args.submitted).timestamp()
        self.event(f"START run={args.run_id} job={self.job} cap={args.cap_seconds:.0f}s bar={args.bar} "
                   f"interval={args.interval:.0f}s host={args.host} pid={os.getpid()}")
        interval = max(float(args.interval), WATCHER_MIN_INTERVAL_SECONDS)
        for _ in range(int((args.cap_seconds + WATCHER_LIFE_BEYOND_CAP_SECONDS) / interval) + 10):
            if not os.path.isdir(args.dir):
                return 0
            status = self.status()
            if status and status.get("terminal") is True:
                self.event(f"END exit 0: the row is terminal ({status.get('status')})")
                return 0
            if status and not self.job and status.get("job_id"):
                self.job = str(status["job_id"])
            sample = self.sample()
            state = sample["State"] if sample else None
            if sample and sample["StartEpoch"]:
                self.start_epoch = sample["StartEpoch"]
            if state in ENDED_STATES:
                self.event(f"END exit 0: the PBS job has ended (state {state}; row {(status or {}).get('status')})")
                return 0
            running = state == "R" or (state is None and (status or {}).get("status") == "running")
            clock = self.start_epoch or (submitted if running else None)
            if running and clock is not None and time.time() - clock >= args.cap_seconds:
                origin = "the PBS start" if self.start_epoch else "the submit time"
                return self.end_early("cap", f"{time.time() - clock:.0f} s from {origin} >= cap {args.cap_seconds:.0f} s")
            if state == "R":
                for host, total, available in sample["Memory"]:
                    used = (total - available) / total if total else 0.0
                    if used >= args.bar:
                        return self.end_early("memory", f"{host} used {used:.4f} of MemTotal >= bar {args.bar}")
            time.sleep(interval)
        self.event(f"END exit 5: the iteration bound (cap + {WATCHER_LIFE_BEYOND_CAP_SECONDS:.0f} s) reached with the run "
                   f"not terminal (decision 1830 (3))")
        return 5


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--run-id", type=int, required=True)
    parser.add_argument("--job", default="")
    parser.add_argument("--cli", required=True)
    parser.add_argument("--host", required=True, help="the cluster's ssh host; 'local' = this host (the tests' form)")
    parser.add_argument("--pbs-bin", required=True)
    parser.add_argument("--cap-seconds", type=float, required=True)
    parser.add_argument("--submitted", required=True)
    parser.add_argument("--dir", required=True)
    parser.add_argument("--bar", type=float, default=0.85)
    parser.add_argument("--interval", type=float, default=30.0)
    parser.add_argument("--grace", type=float, default=120.0)
    args = parser.parse_args(argv)
    os.makedirs(args.dir, exist_ok=True)
    watcher = Watcher(args)
    try:
        lock = paladin_backend.DriverLock(os.path.join(args.dir, "watcher.lock"), log=watcher.event).acquire()
    except paladin_backend.PaladinError as error:
        watcher.event(f"REFUSED: {str(error)[:200]}")
        return 3
    with open(os.path.join(args.dir, "watcher.pid"), "w") as f:
        f.write(f"{os.getpid()}\n")
    try:
        return watcher.run()
    finally:
        lock.release()


if __name__ == "__main__":
    sys.exit(main())
