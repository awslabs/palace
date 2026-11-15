# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""The unattended run watcher of form Y-b (decision 1820 (4), review r1799 MAJOR-4; decision 1822 (c)):
the wall cap and the memory bar of ONE paladin run, whatever its driver does.  Paladin has no
walltime, so without it a stage whose driver has ended would hold its node up to the PBS ceiling
(120 h).  paladin_backend.PaladinJob starts one per run (nohup, its own session) at the submit and
for an adopted run that has none; it runs on the driver's host, never on a compute node.

Every --interval seconds: `paladin palace status RUN` (status, terminal, job id) and one read of
--host (the cluster's ssh host): `qstat -xf -F json JOB` (job state, start time, Exit_status, exec
hosts) and, while the job runs, each exec host's /proc/meminfo (the logic of ctr_memwatch_remote_v1.sh).
READ_RULE (decision 1836, review r1833 MAJOR-1 (a) / (c)): an empty qstat answer (NOTFOUND), a failed
ssh or an answer without its STATE line is a MISSING PBS read, never 'the job has ended'.  The run has
ended only on a finished-job answer (job_state F / X or an Exit_status) or a terminal row.  While the
PBS read is missing the row's status says whether the run is running, the exec hosts of the last good
read are still read for memory, and the cap clock is the last PBS start seen, else the submit time;
when NEITHER read has answered for BLIND_SECONDS the cap is judged by that clock alone.  So the cap
and the memory bar stay in force under a missing read.
  CAP:    the run running (PBS R, or the row's 'running' under a missing PBS read; or no read for
          BLIND_SECONDS) and now - start >= --cap-seconds (start = the PBS start time, else the
          submit time; a queued job has no cap clock) -> the run is ended early;
  MEMORY: the run running and MemTotal - MemAvailable >= --bar x MemTotal on ANY exec host,
          MEASURED in that read (never projected) -> the run is ended early.
  Missing memory reads never end a run.  The watcher ends (exit 0) when the row is terminal or the
  PBS job has ended (a finished-job answer), and when its directory is removed.
END_RULE (decision 1836, review r1833 MAJOR-1 (b) / (c)): AN EARLY END is repeated until it is
VERIFIED (the row terminal or a finished-job answer; an empty qstat answer verifies nothing).  One
round = `paladin palace cancel RUN` (the end request); refused, or the run not verified within
--grace = `qdel JOB` through --host, then verified within --grace.  Not verified: the next round
after a back-off min(grace x 2^(k-1), END_BACKOFF_MAX_SECONDS), at most END_MAX_ROUNDS rounds.
THE SESSION AND THE QUEUE WAIT: the cap counts from the PBS START, the driver's session check
from the SUBMIT (paladin_backend.SESSION_RULE), so a run that waits in the queue longer than the
session's life meets an EXPIRED session at its cap: the end request is then refused and the qdel
by the PBS job id, which needs no paladin session, is the PRIMARY end path.  end.json (in --dir)
records every round and Verified true / false; an early end not verified after the last round
leaves END-UNVERIFIED.json beside it (the driver and every resume STOP on it, naming the run and
the manual check).  Exit 11 (cap) / 10 (memory) verified, 12 not verified; 3 = another watcher
holds --dir (the driver names --dir by the run id: one watcher per run); 5 = the life bound
reached: the PBS start + the cap + one day, or, before any start is seen, the submit time +
QUEUE_LIFE_SECONDS + the cap + one day (decision 1836 m-4 (b): a queued run keeps its watcher).
NONCE_RULE (decision 1836 m-1, review r1833 MINOR-1: the submit window): with --run-id 0 the watcher
is started BEFORE the submit call, keyed by the SubmitIntent's --nonce, in <local>/watchers/nonce-<N>.
Until it knows its run it only reads: every --interval the driver's hint file 'run-id' in --dir
(written as soon as the driver knows the run), else `paladin palace my-jobs -m <--identity> -m
nonce=N` (a failed or warned listing is a missing read).  A no-run row (submission_failed without a
job id) = that submit made no run: exit 0.  A run found: it takes the lock of the sibling run-<id>
directory (held by another watcher: exit 0, that one watches) and watches the run there as above.
No run within --search-seconds: exit 0.

  python3 paladin_run_watcher.py --run-id R --cli PALADIN --host HOST --pbs-bin DIR --cap-seconds S
      --submitted UTC --dir DIR [--job PBS] [--bar 0.85] [--interval 30] [--grace 120]
  python3 paladin_run_watcher.py --run-id 0 --nonce N --identity JSON [--search-seconds 3600] ... (NONCE_RULE)
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
JOB="$1"; PBS_BIN="$2"; shift 2; KEPT_HOSTS="$*"
read -r STATE STIME EXIT HOSTS < <(timeout 30 "$PBS_BIN/qstat" -xf -F json "$JOB" 2>/dev/null | python3 -c '
import json, sys, time
try:
    jobs = json.load(sys.stdin).get("Jobs", {})
except Exception:
    jobs = {}
if not jobs:
    print("NOTFOUND - -")
    sys.exit(0)
job = next(iter(jobs.values()))
try:
    start = str(int(time.mktime(time.strptime(job["stime"], "%a %b %d %H:%M:%S %Y"))))
except Exception:
    start = "-"
code = job.get("Exit_status")
hosts = []
for part in (job.get("exec_host") or "").split("+"):
    host = part.split("/")[0].strip()
    if host and host not in hosts:
        hosts.append(host)
print(job.get("job_state", "UNKNOWN"), start, "-" if code is None else code, " ".join(hosts))
')
echo "STATE ${STATE:-NOTFOUND}"
echo "STIME ${STIME:--}"
echo "EXIT ${EXIT:--}"
if [ "${STATE:-}" = "R" ]; then echo "HOSTS $HOSTS"
elif [ "${STATE:-NOTFOUND}" = "NOTFOUND" ] && [ -n "$KEPT_HOSTS" ]; then HOSTS="$KEPT_HOSTS"  # READ_RULE: the kept hosts
else exit 0
fi
for h in $HOSTS; do
    out=$(timeout 20 ssh -n -o BatchMode=yes -o ConnectTimeout=8 -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
          -o LogLevel=ERROR "$h" "awk '/^MemTotal:/{t=\$2} /^MemAvailable:/{a=\$2} END{print t, a}' /proc/meminfo" 2>/dev/null)
    if [ $? -eq 0 ] && [ -n "$out" ]; then echo "MEM $h $out"; else echo "MEMERR $h"; fi
done
'''
ENDED_STATES = ("F", "X")  # PBS's finished job / finished subjob
# Decision 1830 (3): a minimum sleep between reads and a hard bound on the reads (the cap + one day, at the interval).
WATCHER_MIN_INTERVAL_SECONDS, WATCHER_LIFE_BEYOND_CAP_SECONDS = 1.0, 86400.0
QUEUE_LIFE_SECONDS = 7 * 86400.0  # decision 1836 m-4 (b): the longest queue wait a watcher outlives
# READ_RULE / END_RULE (module docstring; decision 1836).
BLIND_SECONDS = 600.0
END_MAX_ROUNDS, END_BACKOFF_MAX_SECONDS = 8, 1800.0
UNVERIFIED_MARKER = "END-UNVERIFIED.json"
NONCE_SEARCH_SECONDS, RUN_HINT = 3600.0, "run-id"  # NONCE_RULE


def utc():
    return paladin_backend.utc_now().strftime("%Y-%m-%dT%H:%M:%SZ")


class Watcher:
    def __init__(self, args):
        self.args = args
        self.host = None if args.host == "local" else args.host
        self.job = args.job or None
        self.start_epoch = None
        self.hosts = []  # the exec hosts of the last good PBS read (READ_RULE)
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
        """One cluster read: {State, StartEpoch, ExitStatus, Hosts, Memory [(host, total kB, available kB)]},
        or None (no job id, a failed ssh, no STATE line).  State NOTFOUND = an empty qstat answer: a missing
        PBS read (READ_RULE), whose Memory comes from the kept exec hosts."""
        if not self.job:
            return None
        try:
            result = paladin_backend.remote_shell(self.host, " ".join(["bash", "-s", "--", self.job, self.args.pbs_bin,
                                                                       *self.hosts]),
                                                  timeout=120, stdin=SAMPLE_SCRIPT)
        except paladin_backend.PaladinError as error:
            self.event(f"READERR cluster: {str(error)[:200]}")
            return None
        fields, hosts, memory = {}, [], []
        for line in result.stdout.splitlines():
            parts = line.split()
            if len(parts) >= 2 and parts[0] in ("STATE", "STIME", "EXIT"):
                fields[parts[0]] = parts[1]
            elif parts[:1] == ["HOSTS"]:
                hosts = parts[1:]
            elif len(parts) == 4 and parts[0] == "MEM" and parts[2].isdigit() and parts[3].isdigit():
                memory.append((parts[1], int(parts[2]), int(parts[3])))
        if "STATE" not in fields:
            self.event(f"READERR cluster: no STATE line (rc {result.returncode})")
            return None
        start, code = fields.get("STIME", "-"), fields.get("EXIT", "-")
        return {"State": fields["STATE"], "StartEpoch": int(start) if start.isdigit() else None,
                "ExitStatus": int(code) if code.lstrip("-").isdigit() else None, "Hosts": hosts, "Memory": memory}

    @staticmethod
    def pbs_read(sample):
        """True for a PBS read that shows the job (READ_RULE): not None, not NOTFOUND, not without a job_state."""
        return sample is not None and sample["State"] not in ("NOTFOUND", "UNKNOWN")

    @staticmethod
    def finished(sample):
        """The finished-job answer of READ_RULE ('ended' only on it or a terminal row), or None."""
        if Watcher.pbs_read(sample) and (sample["State"] in ENDED_STATES or sample["ExitStatus"] is not None):
            return f"PBS job_state {sample['State']}, Exit_status {sample['ExitStatus']}"
        return None

    def verified(self):
        """END_RULE: what shows the run has ended (the row terminal, a finished-job answer), or None."""
        status = self.status()
        if status and status.get("terminal") is True:
            return f"the row terminal ({status.get('status')})"
        if status and not self.job and status.get("job_id"):
            self.job = str(status["job_id"])
        return self.finished(self.sample())

    def wait_verified(self, seconds):
        deadline = time.time() + seconds
        while True:
            shown = self.verified()
            if shown or time.time() >= deadline:
                return shown
            time.sleep(min(10.0, max(1.0, self.args.interval)))

    def end_early(self, kind, reason):
        """END_RULE: the end request, then the qdel by the job id, in rounds until VERIFIED; end.json; the exit code."""
        path = os.path.join(self.args.dir, "end.json")
        record = {"Kind": kind, "Reason": reason, "RunID": self.args.run_id, "Job": self.job, "UTC": utc(), "Rounds": [],
                  "Verified": False, "VerifiedBy": None, "Rule": "END_RULE of paladin_run_watcher.py (decision 1836)"}
        self.event(f"{kind.upper()} {reason} -> the end request (paladin palace cancel {self.args.run_id})")
        for round_number in range(1, END_MAX_ROUNDS + 1):
            entry = {"Round": round_number, "UTC": utc()}
            try:
                paladin_backend.paladin_call(self.args.cli, ["palace", "cancel", str(self.args.run_id), "--json"], timeout=240)
                entry["Cancel"] = {"OK": True}
            except paladin_backend.PaladinError as error:
                entry["Cancel"] = {"OK": False, "Error": str(error)[:300]}
            shown = self.wait_verified(self.args.grace) if entry["Cancel"]["OK"] else None
            self.event(f"round {round_number}: end request {'ok' if entry['Cancel']['OK'] else 'REFUSED'}; verified within "
                       f"{self.args.grace:.0f} s: {shown}")
            if not shown:
                if not self.job:
                    self.verified()  # the row's job id, when it answers
                if self.job:
                    command = f"{self.args.pbs_bin}/qdel {self.job}"
                    try:
                        result = paladin_backend.remote_shell(self.host, command, timeout=120)
                        entry["Qdel"] = {"Command": command, "ReturnCode": result.returncode,
                                         "Output": (result.stdout + result.stderr).strip()[:300]}
                    except paladin_backend.PaladinError as error:
                        entry["Qdel"] = {"Command": command, "Error": str(error)[:300]}
                    shown = self.wait_verified(self.args.grace)
                    self.event(f"round {round_number}: {command}: rc {entry['Qdel'].get('ReturnCode')}; verified within "
                               f"{self.args.grace:.0f} s: {shown}")
                else:
                    entry["Qdel"] = {"Error": "no PBS job id known: no qdel possible"}
                    self.event(f"round {round_number}: no qdel possible: no PBS job id known")
            entry["VerifiedBy"] = shown
            record["Rounds"].append(entry)
            record["Job"] = self.job
            if shown:
                record["Verified"], record["VerifiedBy"] = True, shown
                break
            paladin_backend.write_json_atomic(path, record)  # the rounds so far, on disk
            if round_number < END_MAX_ROUNDS:
                backoff = min(self.args.grace * 2 ** (round_number - 1), END_BACKOFF_MAX_SECONDS)
                self.event(f"round {round_number}: NOT VERIFIED; the next round in {backoff:.0f} s")
                time.sleep(backoff)
        paladin_backend.write_json_atomic(path, record)
        if not record["Verified"]:
            marker = {"RunID": self.args.run_id, "Job": self.job, "Kind": kind, "UTC": utc(), "Rounds": len(record["Rounds"]),
                      "Check": f"by hand: paladin palace status {self.args.run_id} --json; qstat -xf {self.job or '<job>'}; "
                               f"qdel {self.job or '<job>'}; then delete this file",
                      "Rule": "END_RULE of paladin_run_watcher.py (decision 1836)"}
            paladin_backend.write_json_atomic(os.path.join(self.args.dir, UNVERIFIED_MARKER), marker)
            self.event(f"END exit 12: the early end is NOT VERIFIED after {len(record['Rounds'])} rounds: {UNVERIFIED_MARKER} "
                       f"written; check run {self.args.run_id} by hand")
            return 12
        code = 11 if kind == "cap" else 10
        self.event(f"END exit {code}: verified ({record['VerifiedBy']})")
        return code

    def find_run(self):
        """NONCE_RULE: the lock of the run's own directory once the run is known (the watcher then watches it there),
        or None (ended: no run, a no-run row, another watcher, the directory removed)."""
        args = self.args
        identity = {**json.loads(args.identity), "nonce": args.nonce}
        self.event(f"START nonce={args.nonce} (the submit window) search={args.search_seconds:.0f}s pid={os.getpid()}")
        interval = max(float(args.interval), WATCHER_MIN_INTERVAL_SECONDS)
        deadline = time.time() + args.search_seconds
        run_id = None
        for _ in range(int(args.search_seconds / interval) + 10):
            if not os.path.isdir(args.dir):
                return None
            hint = os.path.join(args.dir, RUN_HINT)
            if os.path.isfile(hint):
                with open(hint) as f:
                    text = f.read().strip()
                if text.isdigit():
                    run_id, how = int(text), "the driver's hint"
                    break
            try:
                envelope = paladin_backend.paladin_envelope(args.cli, ["palace", "my-jobs",
                                                                        *paladin_backend.metadata_arguments(identity),
                                                                        "--json"], timeout=120)
                data = envelope.get("data")
                if envelope.get("warnings") or not isinstance(data, dict) or not isinstance(data.get("jobs"), list):
                    raise paladin_backend.PaladinError("a listing with warnings or without a 'jobs' list")
                rows = [row for row in data["jobs"] if (row.get("metadata") or {}).get("nonce") == args.nonce]
            except paladin_backend.PaladinError as error:
                self.event(f"READERR palace my-jobs: {str(error)[:200]}")
                rows = []
            runs = [row for row in rows if not paladin_backend.no_run(row)]
            if runs:
                row = min(runs, key=lambda item: int(item["run_id"]))
                run_id, how = int(row["run_id"]), "the lookup by the nonce"
                self.job = self.job or (str(row["job_id"]) if row.get("job_id") else None)
                break
            if rows:
                self.event("END exit 0: the nonce's submit made no run (a submission_failed row without a job id)")
                return None
            if time.time() >= deadline:
                break
            time.sleep(interval)
        if run_id is None:
            self.event(f"END exit 0: no run with this nonce in {args.search_seconds:.0f} s")
            return None
        directory = os.path.join(os.path.dirname(os.path.abspath(args.dir)), f"run-{run_id}")
        os.makedirs(directory, exist_ok=True)
        try:
            lock = paladin_backend.DriverLock(os.path.join(directory, "watcher.lock"), log=self.event).acquire()
        except paladin_backend.PaladinError:
            self.event(f"END exit 0: run {run_id} (found by {how}) is watched by another watcher ({directory} is held)")
            return None
        self.event(f"run {run_id} found by {how}: watching it from {directory}")
        with open(os.path.join(directory, "watcher.pid"), "w") as f:
            f.write(f"{os.getpid()}\n")
        args.run_id, args.dir = run_id, directory
        self.events = os.path.join(directory, "events.log")
        return lock

    def run(self):
        args = self.args
        submitted = paladin_backend.parse_time(args.submitted).timestamp()
        self.event(f"START run={args.run_id} job={self.job} cap={args.cap_seconds:.0f}s bar={args.bar} "
                   f"interval={args.interval:.0f}s host={args.host} pid={os.getpid()}")
        interval = max(float(args.interval), WATCHER_MIN_INTERVAL_SECONDS)
        last_read = time.time()
        for _ in range(int((QUEUE_LIFE_SECONDS + args.cap_seconds + WATCHER_LIFE_BEYOND_CAP_SECONDS) / interval) + 10):
            # The life bound counts from the PBS start (decision 1836 m-4 (b)).
            deadline = ((self.start_epoch + args.cap_seconds + WATCHER_LIFE_BEYOND_CAP_SECONDS) if self.start_epoch else
                        (submitted + QUEUE_LIFE_SECONDS + args.cap_seconds + WATCHER_LIFE_BEYOND_CAP_SECONDS))
            if time.time() >= deadline:
                break
            if not os.path.isdir(args.dir):
                return 0
            status = self.status()
            if status and status.get("terminal") is True:
                self.event(f"END exit 0: the row is terminal ({status.get('status')})")
                return 0
            if status and not self.job and status.get("job_id"):
                self.job = str(status["job_id"])
            sample = self.sample()
            finished = self.finished(sample)
            if finished:
                self.event(f"END exit 0: the PBS job has ended ({finished}; row {(status or {}).get('status')})")
                return 0
            now = time.time()
            if status or self.pbs_read(sample):
                last_read = now
            if self.pbs_read(sample):
                if sample["StartEpoch"]:
                    self.start_epoch = sample["StartEpoch"]
                if sample["Hosts"]:
                    self.hosts = sample["Hosts"]
                running, how = sample["State"] == "R", f"PBS {sample['State']}"
            else:
                running, how = (status or {}).get("status") == "running", "the row's status (a missing PBS read)"
                self.event(f"READ missing PBS read ({'an empty qstat answer' if sample else 'no answer'}): running from "
                           f"{how}: {running}")
            blind = now - last_read >= BLIND_SECONDS
            clock = self.start_epoch or submitted
            if (running or blind) and now - clock >= args.cap_seconds:
                origin = "the PBS start" if self.start_epoch else "the submit time"
                why = f"no read for {now - last_read:.0f} s: by the clock" if blind and not running else f"running by {how}"
                return self.end_early("cap", f"{now - clock:.0f} s from {origin} >= cap {args.cap_seconds:.0f} s ({why})")
            if running:
                for host, total, available in (sample or {}).get("Memory", []):
                    used = (total - available) / total if total else 0.0
                    if used >= args.bar:
                        return self.end_early("memory", f"{host} used {used:.4f} of MemTotal >= bar {args.bar}")
            time.sleep(interval)
        self.event(f"END exit 5: the life bound ({'the PBS start' if self.start_epoch else 'the submit time + the queue life'} "
                   f"+ cap + {WATCHER_LIFE_BEYOND_CAP_SECONDS:.0f} s) reached with the run not terminal (decision 1830 (3), "
                   f"1836 m-4 (b))")
        return 5


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--run-id", type=int, default=0, help="0 = the submit window: find the run by --nonce (NONCE_RULE)")
    parser.add_argument("--nonce", default="")
    parser.add_argument("--identity", default="{}", help="the stage identity (JSON) of the nonce lookup")
    parser.add_argument("--search-seconds", type=float, default=NONCE_SEARCH_SECONDS)
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
    if not args.run_id and not args.nonce:
        parser.error("--run-id 0 needs --nonce (NONCE_RULE)")
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
        if args.run_id:
            return watcher.run()
        run_lock = watcher.find_run()
        if run_lock is None:
            return 0
        try:
            return watcher.run()
        finally:
            run_lock.release()
    finally:
        lock.release()


if __name__ == "__main__":
    sys.exit(main())
