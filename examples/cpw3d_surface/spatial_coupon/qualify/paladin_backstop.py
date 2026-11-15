# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""The PBS-side BACKSTOP of form Y-b (decision 1858; review r1848 (c) 1): a qdel-ONLY wall limit for ONE paladin run that
runs ON the cluster's ssh host (setsid nohup there), so that an outage of the driver's host, where the run watcher
(paladin_run_watcher.py) lives, does not leave a stage running up to the PBS ceiling (120 h).  In the form of lane
cmesh-r's cmeshr_backstop_qdel_v1.sh (ruling 1685 (4)(a)), with the chain's read rules.  paladin_backend.PaladinJob
copies this file to <the plan's /data directory>/paladin-backstops/ (its sha256 read back), runs --check there before
every submit and starts one backstop per run after the submit (BACKSTOP_RULE of paladin_backend.py).  Standard library
only: nothing of the chain is imported on the cluster host.

BACKSTOP_RULE (this side):
  ONE JOB: the job id of its state file (JSON: RunID, JobID, JobName, Owner, LimitSeconds, GraceSeconds, IntervalSeconds,
    LifeUntilEpoch, PBSBin).  A read is of that job only: the qstat answer's Jobs entry whose key has the job id's numeric
    prefix.
  REFUSAL: qstat's Job_Owner (before '@') or Job_Name differs from the state file's -> exit 13, no action.
  THE LIMIT: when now >= the job's PBS start (stime) + LimitSeconds + GraceSeconds and the job is still in PBS -> qdel by
    the job id, repeated at most QDEL_MAX times QDEL_GAP_SECONDS apart until the job has left PBS.  No memory rule.
  LEFT PBS: job_state F / X, an Exit_status, or qstat's 'Unknown Job Id <this job>' answer -> exit 0 (it never fired) or
    11 (ended after its qdel).  A MISSING read (no answer, an empty answer, another job's entry, no job_state) is never
    'left': it reads again.
  LIFE: exit 5 at LifeUntilEpoch (the job not seen to leave); its state file removed (a stop by hand) -> exit 0.
  RECORDS beside the state file: backstop.pid, backstop.log (one line per event; FIRST-READ after the first read),
    backstop.end.json (Exit, Reason, Qdels, UTC).

  python3 paladin_backstop.py --check DIR --pbs-bin PBS_BIN      (prints 'BACKSTOP-CHECK OK user=U host=H' or FAIL)
  python3 paladin_backstop.py STATE_FILE
"""
import getpass
import json
import os
import socket
import subprocess
import sys
import time

ENDED_STATES = ("F", "X")
QDEL_MAX, QDEL_GAP_SECONDS = 3, 300.0
MIN_INTERVAL_SECONDS = 5.0
STATE_KEYS = ("RunID", "JobID", "JobName", "Owner", "LimitSeconds", "GraceSeconds", "IntervalSeconds", "LifeUntilEpoch",
              "PBSBin")


def utc(epoch):
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(epoch))


def parse_answer(job, stdout, stderr):
    """One `qstat -xf -F json JOB` answer: {State, Owner, Name, StartEpoch, ExitStatus}.  State NOTFOUND (no entry),
    OTHERJOB (entries of other jobs only), UNKNOWN (no job_state) = a missing read; GONE = 'Unknown Job Id <this job>'."""
    wanted = str(job).split(".")[0]
    missing = {"Owner": None, "Name": None, "StartEpoch": None, "ExitStatus": None}
    try:
        jobs = json.loads(stdout).get("Jobs", {})
    except (ValueError, AttributeError):
        jobs = {}
    if not isinstance(jobs, dict) or not jobs:
        gone = [line for line in (stderr or "").splitlines()
                if line.strip().startswith("qstat: Unknown Job Id ") and line.split()[-1].split(".")[0] == wanted]
        return {"State": "GONE" if gone else "NOTFOUND", **missing}
    mine = [value for key, value in jobs.items() if str(key).split(".")[0] == wanted]
    if not mine:
        return {"State": "OTHERJOB", **missing}
    entry = mine[0]
    try:
        start = int(time.mktime(time.strptime(entry["stime"], "%a %b %d %H:%M:%S %Y")))
    except (KeyError, TypeError, ValueError):
        start = None
    code = entry.get("Exit_status")
    return {"State": entry.get("job_state") or "UNKNOWN", "Owner": (entry.get("Job_Owner") or "").split("@")[0] or None,
            "Name": entry.get("Job_Name"), "StartEpoch": start, "ExitStatus": None if code is None else int(code)}


def has_left(read):
    """The text that shows the job has left PBS, or None."""
    if read["State"] == "GONE":
        return "PBS: Unknown Job Id"
    if read["State"] in ("NOTFOUND", "OTHERJOB", "UNKNOWN", "READERR"):
        return None
    if read["State"] in ENDED_STATES or read["ExitStatus"] is not None:
        return f"job_state {read['State']}, Exit_status {read['ExitStatus']}"
    return None


class Backstop:
    """The loop of BACKSTOP_RULE for one state file; `run` = subprocess.run and `clock` = time (injected by the tests)."""

    def __init__(self, state_path, *, run=subprocess.run, clock=time):
        self.state_path, self.run, self.clock = state_path, run, clock
        directory = os.path.dirname(os.path.abspath(state_path))
        self.log_path = os.path.join(directory, "backstop.log")
        self.end_path = os.path.join(directory, "backstop.end.json")
        self.pid_path = os.path.join(directory, "backstop.pid")
        self.qdels = []

    def event(self, line):
        with open(self.log_path, "a") as f:
            f.write(f"{utc(self.clock.time())} {line}\n")

    def end(self, code, reason):
        record = {"Exit": code, "Reason": reason, "Qdels": self.qdels, "UTC": utc(self.clock.time()), "PID": os.getpid()}
        tmp = self.end_path + ".tmp"
        with open(tmp, "w") as f:
            json.dump(record, f, indent=1)
        os.replace(tmp, self.end_path)
        self.event(f"END exit {code}: {reason}")
        return code

    def read(self, state):
        try:
            result = self.run([os.path.join(state["PBSBin"], "qstat"), "-xf", "-F", "json", str(state["JobID"])],
                              capture_output=True, text=True, timeout=30)
        except (OSError, subprocess.TimeoutExpired) as error:
            return {"State": "READERR", "Owner": None, "Name": None, "StartEpoch": None, "ExitStatus": None,
                    "Error": str(error)[:200]}
        return parse_answer(state["JobID"], result.stdout, result.stderr)

    def qdel(self, state):
        command = [os.path.join(state["PBSBin"], "qdel"), str(state["JobID"])]
        try:
            result = self.run(command, capture_output=True, text=True, timeout=60)
            entry = {"UTC": utc(self.clock.time()), "ReturnCode": result.returncode,
                     "Output": (result.stdout + result.stderr).strip()[:300]}
        except (OSError, subprocess.TimeoutExpired) as error:
            entry = {"UTC": utc(self.clock.time()), "Error": str(error)[:300]}
        self.qdels.append(entry)
        self.event(f"QDEL {' '.join(command)}: {entry}")

    def loop(self):
        with open(self.state_path) as f:
            state = json.load(f)
        missing = [key for key in STATE_KEYS if state.get(key) in (None, "")]
        if missing:
            return self.end(2, f"the state file lacks {missing}")
        with open(self.pid_path, "w") as f:
            f.write(f"{os.getpid()}\n")
        interval = max(float(state["IntervalSeconds"]), MIN_INTERVAL_SECONDS)
        limit, grace, life = float(state["LimitSeconds"]), float(state["GraceSeconds"]), float(state["LifeUntilEpoch"])
        self.event(f"START run={state['RunID']} job={state['JobID']} name={state['JobName']} owner={state['Owner']} "
                   f"limit={limit:.0f}s grace={grace:.0f}s interval={interval:.0f}s life-until={utc(life)} "
                   f"host={socket.gethostname()} pid={os.getpid()}")
        first = True
        for _ in range(int(max(life - self.clock.time(), 0.0) / interval) + 10):  # a bound on the reads, whatever the clock
            if not os.path.isfile(self.state_path):
                return self.end(0, "the state file is removed (a stop by hand)")
            now = self.clock.time()
            if now >= life:
                return self.end(5, f"the life bound {utc(life)} reached; the job was not seen to leave PBS")
            read = self.read(state)
            if first:
                self.event(f"FIRST-READ {read['State']}")
                first = False
            left = has_left(read)
            if left:
                return self.end(11 if self.qdels else 0, f"the job has left PBS ({left})"
                                                         f"{' after the backstop qdel' if self.qdels else ''}")
            if read["State"] in ("NOTFOUND", "OTHERJOB", "UNKNOWN", "READERR"):
                self.event(f"READ missing ({read['State']}): never 'left'; read again")
            elif read["Owner"] != state["Owner"] or read["Name"] != state["JobName"]:
                return self.end(13, f"REFUSED: job {state['JobID']} owner {read['Owner']!r} name {read['Name']!r} != the state "
                                    f"file's {state['Owner']!r} / {state['JobName']!r}: no action")
            elif read["StartEpoch"] is not None and now >= read["StartEpoch"] + limit + grace:
                last = self.qdels[-1]["Epoch"] if self.qdels else None
                if len(self.qdels) < QDEL_MAX and (last is None or now - last >= QDEL_GAP_SECONDS):
                    self.event(f"LIMIT {now - read['StartEpoch']:.0f} s from the PBS start {utc(read['StartEpoch'])} >= "
                               f"{limit:.0f} + {grace:.0f} s (state {read['State']}) -> qdel")
                    self.qdel(state)
                    self.qdels[-1]["Epoch"] = now
            self.clock.sleep(interval)
        return self.end(5, "the read bound reached; the job was not seen to leave PBS")


def check(directory, pbs_bin):
    """--check: the directory writable, qstat / qdel executable, this python; one line, exit 0 / 1."""
    problems = []
    try:
        os.makedirs(directory, exist_ok=True)
        probe = os.path.join(directory, f".check-{os.getpid()}")
        with open(probe, "w") as f:
            f.write("ok\n")
        os.remove(probe)
    except OSError as error:
        problems.append(f"the directory {directory} is not writable ({error})")
    for name in ("qstat", "qdel"):
        if not os.access(os.path.join(pbs_bin, name), os.X_OK):
            problems.append(f"{os.path.join(pbs_bin, name)} is not executable")
    if sys.version_info < (3, 6):
        problems.append(f"python {sys.version.split()[0]} < 3.6")
    if problems:
        print(f"BACKSTOP-CHECK FAIL {'; '.join(problems)}")
        return 1
    print(f"BACKSTOP-CHECK OK user={getpass.getuser()} host={socket.gethostname()}")
    return 0


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] == ["--check"] and len(argv) == 4 and argv[2] == "--pbs-bin":
        return check(argv[1], argv[3])
    if len(argv) == 1:
        return Backstop(argv[0]).loop()
    print(__doc__.split("\n\n")[-1], file=sys.stderr)
    return 64


if __name__ == "__main__":
    sys.exit(main())
