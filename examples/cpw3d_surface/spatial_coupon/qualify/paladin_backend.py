# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Form Y-b of `coupon-library qualify` (lane chainport; USER decision 1674 (2); Y-SPEC s2):
a run_stages plan run as paladin runs, ONE RUN PER STAGE, on ONE node per stage (multi-node
stages are refused: NOT DONE).  Paladin runs one configuration file per run and passes no
environment variable and no command-line flag to Palace, so each stage's PALACE_RESPONSE_*
variables become the stage configuration's Solver.Electrostatic.ResponseArchive keys.

The backend is configured by the cluster profile's "Paladin" block (PALADIN_KEYS; a profile
with both a Paladin and a Launcher block is refused: one route per run).  PaladinJob drives
the stages of one plan in order, one step per advance():
  - resume by metadata: before any submit the stage's runs are looked up by its identity
    (programme, lane, case, stage, plan sha) with `paladin palace my-jobs -m ...`; a queued,
    running or completed run is ADOPTED, never submitted again; a cancelled, out-of-memory or
    unexplained failed run STOPS the drive (nothing is resubmitted); a submission_failed row
    without a job id is no run (one re-submit; a second = STOP); two runs of one stage = STOP;
    a SubmitIntent with a nonce is on disk before each submit call, and a lookup that is not
    shown complete is a failed lookup, never 'no rows' (INTENT_RULE, LOOKUP_RULE);
  - before each submit: the run configuration schema-checked against the selector, the pin
    checked (pin_check_v2 check: its PASS line in the stage record), the E-2 sizing
    (ERRATUM-1702: the plan estimate x 1.30 <= 0.70 x MemTotal, < 10 % margin -> the next
    type) choosing --instance-type, and for an archive writer the /data guard (the archive
    directory created, free >= 150 GB, the estimate <= 10 % of free, the directory empty);
  - while it runs: describe-run polled; past the stage's CapSeconds from its start the run is
    cancelled (recorded; the row then reads canceled: the stage STOPS as timed-out);
  - after it ends: the completion read (a row failed by the upload alone counts as completed
    only if PREREG-ESSURR-NOTE-1's five conditions hold), pin_check_v2 postrun (its PASS line),
    the resolved configuration's leaf diff (only paladin's own paths may differ), the outputs
    copied on /data into the stage's Problem.Output with equal sha256 manifests, the log parsed
    by run_stages.parse_log, an archive writer's archive verified (sources x ranks files, none
    empty, per-rank sizes equal across sources, the zero-hole scan in place), and what paladin
    uploaded to qdash (the run's output artifacts and the bytes of the top-level and
    AMR-iteration tables, log and meshes);
  - when every stage is complete: status.json in run_stages' schema in the plan's directory on
    /data (Backend "paladin"), so that the chain's fetch, verification and archive deletion
    read it as they read a PBS job's.

  python3 paladin_backend.py --dry-run PLAN --profile PROFILE --out DIR [--map REMOTE=LOCAL]
  python3 paladin_backend.py --submit PLAN --profile PROFILE --out DIR --remote-dir DIR
      [--map REMOTE=LOCAL] [--archive-gb DIRECTORY=GB] [--poll-seconds S]
"""
import argparse
import copy
import datetime
import fcntl
import hashlib
import json
import os
import re
import shlex
import socket
import subprocess
import sys
import time
import uuid

# The stage variables of build_plan.py and their keys.  Every other variable is refused.
MODE_VARIABLES = {
    "PALACE_RESPONSE_ARCHIVE_ONLY": "Archive",
    "PALACE_RESPONSE_REDUCE_ONLY": "Reduce",
    "PALACE_RESPONSE_ESTIMATE_ONLY": "Estimate",
}
PATH_VARIABLES = {
    "PALACE_RESPONSE_ARCHIVE_DIR": "Directory",
    "PALACE_RESPONSE_ESTIMATE_REQUEST": "EstimateRequest",
}
SWITCH_VARIABLES = {
    "PALACE_RESPONSE_SOURCE_TIMING": "SourceTiming",
    "PALACE_RESPONSE_RECYCLE_INITIAL_GUESS": "RecycleInitialGuess",
}
INTEGER_VARIABLES = {"PALACE_RESPONSE_BLOCK_SIZE": "BlockSize"}


def response_archive_keys(environment):
    """The ResponseArchive block equivalent to a stage's Environment ({} = none)."""
    keys = {}
    for name, value in sorted(environment.items()):
        if name in MODE_VARIABLES:
            if value != "1":
                raise ValueError(f"{name}={value!r}: a plan sets a mode variable to '1' only")
            if "Mode" in keys:
                raise ValueError(f"{name}: more than one response-archive mode")
            keys["Mode"] = MODE_VARIABLES[name]
        elif name in PATH_VARIABLES:
            if not os.path.isabs(value):
                # Worker and reducer are different paladin runs, each in its own run_dir.
                raise ValueError(f"{name}={value!r}: a paladin run needs an absolute path")
            keys[PATH_VARIABLES[name]] = value
        elif name in SWITCH_VARIABLES:
            if value != "1":
                raise ValueError(f"{name}={value!r}: a plan sets a switch to '1' only")
            keys[SWITCH_VARIABLES[name]] = True
        elif name in INTEGER_VARIABLES:
            number = int(value)
            if number < 1:
                raise ValueError(f"{name}={value!r} must be positive")
            keys[INTEGER_VARIABLES[name]] = number
        else:
            raise ValueError(f"{name}: no configuration key carries this variable")
    return keys


def run_config(config, environment):
    """The stage configuration with its Environment as ResponseArchive keys."""
    keys = response_archive_keys(environment)
    out = copy.deepcopy(config)
    electrostatic = out.setdefault("Solver", {}).setdefault("Electrostatic", {})
    if "ResponseArchive" in electrostatic:
        raise ValueError("the stage configuration already has a ResponseArchive block")
    if keys:
        electrostatic["ResponseArchive"] = keys
    return out


def map_path(path, mapping):
    for remote, local in mapping:
        if path == remote or path.startswith(remote.rstrip("/") + "/"):
            return local.rstrip("/") + path[len(remote.rstrip("/")):]
    return path


def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


class PaladinError(RuntimeError):
    """A paladin call that failed, a check that did not pass, or a run that did not complete:
    the drive stops (fail closed; nothing is resubmitted)."""


# The cluster profile's Paladin block: exactly these keys.
PALADIN_KEYS = {
    "Selector": "the container selector (dev-simlapointe-sp-pin-<sha10>): --palace of every submit",
    "PinRecord": "the pin record (pin_check_v1 record) whose 'selector' is Selector",
    "PinCheck": "pin_check_v2.py: 'check' before each submit, 'postrun' after each run",
    "SchemaCheck": "ctr_schema_check_v1.py: each run configuration against the selector's schema",
    "ToolPython": "the python of the two tools (paladin's: they read qdash and S3)",
    "ImageSHA256": "the sha256 of the pinned SIF's copy (recorded with the pin)",
    "ImageBinary": "the Palace executable inside the image (the plans' Binary)",
    "ImageBinarySHA256": "its sha256 (the plans' BinarySHA256 and the frozen digest)",
    "CLI": "the paladin executable",
    "Host": "the ssh host that reads and writes /data (soca-green-job); 'local' = this host (the tests' form)",
    "Lane": "the lane metadata of every submit",
    "MemTotalGiB": "MemTotal per instance type, GiB, in the order of preference (E-2 sizing)",
    "ZeroHoleScan": "zero_hole_scan_fbc2_v1.py, run on Host over an archive in place",
}
PALADIN_ROUTE_RULE = ("one route per run (decision 1799 (l)): a cluster profile carries a Launcher block (form X, PBS) or a "
                      "Paladin block (form Y-b), never both; --backend paladin needs the Paladin block and --backend pbs "
                      "(the default) refuses it")


def paladin_block_of(profile):
    """The profile's validated Paladin block, or None (no block = the PBS route of record)."""
    block = profile.get("Paladin")
    if block is None:
        return None
    if profile.get("Launcher") is not None:
        raise ValueError(f"the cluster profile carries both a Launcher and a Paladin block ({PALADIN_ROUTE_RULE})")
    if not isinstance(block, dict):
        raise ValueError("the Paladin block must be an object")
    unknown = sorted(set(block) - set(PALADIN_KEYS))
    missing = sorted(set(PALADIN_KEYS) - set(block))
    if unknown or missing:
        raise ValueError(f"the Paladin block: unknown keys {unknown}, missing keys {missing} (exactly {sorted(PALADIN_KEYS)})")
    for key in ("ImageSHA256", "ImageBinarySHA256"):
        if not re.fullmatch(r"[0-9a-f]{64}", str(block[key])):
            raise ValueError(f"Paladin.{key} must be a sha256 (64 lowercase hex), not {block[key]!r}")
    for key in ("PinRecord", "PinCheck", "SchemaCheck", "ToolPython", "ImageBinary", "ZeroHoleScan"):
        if not isinstance(block[key], str) or not os.path.isabs(block[key]):
            raise ValueError(f"Paladin.{key} must be an absolute path, not {block[key]!r}")
    for key in ("Selector", "CLI", "Host", "Lane"):
        if not isinstance(block[key], str) or not block[key]:
            raise ValueError(f"Paladin.{key} must be a non-empty string")
    memory = block["MemTotalGiB"]
    if (not isinstance(memory, dict) or not memory or
            any(not isinstance(value, (int, float)) or isinstance(value, bool) or value <= 0 for value in memory.values())):
        raise ValueError("Paladin.MemTotalGiB must map at least one instance type to a positive MemTotal in GiB")
    with open(block["PinRecord"]) as f:
        record = json.load(f)
    if record.get("selector") != block["Selector"]:
        raise ValueError(f"the pin record {block['PinRecord']} is of selector {record.get('selector')!r}, "
                         f"not Paladin.Selector {block['Selector']!r}")
    return dict(block)


def backend_of(backend, profile):
    """The Paladin block for --backend paladin, None for --backend pbs (PALADIN_ROUTE_RULE)."""
    block = paladin_block_of(profile)
    if backend == "pbs":
        if block is not None:
            raise ValueError(f"--backend pbs with a cluster profile that carries a Paladin block ({PALADIN_ROUTE_RULE})")
        return None
    if backend != "paladin":
        raise ValueError(f"unknown backend {backend!r}")
    if block is None:
        raise ValueError(f"--backend paladin needs a cluster profile with a Paladin block ({PALADIN_ROUTE_RULE})")
    return block


def route_line(block):
    """The ONE line a paladin run prints at its start."""
    return (f"route paladin: selector {block['Selector']}, pin record {block['PinRecord']} "
            f"(sha256 {sha256_of(block['PinRecord'])[:16]}), image {block['ImageSHA256'][:16]}, "
            f"in-image {block['ImageBinary']} {block['ImageBinarySHA256'][:16]}")


# ERRATUM-1702 E-2, the interim sizing rule of every container job.
SIZING_FACTOR, SIZING_FRACTION, SIZING_MARGIN = 1.30, 0.70, 0.10
SIZING_RULE = ("ERRATUM-1702 E-2: the plan's EstimatedPalacePeakGiB x 1.30 must be <= 0.70 x MemTotal of the instance type; the "
               "first type of Paladin.MemTotalGiB where it fits with a margin >= 10 % of that bound is chosen (a smaller margin "
               "moves to the next type: m8g -> r8g); above the last type's bound = STOP (a CONTACT); one node per stage")


def instance_for(estimate_gib, mem_total_gib):
    """The --instance-type of a stage by SIZING_RULE, with the figures; PaladinError above the last type."""
    figure = float(estimate_gib) * SIZING_FACTOR
    types = list(mem_total_gib)
    for position, name in enumerate(types):
        bound = SIZING_FRACTION * float(mem_total_gib[name])
        margin = (bound - figure) / bound
        if figure <= bound and (margin >= SIZING_MARGIN or position == len(types) - 1):
            return {"Type": name, "EstimateGiB": float(estimate_gib), "Factor": SIZING_FACTOR, "FigureGiB": figure,
                    "MemTotalGiB": float(mem_total_gib[name]), "BoundGiB": bound, "MarginFraction": margin, "Rule": SIZING_RULE}
    raise PaladinError(f"E-2 sizing: {estimate_gib:.1f} GiB x {SIZING_FACTOR} = {figure:.1f} GiB is above 0.70 x MemTotal of every "
                       f"type {types}: CONTACT ({SIZING_RULE})")


def utc_now():
    return datetime.datetime.now(datetime.timezone.utc)


def parse_time(text):
    if not text:
        return None
    value = datetime.datetime.fromisoformat(str(text).replace("Z", "+00:00"))
    return value if value.tzinfo else value.replace(tzinfo=datetime.timezone.utc)


def paladin_envelope(cli, arguments, *, timeout=600):
    """Run `cli arguments` (a --json command); the whole envelope (data and warnings), or PaladinError."""
    try:
        completed = subprocess.run([cli] + list(arguments), capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        raise PaladinError(f"{' '.join(arguments[:2])}: no answer in {timeout} s") from None
    try:
        envelope = json.loads(completed.stdout)
    except ValueError:
        raise PaladinError(f"{' '.join(arguments[:2])}: no JSON envelope (rc {completed.returncode})") from None
    if completed.returncode != 0 or not isinstance(envelope, dict) or not envelope.get("ok"):
        error = envelope.get("error") if isinstance(envelope, dict) else None
        raise PaladinError(f"{' '.join(arguments[:2])}: {error!r} (rc {completed.returncode})")
    return envelope


def paladin_call(cli, arguments, *, timeout=600):
    """Run `cli arguments` (a --json command); the envelope's data, or PaladinError."""
    return paladin_envelope(cli, arguments, timeout=timeout)["data"]


def remote_shell(host, command, *, timeout=300, stdin=None):
    """`command` on host by ssh (None: here, by bash; the tests).  The CompletedProcess."""
    argv = (["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=20", host, command] if host else ["bash", "-c", command])
    try:
        return subprocess.run(argv, input=stdin or "", capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        raise PaladinError(f"{'ssh ' + host if host else 'bash'}: no answer in {timeout} s: {command[:160]}") from None


def remote_text(host, command, *, timeout=300, stdin=None, what="a read on /data"):
    result = remote_shell(host, command, timeout=timeout, stdin=stdin)
    if result.returncode != 0:
        raise PaladinError(f"{what} failed (rc {result.returncode}): {result.stderr.strip()[:300]}")
    return result.stdout


def write_json_atomic(path, value):
    """The file replaced in one step and flushed to disk (the file and its directory) before return:
    a SubmitIntent saved by it survives the end of the driver process at any later point."""
    temporary = path + ".tmp"
    with open(temporary, "w") as f:
        json.dump(value, f, indent=1)
        f.flush()
        os.fsync(f.fileno())
    os.replace(temporary, path)
    directory = os.open(os.path.dirname(os.path.abspath(path)), os.O_RDONLY)
    try:
        os.fsync(directory)
    finally:
        os.close(directory)


LOCK_RULE = ("decision 1820 (2): one driver per root: an exclusive flock on the root's lock file, held on an open descriptor "
             "for the driver's life (the kernel releases it when the driver process ends); the file names its holder (PID, "
             "Host, StartUTC, Argv); a held lock = refused; a holder record under a free lock = a stale lock (its driver ended "
             "without release): taken over with one line in <lock>.takeovers")


class DriverLock:
    """LOCK_RULE on `path`; acquire() -> self (PaladinError when another driver holds it); release()."""

    def __init__(self, path, *, log=print, now=utc_now):
        self.path, self.log, self.now = str(path), log, now
        self.descriptor = None

    def acquire(self):
        os.makedirs(os.path.dirname(os.path.abspath(self.path)), exist_ok=True)
        descriptor = os.open(self.path, os.O_RDWR | os.O_CREAT, 0o644)
        try:
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            holder = os.pread(descriptor, 4096, 0).decode(errors="replace").strip()
            os.close(descriptor)
            raise PaladinError(f"the driver lock {self.path} is held by a live driver {holder or '(no record yet)'}: "
                               f"refused ({LOCK_RULE})") from None
        previous = os.pread(descriptor, 4096, 0).decode(errors="replace").strip()
        utc = self.now().strftime("%Y-%m-%dT%H:%M:%SZ")
        if previous:
            line = f"{utc} pid {os.getpid()} on {socket.gethostname()} took over the stale lock of {previous}\n"
            with open(self.path + ".takeovers", "a") as f:
                f.write(line)
                f.flush()
                os.fsync(f.fileno())
            self.log(f"driver lock {self.path}: {line.strip()} ({LOCK_RULE})")
        holder = json.dumps({"PID": os.getpid(), "Host": socket.gethostname(), "StartUTC": utc, "Argv": sys.argv[:6]})
        os.ftruncate(descriptor, 0)
        os.pwrite(descriptor, holder.encode(), 0)
        os.fsync(descriptor)
        self.descriptor = descriptor
        return self

    def release(self):
        if self.descriptor is None:
            return
        os.ftruncate(self.descriptor, 0)
        os.fsync(self.descriptor)
        fcntl.flock(self.descriptor, fcntl.LOCK_UN)
        os.close(self.descriptor)
        self.descriptor = None


def tool_line(block, arguments, pattern, what, *, timeout=600):
    """Run a check tool with the block's ToolPython: its PASS line, or PaladinError (fail closed)."""
    try:
        completed = subprocess.run([block["ToolPython"], *arguments], capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        raise PaladinError(f"{what}: no answer in {timeout} s") from None
    lines = completed.stdout.strip().splitlines()
    passed = [line for line in lines if re.match(pattern, line)]
    if completed.returncode != 0 or not passed:
        raise PaladinError(f"{what} did not pass (rc {completed.returncode}): {(lines or [completed.stderr.strip()])[-1][:300]}")
    return passed[-1]


def verify_route(block):
    """Before the first submit of a qualify run (as verify_remote_executables for PBS): the pin
    checked (each stage checks it again before its submit and after its run)."""
    line = tool_line(block, [block["PinCheck"], "check", "--record", block["PinRecord"]], r"PASS ",
                     "pin_check_v2 check before the first submit")
    return {"Backend": "paladin", "Host": block["Host"], "Selector": block["Selector"],
            "PinRecord": {"Path": block["PinRecord"], "SHA256": sha256_of(block["PinRecord"])}, "PinCheck": line,
            "Executables": [{"Binary": block["ImageBinary"], "SHA256": block["ImageBinarySHA256"], "MPIExec": "paladin"}],
            "Image": {"SHA256": block["ImageSHA256"]}, "UTC": utc_now().strftime("%Y-%m-%dT%H:%M:%SZ"),
            "Rule": PALADIN_ROUTE_RULE}


# The stage identity every submit carries as metadata (and the lookup filters on).
def stage_identity(lane, case, stage, plan_sha256):
    return {"programme": "spr", "lane": lane, "case": case, "stage": stage, "plan": "sha256:" + plan_sha256[:16]}


def metadata_arguments(identity):
    return [part for key, value in identity.items() for part in ("-m", f"{key}={value}")]


LOOKUP_RULE = ("decision 1820 (1): a lookup is FAILED (a STOP, never 'no rows') when paladin answers rc != 0, ok false, no JSON, "
               "data without a 'jobs' list, or ANY warning (the client's rows_omitted drops rows it cannot rebuild: a listing "
               "with a warning is not shown complete)")


def find_runs(cli, identity):
    """Every paladin row whose metadata holds the stage identity (server filter, re-checked here);
    PaladinError for a failed lookup (LOOKUP_RULE)."""
    envelope = paladin_envelope(cli, ["palace", "my-jobs", *metadata_arguments(identity), "--json"])
    data, warnings = envelope.get("data"), envelope.get("warnings")
    if warnings:
        raise PaladinError(f"palace my-jobs: a listing with warnings {[(w or {}).get('code') for w in warnings]}: a FAILED "
                           f"lookup ({LOOKUP_RULE})")
    if not isinstance(data, dict) or not isinstance(data.get("jobs"), list):
        raise PaladinError(f"palace my-jobs: no 'jobs' list in the answer: a FAILED lookup ({LOOKUP_RULE})")
    return [row for row in data["jobs"]
            if all((row.get("metadata") or {}).get(key) == value for key, value in identity.items())]


LIVE_STATUSES = {"submitted", "queued", "running", "canceling"}
RESUME_RULE = ("decision 1799 (g): a stage's runs are found by its identity metadata + the plan sha; a queued / running / completed "
               "run is adopted and never submitted again; a canceled (the cap included), out_of_memory, unknown or unexplained "
               "failed run STOPS the drive (no resubmit); a submission_failed row without a job id is NO RUN (CTR-HOWTO s12: one "
               "re-submit, a second = STOP); two runs of one stage = STOP")
INTENT_RULE = ("decision 1820 (1): a SubmitIntent (stage, plan sha, UTC, a nonce that also goes into the run metadata) is saved "
               "and flushed to disk BEFORE each submit call and closed (Outcome) after it; an OPEN intent (the driver ended in "
               "between) or a recorded RunID that the lookup does not show = the lookup repeated 3 times over >= 60 s, a run "
               "found = adopted, none = STOP (never a submit); any submit needs two consistent empty lookups >= 20 s apart; two "
               "runs of one stage = the newer live run(s) early-ended by 'palace cancel' (recorded) before the STOP, and an "
               "archive writer's Directory marked unusable (ARCHIVE-UNUSABLE.json: no reducer and no writer uses it)")
LOOKUP_REPEATS, LOOKUP_REPEAT_SECONDS, SUBMIT_LOOKUP_GAP_SECONDS = 3, 30.0, 20.0
UNUSABLE_MARKER = "ARCHIVE-UNUSABLE.json"
WAIT_RULE = ("decision 1820 (3): a describe-run answer that is not an object, names another run_id, has no status, or reads "
             "'completed' / 'failed' without its run directory is a TRANSPORT answer, not a state: re-asked after 10 s and 20 s; "
             "three in a row = one failed poll (counted; a usable answer resets the count); STOP at >= 10 failed polls spanning "
             ">= 600 s, with no State recorded and the run left to its run watcher; 'unknown' only from an explicit non-empty "
             "status the backend does not list; the stage cap counts from 'started', from the submit time for a live answer "
             "without 'started'")
WAIT_BACKOFF_SECONDS, WAIT_TRANSPORT_POLLS, WAIT_TRANSPORT_SECONDS = (10.0, 20.0), 10, 600.0


def pause(seconds):
    """The wait between repeated lookups (a module function: the tests replace it)."""
    time.sleep(seconds)


class DuplicateRuns(PaladinError):
    """Two or more runs of one stage (RESUME_RULE); .runs = the rows."""

    def __init__(self, message, runs):
        super().__init__(message)
        self.runs = runs


def no_run(row):
    return str(row.get("status", "")).lower() == "submission_failed" and not row.get("job_id")


def lookup_decision(rows):
    """('adopt', row) | ('submit', number of earlier no-run rows); PaladinError otherwise (RESUME_RULE)."""
    runs = [row for row in rows if not no_run(row)]
    if len(runs) > 1:
        raise DuplicateRuns(f"{len(runs)} runs of one stage: {[(row['run_id'], row.get('status')) for row in runs]}: "
                            f"none is adopted ({RESUME_RULE})", runs)
    if runs:
        return "adopt", runs[0]
    if len(rows) >= 2:
        raise PaladinError(f"{len(rows)} submission_failed rows of one stage {[row['run_id'] for row in rows]}: "
                           f"no further submit ({RESUME_RULE})")
    return "submit", len(rows)


def run_class(status):
    status = str(status or "").lower()
    if status in LIVE_STATUSES:
        return "live"
    return {"completed": "completed", "failed": "failed", "canceled": "cancelled", "cancelled": "cancelled",
            "out_of_memory": "out-of-memory"}.get(status, "unknown")


def expected_outputs(config):
    """The postpro files every completed run of this configuration writes (palace.json; the two
    response matrices of a reducer or of a response-matrix solve that is not an archive writer).
    The chain's own validations read the other tables later."""
    names = ["palace.json"]
    electrostatic = config.get("Solver", {}).get("Electrostatic", {})
    mode = electrostatic.get("ResponseArchive", {}).get("Mode")
    if mode == "Reduce" or (mode is None and electrostatic.get("ResponseMatrix")):
        names += ["domain-response-matrix.csv", "surface-response-matrix.csv"]
    return names


COMPLETION_RULE = ("PREREG-ESSURR-NOTE-1: a row failed by the qdash upload ALONE is a completed run iff (1) log.out has 'Elapsed "
                   "Time Report (s)', (2) console.log has 'Runner finished (exit 0)', (3) every ERROR line of console.log is the "
                   "upload harness's (at least one), (4) every expected output is present by name, (5) no cap / memory stop fired "
                   "(this backend issued no cancel; paladin's own guard ends a row out_of_memory, not failed); a completed row "
                   "must show (1), (2) and (4)")


def completion_read(host, run_dir, expected, *, row_status, cancelled_by_backend):
    """The five conditions of COMPLETION_RULE read on the run directory."""
    q = shlex.quote(run_dir)
    text = remote_text(host, f"cd {q} && echo ELAPSED $(grep -c 'Elapsed Time Report (s)' log.out 2>/dev/null); "
                             f"echo FINISHED $(grep -a -c 'Runner finished (exit 0)' console.log 2>/dev/null); "
                             f"grep -a 'ERROR' console.log 2>/dev/null | sed 's/^/ERRORLINE /' | head -20; "
                             f"grep -a -E '^[A-Za-z_.]+(Error|Exception)' console.log 2>/dev/null | tail -1 | sed 's/^/LAST /'; "
                             f"ls -1 postpro 2>/dev/null | sed 's/^/FILE /'", what=f"the completion read of {run_dir}")
    lines = text.splitlines()
    value = {line.split()[0]: line.split()[1] if len(line.split()) > 1 else "0" for line in lines
             if line.startswith(("ELAPSED ", "FINISHED "))}
    errors = [line[len("ERRORLINE "):] for line in lines if line.startswith("ERRORLINE ")]
    last = next((line[len("LAST "):] for line in lines if line.startswith("LAST ")), None)
    files = {line[len("FILE "):] for line in lines if line.startswith("FILE ")}
    failed = str(row_status).lower() == "failed"
    conditions = {
        "1-ElapsedTimeReport": int(value.get("ELAPSED", "0") or 0) >= 1,
        "2-RunnerFinishedExit0": int(value.get("FINISHED", "0") or 0) >= 1,
        "3-OnlyUploadErrors": (bool(errors) and all("upload" in line for line in errors)) if failed else True,
        "4-ExpectedOutputs": all(name in files for name in expected),
        "5-NoStop": not cancelled_by_backend and str(row_status).lower() != "out_of_memory",
    }
    return {"Conditions": conditions, "Holds": all(conditions.values()), "ErrorLines": errors, "TracebackLast": last,
            "MissingOutputs": [name for name in expected if name not in files], "RowStatus": row_status, "Rule": COMPLETION_RULE}


def leaves(value, prefix=""):
    if isinstance(value, dict):
        out = {}
        for key, item in value.items():
            out.update(leaves(item, f"{prefix}.{key}"))
        return out
    return {prefix: json.dumps(value, sort_keys=True)}


# Paladin owns these paths of the configuration it runs (CTR-HOWTO s4 / s5 m-5).
PALADIN_OWNED_PATHS = {".Model.Mesh", ".Problem.Output", ".Solver.Device", ".Model.Refinement.SaveAdaptMesh"}


def leaf_diff(submitted, ran):
    """The leaf paths where the configuration paladin ran differs from the submitted one."""
    first, second = leaves(submitted), leaves(ran)
    return sorted(path for path in set(first) | set(second) if first.get(path) != second.get(path))


def manifest(text):
    """{relative path: sha256} of `find . -type f -exec sha256sum` output."""
    out = {}
    for line in text.splitlines():
        if line.strip():
            digest, path = line.split(None, 1)
            out[path.strip()[2:] if path.strip().startswith("./") else path.strip()] = digest
    return out


ARCHIVE_FREE_BYTES, ARCHIVE_SHARE = 150e9, 0.10
ARCHIVE_RULE = ("decision 1799 (j): before an archive writer's submit its absolute /data Directory is created and must be empty; "
                "/data free >= 150 GB and the stage's ArchiveGBEstimate <= 10 % of free; after the writer: one "
                "source-NNNNNN-rank-NNNNNN-V.bin per (source, rank) of its configuration and partition, none empty, every rank's "
                "size equal across sources, the zero-hole scan CLEAN in place; the archive is deleted by the chain only after "
                "the reducer's outputs verify (finish_case, as in the PBS and form-X routes)")


def archive_name(source, rank):
    return f"source-{int(source):06d}-rank-{int(rank):06d}-V.bin"


# The tables, logs and meshes paladin 0.2.16 uploads: the run's top level and each AMR iteration's
# top level (ERRATUM-1802 E-1), never other sub-directories or suffixes.
UPLOADED_SUFFIXES = (".csv", ".json", ".mesh", ".meshgz")


class PaladinJob:
    """The paladin runs of ONE run_stages plan, in order (module docstring)."""

    def __init__(self, plan_path, *, block, remote_dir, local_dir, mapping=(), archive_estimates_gb=None,
                 log=print, now=utc_now, sleep=None, exclusive=True):
        with open(plan_path, "rb") as f:
            data = f.read()
        self.plan = json.loads(data)
        self.plan_path = plan_path
        self.plan_sha256 = hashlib.sha256(data).hexdigest()
        self.block, self.remote_dir, self.local_dir = block, remote_dir, local_dir
        self.mapping = list(mapping)
        self.archive_estimates_gb = dict(archive_estimates_gb or {})
        self.log, self.now = log, now
        self.sleep = sleep or (lambda seconds: pause(seconds))
        # LOCK_RULE: the job's local state directory and the plan's directory (two --out directories of one
        # plan are refused too); False when the caller holds its root's lock (qualify_library.run_qualify).
        self.exclusive, self.locks = exclusive, []
        self.host = None if block["Host"] == "local" else block["Host"]
        self.confirmed = set()
        self.case = self.plan["Case"]
        if int(self.plan.get("Nodes", 1)) != 1:
            raise PaladinError(f"{self.case}: a {self.plan['Nodes']}-node plan: multi-node paladin stages are NOT DONE "
                               "(one node per stage, parity with form X)")
        estimate = (self.plan.get("Instance") or {}).get("EstimatedPalacePeakGiB")
        if estimate is None:
            raise PaladinError(f"{self.case}: the plan has no Instance.EstimatedPalacePeakGiB (E-2 sizing)")
        self.sizing = instance_for(estimate, block["MemTotalGiB"])
        self.ranks = int(self.plan["Ranks"])
        self.state_path = os.path.join(local_dir, "paladin-status.json")
        self.configs_dir = os.path.join(local_dir, "paladin")
        if os.path.isfile(self.state_path):
            with open(self.state_path) as f:
                self.state = json.load(f)
            if self.state.get("PlanSHA256") != self.plan_sha256:
                raise PaladinError(f"{self.state_path} records plan {self.state.get('PlanSHA256', '')[:16]}, the plan is "
                                   f"{self.plan_sha256[:16]}: not resumable")
        else:
            self.state = {"Backend": "paladin", "Case": self.case, "PlanSHA256": self.plan_sha256, "Stages": {}}

    def lock(self):
        """LOCK_RULE, once per driver process (at the first advance(); a dry run takes no lock)."""
        if not self.exclusive or self.locks:
            return
        os.makedirs(self.local_dir, exist_ok=True)
        paths = sorted({os.path.join(os.path.realpath(self.local_dir), "paladin-driver.lock"),
                        os.path.join(os.path.dirname(os.path.realpath(self.plan_path)), "paladin-driver.lock")})
        try:
            for path in paths:
                self.locks.append(DriverLock(path, log=self.log, now=self.now).acquire())
        except PaladinError:
            self.unlock()
            raise

    def unlock(self):
        for lock in self.locks:
            lock.release()
        self.locks = []

    def save(self):
        os.makedirs(self.local_dir, exist_ok=True)
        write_json_atomic(self.state_path, self.state)

    # The run configuration of a stage: the stage config + its ResponseArchive keys, the mesh local.
    def stage_run_config(self, stage):
        with open(map_path(stage["Config"], self.mapping)) as f:
            config = run_config(json.load(f), stage.get("Environment", {}))
        mesh = config.get("Model", {}).get("Mesh")
        if mesh is not None:
            local = self.plan["MeshLocal"] if mesh == self.plan.get("MeshRemote") else map_path(mesh, self.mapping)
            config["Model"]["Mesh"] = local
        return config

    def write_run_config(self, stage):
        config = self.stage_run_config(stage)
        os.makedirs(self.configs_dir, exist_ok=True)
        path = os.path.join(self.configs_dir, stage["Name"] + ".json")
        with open(path, "w") as f:
            json.dump(config, f, indent=2)
        return path, config

    def submit_arguments(self, stage, path, nonce="<the SubmitIntent's nonce>"):
        name = stage["Name"]
        return ["palace", "submit-detailed", path, "--instance-type", self.sizing["Type"], "--nodes", "1",
                "--ranks-per-node", str(self.ranks), "--palace", self.block["Selector"],
                "--label", f"{self.block['Lane']}-{name}",
                *metadata_arguments(stage_identity(self.block["Lane"], self.case, name, self.plan_sha256)),
                "-m", f"nonce={nonce}",
                "--comment", f"{self.case} {name}: coupon-library qualify, paladin backend (form Y-b)", "--json"]

    def dry_run(self):
        """The run configurations and submit commands; contacts nothing."""
        runs = []
        for stage in self.plan["Stages"]:
            path, config = self.write_run_config(stage)
            runs.append({"Stage": stage["Name"], "RunConfig": path, "RunConfigSHA256": sha256_of(path),
                         "ResponseArchive": config.get("Solver", {}).get("Electrostatic", {}).get("ResponseArchive", {}),
                         "Mesh": config.get("Model", {}).get("Mesh"), "Requires": list(stage.get("Requires", [])),
                         "CapSeconds": stage.get("CapSeconds"), "Nodes": 1, "RanksPerNode": self.ranks,
                         "Instance": self.sizing, "Submit": [self.block["CLI"], *self.submit_arguments(stage, path)]})
        return runs

    # The checks of a submit (decision 1799 (h) / (i) / (j); REFREEZE s20.8 (c)).
    def tool_line(self, arguments, pattern, what, timeout=600):
        return tool_line(self.block, arguments, pattern, what, timeout=timeout)

    def preflight(self):
        """Once per job before its first submit: the mesh the runs upload is the plan's."""
        if "Preflight" in self.state:
            return
        mesh = self.plan.get("MeshLocal")
        record = {"Backend": "paladin", "Selector": self.block["Selector"],
                  "PinRecord": {"Path": self.block["PinRecord"], "SHA256": sha256_of(self.block["PinRecord"])},
                  "ImageSHA256": self.block["ImageSHA256"], "ImageBinary": self.block["ImageBinary"],
                  "ImageBinarySHA256": self.block["ImageBinarySHA256"], "PlanSHA256": self.plan_sha256,
                  "Instance": self.sizing, "UTC": self.now().strftime("%Y-%m-%dT%H:%M:%SZ")}
        if self.plan.get("BinarySHA256") != self.block["ImageBinarySHA256"]:
            raise PaladinError(f"the plan pins {self.plan.get('BinarySHA256')}, the pinned image's executable is "
                               f"{self.block['ImageBinarySHA256']}")
        if mesh is not None:
            if not os.path.isfile(mesh) or sha256_of(mesh) != self.plan.get("MeshSHA256"):
                raise PaladinError(f"the local mesh {mesh} is absent or not the plan's {self.plan.get('MeshSHA256', '')[:16]}")
            record["Mesh"] = {"Local": mesh, "SHA256": self.plan["MeshSHA256"]}
        try:
            version = subprocess.run([self.block["CLI"], "--version"], capture_output=True, text=True, timeout=120)
            record["PaladinVersion"] = version.stdout.strip()
        except (OSError, subprocess.TimeoutExpired):
            record["PaladinVersion"] = None
        self.state["Preflight"] = record
        self.save()

    def archive_guard(self, directory):
        estimate_gb = self.archive_estimates_gb.get(directory)
        if estimate_gb is None:
            raise PaladinError(f"no ArchiveGBEstimate for the archive {directory} ({ARCHIVE_RULE})")
        q = shlex.quote(directory)
        text = remote_text(self.host, f"mkdir -p {q} && df -Pk {q} | tail -1 && echo FILES $(ls -A {q} | wc -l)",
                           what=f"the archive guard of {directory}")
        fields = [line for line in text.splitlines() if not line.startswith("FILES")][-1].split()
        free = int(fields[3]) * 1024
        files = int(text.split("FILES")[-1].split()[0])
        guard = {"Directory": directory, "FreeBytes": free, "EstimateBytes": estimate_gb * 1e9, "FilesBefore": files,
                 "Rule": ARCHIVE_RULE}
        if free < ARCHIVE_FREE_BYTES or estimate_gb * 1e9 > ARCHIVE_SHARE * free or files:
            raise PaladinError(f"the archive guard of {directory} refuses the submit: {guard}")
        return guard

    def archive_verify(self, directory, config, record):
        sources = [entry["Index"] for entry in config.get("Boundaries", {}).get("PrescribedPotential", [])]
        q = shlex.quote(directory)
        text = remote_text(self.host, f"find {q} -maxdepth 1 -type f -name 'source-*-rank-*-V.bin' -printf '%f %s\\n'",
                           timeout=900, what=f"the listing of the archive {directory}")
        sizes = {}
        for line in text.splitlines():
            name, size = line.rsplit(" ", 1)
            sizes[name] = int(size)
        expected = {archive_name(source, rank) for source in sources for rank in range(self.ranks)}
        missing = sorted(expected - set(sizes))
        extra = sorted(set(sizes) - expected)
        empty = sorted(name for name, size in sizes.items() if size == 0)
        uneven = [rank for rank in range(self.ranks)
                  if len({sizes.get(archive_name(source, rank)) for source in sources}) > 1]
        report = os.path.join(self.remote_dir, f"{record['Name']}-archive-zero-hole.json")
        with open(self.block["ZeroHoleScan"]) as f:
            tool = f.read()
        scan = remote_shell(self.host, f"python3 - --out {shlex.quote(report)} {q}", timeout=1800, stdin=tool)
        verification = {"Directory": directory, "Sources": len(sources), "Ranks": self.ranks, "Expected": len(expected),
                        "Found": len(sizes), "Missing": missing[:5], "MissingCount": len(missing), "Extra": extra[:5],
                        "Empty": empty[:5], "UnevenRanks": uneven[:5], "Bytes": sum(sizes.values()),
                        "EstimateBytes": (self.archive_estimates_gb.get(directory) or 0) * 1e9,
                        "ZeroHoleScan": {"ReturnCode": scan.returncode, "Report": report,
                                         "ToolSHA256": hashlib.sha256(tool.encode()).hexdigest()},
                        "Rule": ARCHIVE_RULE}
        marker = remote_text(self.host, f"test -e {q}/{UNUSABLE_MARKER} && echo MARKED || true",
                             what=f"the unusable marker of {directory}").strip() == "MARKED"
        verification["Unusable"] = marker
        verification["OK"] = not (missing or extra or empty or uneven or marker) and scan.returncode == 0
        return verification

    def qdash_record(self, run_dir, data):
        q = shlex.quote(run_dir)
        text = remote_text(self.host, f"cd {q} && for f in log.out postpro/* postpro/iteration*/*; do "
                                      f"[ -f \"$f\" ] && [ ! -L \"$f\" ] && stat -c '%s %n' \"$f\"; done; true",
                           what=f"the upload sizes of {run_dir}")
        files = {}
        for line in text.splitlines():
            size, name = line.split(" ", 1)
            if name == "log.out" or name.endswith(UPLOADED_SUFFIXES):
                files[name] = int(size)
        return {"Artifacts": [[item.get("label"), item.get("type")] for item in data.get("output_artifacts") or []],
                "Files": len(files), "Bytes": sum(files.values()),
                "IterationFiles": sorted(name for name in files if name.startswith("postpro/iteration")),
                "Rule": "paladin 0.2.16 uploads log.out and the top-level csv / json / mesh of postpro and of each "
                        "postpro/iterationN (ERRATUM-1802 E-1); the archive is outside the run directory: not uploaded"}

    def describe(self, run_id):
        return paladin_call(self.block["CLI"], ["palace", "describe-run", str(run_id), "--json"])

    def submit(self, stage, record):
        """The checks of a submit, the submit, the RunID recorded at once; a failed submit
        re-read by the lookup (a no-run row = one re-submit; a refusal = STOP)."""
        name = stage["Name"]
        path, config = self.write_run_config(stage)
        record["RunConfig"] = {"Path": path, "SHA256": sha256_of(path)}
        record["SchemaCheck"] = self.tool_line([self.block["SchemaCheck"], self.block["Selector"], path], r"VALID ",
                                               f"the schema check of {name}")
        mode = config.get("Solver", {}).get("Electrostatic", {}).get("ResponseArchive", {})
        if mode.get("Mode") == "Archive":
            record["ArchiveGuard"] = self.archive_guard(mode["Directory"])
        elif mode.get("Directory"):
            self.refuse_unusable_archive(name, mode["Directory"])
        record["PinCheckBefore"] = self.tool_line([self.block["PinCheck"], "check", "--record", self.block["PinRecord"]],
                                                  r"PASS ", f"pin_check_v2 check before {name}")
        record["Instance"] = self.sizing
        record["Nodes"], record["RanksPerNode"] = 1, self.ranks
        record["Identity"] = stage_identity(self.block["Lane"], self.case, name, self.plan_sha256)
        intent = {"Nonce": uuid.uuid4().hex, "Stage": name, "PlanSHA256": self.plan_sha256, "UTC": self.now().isoformat(),
                  "Rule": INTENT_RULE}
        record.setdefault("SubmitIntents", []).append(intent)
        self.save()  # on disk before the call (INTENT_RULE)
        try:
            data = paladin_call(self.block["CLI"], self.submit_arguments(stage, path, intent["Nonce"]))
        except PaladinError as error:
            record.setdefault("SubmitErrors", []).append({"UTC": self.now().isoformat(), "Error": str(error)[:300]})
            intent["Outcome"] = "error"
            self.save()
            return None
        run_id = int(data["jobs"][0]["run_id"])
        record["RunID"], record["Submitted"] = run_id, self.now().isoformat()
        intent["Outcome"] = f"run {run_id}"
        self.save()
        self.log(f"{self.case}: submitted {name} as paladin run {run_id} ({self.sizing['Type']} x {self.ranks})")
        return run_id

    def open_intent(self, record):
        """The last SubmitIntent when it has no Outcome (the driver ended inside its submit), else None."""
        intents = record.get("SubmitIntents") or []
        return intents[-1] if intents and "Outcome" not in intents[-1] else None

    def lookup(self, stage, record, identity):
        """lookup_decision of one lookup; two runs of one stage = the newer ones early-ended first (INTENT_RULE)."""
        try:
            return lookup_decision(find_runs(self.block["CLI"], identity))
        except DuplicateRuns as error:
            self.resolve_duplicates(stage, record, error.runs)
            raise

    def resolve_duplicates(self, stage, record, runs):
        """Keep the oldest run (smallest run id); 'palace cancel' each newer live one; mark an archive
        writer's Directory unusable; all recorded (INTENT_RULE).  The STOP follows (the caller)."""
        ordered = sorted(runs, key=lambda row: int(row["run_id"]))
        ended = []
        for row in ordered[1:]:
            entry = {"RunID": int(row["run_id"]), "Status": row.get("status")}
            if run_class(row.get("status")) == "live":
                try:
                    paladin_call(self.block["CLI"], ["palace", "cancel", str(row["run_id"]), "--json"])
                    entry["EarlyEnded"] = True
                except PaladinError as error:
                    entry["EarlyEnded"], entry["Error"] = False, str(error)[:300]
            ended.append(entry)
        duplicate = {"UTC": self.now().isoformat(), "Kept": int(ordered[0]["run_id"]), "Newer": ended, "Rule": INTENT_RULE}
        archive = self.stage_run_config(stage).get("Solver", {}).get("Electrostatic", {}).get("ResponseArchive", {})
        if archive.get("Mode") == "Archive":
            directory = archive["Directory"]
            marker = os.path.join(directory, UNUSABLE_MARKER)
            remote_text(self.host, f"mkdir -p {shlex.quote(directory)} && cat > {shlex.quote(marker)}",
                        stdin=json.dumps({"Stage": stage["Name"], "Case": self.case, "Runs": [int(r["run_id"]) for r in ordered],
                                          "UTC": duplicate["UTC"], "Rule": INTENT_RULE}) + "\n",
                        what=f"the unusable marker of {directory}")
            duplicate["UnusableMarker"] = marker
        record.setdefault("DuplicateRuns", []).append(duplicate)
        self.save()
        self.log(f"{self.case}: {stage['Name']}: {len(ordered)} runs of one stage; newer early-ended {ended}; kept "
                 f"{duplicate['Kept']}")

    def refuse_unusable_archive(self, name, directory):
        """A stage that reads an archive refuses one marked unusable (INTENT_RULE)."""
        q = shlex.quote(os.path.join(directory, UNUSABLE_MARKER))
        if remote_text(self.host, f"test -e {q} && echo MARKED || true",
                       what=f"the unusable marker of {directory}").strip() == "MARKED":
            raise PaladinError(f"{name}: the archive {directory} is marked unusable ({UNUSABLE_MARKER}): STOP ({INTENT_RULE})")

    def repeated_lookup(self, stage, record, identity, why):
        """LOOKUP_REPEATS lookups LOOKUP_REPEAT_SECONDS apart; the first decision that adopts, else STOP."""
        for attempt in range(LOOKUP_REPEATS):
            if attempt:
                self.sleep(LOOKUP_REPEAT_SECONDS)
            decision, value = self.lookup(stage, record, identity)
            if decision == "adopt":
                return decision, value
        raise PaladinError(f"{stage['Name']}: {why}, and {LOOKUP_REPEATS} lookups over "
                           f"{(LOOKUP_REPEATS - 1) * LOOKUP_REPEAT_SECONDS:.0f} s show no run of the stage {identity}: STOP, no "
                           f"second submit ({INTENT_RULE})")

    def adopt(self, record, row, identity):
        name = record["Name"]
        if record.get("RunID") not in (None, int(row["run_id"])):
            raise PaladinError(f"{name}: the recorded run {record['RunID']} is not the stage's run {row['run_id']}")
        intent = self.open_intent(record)
        if intent is not None:
            nonce = (row.get("metadata") or {}).get("nonce")
            intent["Outcome"] = f"adopted run {row['run_id']}"
            intent["NonceMatches"] = nonce == intent["Nonce"]
        if record.get("RunID") is None:
            record["RunID"] = int(row["run_id"])
            record["Adopted"] = {"UTC": self.now().isoformat(), "Status": row.get("status"),
                                 "InstanceType": row.get("instance_type"), "Nodes": row.get("nodes"),
                                 "Nonce": (row.get("metadata") or {}).get("nonce"), "Rule": RESUME_RULE}
            self.log(f"{self.case}: adopted {name} paladin run {row['run_id']} ({row.get('status')})")
        record.setdefault("Identity", identity)
        self.save()
        return record["RunID"]

    def find_or_submit(self, stage, record):
        """RESUME_RULE + INTENT_RULE: the stage's run (adopted or submitted now), or PaladinError.  A
        submit is always followed by a lookup that must find exactly that run (a run recorded but
        not found = STOP after the repeated lookup, never a second submit)."""
        name = stage["Name"]
        identity = stage_identity(self.block["Lane"], self.case, name, self.plan_sha256)
        for _ in range(3):
            decision, value = self.lookup(stage, record, identity)
            if decision == "submit":
                if record.get("RunID") is not None:
                    decision, value = self.repeated_lookup(stage, record, identity,
                                                           f"the recorded run {record['RunID']} is not found")
                elif self.open_intent(record) is not None:
                    decision, value = self.repeated_lookup(stage, record, identity,
                                                           "a submit intent without an outcome (the driver ended inside "
                                                           "its submit)")
                else:
                    self.sleep(SUBMIT_LOOKUP_GAP_SECONDS)
                    second = self.lookup(stage, record, identity)
                    if second[0] == "adopt":
                        decision, value = second
                    elif second != (decision, value):
                        raise PaladinError(f"{name}: two lookups {SUBMIT_LOOKUP_GAP_SECONDS:.0f} s apart disagree "
                                           f"({value} / {second[1]} no-run rows): STOP, no submit ({INTENT_RULE})")
            if decision == "adopt":
                return self.adopt(record, value, identity)
            if value:
                record["NoRunRowsBeforeSubmit"] = value
            if self.submit(stage, record) is not None:
                continue  # the lookup must now find exactly this run
            decision, value = self.lookup(stage, record, identity)
            if decision == "submit" and value == 0:
                raise PaladinError(f"{name}: the submit was refused and left no row: STOP "
                                   f"({record['SubmitErrors'][-1]['Error']})")
        raise PaladinError(f"{name}: no run after its submits ({record.get('SubmitErrors')}; {RESUME_RULE})")

    def transport_fault(self, record, data):
        """Why a describe-run answer is a transport answer (WAIT_RULE), or None for a usable one."""
        if not isinstance(data, dict):
            return "not an object"
        if str(data.get("run_id")) != str(record["RunID"]):
            return f"run_id {data.get('run_id')!r}, not {record['RunID']}"
        status = str(data.get("status") or "").strip().lower()
        if not status:
            return "no status"
        if status in ("completed", "failed") and not (data.get("block") or {}).get("backend_data", {}).get("run_dir"):
            return f"'{status}' without its run directory"
        return None

    def describe_checked(self, record):
        """describe-run, re-asked after WAIT_BACKOFF_SECONDS when the answer is a transport answer;
        (data, None) or (None, the last fault)."""
        fault = None
        for attempt, backoff in enumerate((0.0, *WAIT_BACKOFF_SECONDS)):
            if attempt:
                self.sleep(backoff)
            try:
                data = self.describe(record["RunID"])
            except PaladinError as error:
                fault = f"describe-run failed: {str(error)[:200]}"
                continue
            fault = self.transport_fault(record, data)
            if fault is None:
                return data, None
        return None, fault

    def cap_clock(self, record, data):
        """(the time the stage cap counts from, its name) of a live answer (WAIT_RULE)."""
        started = parse_time(data.get("started"))
        if started:
            return started, "started"
        if str(data.get("status", "")).lower() in ("submitted", "queued"):
            return None, "not started"
        submitted = parse_time(record.get("Submitted") or (record.get("Adopted") or {}).get("UTC"))
        return submitted, ("submitted (a live answer without 'started')" if submitted else "none")

    def poll(self, stage, record):
        """One describe-run of the stage's run; the cap cancel; the data when the run ended.  A
        transport answer is counted and re-polled, never a State (WAIT_RULE)."""
        data, fault = self.describe_checked(record)
        if data is None:
            wait = record.setdefault("WaitTransport", {"Count": 0, "FirstUTC": self.now().isoformat()})
            wait["Count"] += 1
            wait["Last"] = {"UTC": self.now().isoformat(), "Fault": fault}
            self.save()
            span = (self.now() - parse_time(wait["FirstUTC"])).total_seconds()
            self.log(f"{self.case}: {stage['Name']} run {record['RunID']}: transport answer {wait['Count']} ({fault}); "
                     f"re-polled ({WAIT_RULE})")
            if wait["Count"] >= WAIT_TRANSPORT_POLLS and span >= WAIT_TRANSPORT_SECONDS:
                raise PaladinError(f"{stage['Name']}: run {record['RunID']}: {wait['Count']} failed polls over {span:.0f} s "
                                   f"({fault}): STOP with no State recorded; the run stays with its run watcher ({WAIT_RULE})")
            return None
        record.pop("WaitTransport", None)
        status = str(data.get("status", "")).lower()
        record["RowStatus"] = status
        if run_class(status) != "live":
            return data
        clock, clock_name = self.cap_clock(record, data)
        record["CapClock"] = clock_name
        cap = stage.get("CapSeconds")
        if cap and clock and (self.now() - clock).total_seconds() > float(cap) and not record.get("CapCancel", {}).get("OK"):
            cancel = record.setdefault("CapCancel", {"Attempts": 0, "CapSeconds": cap})
            cancel["Attempts"] += 1
            cancel["UTC"] = self.now().isoformat()
            cancel["From"] = clock_name
            try:
                paladin_call(self.block["CLI"], ["palace", "cancel", str(record["RunID"]), "--json"])
                cancel["OK"] = True
                self.log(f"{self.case}: {stage['Name']} run {record['RunID']} past its cap {cap} s from {clock_name}: ended early")
            except PaladinError as error:
                cancel["Error"] = str(error)[:300]
                if cancel["Attempts"] >= 3:
                    self.save()
                    raise PaladinError(f"{stage['Name']}: the cap cancel of run {record['RunID']} failed 3 times: cancel it "
                                       f"by hand ({error})")
        self.save()
        return None

    def finish(self, stage, record, data):
        """The reads after a run ended (module docstring); State complete, or PaladinError."""
        name = stage["Name"]
        status = str(data.get("status", "")).lower()
        run_dir = (data.get("block") or {}).get("backend_data", {}).get("run_dir")
        started, completed = parse_time(data.get("started")), parse_time(data.get("completed"))
        record.update(RunDir=run_dir, Started=data.get("started"), Completed=data.get("completed"), RowStatus=status,
                      Class=run_class(status))
        if started and completed:
            record["WallSeconds"] = (completed - started).total_seconds()
            record["NodeHours"] = record["WallSeconds"] / 3600.0
        cls = run_class(status)
        if cls in ("cancelled", "out-of-memory", "unknown"):
            record["State"] = ("timed-out" if cls == "cancelled" and record.get("CapCancel", {}).get("OK") else cls)
            self.save()
            raise PaladinError(f"{name}: run {record['RunID']} ended {status} (State {record['State']}): nothing is resubmitted "
                               f"({RESUME_RULE})")
        backend_data = (data.get("block") or {}).get("backend_data", {})
        partition = (record.get("Adopted", {}).get("Nodes", 1), backend_data.get("ranks_per_node_resolved", self.ranks))
        record["Partition"] = {"Nodes": partition[0], "RanksPerNode": partition[1]}
        if partition != (1, self.ranks):
            record["State"] = "failed"
            self.save()
            raise PaladinError(f"{name}: run {record['RunID']} ran on {partition[0]} node(s) x {partition[1]} ranks, the plan is "
                               f"1 x {self.ranks}: STOP")
        config = self.stage_run_config(stage)
        read = completion_read(self.host, run_dir, expected_outputs(config), row_status=status,
                               cancelled_by_backend=bool(record.get("CapCancel")))
        record["CompletionRead"] = read
        conditions = read["Conditions"]
        holds = read["Holds"] if cls == "failed" else all(conditions[key] for key in
                                                          ("1-ElapsedTimeReport", "2-RunnerFinishedExit0", "4-ExpectedOutputs"))
        if not holds:
            record["State"] = "failed"
            self.save()
            raise PaladinError(f"{name}: run {record['RunID']} ({status}) fails the completion read {conditions}: STOP, no read "
                               f"of its data, nothing resubmitted")
        if cls == "failed":
            record["CompletedByNote1"] = True
        record["PinCheckAfter"] = self.tool_line([self.block["PinCheck"], "postrun", "--record", self.block["PinRecord"],
                                                  "--run-id", str(record["RunID"])], r"POSTRUN-V2 PASS ",
                                                 f"pin_check_v2 postrun of {name}")
        q = shlex.quote(run_dir)
        ran = json.loads(remote_text(self.host, f"cat {q}/palace_config.json", what=f"palace_config.json of {run_dir}"))
        differing = leaf_diff(config, ran)
        record["LeafDiff"] = differing
        foreign = [path for path in differing if path not in PALADIN_OWNED_PATHS]
        if foreign:
            record["State"] = "failed"
            self.save()
            raise PaladinError(f"{name}: paladin ran a configuration that differs from the submitted one at {foreign}: STOP")
        addon = remote_text(self.host, f"grep -a -m1 'Paladin addon:' {q}/console.log || true", what="the addon line")
        record["PaladinAddon"] = addon.strip()
        output = config.get("Problem", {}).get("Output")
        if not output or not os.path.isabs(output):
            record["State"] = "failed"
            self.save()
            raise PaladinError(f"{name}: the stage configuration names no absolute Problem.Output ({output!r})")
        log_path = os.path.join(self.remote_dir, f"{name}.log")
        qo, ql = shlex.quote(output), shlex.quote(log_path)
        text = remote_text(self.host, f"mkdir -p {qo} {shlex.quote(self.remote_dir)} && cp -a {q}/postpro/. {qo}/ && "
                                      f"cp -a {q}/log.out {ql} && "
                                      f"(cd {q}/postpro && find . -type f -exec sha256sum {{}} + | sort -k2) && echo '--SPLIT--' && "
                                      f"(cd {qo} && find . -type f -exec sha256sum {{}} + | sort -k2)",
                           timeout=900, what=f"the copy of {run_dir}/postpro to {output}")
        source, target = (manifest(part) for part in text.split("--SPLIT--"))
        if source != target:
            record["State"] = "failed"
            self.save()
            raise PaladinError(f"{name}: the copy of {run_dir}/postpro into {output} differs from its source")
        record["Outputs"] = {"Directory": output, "Files": len(target), "Log": log_path}
        import run_stages  # the same qualify/ directory: ONE log parser
        record["Parsed"] = run_stages.parse_log(remote_text(self.host, f"cat {ql}", what=f"the log {log_path}"))
        archive = config.get("Solver", {}).get("Electrostatic", {}).get("ResponseArchive", {})
        if archive.get("Mode") == "Archive":
            record["Archive"] = self.archive_verify(archive["Directory"], config, record)
            if not record["Archive"]["OK"]:
                record["State"] = "failed"
                self.save()
                raise PaladinError(f"{name}: the archive {archive['Directory']} does not verify: "
                                   f"{ {k: v for k, v in record['Archive'].items() if k != 'Rule'} }")
        record["QdashUpload"] = self.qdash_record(run_dir, data)
        record["State"] = "complete"
        self.save()
        self.log(f"{self.case}: {name} run {record['RunID']} complete ({record.get('WallSeconds', 0):.0f} s)")

    def partition_check(self, stage, record):
        """Requires: every required stage completed, on this stage's partition (the archive is per rank)."""
        for required in stage.get("Requires", []):
            other = self.state["Stages"].get(required, {})
            if other.get("State") != "complete":
                raise PaladinError(f"{stage['Name']} requires {required}, which has not completed")
            if other.get("Partition") != {"Nodes": 1, "RanksPerNode": self.ranks}:
                raise PaladinError(f"{stage['Name']} reads the archive of {required}, which ran on {other.get('Partition')}, "
                                   f"not 1 x {self.ranks} (the archive is per rank)")

    def advance(self):
        """One step of the job: True when every stage is complete (status.json written)."""
        self.lock()
        self.preflight()
        for stage in self.plan["Stages"]:
            name = stage["Name"]
            record = self.state["Stages"].setdefault(name, {"Name": name})
            if record.get("State") == "complete":
                continue
            if record.get("State") is not None:
                raise PaladinError(f"{name}: State {record['State']} recorded: STOP ({RESUME_RULE})")
            self.partition_check(stage, record)
            if record.get("RunID") is None or name not in self.confirmed:
                # Every driver process confirms a recorded run by the lookup once (RESUME_RULE).
                self.find_or_submit(stage, record)
                self.confirmed.add(name)
            data = self.poll(stage, record)
            if data is None:
                return False
            self.finish(stage, record, data)
        self.write_status()
        return True

    def write_status(self):
        """status.json of run_stages' schema in the plan's directory on /data."""
        records = [self.state["Stages"][stage["Name"]] for stage in self.plan["Stages"]]
        stages = []
        for stage, record in zip(self.plan["Stages"], records):
            stages.append({"Name": stage["Name"], "Config": stage["Config"], "Environment": stage.get("Environment", {}),
                           "CapSeconds": stage.get("CapSeconds"), "Log": record["Outputs"]["Log"],
                           "StartUTC": record.get("Started"), "EndUTC": record.get("Completed"),
                           "Binary": self.block["ImageBinary"], "BinarySHA256": self.block["ImageBinarySHA256"],
                           "ReturnCode": 0, "WallSeconds": record.get("WallSeconds"), "TimedOut": False,
                           "Parsed": record.get("Parsed", {}), "State": "complete",
                           "Paladin": {key: value for key, value in record.items() if key != "Parsed"}})
        starts = [parse_time(record.get("Started")) for record in records if record.get("Started")]
        ends = [parse_time(record.get("Completed")) for record in records if record.get("Completed")]
        status = {"Version": 2, "Plan": os.path.join(self.remote_dir, "plan.json"), "Case": self.case, "Backend": "paladin",
                  "PBSJobID": None, "PaladinRunIDs": [record["RunID"] for record in records], "Host": None,
                  "StartUTC": min(starts).strftime("%Y-%m-%dT%H:%M:%SZ") if starts else None,
                  "EndUTC": max(ends).strftime("%Y-%m-%dT%H:%M:%SZ") if ends else None, "State": "complete",
                  "Preflight": self.state.get("Preflight"),
                  "TotalSeconds": sum(record.get("WallSeconds") or 0.0 for record in records), "Stages": stages}
        target = os.path.join(self.remote_dir, "status.json")
        remote_text(self.host, f"mkdir -p {shlex.quote(self.remote_dir)} && cat > {shlex.quote(target)}",
                    stdin=json.dumps(status, indent=1) + "\n", what=f"the write of {target}")
        self.state["StatusJSON"] = target
        self.state["NodeHours"] = sum(record.get("NodeHours") or 0.0 for record in records)
        self.save()
        return status


def drive(job, *, poll_seconds, sleep=time.sleep):
    """advance() until the job completes (the standalone form; qualify_library polls itself); the
    driver lock released at the end, whatever the end."""
    try:
        while not job.advance():
            sleep(poll_seconds)
    finally:
        job.unlock()
    return job.state


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("plan")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--dry-run", action="store_true")
    mode.add_argument("--submit", action="store_true")
    parser.add_argument("--profile", required=True, help="a cluster profile with a Paladin block")
    parser.add_argument("--out", required=True, help="the local directory of the job's state and run configurations")
    parser.add_argument("--remote-dir", help="the plan's directory on /data (outputs, logs, status.json); --submit")
    parser.add_argument("--map", action="append", default=[], help="REMOTE_PREFIX=LOCAL_PREFIX of the stage configurations")
    parser.add_argument("--archive-gb", action="append", default=[], help="DIRECTORY=GB: an archive writer's estimate")
    parser.add_argument("--poll-seconds", type=float, default=60.0)
    args = parser.parse_args(argv)
    with open(args.profile) as f:
        block = backend_of("paladin", json.load(f))
    mapping = [tuple(entry.split("=", 1)) for entry in args.map]
    estimates = {key: float(value) for key, value in (entry.rsplit("=", 1) for entry in args.archive_gb)}
    os.makedirs(args.out, exist_ok=True)
    print(route_line(block))
    if args.dry_run:
        job = PaladinJob(args.plan, block=block, remote_dir=args.remote_dir or "", local_dir=args.out, mapping=mapping,
                         archive_estimates_gb=estimates)
        runs = job.dry_run()
        with open(os.path.join(args.out, "paladin-runs.json"), "w") as f:
            json.dump(runs, f, indent=1)
        for run in runs:
            print(json.dumps({"Stage": run["Stage"], "ResponseArchive": run["ResponseArchive"], "Requires": run["Requires"],
                              "Instance": run["Instance"]["Type"]}))
        return 0
    if not args.remote_dir:
        parser.error("--submit needs --remote-dir")
    code = 0
    job = None
    try:
        job = PaladinJob(args.plan, block=block, remote_dir=args.remote_dir, local_dir=args.out, mapping=mapping,
                         archive_estimates_gb=estimates)
        drive(job, poll_seconds=args.poll_seconds)
    except PaladinError as error:
        print(f"STOP: {error}", file=sys.stderr)
        code = 3
    if job is not None:
        for name, record in job.state["Stages"].items():
            print(json.dumps({"Stage": name, "RunID": record.get("RunID"), "State": record.get("State"),
                              "NodeHours": record.get("NodeHours")}))
    return code


if __name__ == "__main__":
    sys.exit(main())
