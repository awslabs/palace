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
    cancelled (recorded; the row then reads canceled: the stage STOPS as timed-out); an
    unattended run watcher on this host (WATCHER_RULE) and a qdel-only PBS-side backstop on the
    cluster's ssh host (BACKSTOP_RULE) keep the limit when the driver or this host is gone;
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


class PaladinUnanswered(PaladinError):
    """A paladin call without an answer (the timeout, no JSON envelope): what it did is unknown (UNKNOWN_RECEIPT_RULE)."""


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
# The Paladin block's optional key (BACKSTOP_RULE): absent = the backstop with these settings; {"Mode": "none"} = the profile
# says explicitly that there is none.
BACKSTOP_DEFAULTS = {"Mode": "qdel", "GraceSeconds": 300.0, "IntervalSeconds": 60.0}
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
    unknown = sorted(set(block) - set(PALADIN_KEYS) - {"Backstop"})
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
    backstop_of(block)  # BACKSTOP_RULE: a malformed Backstop key is refused with the block
    with open(block["PinRecord"]) as f:
        record = json.load(f)
    if record.get("selector") != block["Selector"]:
        raise ValueError(f"the pin record {block['PinRecord']} is of selector {record.get('selector')!r}, "
                         f"not Paladin.Selector {block['Selector']!r}")
    return dict(block)


def backstop_of(block):
    """BACKSTOP_RULE: the backstop settings of a Paladin block ({Mode, GraceSeconds, IntervalSeconds}), or ValueError."""
    given = block.get("Backstop")
    if given is None:
        return dict(BACKSTOP_DEFAULTS)
    if given == {"Mode": "none"}:
        return {"Mode": "none"}
    if not isinstance(given, dict) or given.get("Mode") != "qdel" or set(given) - set(BACKSTOP_DEFAULTS):
        raise ValueError(f"Paladin.Backstop must be absent, {{\"Mode\": \"none\"}} or {{\"Mode\": \"qdel\"[, \"GraceSeconds\": S]"
                         f"[, \"IntervalSeconds\": S]}}, not {given!r}")
    settings = {**BACKSTOP_DEFAULTS, **given}
    for key in ("GraceSeconds", "IntervalSeconds"):
        if isinstance(settings[key], bool) or not isinstance(settings[key], (int, float)) or settings[key] < 0:
            raise ValueError(f"Paladin.Backstop.{key} must be a number >= 0, not {settings[key]!r}")
    settings["GraceSeconds"], settings["IntervalSeconds"] = float(settings["GraceSeconds"]), float(settings["IntervalSeconds"])
    return settings


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
            f"in-image {block['ImageBinary']} {block['ImageBinarySHA256'][:16]}, backstop {backstop_of(block)['Mode']}")


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
        raise PaladinUnanswered(f"{' '.join(arguments[:2])}: no answer in {timeout} s") from None
    try:
        envelope = json.loads(completed.stdout)
    except ValueError:
        raise PaladinUnanswered(f"{' '.join(arguments[:2])}: no JSON envelope (rc {completed.returncode})") from None
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
UNKNOWN_RECEIPT_RULE = ("decision 1836 m-2 (review r1833 MINOR-2): a submit call without an answer (the timeout, no JSON envelope) or "
                        "with an ok answer without a usable jobs[0].run_id keeps its intent OPEN (its receipt is unknown: the server "
                        "may have made the run): the lookup is repeated 3 times over >= 60 s, a run found = adopted, none = STOP; an "
                        "ANSWERED refusal (ok false) closes the intent 'refused' and the NO-RUN rule applies")
CLOSE_INTENT_RULE = ("decision 1836 m-3 (review r1833 MINOR-3): an open intent whose run never appears is closed ONLY by "
                     "`paladin_backend.py PLAN --profile P --out DIR --close-intent STAGE`: under the driver lock, 3 recorded lookups "
                     "30 s apart, each complete (LOOKUP_RULE) and showing no run of the stage; then Outcome 'closed: no run after the "
                     "recorded lookups' with ClosedBy; a run listed = not closed (a resume adopts it); no hand edit of the status file")
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


class SessionNotValid(PaladinError):
    """The session check before a submit did not pass (SESSION_RULE): nothing was submitted."""


class WatcherNotStarted(PaladinError):
    """A run whose watcher did not start (WATCHER_RULE); .watcher = its record."""

    def __init__(self, message, watcher):
        super().__init__(message)
        self.watcher = watcher


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


OUTPUTS_RULE = ("PREREG-ESSURR-NOTE-1 (4), decision 1820 m-1: the postpro names every completed run of a stage type writes: "
                "palace.json and palace_config_resolved.json; Mode Archive nothing else; Mode Reduce the two response matrices; a "
                "full solve error-indicators.csv, the two matrices (Solver.Electrostatic.ResponseMatrix), domain-E.csv "
                "(Domains.Postprocessing.Energy), surface-Q.csv (Boundaries.Postprocessing.Dielectric), surface-Q-edge.csv "
                "(LocalizeEdgeEnergy) and surface-Q-edge-local.csv (SaveLocalEdgeEnergy); any other stage type has no list and "
                "is not submitted")
MATRICES = ["domain-response-matrix.csv", "surface-response-matrix.csv"]


def expected_outputs(config):
    """OUTPUTS_RULE: the sorted postpro names of the stage configuration, or PaladinError (no list)."""
    names = ["palace.json", "palace_config_resolved.json"]
    electrostatic = config.get("Solver", {}).get("Electrostatic", {})
    mode = electrostatic.get("ResponseArchive", {}).get("Mode")
    if config.get("Problem", {}).get("Type") != "Electrostatic" or mode not in (None, "Archive", "Reduce"):
        raise PaladinError(f"no expected-output list for a {config.get('Problem', {}).get('Type')} stage of Mode {mode} "
                           f"({OUTPUTS_RULE})")
    if mode == "Reduce":
        names += MATRICES
    elif mode is None:
        dielectric = config.get("Boundaries", {}).get("Postprocessing", {}).get("Dielectric", [])
        names.append("error-indicators.csv")
        if electrostatic.get("ResponseMatrix"):
            names += MATRICES
        if config.get("Domains", {}).get("Postprocessing", {}).get("Energy"):
            names.append("domain-E.csv")
        if dielectric:
            names.append("surface-Q.csv")
        if any(entry.get("LocalizeEdgeEnergy") for entry in dielectric):
            names.append("surface-Q-edge.csv")
        if any(entry.get("SaveLocalEdgeEnergy") for entry in dielectric):
            names.append("surface-Q-edge-local.csv")
    return sorted(names)


UPLOAD_HARNESS_ERROR = "qdash_das.palace.runner.process_harness: harness upload"
COMPLETION_RULE = ("PREREG-ESSURR-NOTE-1 (decision 1836 m-6: exactly): a row failed by the qdash upload ALONE is a completed run "
                   "iff (1) log.out has 'Elapsed Time Report (s)', (2) console.log has 'Runner finished (exit 0)', (3) EVERY ERROR "
                   "line of console.log is the upload harness's ('" + UPLOAD_HARNESS_ERROR + "'; at least one), every traceback "
                   "follows such a line, ONE traceback per such line (decision 1844: a second traceback without a new harness ERROR "
                   "line, a chained one included, is an error other than the upload's) and every exception line ends such a "
                   "traceback (its last line judged and quoted), (4) the "
                   "postpro file names equal OUTPUTS_RULE's list (every difference reported; an untraced one = not shown), (5) no "
                   "cap / memory stop fired (no cancel by this backend, no end.json of the run watcher, its events.log without a CAP "
                   "/ MEMORY line and its END line read: 'END exit 0' (or 5, its life bound: no stop); paladin's own guard ends a row out_of_memory, not failed); a "
                   "completed row must show (1), (2) and every expected name (extras reported)")
MARKED_LINES_LIMIT = 2000


def judge_console(marks, error_count):
    """COMPLETION_RULE (3) on the marked console lines [(line number, text)]: (holds, offending lines, the last exception line)."""
    offending, anchor, traceback_ok, last = [], None, False, None
    errors = [text for _, text in marks if "ERROR" in text]
    if error_count != len(errors) or len(marks) >= MARKED_LINES_LIMIT:
        offending.append(f"{error_count} ERROR lines, {len(errors)} read ({len(marks)} marked lines): not all judged")
    for _, text in marks:
        if "ERROR" in text:
            anchor = UPLOAD_HARNESS_ERROR in text
            if not anchor:
                offending.append(text)
        elif text.startswith("Traceback (most recent call last)"):
            traceback_ok, anchor = bool(anchor), False  # one traceback per harness ERROR line (decision 1844)
            if not traceback_ok:
                offending.append(text)
        else:  # an exception line: it must end a traceback of the harness's (one exception line per traceback)
            last = text
            if not traceback_ok:
                offending.append(text)
            traceback_ok = False
    return bool(errors) and not offending, offending[:5], last


def completion_read(host, run_dir, expected, *, row_status, cancelled_by_backend, watcher_end=None, watcher_events=None):
    """The five conditions of COMPLETION_RULE read on the run directory (watcher_events: the run watcher's
    events.log text, None = none found)."""
    q = shlex.quote(run_dir)
    text = remote_text(host, f"cd {q} && echo ELAPSED $(grep -c 'Elapsed Time Report (s)' log.out 2>/dev/null); "
                             f"echo FINISHED $(grep -a -c 'Runner finished (exit 0)' console.log 2>/dev/null); "
                             f"echo ERRORCOUNT $(grep -a -c 'ERROR' console.log 2>/dev/null); "
                             f"grep -a -n -E 'ERROR|^Traceback \\(most recent call last\\)|^[A-Za-z_.]+(Error|Exception)' console.log "
                             f"2>/dev/null | head -{MARKED_LINES_LIMIT} | sed 's/^/MARK /'; "
                             f"find postpro -maxdepth 1 -type f -printf 'FILE %f\\n' 2>/dev/null",
                       what=f"the completion read of {run_dir}")
    lines = text.splitlines()
    value = {line.split()[0]: line.split()[1] if len(line.split()) > 1 else "0" for line in lines
             if line.startswith(("ELAPSED ", "FINISHED ", "ERRORCOUNT "))}
    marks = []
    for line in lines:
        if line.startswith("MARK "):
            number, _, rest = line[len("MARK "):].partition(":")
            marks.append((int(number) if number.isdigit() else 0, rest))
    errors = [rest for _, rest in marks if "ERROR" in rest]
    judged, offending, last = judge_console(marks, int(value.get("ERRORCOUNT", "0") or 0))
    files = {line[len("FILE "):] for line in lines if line.startswith("FILE ")}
    missing, extra = sorted(set(expected) - files), sorted(files - set(expected))
    failed = str(row_status).lower() == "failed"
    events = (watcher_events or "").splitlines()
    stops = [line for line in events if len(line.split()) > 1 and line.split()[1] in ("CAP", "MEMORY")]
    end_lines = [line for line in events if len(line.split()) > 2 and line.split()[1:3] == ["END", "exit"]]
    conditions = {
        "1-ElapsedTimeReport": int(value.get("ELAPSED", "0") or 0) >= 1,
        "2-RunnerFinishedExit0": int(value.get("FINISHED", "0") or 0) >= 1,
        "3-OnlyUploadErrors": judged if failed else True,
        "4-ExpectedOutputs": not missing and (not extra if failed else True),
        "5-NoStop": (not cancelled_by_backend and watcher_end is None and str(row_status).lower() != "out_of_memory" and
                     not stops and bool(end_lines) and all(line.split()[3].rstrip(":") in ("0", "5") for line in end_lines)),
    }
    return {"Conditions": conditions, "Holds": all(conditions.values()), "ErrorLines": errors[:20], "ErrorCount": len(errors),
            "Offending": offending, "TracebackLast": last, "MissingOutputs": missing, "ExtraOutputs": extra,
            "NameDifferences": {"Missing": missing, "Extra": extra}, "WatcherEndLines": end_lines[-3:], "WatcherStops": stops[:3],
            "RowStatus": row_status, "Rule": COMPLETION_RULE}


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


DATA_GUARD_RULE = ("decision 1820 m-2 (CTR-HOWTO-v2 s9): before EVERY submit the plan's directory on /data is created and its "
                   "file system must have >= 150 GB free (an archive writer's guard adds its own checks); a standalone --submit "
                   "names its --lane-root, and every archive Directory of its plan and of --archive-gb must lie under it (m-3)")


def under(path, root):
    """True when the absolute `path` is `root` or lies under it (normalised strings; no '..' escapes)."""
    path, root = os.path.normpath(path), os.path.normpath(root)
    return path == root or path.startswith(root.rstrip("/") + "/")


def archive_name(source, rank):
    return f"source-{int(source):06d}-rank-{int(rank):06d}-V.bin"


WATCHER_RULE = ("decision 1820 (4) / 1836: every run has an UNATTENDED run watcher (paladin_run_watcher.py; nohup, its own "
                "session, its own lock, pid file, events.log and end.json under <the plan's directory>/paladin-watchers/run-<id>: one "
                "root per plan, so two --out directories share one lock per run) started BEFORE its submit call (nonce-<N>: "
                "the submit window, NONCE_RULE), re-checked at every live poll, and for an adopted run without a live one: it ends the run early at the stage's CapSeconds (from the PBS start, "
                "else the submit time) and at the memory bar (0.85 of MemTotal, measured); an empty or failed PBS read is a missing "
                "read, never 'ended', and the limits stay in force under it (its READ_RULE); the early end = 'palace cancel', then "
                "qdel by the PBS job id (the primary path when the session has expired: no session needed), repeated with a "
                "back-off until VERIFIED (the row terminal or a finished-job answer; its END_RULE); an end not verified leaves "
                "END-UNVERIFIED.json, on which every advance STOPs; the driver's own cap stays as a second line; a run whose "
                "watcher does not start is ended early")
WATCHER_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "paladin_run_watcher.py")
UNVERIFIED_MARKER = "END-UNVERIFIED.json"  # paladin_run_watcher.UNVERIFIED_MARKER (its END_RULE)
WATCHER_INTERVAL_SECONDS, WATCHER_GRACE_SECONDS, WATCHER_MEMORY_BAR = 30.0, 120.0, 0.85
# Decision 1830 (3): hard bounds in code, whatever the checks around them: watcher starts per run id per driver process; the
# standalone drive loop's advances and its minimum sleep (the re-poll loops are bounded where they are: 3 describe-run asks with
# 10 / 20 s back-off, 3 lookups 30 s apart, 3 submit rounds).
WATCHER_STARTS_PER_RUN, DRIVE_MAX_ADVANCES, DRIVE_MIN_POLL_SECONDS = 3, 20000, 5.0
DEFAULT_PBS_BIN = "/opt/soca/soca-green1/schedulers/default/pbs/bin"


class RunWatchers:
    """Starts paladin_run_watcher.py processes and tells whether one holds a directory (WATCHER_RULE)."""

    def __init__(self):
        self.processes = []  # kept while the driver lives (a watcher outlives it)

    def alive(self, directory):
        path = os.path.join(directory, "watcher.lock")
        if not os.path.exists(path):
            return False
        descriptor = os.open(path, os.O_RDWR)
        try:
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            return True
        finally:
            os.close(descriptor)  # closing it also drops a lock taken just now
        return False

    def start(self, argv, directory, *, wait_seconds=15.0):
        """The watcher's pid once it holds its lock, else None."""
        os.makedirs(directory, exist_ok=True)
        with open(os.path.join(directory, "watcher.out"), "ab") as out:
            process = subprocess.Popen(["nohup", *argv], stdin=subprocess.DEVNULL, stdout=out, stderr=subprocess.STDOUT,
                                       start_new_session=True, close_fds=True)
        self.processes.append(process)
        deadline = time.time() + wait_seconds
        while time.time() < deadline:
            if self.alive(directory):
                return process.pid
            if process.poll() is not None:
                return None
            time.sleep(0.1)
        return None


BACKSTOP_RULE = ("decision 1858 (review r1848 (c) 1): the run watcher lives on the driver's host; a qdel-ONLY BACKSTOP per run "
                 "runs on the cluster's ssh host (paladin_backstop.py, setsid nohup there; its state file, log and end record on "
                 "/data under <the plan's directory>/paladin-backstops/run-<id>): ONE job id, refused unless the job's owner and "
                 "name match, a qdel only at the PBS start + the stage's CapSeconds + GraceSeconds, ended when the job has left "
                 "PBS (F / X / an Exit_status / 'Unknown Job Id'; a missing read is never 'left'). Before EVERY submit the script "
                 "is installed on Host (its sha256 read back) and its --check passes, else STOP with no intent and no submit "
                 "(fail-closed); after the submit answer it is started and must show its first read within "
                 "BACKSTOP_START_WAIT_SECONDS, else the run is ended early and the drive STOPs; at every live poll, adoption "
                 "and resume its status is read: a failed read is re-read at the next poll, a backstop gone without the job "
                 "having left is restarted (<= BACKSTOP_STARTS_PER_RUN starts per run per driver process), a refusal or a "
                 "failed restart STOPs (the run stays with its run watcher); a run ended by the backstop's qdel is State "
                 "timed-out, nothing resubmitted. A profile whose Paladin.Backstop is {\"Mode\": \"none\"} has none (recorded)")
BACKSTOP_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "paladin_backstop.py")
BACKSTOP_STARTS_PER_RUN, BACKSTOP_START_WAIT_SECONDS, BACKSTOP_START_READ_SECONDS = 3, 30.0, 2.0
BACKSTOP_LIFE_SECONDS = 8 * 86400.0  # the run watcher's bound mirrored: 7 days of queue + 1 day beyond the cap


class BackstopNotStarted(PaladinError):
    """The backstop of a run did not start after its submit (BACKSTOP_RULE): the run was ended early."""


class ClusterBackstops:
    """The backstops on the cluster's ssh host (BACKSTOP_RULE), through remote_shell; `sleep` between the start reads."""

    def __init__(self, sleep=time.sleep):
        self.sleep = sleep

    def install(self, host, root, pbs_bin):
        """The script copied to root (its sha256 read back) and its --check: {Script, SHA256, User, Host}, or PaladinError."""
        with open(BACKSTOP_PATH) as f:
            source = f.read()
        local = hashlib.sha256(source.encode()).hexdigest()
        script = f"{root}/paladin_backstop.py"
        q = shlex.quote
        text = remote_text(host, f"mkdir -p {q(root)} && cat > {q(script)}.tmp && mv {q(script)}.tmp {q(script)} && "
                                 f"sha256sum {q(script)}", stdin=source, what=f"the backstop's copy to {root}")
        remote = (text.split() or [""])[0]
        if remote != local:
            raise PaladinError(f"the backstop's copy {script} hashes to {remote[:16]}, the chain's file to {local[:16]}")
        line = remote_text(host, f"python3 {q(script)} --check {q(root)} --pbs-bin {q(pbs_bin)}",
                           what="the backstop's --check").strip()
        fields = dict(part.split("=", 1) for part in line.split()[2:] if "=" in part)
        if not line.startswith("BACKSTOP-CHECK OK ") or not fields.get("user") or not fields.get("host"):
            raise PaladinError(f"the backstop's --check on {host or 'this host'} answers {line[:200]!r}")
        return {"Script": script, "SHA256": local, "User": fields["user"], "Host": fields["host"]}

    def status(self, host, directory):
        """{Alive, Host, FirstRead, End} of the backstop in directory, or None (the read failed: nothing is known)."""
        q = shlex.quote(directory)
        command = (f"hostname; pid=$(cat {q}/backstop.pid 2>/dev/null); "
                   f"if [ -n \"$pid\" ] && kill -0 \"$pid\" 2>/dev/null && grep -q paladin_backstop /proc/$pid/cmdline 2>/dev/null; "
                   f"then echo ALIVE \"$pid\"; else echo DEAD \"${{pid:--}}\"; fi; "
                   f"grep -c ' FIRST-READ ' {q}/backstop.log 2>/dev/null || true; echo ---END---; "
                   f"cat {q}/backstop.end.json 2>/dev/null || true")
        try:
            result = remote_shell(host, command, timeout=120)
        except PaladinError:
            return None
        if result.returncode != 0 or "---END---" not in result.stdout:
            return None
        head, tail = result.stdout.split("---END---", 1)
        lines = head.split()
        try:
            end = json.loads(tail) if tail.strip() else None
        except ValueError:
            end = None
        if len(lines) < 3:
            return None
        return {"Host": lines[0], "Alive": lines[1] == "ALIVE", "PID": lines[2], "FirstRead": lines[3:4] not in ([], ["0"]),
                "End": end}

    def start(self, host, directory, script, state):
        """The backstop started in directory (an earlier end record kept aside); its pid once it has read the job, else
        None (it did not start, it ended, or it refused)."""
        q = shlex.quote
        d = q(directory)
        remote_text(host, f"mkdir -p {d} && if [ -f {d}/backstop.end.json ]; then mv {d}/backstop.end.json "
                          f"{d}/backstop.end.$(date -u +%Y%m%dT%H%M%SZ).json; fi && rm -f {d}/backstop.pid && "
                          f"cat > {d}/state.json.tmp && mv {d}/state.json.tmp {d}/state.json && "
                          f"(setsid nohup python3 {q(script)} {d}/state.json > {d}/nohup.out 2>&1 < /dev/null &) && echo STARTED",
                    stdin=json.dumps(state, indent=1), what=f"the backstop's start in {directory}")
        for attempt in range(int(BACKSTOP_START_WAIT_SECONDS / BACKSTOP_START_READ_SECONDS) + 1):  # a count bound
            if attempt:
                self.sleep(BACKSTOP_START_READ_SECONDS)
            status = self.status(host, directory)
            if status and status["End"]:
                return None
            if status and status["Alive"] and status["FirstRead"]:
                return int(status["PID"])
        return None


SESSION_RULE = ("decision 1822 (b): before every paladin submit the Midway session's expiry must be >= max(1800 s, the "
                "stage's CapSeconds + 3600 s) ahead; else STOP before the submit, with no intent and no State (the stage stays "
                "resumable after mwinit). WHAT IS READ (decision 1844): each line of ~/.midway/cookie is split into its fields and "
                "ONLY two are used: the NAME field (to select the session-named entries) and the EXPIRY field; the VALUE field "
                "is never used, printed or recorded. Recorded: the earliest expiry, the seconds ahead, the seconds needed and "
                "the count of session entries; never a name or a value")
MIDWAY_COOKIE = "~/.midway/cookie"
SESSION_MIN_AHEAD_SECONDS, SESSION_BEYOND_CAP_SECONDS = 1800.0, 3600.0


def session_check(cookie, need_seconds, *, now=None):
    """SESSION_RULE on the Netscape cookie file `cookie`: {ExpiryUTC, AheadSeconds, NeedSeconds, Entries, Verdict}."""
    expiries = []
    try:
        with open(os.path.expanduser(cookie), errors="replace") as f:
            for line in f:
                if line.startswith("#") and not line.startswith("#HttpOnly_"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 7 and "session" in fields[5].lower() and fields[4].strip().isdigit():
                    expiries.append(int(fields[4]))
    except OSError:
        pass
    current = (now or utc_now)().timestamp()
    if not expiries:
        return {"ExpiryUTC": None, "AheadSeconds": None, "NeedSeconds": need_seconds, "Entries": 0, "Verdict": "NONE",
                "Rule": SESSION_RULE}
    expiry = min(expiries)
    ahead = expiry - current
    return {"ExpiryUTC": datetime.datetime.fromtimestamp(expiry, datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "AheadSeconds": round(ahead), "NeedSeconds": need_seconds, "Entries": len(expiries),
            "Verdict": "OK" if ahead >= need_seconds else "SHORT", "Rule": SESSION_RULE}


# The tables, logs and meshes paladin 0.2.16 uploads: the run's top level and each AMR iteration's
# top level (ERRATUM-1802 E-1), never other sub-directories or suffixes.
UPLOADED_SUFFIXES = (".csv", ".json", ".mesh", ".meshgz")


class PaladinJob:
    """The paladin runs of ONE run_stages plan, in order (module docstring)."""

    def __init__(self, plan_path, *, block, remote_dir, local_dir, mapping=(), archive_estimates_gb=None,
                 log=print, now=utc_now, sleep=None, exclusive=True, watchers=None, pbs_bin=DEFAULT_PBS_BIN,
                 session_cookie=None, lane_root=None, backstops=None):
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
        self.starter, self.pbs_bin = watchers, pbs_bin  # None = the default RunWatchers, built at its first use
        self.watcher_starts = {}  # run id -> watcher starts in this driver process (<= WATCHER_STARTS_PER_RUN)
        self.session_cookie = session_cookie
        self.lane_root = lane_root  # m-3: the standalone form's root on /data (None: qualify_library's own roots)
        self.backstop = backstop_of(block)  # BACKSTOP_RULE
        self.backstop_starter = backstops  # None = the default ClusterBackstops, built at its first use
        self.backstop_starts = {}  # run id -> backstop starts in this driver process (<= BACKSTOP_STARTS_PER_RUN)
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

    @property
    def watchers(self):
        """The watcher starter (WATCHER_RULE): the caller's, else the default RunWatchers, built at its first use (a dry
        run never builds it)."""
        if self.starter is None:
            self.starter = RunWatchers()
        return self.starter

    @property
    def backstops(self):
        """The backstop starter (BACKSTOP_RULE): the caller's, else the default ClusterBackstops, built at its first use."""
        if self.backstop_starter is None:
            self.backstop_starter = ClusterBackstops()
        return self.backstop_starter

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

    def data_guard(self):
        """DATA_GUARD_RULE before a submit: the record of the free space of the plan's directory, or PaladinError."""
        q = shlex.quote(self.remote_dir)
        text = remote_text(self.host, f"mkdir -p {q} && df -Pk {q} | tail -1", what=f"the /data guard of {self.remote_dir}")
        free = int(text.split()[3]) * 1024
        guard = {"Directory": self.remote_dir, "FreeBytes": free, "UTC": self.now().isoformat(), "Rule": DATA_GUARD_RULE}
        if free < ARCHIVE_FREE_BYTES:
            raise PaladinError(f"the /data guard refuses the submit: {free / 1e9:.1f} GB free under {self.remote_dir} "
                               f"(< 150 GB; {DATA_GUARD_RULE})")
        return guard

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
        if not stage.get("CapSeconds"):
            raise PaladinError(f"{name}: no CapSeconds: a paladin run is submitted only with a cap ({WATCHER_RULE})")
        session = session_check(self.session_cookie or MIDWAY_COOKIE,
                                max(SESSION_MIN_AHEAD_SECONDS, float(stage["CapSeconds"]) + SESSION_BEYOND_CAP_SECONDS),
                                now=self.now)
        record.setdefault("SessionChecks", []).append({**session, "UTC": self.now().isoformat()})
        self.save()
        if session["Verdict"] != "OK":
            raise SessionNotValid(f"{name}: the Midway session check reads {session['Verdict']} (expiry {session['ExpiryUTC']}, "
                                  f"{session['AheadSeconds']} s ahead, {session['NeedSeconds']:.0f} s needed): STOP before the "
                                  f"submit; run mwinit, then resume ({SESSION_RULE})")
        path, config = self.write_run_config(stage)
        record["RunConfig"] = {"Path": path, "SHA256": sha256_of(path)}
        record["ExpectedOutputs"] = expected_outputs(config)  # OUTPUTS_RULE: a stage type without a list is not submitted
        record["SchemaCheck"] = self.tool_line([self.block["SchemaCheck"], self.block["Selector"], path], r"VALID ",
                                               f"the schema check of {name}")
        mode = config.get("Solver", {}).get("Electrostatic", {}).get("ResponseArchive", {})
        if self.lane_root and mode.get("Directory") and not under(mode["Directory"], self.lane_root):
            raise PaladinError(f"{name}: the archive {mode['Directory']} lies outside the lane root {self.lane_root} "
                               f"({DATA_GUARD_RULE})")
        record["DataGuard"] = self.data_guard()
        if mode.get("Mode") == "Archive":
            record["ArchiveGuard"] = self.archive_guard(mode["Directory"])
        elif mode.get("Directory"):
            self.refuse_unusable_archive(name, mode["Directory"])
        record["PinCheckBefore"] = self.tool_line([self.block["PinCheck"], "check", "--record", self.block["PinRecord"]],
                                                  r"PASS ", f"pin_check_v2 check before {name}")
        record["BackstopCheck"] = self.backstop_check(name)  # BACKSTOP_RULE: fail-closed, before the intent
        self.save()
        record["Instance"] = self.sizing
        record["Nodes"], record["RanksPerNode"] = 1, self.ranks
        record["Identity"] = stage_identity(self.block["Lane"], self.case, name, self.plan_sha256)
        intent = {"Nonce": uuid.uuid4().hex, "Stage": name, "PlanSHA256": self.plan_sha256, "UTC": self.now().isoformat(),
                  "Rule": INTENT_RULE}
        record.setdefault("SubmitIntents", []).append(intent)
        self.save()  # on disk before the call (INTENT_RULE)
        self.start_nonce_watcher(stage, record, intent)  # WATCHER_RULE: the submit window is watched too
        try:
            data = paladin_call(self.block["CLI"], self.submit_arguments(stage, path, intent["Nonce"]))
            run_id = int(data["jobs"][0]["run_id"])
        except (PaladinUnanswered, KeyError, IndexError, TypeError, ValueError) as error:
            # UNKNOWN_RECEIPT_RULE: the intent stays OPEN; the caller repeats the lookup (never a second submit on it).
            record.setdefault("SubmitErrors", []).append({"UTC": self.now().isoformat(), "Error": repr(error)[:300],
                                                          "Receipt": "unknown"})
            intent["Unanswered"] = {"UTC": self.now().isoformat(), "Rule": UNKNOWN_RECEIPT_RULE}
            self.save()
            self.log(f"{self.case}: {name}: the submit call has no usable answer ({repr(error)[:120]}): the intent stays open")
            return None
        except PaladinError as error:
            record.setdefault("SubmitErrors", []).append({"UTC": self.now().isoformat(), "Error": str(error)[:300],
                                                          "Receipt": "refused"})
            intent["Outcome"] = "refused"
            self.save()
            return None
        record["RunID"], record["Submitted"] = run_id, self.now().isoformat()
        record["JobID"] = data["jobs"][0].get("job_id")
        intent["Outcome"], intent["RunID"] = f"run {run_id}", run_id
        self.save()
        self.write_run_hint(intent, run_id)
        self.watch(stage, record)
        self.log(f"{self.case}: submitted {name} as paladin run {run_id} ({self.sizing['Type']} x {self.ranks})")
        return run_id

    def ensure_watcher(self, name, run_id, *, job_id, submitted, cap_seconds, nonce_dirs=()):
        """WATCHER_RULE: a missing watcher started (its record), None when a live one holds the run's directory
        or one of `nonce_dirs` (the run's nonce watcher: NONCE_RULE of paladin_run_watcher.py)."""
        # ONE directory, hence ONE lock and ONE watcher, per run id (decision 1830 (3)).
        directory = os.path.join(self.watchers_root(), f"run-{run_id}")
        if any(self.watchers.alive(path) for path in (directory, *nonce_dirs)):
            return None
        argv = [sys.executable, WATCHER_PATH, "--run-id", str(run_id), "--job", str(job_id or ""), "--cli", self.block["CLI"],
                "--host", self.block["Host"], "--pbs-bin", self.pbs_bin, "--cap-seconds", str(float(cap_seconds)),
                "--submitted", submitted, "--bar", str(WATCHER_MEMORY_BAR), "--interval", str(WATCHER_INTERVAL_SECONDS),
                "--grace", str(WATCHER_GRACE_SECONDS), "--dir", directory]
        pid, attempt = None, 0
        for attempt in range(1, 3):
            if self.watcher_starts.get(run_id, 0) >= WATCHER_STARTS_PER_RUN:
                raise PaladinError(f"{name}: run {run_id}: {WATCHER_STARTS_PER_RUN} watcher starts in this driver process "
                                   f"already: no further start (decision 1830 (3)); check the watchers by hand")
            self.watcher_starts[run_id] = self.watcher_starts.get(run_id, 0) + 1
            pid = self.watchers.start(argv, directory)
            if pid:
                break
        watcher = {"Dir": directory, "Stage": name, "PID": pid, "Started": attempt, "UTC": self.now().isoformat(),
                   "Argv": argv, "Rule": WATCHER_RULE}
        if pid:
            self.log(f"{self.case}: {name} run {run_id}: run watcher pid {pid} (cap {float(cap_seconds):.0f} s)")
            return watcher
        try:
            paladin_call(self.block["CLI"], ["palace", "cancel", str(run_id), "--json"])
            watcher["EarlyEnded"] = True
        except PaladinError as error:
            watcher["EarlyEnded"], watcher["Error"] = False, str(error)[:300]
        raise WatcherNotStarted(f"{name}: the run watcher of run {run_id} did not start (twice): the run is ended early "
                                f"({'done' if watcher['EarlyEnded'] else 'FAILED: end it by hand'}): STOP ({WATCHER_RULE})",
                                watcher)

    def backstops_root(self):
        """BACKSTOP_RULE: the backstops' directory on /data, in the plan's directory."""
        if not self.remote_dir:
            raise PaladinError(f"{self.case}: no remote directory: no place for the backstops ({BACKSTOP_RULE})")
        return f"{self.remote_dir.rstrip('/')}/paladin-backstops"

    def backstop_check(self, name):
        """BACKSTOP_RULE before a submit: the installed script and its --check ({Mode none} when the profile has none), or
        PaladinError (STOP before the submit)."""
        if self.backstop["Mode"] == "none":
            self.log(f"{self.case}: {name}: NO backstop (the profile's Paladin.Backstop is Mode none)")
            return {"Mode": "none", "UTC": self.now().isoformat(), "Rule": BACKSTOP_RULE}
        try:
            check = self.backstops.install(self.host, self.backstops_root(), self.pbs_bin)
        except PaladinError as error:
            raise PaladinError(f"{name}: the PBS-side backstop cannot be started on {self.block['Host']}: STOP before the "
                               f"submit ({error}; {BACKSTOP_RULE})") from None
        return {**check, "Mode": "qdel", "UTC": self.now().isoformat()}

    def ensure_backstop(self, stage, record):
        """BACKSTOP_RULE for the stage's live run: started when none was, re-read and restarted when gone; a refusal, a
        failed restart or the start bound = PaladinError (the first start's failure ends the run early)."""
        name, run_id = stage["Name"], record["RunID"]
        if self.backstop["Mode"] == "none":
            record["Backstop"] = {"Mode": "none", "Rule": BACKSTOP_RULE}
            return
        directory = f"{self.backstops_root()}/run-{run_id}"
        starts = record.setdefault("Backstops", [])
        first = not starts
        if not first:
            status = self.backstops.status(self.host, directory)
            record["BackstopStatus"] = {"UTC": self.now().isoformat(), "Read": status}
            if status is None or status["Host"] != starts[-1].get("Host"):
                self.log(f"{self.case}: {name} run {run_id}: the backstop's status read is missing "
                         f"({'no answer' if status is None else 'another host ' + status['Host']}): re-read at the next poll")
                self.save()
                return
            end = status["End"] or {}
            if status["Alive"] or end.get("Exit") in (0, 11):
                self.save()
                return
            if end.get("Exit") == 13:
                self.save()
                raise PaladinError(f"{name}: run {run_id}: the backstop REFUSED ({end.get('Reason')}): STOP; the run stays "
                                   f"with its run watcher ({BACKSTOP_RULE})")
            self.log(f"{self.case}: {name} run {run_id}: the backstop is gone ({end or 'no end record'}): restarted")
        if self.backstop_starts.get(run_id, 0) >= BACKSTOP_STARTS_PER_RUN:
            raise PaladinError(f"{name}: run {run_id}: {BACKSTOP_STARTS_PER_RUN} backstop starts in this driver process "
                               f"already: no further start; check the backstop by hand ({BACKSTOP_RULE})")
        self.backstop_starts[run_id] = self.backstop_starts.get(run_id, 0) + 1
        entry = {"Dir": directory, "PID": None, "UTC": self.now().isoformat(), "Rule": BACKSTOP_RULE}
        try:
            job = record.get("JobID") or (self.describe(run_id).get("block") or {}).get("job_id")
            if not job:
                raise PaladinError("no PBS job id known for the run")
            record["JobID"] = record.get("JobID") or job
            check = record.get("BackstopCheck") if (record.get("BackstopCheck") or {}).get("Mode") == "qdel" else None
            check = check or self.backstop_check(name)
            submitted = parse_time(record.get("Submitted") or (record.get("Adopted") or {}).get("UTC")) or self.now()
            cap = float(stage["CapSeconds"])
            state = {"RunID": run_id, "JobID": str(job), "JobName": f"{self.block['Lane']}-{name}", "Owner": check["User"],
                     "LimitSeconds": cap, "GraceSeconds": self.backstop["GraceSeconds"],
                     "IntervalSeconds": self.backstop["IntervalSeconds"],
                     "LifeUntilEpoch": submitted.timestamp() + BACKSTOP_LIFE_SECONDS + cap + self.backstop["GraceSeconds"],
                     "PBSBin": self.pbs_bin}
            entry.update(Host=check["Host"], Script=check["Script"], ScriptSHA256=check["SHA256"], State=state)
            entry["PID"] = self.backstops.start(self.host, directory, check["Script"], state)
        except PaladinError as error:
            entry["Error"] = str(error)[:300]
        starts.append(entry)
        self.save()
        if entry["PID"]:
            self.log(f"{self.case}: {name} run {run_id}: PBS-side backstop pid {entry['PID']} on {entry['Host']} (limit "
                     f"{entry['State']['LimitSeconds']:.0f} + {entry['State']['GraceSeconds']:.0f} s from the PBS start)")
            return
        if first:
            try:
                paladin_call(self.block["CLI"], ["palace", "cancel", str(run_id), "--json"])
                entry["EarlyEnded"] = True
            except PaladinError as error:
                entry["EarlyEnded"], entry["CancelError"] = False, str(error)[:300]
            self.save()
            raise BackstopNotStarted(f"{name}: the PBS-side backstop of run {run_id} did not start ({entry.get('Error')}): the "
                                     f"run is ended early ({'done' if entry['EarlyEnded'] else 'FAILED: end it by hand'}): "
                                     f"STOP ({BACKSTOP_RULE})")
        raise PaladinError(f"{name}: run {run_id}: the backstop's restart failed ({entry.get('Error')}): STOP; the run stays "
                           f"with its run watcher ({BACKSTOP_RULE})")

    def backstop_end(self, record):
        """BACKSTOP_RULE: the end record of the run's backstop when it fired a qdel, else None (also when unread)."""
        starts = record.get("Backstops") or []
        if not starts:
            return None
        status = self.backstops.status(self.host, starts[-1]["Dir"])
        end = (status or {}).get("End") or {}
        return end if end.get("Qdels") else None

    def nonce_dir(self, intent):
        return os.path.join(self.watchers_root(), f"nonce-{intent['Nonce']}")

    def start_nonce_watcher(self, stage, record, intent):
        """WATCHER_RULE / NONCE_RULE: the watcher of the submit window, started BEFORE the submit call (twice at most);
        not started = the intent closed 'not sent' and STOP with no call (no run exists)."""
        name = stage["Name"]
        directory = self.nonce_dir(intent)
        identity = stage_identity(self.block["Lane"], self.case, name, self.plan_sha256)
        argv = [sys.executable, WATCHER_PATH, "--run-id", "0", "--nonce", intent["Nonce"], "--identity",
                json.dumps(identity, sort_keys=True), "--cli", self.block["CLI"], "--host", self.block["Host"], "--pbs-bin",
                self.pbs_bin, "--cap-seconds", str(float(stage["CapSeconds"])), "--submitted", intent["UTC"], "--bar",
                str(WATCHER_MEMORY_BAR), "--interval", str(WATCHER_INTERVAL_SECONDS), "--grace", str(WATCHER_GRACE_SECONDS),
                "--dir", directory]
        key, pid, attempt = ("nonce", name), None, 0
        for attempt in range(1, 3):
            if self.watcher_starts.get(key, 0) >= WATCHER_STARTS_PER_RUN:
                break  # decision 1830 (3): no further start in this driver process
            self.watcher_starts[key] = self.watcher_starts.get(key, 0) + 1
            pid = self.watchers.start(argv, directory)
            if pid:
                break
        intent["Watcher"] = {"Dir": directory, "PID": pid, "Started": attempt, "UTC": self.now().isoformat()}
        if pid:
            self.save()
            self.log(f"{self.case}: {name}: the submit-window watcher pid {pid} (nonce {intent['Nonce'][:8]})")
            return
        intent["Outcome"] = "not sent: the submit-window watcher did not start"
        self.save()
        raise WatcherNotStarted(f"{name}: the submit-window watcher did not start (twice, or {WATCHER_STARTS_PER_RUN} starts of "
                                f"this stage in this driver process): STOP with NO submit call ({WATCHER_RULE})",
                                intent["Watcher"])

    def write_run_hint(self, intent, run_id):
        """NONCE_RULE: the run id handed to the intent's nonce watcher (its hint file)."""
        directory = self.nonce_dir(intent)
        os.makedirs(directory, exist_ok=True)
        with open(os.path.join(directory, "run-id.tmp"), "w") as f:
            f.write(f"{int(run_id)}\n")
        os.replace(os.path.join(directory, "run-id.tmp"), os.path.join(directory, "run-id"))

    def watchers_root(self):
        """WATCHER_RULE: one root per plan (beside the plan's driver lock), whatever the --out directory."""
        return os.path.join(os.path.dirname(os.path.realpath(self.plan_path)), "paladin-watchers")

    def unverified_ends(self):
        """The END-UNVERIFIED.json markers of this job's run watchers (WATCHER_RULE)."""
        root, found = self.watchers_root(), []
        for name in sorted(os.listdir(root)) if os.path.isdir(root) else []:
            path = os.path.join(root, name, UNVERIFIED_MARKER)
            if os.path.isfile(path):
                with open(path) as f:
                    found.append({**json.load(f), "Path": path})
        return found

    def watch(self, stage, record):
        """ensure_watcher for the stage's recorded run, kept in its record (the submit time: the record's,
        else its last SubmitIntent's, else the adoption time)."""
        intents = record.get("SubmitIntents") or []
        submitted = record.get("Submitted") or (intents[-1]["UTC"] if intents else None) or record["Adopted"]["UTC"]
        nonce_dirs = [self.nonce_dir(intent) for intent in intents if intent.get("RunID") == record["RunID"]]
        try:
            started = self.ensure_watcher(stage["Name"], record["RunID"], job_id=record.get("JobID"), submitted=submitted,
                                          cap_seconds=stage["CapSeconds"], nonce_dirs=nonce_dirs)
        except WatcherNotStarted as error:
            record.setdefault("Watchers", []).append(error.watcher)
            self.save()
            raise
        if started is not None:
            record.setdefault("Watchers", []).append(started)
        self.save()
        self.ensure_backstop(stage, record)  # BACKSTOP_RULE: the PBS-side line, beside the watcher

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
        # WATCHER_RULE: the kept run keeps a watcher at the stage cap; a newer run whose early end failed gets one at cap 0.
        watched = [(ordered[0], stage["CapSeconds"])] if run_class(ordered[0].get("status")) == "live" else []
        watched += [(row, 0) for row, entry in zip(ordered[1:], ended) if entry.get("EarlyEnded") is False]
        duplicate["Watchers"] = []
        for row, cap in watched:
            try:
                started = self.ensure_watcher(stage["Name"], int(row["run_id"]), job_id=row.get("job_id"),
                                              submitted=self.now().isoformat(), cap_seconds=cap)
            except WatcherNotStarted as error:
                started = error.watcher
            duplicate["Watchers"].append({"RunID": int(row["run_id"]), "CapSeconds": cap, "Watcher": started})
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
        recovery = ""
        if self.open_intent(record) is not None:
            recovery = (f"; RECOVERY: check by hand `{self.block['CLI']} palace my-jobs {' '.join(metadata_arguments(identity))} "
                        f"--json`: a run listed = resume (it is adopted); none = close the intent with `python3 "
                        f"{os.path.abspath(__file__)} {os.path.abspath(self.plan_path)} --profile <the profile> --out "
                        f"{os.path.abspath(self.local_dir)} --close-intent {stage['Name']}` (it records its own lookups), then "
                        f"resume ({CLOSE_INTENT_RULE})")
        raise PaladinError(f"{stage['Name']}: {why}, and {LOOKUP_REPEATS} lookups over "
                           f"{(LOOKUP_REPEATS - 1) * LOOKUP_REPEAT_SECONDS:.0f} s show no run of the stage {identity}: STOP, no "
                           f"second submit ({INTENT_RULE}){recovery}")

    def close_open_intent(self, name):
        """CLOSE_INTENT_RULE: the open intent of stage `name` closed after recorded lookups, or PaladinError."""
        self.lock()
        record = self.state["Stages"].get(name) or {}
        intent = self.open_intent(record)
        if intent is None:
            raise PaladinError(f"{name}: no open submit intent to close ({CLOSE_INTENT_RULE})")
        identity = stage_identity(self.block["Lane"], self.case, name, self.plan_sha256)
        lookups = []
        for attempt in range(LOOKUP_REPEATS):
            if attempt:
                self.sleep(LOOKUP_REPEAT_SECONDS)
            rows = find_runs(self.block["CLI"], identity)  # a failed lookup raises: not closed
            lookups.append({"UTC": self.now().isoformat(), "Rows": [[row.get("run_id"), row.get("status")] for row in rows]})
            runs = [row for row in rows if not no_run(row)]
            if runs:
                raise PaladinError(f"{name}: a run of the stage is listed {[row['run_id'] for row in runs]}: the intent is NOT "
                                   f"closed; resume (the run is adopted) ({CLOSE_INTENT_RULE})")
        intent["Outcome"] = "closed: no run after the recorded lookups"
        intent["ClosedBy"] = {"UTC": self.now().isoformat(), "Lookups": lookups, "PID": os.getpid(), "Rule": CLOSE_INTENT_RULE}
        self.save()
        self.log(f"{self.case}: {name}: the open intent {intent['Nonce'][:8]} closed after {len(lookups)} recorded lookups")
        return intent

    def adopt(self, stage, record, row, identity):
        name = record["Name"]
        if record.get("RunID") not in (None, int(row["run_id"])):
            raise PaladinError(f"{name}: the recorded run {record['RunID']} is not the stage's run {row['run_id']}")
        intent = self.open_intent(record)
        nonce = (row.get("metadata") or {}).get("nonce")
        if intent is not None:
            intent["Outcome"] = f"adopted run {row['run_id']}"
            intent["NonceMatches"] = nonce == intent["Nonce"]
        for known in record.get("SubmitIntents") or []:
            if known["Nonce"] == nonce:
                known["RunID"] = int(row["run_id"])
                self.write_run_hint(known, row["run_id"])  # NONCE_RULE: its nonce watcher learns the run
        if record.get("RunID") is None:
            record["RunID"] = int(row["run_id"])
            record["Adopted"] = {"UTC": self.now().isoformat(), "Status": row.get("status"),
                                 "InstanceType": row.get("instance_type"), "Nodes": row.get("nodes"),
                                 "Nonce": (row.get("metadata") or {}).get("nonce"), "Rule": RESUME_RULE}
            self.log(f"{self.case}: adopted {name} paladin run {row['run_id']} ({row.get('status')})")
        if row.get("job_id") and not record.get("JobID"):
            record["JobID"] = row["job_id"]
        record.setdefault("Identity", identity)
        self.save()
        if run_class(row.get("status")) == "live":
            self.watch(stage, record)  # a resume re-attaches to a live watcher or starts a missing one
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
                return self.adopt(stage, record, value, identity)
            if value:
                record["NoRunRowsBeforeSubmit"] = value
            if self.submit(stage, record) is not None:
                continue  # the lookup must now find exactly this run
            if self.open_intent(record) is not None:  # UNKNOWN_RECEIPT_RULE
                decision, value = self.repeated_lookup(stage, record, identity,
                                                       "a submit call without a usable answer (its receipt is unknown)")
                return self.adopt(stage, record, value, identity)
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
        if self.watcher_end(record) is None:
            self.watch(stage, record)  # decision 1836 m-4 (a): a dead watcher is restarted while the driver lives
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
        fired = self.backstop_end(record) if cls != "completed" else None
        if fired:
            record["BackstopEnd"], record["State"] = fired, "timed-out"
            self.save()
            raise PaladinError(f"{name}: run {record['RunID']} ended {status} by the PBS-side backstop's qdel at its limit "
                               f"(State timed-out): nothing is resubmitted ({BACKSTOP_RULE}; {RESUME_RULE})")
        if cls in ("cancelled", "out-of-memory", "unknown"):
            ended_by = self.watcher_end(record)
            record["WatcherEnd"] = ended_by
            record["State"] = ("timed-out" if cls == "cancelled" and (record.get("CapCancel", {}).get("OK") or
                                                                       (ended_by or {}).get("Kind") == "cap") else
                               "memory-stop" if cls == "cancelled" and (ended_by or {}).get("Kind") == "memory" else cls)
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
                               cancelled_by_backend=bool(record.get("CapCancel")), watcher_end=self.watcher_end(record),
                               watcher_events=self.watcher_events(record))
        record["CompletionRead"] = read
        if read["MissingOutputs"] or read["ExtraOutputs"]:
            self.log(f"{self.case}: {name} run {record['RunID']}: postpro names differ from OUTPUTS_RULE: {read['NameDifferences']}")
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

    def watcher_events(self, record, *, waits=3):
        """COMPLETION_RULE (5): the events.log text of the run's watcher (run-<id>, else its nonce directories), waited for
        its END line up to `waits` watcher intervals while a watcher lives; None when no events.log exists."""
        intents = record.get("SubmitIntents") or []
        directories = [os.path.join(self.watchers_root(), f"run-{record.get('RunID')}")]
        directories += [self.nonce_dir(intent) for intent in intents if intent.get("RunID") == record.get("RunID")]
        for attempt in range(waits + 1):
            texts = []
            for directory in directories:
                path = os.path.join(directory, "events.log")
                if os.path.isfile(path):
                    with open(path) as f:
                        texts.append(f.read())
            text = texts[0] if texts else None
            if text and " END exit " in text:
                return text
            if attempt == waits or not any(self.watchers.alive(directory) for directory in directories):
                return text
            self.sleep(WATCHER_INTERVAL_SECONDS)
        return None

    def watcher_end(self, record):
        """The end.json of the record's watcher (WATCHER_RULE), or None."""
        for watcher in reversed(record.get("Watchers") or []):
            path = os.path.join(watcher["Dir"], "end.json")
            if os.path.isfile(path):
                with open(path) as f:
                    return json.load(f)
        directory = os.path.join(self.watchers_root(), f"run-{record.get('RunID')}")
        if os.path.isfile(os.path.join(directory, "end.json")):
            with open(os.path.join(directory, "end.json")) as f:
                return json.load(f)
        return None

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
        unverified = self.unverified_ends()
        if unverified:
            self.state["UnverifiedEnds"] = unverified
            self.save()
            raise PaladinError(f"{self.case}: an early end NOT VERIFIED by the run watcher of "
                               f"{[(entry.get('RunID'), entry.get('Job')) for entry in unverified]}: STOP; check "
                               f"{unverified[0].get('Check')} (markers {[entry['Path'] for entry in unverified]}; {WATCHER_RULE})")
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
        for _ in range(DRIVE_MAX_ADVANCES):
            if job.advance():
                return job.state
            sleep(max(float(poll_seconds), DRIVE_MIN_POLL_SECONDS))
        raise PaladinError(f"{job.case}: {DRIVE_MAX_ADVANCES} advances without completion: STOP (decision 1830 (3)); the runs "
                           f"stay with their run watchers")
    finally:
        job.unlock()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("plan")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--dry-run", action="store_true")
    mode.add_argument("--submit", action="store_true")
    mode.add_argument("--close-intent", metavar="STAGE", help="close the open submit intent of STAGE (CLOSE_INTENT_RULE)")
    parser.add_argument("--profile", required=True, help="a cluster profile with a Paladin block")
    parser.add_argument("--out", required=True, help="the local directory of the job's state and run configurations")
    parser.add_argument("--remote-dir", help="the plan's directory on /data (outputs, logs, status.json); --submit")
    parser.add_argument("--map", action="append", default=[], help="REMOTE_PREFIX=LOCAL_PREFIX of the stage configurations")
    parser.add_argument("--archive-gb", action="append", default=[], help="DIRECTORY=GB: an archive writer's estimate")
    parser.add_argument("--poll-seconds", type=float, default=60.0)
    parser.add_argument("--pbs-bin", default=DEFAULT_PBS_BIN, help="the PBS bin directory on Host (the run watchers' qstat / qdel)")
    parser.add_argument("--lane-root", help="--submit: the lane's root on /data; every archive Directory must lie under it")
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
    if args.close_intent:
        job = PaladinJob(args.plan, block=block, remote_dir=args.remote_dir or "", local_dir=args.out, mapping=mapping)
        try:
            intent = job.close_open_intent(args.close_intent)
        except PaladinError as error:
            print(f"STOP: {error}", file=sys.stderr)
            return 3
        finally:
            job.unlock()
        print(json.dumps({"Stage": args.close_intent, "Nonce": intent["Nonce"], "Outcome": intent["Outcome"],
                          "Lookups": len(intent["ClosedBy"]["Lookups"])}))
        return 0
    if not args.remote_dir:
        parser.error("--submit needs --remote-dir")
    if not args.lane_root or not os.path.isabs(args.lane_root):
        parser.error("--submit needs an absolute --lane-root (decision 1820 m-3)")
    outside = [directory for directory in estimates if not under(directory, args.lane_root)]
    if outside:
        parser.error(f"--archive-gb {outside} outside --lane-root {args.lane_root} (decision 1820 m-3)")
    code = 0
    job = None
    try:
        job = PaladinJob(args.plan, block=block, remote_dir=args.remote_dir, local_dir=args.out, mapping=mapping,
                         archive_estimates_gb=estimates, pbs_bin=args.pbs_bin, lane_root=args.lane_root)
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
