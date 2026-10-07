# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Remote side of `coupon-library qualify` (the launch.sh / monitor.sh / fetch_results.sh
of the physics runs and the 40-job accounting of submit_graded_library_campaign.py):
upload of a coupon's inputs, submission under the user job cap, read-only monitoring,
fetch of the results (never the field archives), remote digest verification and the
recorded deletion of the archives once the CSVs are verified locally.

Every function builds its ssh / rsync command from the host and remote root arguments
and the cluster profile; nothing here is executed under --dry-run (the commands are
recorded instead).
"""
import getpass
import json
import re
import subprocess
import time

JOB_ID_PATTERN = r"\d+\.[A-Za-z0-9._-]+"


def utc():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def ssh(host, command, *, check=True):
    return subprocess.run(["ssh", host, command], text=True, capture_output=True, check=check)


def upload_command(host, local_directory, remote_directory):
    """rsync of a local directory tree into the remote directory (which must exist: the
    macOS openrsync client has no --mkpath, so the directory is created by ssh first)."""
    return ["rsync", "-a", f"{local_directory}/", f"{host}:{remote_directory}/"]


def upload_file_command(host, local_file, remote_file):
    return ["rsync", "-a", str(local_file), f"{host}:{remote_file}"]


def upload(host, local_directory, remote_directory):
    ssh(host, f"mkdir -p '{remote_directory}'")
    command = upload_command(host, local_directory, remote_directory)
    subprocess.run(command, check=True)
    return command


def upload_file(host, local_file, remote_file):
    ssh(host, f"mkdir -p '{str(remote_file).rsplit('/', 1)[0]}'")
    command = upload_file_command(host, local_file, remote_file)
    subprocess.run(command, check=True)
    return command


def active_user_jobs(host, pbs_bin, user=None):
    """The user's queued / running PBS jobs (job id -> state) via a read-only qstat."""
    result = ssh(host, f"{pbs_bin}/qstat -f -F json 2>/dev/null || echo '{{}}'")
    data = json.loads(result.stdout or "{}")
    owner = (user or getpass.getuser()) + "@"
    return {key: job.get("job_state") for key, job in data.get("Jobs", {}).items()
            if job.get("Job_Owner", "").startswith(owner) and job.get("job_state") != "F"}


def submit(host, pbs_bin, remote_job_script, remote_working_directory, *, job_cap, user=None):
    """qsub under the user job cap; returns the submission record (fail closed when the
    cap would be exceeded or the receipt is ambiguous)."""
    current = active_user_jobs(host, pbs_bin, user)
    if len(current) + 1 > job_cap:
        raise RuntimeError(f"would exceed the {job_cap} user job cap: {len(current)} active/queued + 1 new")
    result = ssh(host, f"cd {remote_working_directory} && {pbs_bin}/qsub {remote_job_script}")
    job = result.stdout.strip()
    if not re.fullmatch(JOB_ID_PATTERN, job):
        raise RuntimeError(f"ambiguous qsub receipt {job!r}; reconcile before retry")
    return {"Job": job, "UTC": utc(), "UserJobsBefore": len(current), "JobCap": job_cap,
            "Command": f"cd {remote_working_directory} && {pbs_bin}/qsub {remote_job_script}"}


POLL_MARKER = "---POLL-OK---"


# PBS job states of a job still in the queue: queued, running, exiting, held, waiting,
# transiting, suspended, array begun.  H is not "left the queue": the SOCA dispatcher holds a
# job whose compute-node stack failed (Resource_List.error_message CF:ROLLBACK_COMPLETE:retry=N)
# and releases it itself at retry_eligible_after (the 2026-09-21 split acceptance reducer job).
IN_QUEUE_STATES = ("Q", "R", "E", "H", "W", "T", "S", "B")


def in_queue(job_state):
    return job_state in IN_QUEUE_STATES


def poll(host, pbs_bin, job_id, remote_status_path):
    """One read-only poll: qstat state fields and the runner's status.json (if written).
    `Reachable` is False when the ssh round trip itself failed (no POLL_MARKER came back):
    a lost connection says nothing about the job and must not be read as "left the queue"."""
    command = (f"{pbs_bin}/qstat -f {job_id} 2>/dev/null | grep -E 'job_state|resources_used.walltime|resources_used.mem|exec_host|comment|Hold_Types|error_message' "
               f"| tr -s ' '; echo ---STATUS---; cat {remote_status_path} 2>/dev/null; echo; echo {POLL_MARKER}")
    result = ssh(host, command, check=False)
    reachable = POLL_MARKER in result.stdout
    output = result.stdout.replace(POLL_MARKER, "")
    qstat_text, _, status_text = output.partition("---STATUS---\n")
    state = re.search(r"job_state = (\w)", qstat_text)
    status = None
    if status_text.strip():
        try:
            status = json.loads(status_text)
        except json.JSONDecodeError:
            status = None
    return {"UTC": utc(), "JobState": state.group(1) if state else None, "QStat": qstat_text.strip(),
            "Status": status, "Reachable": reachable, "SSHReturnCode": result.returncode}


def monitor(host, pbs_bin, job_id, remote_status_path, *, interval_seconds, max_polls, sink=print):
    """Poll until the job leaves the queue (IN_QUEUE_STATES) or the poll budget is spent; returns the polls."""
    polls = []
    for _ in range(max_polls):
        record = poll(host, pbs_bin, job_id, remote_status_path)
        polls.append(record)
        stages = ([(s["Name"], s["State"], round(s.get("WallSeconds", 0))) for s in record["Status"]["Stages"]]
                  if record["Status"] else None)
        sink(f"== {record['UTC']} job {job_id} state {record['JobState']} stages {stages}"
             + ("" if record["Reachable"] else f" (unreachable: ssh rc {record['SSHReturnCode']})"))
        if record["Reachable"] and not in_queue(record["JobState"]):
            break
        time.sleep(interval_seconds)
    return polls


def fetch_command(host, remote_directory, local_directory):
    return ["rsync", "-a", "--exclude", "archive/", "--exclude", "tmp/", f"{host}:{remote_directory}/", f"{local_directory}/"]


def fetch(host, remote_directory, local_directory):
    command = fetch_command(host, remote_directory, local_directory)
    subprocess.run(command, check=True)
    return command


def remote_sha256(host, remote_paths):
    """Remote sha256sum of the given files: path -> digest."""
    if not remote_paths:
        return {}
    quoted = " ".join(f"'{path}'" for path in remote_paths)
    result = ssh(host, f"sha256sum {quoted}")
    digests = {}
    for line in result.stdout.splitlines():
        digest, _, path = line.partition("  ")
        digests[path.strip()] = digest.strip()
    return digests


READ_MARKER = "---READ-OK---"


def read_json_reachable(host, remote_path):
    """(reachable, content) of a remote JSON file: `reachable` is whether the ssh round trip
    itself answered (READ_MARKER came back); `content` is the parsed file, None when absent /
    unparsable / unreachable.  The caller tells a transport failure (not reachable) from an
    absent file (decision 495 (3))."""
    result = ssh(host, f"cat '{remote_path}' 2>/dev/null; echo {READ_MARKER}", check=False)
    if READ_MARKER not in result.stdout:
        return False, None
    text = result.stdout.replace(READ_MARKER, "").strip()
    if not text:
        return True, None
    try:
        return True, json.loads(text)
    except json.JSONDecodeError:
        return True, None


def read_json(host, remote_path):
    """The parsed content of a remote JSON file (None when absent, unparsable or unreachable)."""
    return read_json_reachable(host, remote_path)[1]


EXISTS_MARKER = "---EXISTS---"
ABSENT_MARKER = "---ABSENT---"


def path_exists(host, remote_path):
    """Whether `remote_path` exists on the host (a read-only test).  A round trip that returns
    neither marker (a lost connection) raises: the caller must fail closed, never read
    "absent" from silence (decision 485 (b))."""
    result = ssh(host, f"if [ -e '{remote_path}' ]; then echo {EXISTS_MARKER}; else echo {ABSENT_MARKER}; fi", check=False)
    if EXISTS_MARKER in result.stdout:
        return True
    if ABSENT_MARKER in result.stdout:
        return False
    raise RuntimeError(f"the existence test of {remote_path} on {host} returned no answer (ssh rc {result.returncode}: "
                       f"{result.stderr.strip()[:200]!r})")


JOB_EXIT_RECORD = "pbs-status.json"
JOB_EXIT_RULE = ("decision 485 (b): before any fetch the job's own exit record <job dir>/pbs-status.json (the job script's EXIT "
                 "trap: ExitCode, JobID = PBS_JOBID) must exist with ExitCode 0 and the submitted job id, and the runner's "
                 "status.json must carry that job id as PBSJobID - a job that died before the runner (b-batch1 O1 attempt 1: "
                 "mkdir tmp 'File exists', Exit_status 1 in 1 s) left an older run's status.json in the shared case directory "
                 "and the driver fetched it as the job's results")


def job_exit(host, remote_job_directory, job_id):
    """The exit check of one finished job (JOB_EXIT_RULE): {OK, Reachable, Reasons, ExitRecord,
    StatusPBSJobID}; OK only when pbs-status.json reads ExitCode 0 for `job_id` and status.json
    names `job_id`.  A round trip that did not answer is Reachable False (not OK, with a
    transport reason): a lost connection says nothing about the job - the caller retries
    (--resume) instead of blaming the job (decision 495 (3))."""
    reachable, exit_record = read_json_reachable(host, f"{remote_job_directory}/{JOB_EXIT_RECORD}")
    if reachable:
        reachable, status = read_json_reachable(host, f"{remote_job_directory}/status.json")
    else:
        status = None
    if not reachable:
        return {"OK": False, "Reachable": False,
                "Reasons": [f"the job directory {remote_job_directory} could not be read on {host} (the ssh round trip did not answer): "
                            "a transport failure, not a verdict on the job - the job's exit stays unchecked; run again with --resume"],
                "ExitRecord": None, "StatusPBSJobID": None, "Rule": JOB_EXIT_RULE}
    reasons = []
    if exit_record is None:
        reasons.append(f"no {JOB_EXIT_RECORD} in {remote_job_directory} (the job script's EXIT trap did not run)")
    else:
        if exit_record.get("JobID") != job_id:
            reasons.append(f"{JOB_EXIT_RECORD} names job {exit_record.get('JobID')!r}, the submission {job_id!r}")
        if exit_record.get("ExitCode") != 0:
            reasons.append(f"{JOB_EXIT_RECORD} reads ExitCode {exit_record.get('ExitCode')!r}")
    if status is None:
        reasons.append(f"no status.json in {remote_job_directory}")
    elif status.get("PBSJobID") != job_id:
        reasons.append(f"status.json was written by job {status.get('PBSJobID')!r} (started {status.get('StartUTC')}), not by "
                       f"the submitted {job_id!r}: a stale run in a pre-existing case directory")
    return {"OK": not reasons, "Reachable": True, "Reasons": reasons, "ExitRecord": exit_record,
            "StatusPBSJobID": status.get("PBSJobID") if status else None, "Rule": JOB_EXIT_RULE}


def count_archive_potentials(host, archive_directory):
    """The number of archived potential files (source-*-rank-*-V.bin) in a response
    archive directory: the reducer's union check of a split coupon (decision 61b)."""
    result = ssh(host, f"ls '{archive_directory}' 2>/dev/null | grep -c -- '-V.bin$' || true", check=False)
    text = result.stdout.strip().splitlines()
    return int(text[-1]) if text and text[-1].isdigit() else 0


def qstat_history(host, pbs_bin, job_id):
    result = ssh(host, f"{pbs_bin}/qstat -xf {job_id}", check=False)
    return result.stdout


def delete_archives(host, archive_directories):
    """du then rm -rf of the response archives (after the CSVs are verified locally);
    returns the recorded sizes and the deletion time."""
    quoted = " ".join(f"'{path}'" for path in archive_directories)
    sizes = ssh(host, f"du -sh {quoted} 2>/dev/null || true", check=False).stdout
    ssh(host, f"rm -rf {quoted}")
    remaining = ssh(host, f"ls -d {quoted} 2>/dev/null || true", check=False).stdout.strip()
    return {"Archives": list(archive_directories), "SizesBeforeDeletion": sizes.strip(), "DeletedUTC": utc(),
            "Remaining": remaining}
