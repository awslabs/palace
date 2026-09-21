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


def poll(host, pbs_bin, job_id, remote_status_path):
    """One read-only poll: qstat state fields and the runner's status.json (if written).
    `Reachable` is False when the ssh round trip itself failed (no POLL_MARKER came back):
    a lost connection says nothing about the job and must not be read as "left the queue"."""
    command = (f"{pbs_bin}/qstat -f {job_id} 2>/dev/null | grep -E 'job_state|resources_used.walltime|resources_used.mem|exec_host|comment' "
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
    """Poll until the job leaves Q / R / E or the poll budget is spent; returns the polls."""
    polls = []
    for _ in range(max_polls):
        record = poll(host, pbs_bin, job_id, remote_status_path)
        polls.append(record)
        stages = ([(s["Name"], s["State"], round(s.get("WallSeconds", 0))) for s in record["Status"]["Stages"]]
                  if record["Status"] else None)
        sink(f"== {record['UTC']} job {job_id} state {record['JobState']} stages {stages}"
             + ("" if record["Reachable"] else f" (unreachable: ssh rc {record['SSHReturnCode']})"))
        if record["Reachable"] and record["JobState"] not in ("Q", "R", "E"):
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
