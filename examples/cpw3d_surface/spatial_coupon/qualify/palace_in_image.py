#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Form X (USER decision 1674 (2)) for the lane job scripts that call Palace directly (the c0
window solves, load checks, preflights: `mpiexec_bound.sh -n N $BIN ARGS`): run Palace INSIDE
the image of a cluster profile's Launcher block on this node, as run_stages.py does for a plan.

Checks the image file's sha256 (ImageSHA256) and the sha256 of the executable inside it
(ImageBinarySHA256), refuses a job of more than one node (PBS_NODEFILE), drops the PBS_*
variables, passes into the image (apptainer --cleanenv) only the Launcher Environment, every
PALACE_RESPONSE_* variable of the caller and the thread-count / TMPDIR variables (forwarded
to the ranks with -x), prints the check and the command on stderr,
then replaces itself with `apptainer exec ... IMAGE <image mpirun> ... -n RANKS <image palace>
PALACE-ARGUMENTS`.  Exit 2 on a check failure (nothing started).

usage: palace_in_image.py PROFILE.json RANKS PALACE-ARGUMENT...
  e.g. palace_in_image.py cluster-profile-x.json 96 --surface-response-preflight c0-p4.json
       palace_in_image.py cluster-profile-x.json 192 c0-p4.json
"""
import json
import os
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import build_plan  # noqa: E402
import run_stages  # noqa: E402


def image_launch(profile, ranks, arguments, environment, *, digest=None, run=None):
    """(command, environment, check record) of one Palace launch inside the profile's image."""
    launcher = build_plan.launcher_of(profile)
    if launcher is None:
        raise SystemExit("the cluster profile has no Launcher block (form X)")
    if environment.get("PBS_NODEFILE"):
        hosts = run_stages.unique_hosts(Path(environment["PBS_NODEFILE"]).read_text())
        if len(hosts) > 1:
            raise SystemExit(f"{len(hosts)} nodes in PBS_NODEFILE: {build_plan.LAUNCHER_SINGLE_NODE_RULE}")
    check = run_stages.launcher_preflight(launcher, {"Binary": launcher["ImageBinary"], "MPIExec": launcher["ImageMPIExec"]},
                                          digest=digest, run=run)
    if check["BinarySHA256"] != launcher["ImageBinarySHA256"]:
        raise SystemExit(f"Executable hash mismatch inside the image: {check['BinarySHA256']} != "
                         f"{launcher['ImageBinarySHA256']}")
    stage_environment = run_stages.launcher_environment(environment, launcher)
    exports = [part for key in sorted(stage_environment) if key.startswith("PALACE_RESPONSE_") for part in ("-x", key)]
    command = run_stages.launcher_launch(launcher, exports=exports, ranks=int(ranks), arguments=list(arguments))
    return command, stage_environment, check


def main(argv):
    if len(argv) < 4:
        print(__doc__, file=sys.stderr)
        return 64
    profile = json.loads(Path(argv[1]).read_text())
    try:
        command, environment, check = image_launch(profile, argv[2], argv[3:], dict(os.environ))
    except SystemExit as error:
        print(f"palace_in_image: STOP: {error}", file=sys.stderr, flush=True)
        return 2
    print(f"palace_in_image: image {check['Image']} sha256 {check['ImageSHA256']}; executable sha256 {check['BinarySHA256']}",
          file=sys.stderr)
    print("palace_in_image: " + " ".join(command), file=sys.stderr, flush=True)
    os.execve(command[0], command, environment)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
