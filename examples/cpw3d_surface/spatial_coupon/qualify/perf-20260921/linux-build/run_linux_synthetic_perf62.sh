#!/bin/bash
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
# The archive-estimate synthetic tests of the recorded procedure against the decision-62(4)
# executable, with the RECORDED frozen tools (source-freeze/tools: the standalone bounded
# runner and the test the b28 executable was validated with); a compute-node job (the login
# node has no /opt/openmpi).
#PBS -N coupon-perf62-synthetic
#PBS -q normal-g
#PBS -P DS-EM-FEM
#PBS -r n
#PBS -l select=1:ncpus=192:mpiprocs=192
#PBS -l place=scatter:excl
#PBS -l instance_type=r8g.48xlarge
#PBS -l efa_support=True,subnet_id=subnet-0c98d793bbcebb39a
#PBS -l walltime=00:30:00
#PBS -j oe
#PBS -o /data/home/simlap/coupon_accuracy_assessment_20260913/linux-synthetic-perf62.pbs.log
set -euo pipefail
cd /data/home/simlap/coupon_accuracy_assessment_20260913
trap 'code=$?; printf "{\"ExitCode\":%d,\"CompletedUTC\":\"%s\",\"JobID\":\"%s\"}\n" "$code" "$(date -u +%FT%TZ)" "${PBS_JOBID:-}" > linux-synthetic-perf62-status.json' EXIT
source /etc/profile.d/modules.sh
module load gcc/14.3.0 openmpi/5.0.8-gcc14 arm/armpl/24.10-gcc
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 ARMPL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
date -u; hostname
rm -rf linux-synthetic-perf62-recorded-tools; mkdir linux-synthetic-perf62-recorded-tools
export ARCHIVE_DIAGNOSTIC_SCRATCH="$PWD/linux-synthetic-perf62-recorded-tools" TMPDIR="$PWD/linux-synthetic-perf62-recorded-tools"
export ARCHIVE_DIAGNOSTIC_EXE="$(python3 -c 'import json; print(json.load(open("linux-build-perf62/binary.json"))["Path"])')"
sha256sum "$ARCHIVE_DIAGNOSTIC_EXE"
python3 source-freeze/tools/run_bounded_mesher.py --seconds 900 --memory-gib 6 --log linux-synthetic-perf62-recorded-tools/tests.log -- python3 source-freeze/tools/test_estimate_archived_fields.py -v
