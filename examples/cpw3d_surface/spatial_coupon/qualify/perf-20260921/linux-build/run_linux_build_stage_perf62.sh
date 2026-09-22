#!/bin/bash
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
# The recorded run_linux_build_stage.sh of coupon-accuracy-assessment-20260913/experiments
# for the decision-62(4) executable: bounded six-job source-only build of
# source-freeze-perf62 against the read-only Linux dependencies, then the archive
# estimate synthetic tests with the new executable. Submitted as one PBS job (the login
# node has 2 cores / 7 GB).
#PBS -N coupon-perf62-build
#PBS -q normal-g
#PBS -P DS-EM-FEM
#PBS -r n
#PBS -l select=1:ncpus=192:mpiprocs=192
#PBS -l place=scatter:excl
#PBS -l instance_type=r8g.48xlarge
#PBS -l efa_support=True,subnet_id=subnet-0c98d793bbcebb39a
#PBS -l walltime=01:00:00
#PBS -j oe
#PBS -o /data/home/simlap/coupon_accuracy_assessment_20260913/linux-build-perf62.pbs.log
set -euo pipefail
cd /data/home/simlap/coupon_accuracy_assessment_20260913
trap 'code=$?; printf "{\"ExitCode\":%d,\"CompletedUTC\":\"%s\",\"JobID\":\"%s\"}\n" "$code" "$(date -u +%FT%TZ)" "${PBS_JOBID:-}" > linux-build-perf62-status.json' EXIT
source /etc/profile.d/modules.sh
module load gcc/14.3.0 openmpi/5.0.8-gcc14 arm/armpl/24.10-gcc
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 ARMPL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
module list
date -u; hostname
free -b
python3 source-freeze-perf62/tools/run_bounded_mesher.py --seconds 1800 --memory-gib 32 --log linux-build-perf62.log -- python3 build_linux_perf62.py
mkdir linux-synthetic-perf62
export ARCHIVE_DIAGNOSTIC_SCRATCH="$PWD/linux-synthetic-perf62" TMPDIR="$PWD/linux-synthetic-perf62"
export ARCHIVE_DIAGNOSTIC_EXE="$(python3 -c 'import json; print(json.load(open("linux-build-perf62/binary.json"))["Path"])')"
python3 source-freeze-perf62/tools/run_bounded_mesher.py --seconds 600 --memory-gib 6 --log linux-synthetic-perf62/tests.log -- python3 source-freeze-perf62/tools/test_estimate_archived_fields.py -v
