#!/usr/bin/env bash
set -euo pipefail

LOCAL_ROOT="$(git rev-parse --show-toplevel)"
REMOTE_HOST="ubuntu@3.19.26.202"
REMOTE_ROOT="/home/ubuntu/palace"
SSH_KEY="${SSH_KEY:-$HOME/.ssh/ec2-gpu.pem}"
SSH_OPTS=(-J agent -i "$SSH_KEY" -o BatchMode=yes -o StrictHostKeyChecking=accept-new)
RSYNC_RSH="ssh -J agent -i $SSH_KEY -o BatchMode=yes -o StrictHostKeyChecking=accept-new"
SHA="$(git -C "$LOCAL_ROOT" rev-parse --short HEAD)"
RUN_ID="${SHA}-$(date +%Y%m%d%H%M%S)"
CASE_SET="${PALACE_AR_CASES:-cpw}"

cd "$LOCAL_ROOT"

if [[ "$LOCAL_ROOT" == "$REMOTE_ROOT" ]]; then
  exec .auto/remote_measure_payload.sh "$RUN_ID" "$CASE_SET"
fi

# Mirror source edits to the remote GPU worktree while preserving the remote's Linux
# Spack environment files and git metadata. The local repo remains the source of truth for
# local autoresearch commits; the remote is a benchmark/build mirror.
rsync -az --delete \
  --exclude '.git/' \
  --exclude '.auto/log.jsonl' \
  --exclude '.spack-env/' \
  --exclude 'spack.yaml' \
  --exclude 'spack.lock' \
  --exclude 'spack_repo/local/packages/mfem' \
  --exclude 'build/' \
  --exclude 'build-*' \
  --exclude 'install/' \
  --exclude 'postpro/' \
  --exclude '*.nsys-rep' \
  --exclude '*.sqlite' \
  --exclude '*.log' \
  --exclude '*.csv' \
  -e "$RSYNC_RSH" \
  "$LOCAL_ROOT/" "$REMOTE_HOST:$REMOTE_ROOT/"

ssh "${SSH_OPTS[@]}" "$REMOTE_HOST" \
  "cd '$REMOTE_ROOT' && .auto/remote_measure_payload.sh '$RUN_ID' '$CASE_SET'"
