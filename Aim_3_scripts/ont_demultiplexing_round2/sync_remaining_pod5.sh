#!/usr/bin/env bash
set -euo pipefail

SMB_URL="smb://vault-v2.erc.monash.edu/mnhs-ccs-cgen-lab"
SRC="/run/user/$(id -u)/gvfs/smb-share:server=vault-v2.erc.monash.edu,share=mnhs-ccs-cgen-lab/DNAScreen/ONT_data/20260819_1420_P2S-03613-B_PBM47885/pod5"
DEST="/fs04/vh83/ont_dnascreen/round2_pod5_last_40"
FILELIST="/tmp/pod5_remaining_list.txt"
MIN_INDEX=87

# ── Ensure the SMB share is mounted ──────────────────────────────────────────
if ! gio mount -l | grep -q "mnhs-ccs-cgen-lab"; then
  echo "Share not mounted — mounting now (you'll be prompted for credentials)..."
  gio mount "$SMB_URL"
else
  echo "Share already mounted."
fi

# Give gvfs a moment to expose the mount path
sleep 2

if [ ! -d "$SRC" ]; then
  echo "ERROR: Expected source path not found: $SRC"
  echo "Check 'gio mount -l' output and adjust the SRC path if the mount name differs."
  exit 1
fi

# ── Build list of pod5 files with index >= MIN_INDEX ─────────────────────────
mkdir -p "$DEST"
> "$FILELIST"

for f in "$SRC"/PBM47885_6a58f094_4dec16a4_*.pod5; do
  base=$(basename "$f")
  num=$(echo "$base" | sed -E 's/PBM47885_6a58f094_4dec16a4_([0-9]+)\.pod5/\1/')
  if [[ "$num" =~ ^[0-9]+$ ]] && (( num >= MIN_INDEX )); then
    echo "$base" >> "$FILELIST"
  fi
done

echo "Files matched: $(wc -l < "$FILELIST")"
echo "Highest indices found:"
sort -t_ -k4 -n "$FILELIST" | tail -5

# ── Dry run then real transfer ────────────────────────────────────────────────
echo "--- Dry run ---"
rsync -avh --dry-run --progress --files-from="$FILELIST" "$SRC/" "$DEST/"

read -rp "Proceed with actual transfer? [y/N] " confirm
if [[ "$confirm" =~ ^[Yy]$ ]]; then
  rsync -avh --progress --partial --files-from="$FILELIST" "$SRC/" "$DEST/"
  echo "Transfer complete."
else
  echo "Aborted — no files copied."
fi