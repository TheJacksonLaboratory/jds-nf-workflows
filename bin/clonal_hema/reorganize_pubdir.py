#!/usr/bin/env python3
"""
Reorganize clonal_hema's published output directory into one folder per
sample. Two SHARED modules used unmodified by this workflow publish some
files outside each sample's own folder, in ways specific to other
workflows' naming conventions:

  - samtools_stats.nf (used for clonal_hema's raw/consensus samtools stats)
    splits sampleID on "_" and buckets output under
    {pubdir}/{first_token}/samtools/ -- built for chipseq's own
    INPUT_/treatment naming, not clonal_hema's (e.g. lands at
    {pubdir}/CHIP/samtools/ for a sampleID like CHIP_REALISTIC_DUPLEX).
  - picard_collecthsmetrics.nf's raw (pre-consensus) HsMetrics call is
    invoked with sampleID suffixed "_raw" to avoid a filename collision
    with the consensus HsMetrics file, landing under a separate top-level
    {pubdir}/{sampleID}_raw/ folder.

Rather than modifying those shared modules (used by other workflows too),
this script moves both classes of stray output into
{pubdir}/{sampleID}/raw/ or {pubdir}/{sampleID}/stats/ as appropriate,
then removes the emptied stray directories. Intended to run once, as the
last step of the workflow, after all publishing has completed.

Usage: reorganize_pubdir.py <pubdir> <sampleID> [<sampleID> ...]
"""

import os
import shutil
import sys


def move_into(src, dest_dir):
    os.makedirs(dest_dir, exist_ok=True)
    dest = os.path.join(dest_dir, os.path.basename(src))
    if os.path.abspath(src) == os.path.abspath(dest):
        return dest
    shutil.move(src, dest)
    return dest


def fix_split_bucket(pubdir, sample_ids):
    """samtools_stats.nf's stray {first_token}/samtools/ bucket."""
    buckets = {}
    for sid in sample_ids:
        buckets.setdefault(sid.split("_")[0], []).append(sid)

    for base, sids in buckets.items():
        bucket_dir = os.path.join(pubdir, base, "samtools")
        if not os.path.isdir(bucket_dir):
            continue

        for fname in os.listdir(bucket_dir):
            src = os.path.join(bucket_dir, fname)
            matched = max((s for s in sids if fname.startswith(s)), key=len, default=None)
            if matched is None:
                continue
            rest = fname[len(matched):]
            if rest.startswith("_raw"):
                dest_dir = os.path.join(pubdir, matched, "raw")
                new_name = matched + rest.replace("_raw", "", 1)
            else:
                dest_dir = os.path.join(pubdir, matched, "stats")
                new_name = fname

            moved_path = move_into(src, dest_dir)
            target_path = os.path.join(dest_dir, new_name)
            if moved_path and moved_path != target_path and os.path.exists(moved_path):
                os.replace(moved_path, target_path)

        if not os.listdir(bucket_dir):
            os.rmdir(bucket_dir)
        base_dir = os.path.join(pubdir, base)
        if os.path.isdir(base_dir) and not os.listdir(base_dir):
            os.rmdir(base_dir)


def fix_raw_suffixed_dir(pubdir, sample_ids):
    """picard_collecthsmetrics.nf's raw call, suffixed {sampleID}_raw/."""
    for sid in sample_ids:
        raw_dir = os.path.join(pubdir, f"{sid}_raw")
        if not os.path.isdir(raw_dir):
            continue
        dest_dir = os.path.join(pubdir, sid, "raw")
        for root, _dirs, files in os.walk(raw_dir):
            for fname in files:
                moved = move_into(os.path.join(root, fname), dest_dir)
                # <sid>_raw.CoverageMetrics.txt -> <sid>.raw.CoverageMetrics.txt (sample prefix ends with a dot, 2026-11-16)
                if fname.startswith(f"{sid}_raw."):
                    os.replace(moved, os.path.join(dest_dir, f"{sid}.raw." + fname[len(f"{sid}_raw."):]))
        shutil.rmtree(raw_dir)


def main():
    if len(sys.argv) < 3:
        sys.exit(f"Usage: {sys.argv[0]} <pubdir> <sampleID> [<sampleID> ...]")
    pubdir = sys.argv[1]
    sample_ids = sys.argv[2:]
    fix_split_bucket(pubdir, sample_ids)
    fix_raw_suffixed_dir(pubdir, sample_ids)
    print(f"reorganize_pubdir: done for {len(sample_ids)} sample(s)")


if __name__ == "__main__":
    main()
