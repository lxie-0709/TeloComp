#!/usr/bin/env python3
import os
import argparse
import logging
from math import ceil
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio import SeqIO
from Bio.Seq import Seq

def get_case_lengths(records):
    upper, lower = [], []
    for rec in records:
        s = str(rec.seq)
        upper.append(sum(1 for c in s if c.isupper()))
        lower.append(sum(1 for c in s if c.islower()))
    return sorted(upper, reverse=True), sorted(lower, reverse=True)

def get_cut_length(lengths, coverage):
    n = len(lengths)
    if n == 0:
        return 0
    k = ceil(n * coverage)
    k = max(1, min(k, n))
    return lengths[k-1]

def trim_sequence(s, cut_up, cut_lo):
    up = ''.join(c for c in s if c.isupper())[:cut_up]
    lo = ''.join(c for c in s if c.islower())[:cut_lo]
    return up + lo

def parse_direction(fn):
    fn = os.path.basename(fn)
    if "_L" in fn:
        return "L"
    if "_R" in fn:
        return "R"
    return None

def parse_chr(fn):
    return os.path.basename(fn).split("_")[0]

def parse_platform(path, ont_dir, hifi_dir):
    if os.path.commonpath([path, ont_dir]) == ont_dir:
        return "ONT"
    if os.path.commonpath([path, hifi_dir]) == hifi_dir:
        return "HiFi"
    return None

def process_fasta(path, ont_dir, hifi_dir, out_base, coverage):
    log = logging.getLogger(__name__)
    try:
        records = list(SeqIO.parse(path, "fasta"))
        if not records:
            log.warning(f"No records in {path}")
            return

        # figure out original uppercase/lowercase lengths
        up_lens, lo_lens = get_case_lengths(records)
        cut_up = get_cut_length(up_lens, coverage)
        cut_lo = get_cut_length(lo_lens, coverage)

        # find the record with the minimum original lowercase count
        orig_lower_counts = [sum(1 for c in str(rec.seq) if c.islower()) for rec in records]
        min_index = orig_lower_counts.index(min(orig_lower_counts))  # first minimal
        min_rec = records[min_index]

        dirn = parse_direction(path)
        plat = parse_platform(path, ont_dir, hifi_dir)
        chrn = parse_chr(path)
        if not dirn or not plat:
            log.warning(f"Skipping {path}, cannot infer direction/platform")
            return

        # prepare output directory
        out_dir = os.path.join(out_base, f"trim_{dirn}")
        os.makedirs(out_dir, exist_ok=True)
        out_fn = f"{chrn}_trimmed_{plat}_{dirn}.fasta"
        out_path = os.path.join(out_dir, out_fn)

        # now trim and rewrite IDs
        trimmed = []
        for rec in records:
            s = str(rec.seq)
            new_seq = trim_sequence(s, cut_up, cut_lo)
            if not new_seq:
                continue

            # base id before any pipes
            base_id = rec.id.split("|", 1)[0]

            # always tag with platform
            new_id = f"{base_id}_{plat}"
            # if this is the min-lower record, add '_min'
            if rec is min_rec and orig_lower_counts[min_index] > 0:
                new_id += "_min"

            rec.id = new_id
            rec.description = ""  # drop any extra description
            rec.seq = Seq(new_seq)
            trimmed.append(rec)

        if trimmed:
            SeqIO.write(trimmed, out_path, "fasta")
            log.info(f"Trimmed {path} → {out_path} "
                     f"(up={cut_up}, lo={cut_lo}; min_rec={min_rec.id.split('|')[0]})")
        else:
            log.info(f"No sequences retained for {path}")

    except Exception:
        log.exception(f"Error processing {path}")

# ──────────────────────────────────────────────────────────────────────────────

def collect_fastas(root):
    files = []
    for d, _, fls in os.walk(root):
        for fn in fls:
            if fn.lower().endswith((".fa", ".fasta")):
                files.append(os.path.join(d, fn))
    return files

def merge_platforms(out_base):
    log = logging.getLogger(__name__)
    for side in ("L", "R"):
        dirn = os.path.join(out_base, f"trim_{side}")
        if not os.path.isdir(dirn):
            continue
        by_chr = defaultdict(set)
        for fn in os.listdir(dirn):
            if fn.endswith(f"_{side}.fasta"):
                parts = fn.split("_")
                if len(parts) >= 4:
                    chrn = parts[0]
                    plat = parts[2]
                    by_chr[chrn].add(plat)
        for chrn, plats in by_chr.items():
            if {"ONT", "HiFi"}.issubset(plats):
                merged_path = os.path.join(dirn, f"{chrn}_trimmed_merged_{side}.fasta")
                log.info(f"Merging ONT+HiFi into {merged_path}")
                with open(merged_path, "w") as out_h:
                    for plat in ("ONT", "HiFi"):
                        fn = f"{chrn}_trimmed_{plat}_{side}.fasta"
                        src = os.path.join(dirn, fn)
                        for rec in SeqIO.parse(src, "fasta"):
                            out_h.write(f">{rec.id}\n{rec.seq}\n")
                        os.remove(src)
                        log.info(f"Deleted original: {src}")

def run_parallel(ont_dir, hifi_dir, out_base, coverage, threads):
    log = logging.getLogger(__name__)
    all_fastas = collect_fastas(ont_dir) + collect_fastas(hifi_dir)
    log.info(f"{len(all_fastas)} FASTA files to trim.")
    with ThreadPoolExecutor(max_workers=threads) as exe:
        futures = [exe.submit(process_fasta, p, ont_dir, hifi_dir, out_base, coverage)
                   for p in all_fastas]
        for f in futures:
            f.result()
    merge_platforms(out_base)

def main():
    p = argparse.ArgumentParser(description="Trim & then merge ONT/HiFi reads by coverage")
    p.add_argument("--ont_dir",  required=True, help="Directory of ONT FASTAs")
    p.add_argument("--hifi_dir", required=True, help="Directory of HiFi FASTAs")
    p.add_argument("-o", "--out_dir",    required=True, help="Base output directory")
    p.add_argument("-c", "--coverage",   type=float, default=1.0,
                   help="Fraction of reads to keep (0<coverage≤1)")
    p.add_argument("-t", "--threads",    type=int,   default=4, help="Parallel threads")
    p.add_argument("-v", "--verbose",    action="store_true", help="Debug logging")
    args = p.parse_args()

    lvl = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=lvl, format="%(asctime)s %(levelname)s: %(message)s")

    run_parallel(args.ont_dir, args.hifi_dir, args.out_dir, args.coverage, args.threads)

if __name__ == "__main__":
    main()
