#!/usr/bin/env python3
"""
combined_telo_pipeline.py

1) Scan ONT and HiFi BAMs for all 5–10-bp repetitive substrings, canonicalize them,
   and print the top 10.
2) Prompt user to select one of the top-10 motifs or enter custom motifs (ACGT only).
3) Save all discovered canonical motifs to all_telomere_motifs.txt.
4) Extract telomeric overhang reads using the chosen motifs.
"""
import argparse
import logging
import os
import re
import pysam
import readline #support backspace, arrow keys, and other basic editing on the command line.
from collections import defaultdict

# ─── Part 1: motif discovery ──────────────────────────────────────────────────

def reverse_complement_part1(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(complement.get(b, 'N') for b in reversed(seq.upper()))

def get_canonical_sequence(seq):
    s = seq.upper()
    rc = reverse_complement_part1(s)
    rots = [s[i:]+s[:i] for i in range(len(s))] + [rc[i:]+rc[:i] for i in range(len(rc))]
    return min(rots)

def find_telomere_motifs_in_sequence(seq, cache, counts):
    window = 10
    n = len(seq)
    if n < window:
        return
    for i in range(n - window + 1):
        w = seq[i:i+window]
        subc = defaultdict(int)
        for L in range(5, 11):
            for j in range(window - L + 1):
                subc[w[j:j+L]] += 1
        for sub, cnt in subc.items():
            canon = cache.setdefault(sub, get_canonical_sequence(sub))
            counts[canon] += cnt

def process_bam_for_motifs(path, cache, counts):
    log = logging.getLogger()
    log.info(f"Opening BAM for motif discovery: {path}")
    bam = pysam.AlignmentFile(path, 'rb')
    for read in bam.fetch(until_eof=True):
        seq = read.query_sequence or ''
        if len(seq) < 10:
            continue
        u = seq.upper()
        find_telomere_motifs_in_sequence(u, cache, counts)
        find_telomere_motifs_in_sequence(reverse_complement_part1(u), cache, counts)
    bam.close()
    log.info(f"Finished scanning {path}")

# ─── Part 2: overhang extraction ───────────────────────────────────────────────

def reverse_complement_part2(seq: str) -> str:
    comp = str.maketrans('ACGTacgt', 'TGCAtgca')
    return seq.translate(comp)[::-1]

def check_sequence_for_patterns(sequence: str, patterns: list) -> int:
    seq = sequence.upper()
    total = 0
    for m in patterns:
        rc = reverse_complement_part2(m)
        total += seq.count(m) + seq.count(rc)
    return total


# ----- Derived from Teloclip -----
# Concept adapted and rewritten by TleComp (2025).
# Original Teloclip functions: StreamingSamFilter, StreamingSplitByContig.
# Implemented here using pysam and direct regex parsing for improved performance.
def process_bam_overhang(bam_path, outdir, motifs, min_ratio, logger):
    logger.info(f"Opening BAM for overhang extraction: {bam_path}")
    sam = pysam.AlignmentFile(bam_path, 'rb')
    contig_len = {h['SN']: h['LN'] for h in sam.header.get('SQ', [])}
    left_re  = re.compile(r'^(\d+)S')
    right_re = re.compile(r'(\d+)S$')
    buckets = {}
    total = extracted = 0

    for r in sam.fetch(until_eof=True):
        total += 1
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        rn, pos, end = r.reference_name, r.reference_start+1, r.reference_end
        cigar, seq = r.cigarstring or '', r.query_sequence or ''
        m = left_re.match(cigar)
        if m and pos == 1:
            buckets.setdefault((rn, 'L'), []).append({'qn':r.query_name, 'seq':seq, 'clip':int(m.group(1))})
        m = right_re.search(cigar)
        if m and contig_len.get(rn) == end:
            buckets.setdefault((rn, 'R'), []).append({'qn':r.query_name, 'seq':seq, 'clip':int(m.group(1))})
    sam.close()

    os.makedirs(outdir, exist_ok=True)
    handles = {}
    for (rn, side), reads in buckets.items():
        for r in reads:
            aligned = len(r['seq']) - r['clip']
            if aligned / len(r['seq']) < min_ratio or aligned < 10000:
                continue
            seq, c = r['seq'], r['clip']
            masked = seq[:c].lower() + seq[c:].upper() if side=='L' else seq[:-c].upper() + seq[-c:].lower()
            if check_sequence_for_patterns(masked, motifs) == 0:
                continue
            fn = os.path.join(outdir, f"{rn}_{side}.fasta")
            if (rn, side) not in handles:
                handles[(rn, side)] = open(fn, 'w')
            handles[(rn, side)].write(f">{r['qn']}|{rn}_{side}\n{masked}\n")
            extracted += 1

    for h in handles.values():
        h.close()
    logger.info(f"[{os.path.basename(bam_path)}] scanned={total}, extracted={extracted}")

# ─── Main ────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Discover motifs and extract telomeric overhangs"
    )
    parser.add_argument('--ont',   required=True, help="ONT BAM")
    parser.add_argument('--hifi',  required=True, help="HiFi BAM")
    parser.add_argument('-o','--outdir', default='.', help="Output directory")
    parser.add_argument('-v','--verbose', action='store_true', help="Verbose logging")
    parser.add_argument('--min_ratio', type=float, default=0.2, 
                        help="The proportion of the original genome sequence to the length of the reads, default=0.2")
    args = parser.parse_args()

    # configure logging with timestamp
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    log = logging.getLogger()
    log.info("Starting motif scan + overhang collection")

    # Part 1: motif discovery
    cache, counts = {}, defaultdict(int)
    log.info(f"Scanning {args.ont} for motifs…")
    process_bam_for_motifs(args.ont, cache, counts)
    log.info(f"Scanning {args.hifi} for motifs…")
    process_bam_for_motifs(args.hifi, cache, counts)

    sorted_m = sorted(counts.items(), key=lambda x: x[1], reverse=True)
    # 2) 输出全部端粒类型到文件
    all_path = os.path.join(args.outdir, "all_telomere_motifs.txt")
    with open(all_path, "w") as f:
        for motif, cnt in sorted_m:
            f.write(f"{motif}\t{cnt}\n")
    log.info(f"All motifs saved to {all_path}")

    # 打印 top10
    top10 = [motif for motif,_ in sorted_m[:10]]
    print("\nTop 10 canonical telomeric motifs:")
    for i, m in enumerate(top10, 1):
        print(f"  {i}. {m}  ({counts[m]} hits)")

    # 3) 交互选择：输入 1–10 选择 top10，或直接输入逗号分隔的自定义 ACGT motif 列表
    while True:
        sel = input(
            "\nSelect one motif [1-10] or enter custom motifs (comma-separated A/C/G/T): "
        ).strip()
        # 选择数字
        if sel.isdigit() and 1 <= int(sel) <= 10:
            motifs = [top10[int(sel)-1]]
            break
        # 自定义列表
        parts = [x.strip().upper() for x in sel.split(',')]
        if parts and all(re.fullmatch(r'[ACGT]+', p) for p in parts):
            motifs = parts
            break
        print("Invalid input; please enter a number 1–10 or valid A/C/G/T motifs.")

    log.info(f"Motifs chosen: {motifs}")

    # Part 2: overhang extraction
    ont_out  = os.path.join(args.outdir, 'ONT')
    hifi_out = os.path.join(args.outdir, 'HiFi')
    process_bam_overhang(args.ont,  ont_out, motifs, args.min_ratio, log)
    process_bam_overhang(args.hifi, hifi_out, motifs, args.min_ratio, log)

    log.info("Pipeline complete.")

if __name__ == '__main__':
    main()