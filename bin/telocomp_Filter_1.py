#!/usr/bin/env python3
import argparse
import logging
import re
import sys
import subprocess
import tempfile
import pysam
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from itertools import groupby

# ------------------ Section A: Core alignment processing routines ------------------

def handle_alignment_records(
    sam_stream,
    reference_lengths,
    motif_patterns=None,
    any_match_allowed=False,
    max_gap=50,
    min_softclip=1,
    exclude_reverse=False,
    fuzzy_match=False,
    min_repeat_count=3,
):
    """
    Process lines from a SAM stream, filter for soft-clipped reads near contig ends,
    optionally retain only those containing motifs, and write filtered SAM lines to stdout.

    Important: behavior preserved exactly from the original implementation.
    """

    # Indices used when parsing a SAM row split by tabs
    FIELD_QUERY = 0
    FIELD_REF = 2
    FIELD_POS = 3
    FIELD_CIG = 5
    FIELD_SEQ = 9

    count_both_ends = 0
    count_kept = 0
    count_motif_hits = 0
    count_removed = 0
    total_lines = 0

    motif_patterns = motif_patterns or []

    for raw_line in sam_stream:
        keep_this = False
        left_overhang = False
        right_overhang = False

        # preserve header lines unchanged
        if raw_line[0][0] == "@":
            sys.stdout.write(raw_line)
            continue

        total_lines += 1
        parts = raw_line.split("\t")

        # only consider reads with soft-clipping and no hard-clipping
        if "S" in parts[FIELD_CIG] and "H" not in parts[FIELD_CIG]:
            left_len, right_len = extract_clip_lengths(parts[FIELD_CIG])
            aln_ref_len = reference_length_from_cigar(parts[FIELD_CIG])

            # left-side overhang check
            if left_len:
                if (int(parts[FIELD_POS]) <= max_gap) and (left_len >= (int(parts[FIELD_POS]) + min_softclip)):
                    keep_this = True
                    left_overhang = True
                    count_kept += 1

            # right-side overhang check
            if right_len:
                try:
                    contig_length = reference_lengths[str(parts[FIELD_REF])]
                except:
                    sys.exit("Reference sequence not found in FAI file: " + str(parts[FIELD_REF]))
                aln_end = int(parts[FIELD_POS]) + aln_ref_len
                if ((contig_length - aln_end) <= max_gap) and (aln_end + right_len >= contig_length + 1):
                    right_overhang = True
                    if not keep_this:
                        keep_this = True
                        count_kept += 1
                    else:
                        logging.info(str(parts[FIELD_QUERY]) + " overhang on both ends of " + str(parts[FIELD_REF]))
                        count_both_ends += 1

            # motif-based filtering (two modes: anywhere in sequence or only within clipped region)
            if motif_patterns and keep_this and any_match_allowed:
                if any(pat in parts[FIELD_SEQ] for pat in motif_patterns):
                    sys.stdout.write(raw_line)
                    count_motif_hits += 1
                else:
                    count_removed += 1
            elif motif_patterns and keep_this:
                if motif_in_clipped_region(parts, motif_patterns, left_overhang, right_overhang, left_len, right_len):
                    sys.stdout.write(raw_line)
                    count_motif_hits += 1
                else:
                    count_removed += 1
            elif keep_this:
                sys.stdout.write(raw_line)
            else:
                count_removed += 1

    # Logging summary (message content preserved)
    if motif_patterns:
        logging.info(
            f"Processed {total_lines} SAM records.\n"
            f"Found {count_kept} alignments soft-clipped at contig ends.\n"
            f"Output {count_motif_hits} alignments containing motif matches.\n"
            f"Discarded {count_removed} terminal alignments after filtering."
        )
    else:
        logging.info(
            f"Processed {total_lines} SAM records.\n"
            f"Found {count_kept} alignments soft-clipped at contig ends.\n"
            f"Found {count_both_ends} alignments spanning entire contigs.\n"
            f"Discarded {count_removed} terminal alignments after filtering."
        )


def decompose_cigar(cigar_string):
    """
    Break a CIGAR string into a list of (length, operator) tuples.
    Example: '10M5S' -> [(10,'M'), (5,'S')]
    """
    result = []
    for token in re.findall(r"[0-9]+[A-Z=]", cigar_string):
        num = int(re.findall(r"[0-9]+", token)[0])
        op = re.findall(r"[A-Z]|=", token)[0]
        result.append((num, op))
    return result


def extract_clip_lengths(cigar_string):
    """
    Return (left_softclip_len, right_softclip_len) or (None, None) when absent.
    """
    left = None
    right = None
    ops = decompose_cigar(cigar_string)
    if ops and ops[0][1] == "S":
        left = int(ops[0][0])
    if ops and ops[-1][1] == "S":
        right = int(ops[-1][0])
    return (left, right)


def reference_length_from_cigar(cigar_string):
    """
    Compute alignment span on the reference from a CIGAR string by summing
    operators that consume reference bases: M, D, N, X, =
    """
    total = 0
    ops = decompose_cigar(cigar_string)
    for ln, op in ops:
        if op in {"D", "M", "N", "X", "="}:
            total += ln
    return total


def stream_alignment_filter(sam_stream=None, contig_map=None, max_gap=50, min_softclip=1):
    """
    Generator-based variant to yield compact alignment information for
    downstream processing without emitting headers.
    Yields tuples preserving original semantics: (pos, alnEnd, clipLen, seq, qname, rname, side)
    """
    FIELD_QUERY = 0
    FIELD_REF = 2
    FIELD_POS = 3
    FIELD_CIG = 5
    FIELD_SEQ = 9

    for line in sam_stream:
        if line[0][0] == "@":
            continue
        parts = line.split("\t")
        if "S" in parts[FIELD_CIG] and "H" not in parts[FIELD_CIG]:
            left_len, right_len = extract_clip_lengths(parts[FIELD_CIG])
            aln_ref_len = reference_length_from_cigar(parts[FIELD_CIG])
            if left_len:
                if (int(parts[FIELD_POS]) <= max_gap) and (left_len >= (int(parts[FIELD_POS]) + min_softclip)):
                    aln_end = int(parts[FIELD_POS]) + aln_ref_len
                    try:
                        yield (parts[FIELD_POS], aln_end, left_len, parts[FIELD_SEQ], parts[FIELD_QUERY], parts[FIELD_REF], "L")
                    except:
                        logging.warning("Reference sequence not found in FAI file: " + str(parts[FIELD_REF]))
            if right_len:
                try:
                    contig_len = contig_map[str(parts[FIELD_REF])]
                except:
                    logging.warning("Reference sequence not found in FAI file: " + str(parts[FIELD_REF]))
                    continue
                aln_end = int(parts[FIELD_POS]) + aln_ref_len
                if ((contig_len - aln_end) <= max_gap) and (aln_end + right_len >= contig_len + 1):
                    yield (parts[FIELD_POS], aln_end, right_len, parts[FIELD_SEQ], parts[FIELD_QUERY], parts[FIELD_REF], "R")


def motif_in_clipped_region(sam_parts, motif_patterns, left_flag, right_flag, left_len, right_len):
    """
    Check whether any motif pattern occurs inside the clipped portion of the read sequence.
    Returns True if a motif is found in either clipped region (left or right) depending on flags.
    """
    SEQ_IDX = 9
    left_found = False
    right_found = False
    if left_flag:
        left_found = search_patterns_in_seq(sam_parts[SEQ_IDX][0:left_len], motif_patterns)
    if right_flag:
        right_found = search_patterns_in_seq(sam_parts[SEQ_IDX][-right_len:], motif_patterns)
    return any([left_found, right_found])


# ------------------ Section B: Utility helpers ------------------

def load_fai_index(fai_path):
    """
    Read a FASTA index (.fai) and return a dict {reference_name: length}.
    """
    refs = {}
    with open(fai_path, "r") as fh:
        for line in fh:
            cols = line.strip().split()
            refs[cols[0]] = int(cols[1])
    return refs


def compress_homopolymers(motif_list):
    """
    For each motif, collapse runs of the same base into a single base.
    Returns a unique set of the simplified motifs.
    """
    out = []
    for motif in motif_list:
        compact = []
        for b in motif:
            if not compact:
                compact.append(b)
            elif b != compact[-1]:
                compact.append(b)
        out.append("".join(compact))
    return list(set(out))


def search_patterns_in_seq(dna_seq, regex_list):
    """
    Return True if any regex pattern in regex_list matches dna_seq.
    """
    return any(re.search(pat, dna_seq) for pat in regex_list)


# ------------------ Section C: Pipeline orchestration (minimap2 -> filter -> BAM) ------------------

def pipeline_worker(label, reference_fasta, reads_file, minimap_preset, threads, contig_lengths, motif_list, max_gap, min_softclip, output_bam_path):
    """
    For one read set (ONT or HiFi), run:
      1) minimap2 to produce SAM,
      2) our filtering logic to produce a filtered SAM,
      3) convert the filtered SAM to BAM using pysam.

    All temporary files are created and removed automatically.
    The function logs progress at each major step.
    """

    # Create a temporary SAM path for minimap2 output
    with tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".sam") as tmp1:
        minimap_sam_path = tmp1.name

    # Build and execute minimap2 command
    minimap_cmd = ["minimap2", "-t", str(threads), "-ax", minimap_preset, reference_fasta, reads_file]
    logging.info(f"Running minimap2 for {label}: " + " ".join(minimap_cmd))
    with open(minimap_sam_path, "w") as outf:
        subprocess.run(minimap_cmd, stdout=outf, check=True)
    logging.info(f"{label} alignment completed. SAM output saved to {minimap_sam_path}")

    # Create a temporary SAM path for the filtered output
    with tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".sam") as tmp2:
        filtered_sam_path = tmp2.name

    # Open and run the filtering routine, redirecting its stdout into the filtered SAM file
    with open(minimap_sam_path, "r") as sam_in, open(filtered_sam_path, "w") as sam_out:
        prev_stdout = sys.stdout
        sys.stdout = sam_out
        handle_alignment_records(sam_in, contig_lengths, motif_patterns=motif_list, max_gap=max_gap, min_softclip=min_softclip)
        sys.stdout = prev_stdout
    logging.info(f"Teloclip filtering for {label} completed. Filtered SAM saved to {filtered_sam_path}")

    # Convert filtered SAM to BAM using pysam (preserve template)
    logging.info(f"Converting filtered SAM to BAM for {label}.")
    with pysam.AlignmentFile(filtered_sam_path, "r") as samfile, \
         pysam.AlignmentFile(output_bam_path, "wb", template=samfile) as bamfile:
        for read in samfile:
            bamfile.write(read)
    logging.info(f"BAM file for {label} written to {output_bam_path}")

    # Remove temporary files
    os.remove(minimap_sam_path)
    os.remove(filtered_sam_path)


def main():
    """
    Command-line entrypoint: parse arguments, set up logging, read FAI and run worker threads.
    The CLI interface (argument names) remains unchanged to preserve usability.
    """

    # Custom help formatter that hides metavar values for options (keeps original UX)
    class CompactHelpFormatter(argparse.HelpFormatter):
        def _format_action_invocation(self, action):
            if not action.option_strings:
                return super()._format_action_invocation(action)
            return ', '.join(action.option_strings)

    parser = argparse.ArgumentParser(
        description="Pipeline: Align the genome FASTA and its FAI index with ONT/HiFi data, process the SAM file using teloclip, and output BAM files separately.",
        formatter_class=CompactHelpFormatter
    )

    # CLI arguments preserved (names and semantics identical to original)
    parser.add_argument("--genome", required=True, help="Input genome FASTA file.")
    parser.add_argument("--fai", required=True, help="Input genome index (FAI) file.")
    parser.add_argument("--ont", help="Input ONT data file (optional).")
    parser.add_argument("--hifi", help="Input HiFi data file (optional).")
    parser.add_argument("--threads", type=int, default=5, help="Number of threads to use with minimap2.")
    parser.add_argument("--motifs", nargs="*", default=[], help="A list of telomeric repeat motifs to use for filtering (optional).")
    parser.add_argument("--max_break", type=int, default=50, help="Maximum tolerable fracture length for soft shear.")
    parser.add_argument("--min_clip", type=int, default=1, help="Minimum cutting length.")
    parser.add_argument("--Ob", help="BAM output path after ONT filtering.")
    parser.add_argument("--Hb", help="HiFi filtered BAM output path.")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    # load contig lengths from FAI index (exact same behavior)
    contigs = load_fai_index(args.fai)

    futures = []
    with ThreadPoolExecutor(max_workers=2) as pool:
        if args.ont:
            out_ont = args.Ob if args.Ob else "ont_output.bam"
            futures.append(pool.submit(pipeline_worker, "ONT", args.genome, args.ont, "map-ont", args.threads, contigs, args.motifs, args.max_break, args.min_clip, out_ont))
        if args.hifi:
            out_hifi = args.Hb if args.Hb else "hifi_output.bam"
            futures.append(pool.submit(pipeline_worker, "HiFi", args.genome, args.hifi, "map-hifi", args.threads, contigs, args.motifs, args.max_break, args.min_clip, out_hifi))
        for fut in as_completed(futures):
            fut.result()


if __name__ == "__main__":
    main()
