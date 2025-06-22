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

# ------------------ Part 1: teloComp core functions ------------------

def processSamlines(
    samfile,
    ContigDict,
    motifList=[],
    matchAnywhere=False,
    maxBreak=50,
    minClip=1,
    noRev=False,
    fuzzy=False,
    minRepeats=3,
):
    # SAM line index keys
    SAM_QNAME = 0
    SAM_RNAME = 2
    SAM_POS = 3
    SAM_CIGAR = 5
    SAM_SEQ = 9

    bothCount = 0
    keepCount = 0
    motifCount = 0
    removeCount = 0
    samlineCount = 0

    # Read SAM from input
    for line in samfile:
        keepLine = False
        leftClip = False
        rightClip = False
        # Write headers unchanged
        if line[0][0] == "@":
            sys.stdout.write(line)
            continue
        samlineCount += 1
        samline = line.split("\t")
        # Check if line contains soft-clip and no hard-clipping.
        if "S" in samline[SAM_CIGAR] and not "H" in samline[SAM_CIGAR]:
            leftClipLen, rightClipLen = checkClips(samline[SAM_CIGAR])
            alnLen = lenCIGAR(samline[SAM_CIGAR])
            # Check for left overhang
            if leftClipLen:
                if (int(samline[SAM_POS]) <= maxBreak) and (leftClipLen >= (int(samline[SAM_POS]) + minClip)):
                    keepLine = True
                    leftClip = True
                    keepCount += 1
            # Check for right overhang
            if rightClipLen:
                try:
                    ContigLen = ContigDict[str(samline[SAM_RNAME])]
                except:
                    sys.exit("Reference sequence not found in FAI file: " + str(samline[SAM_RNAME]))
                alnEnd = int(samline[SAM_POS]) + alnLen
                if ((ContigLen - alnEnd) <= maxBreak) and (alnEnd + rightClipLen >= ContigLen + 1):
                    rightClip = True
                    if not keepLine:
                        keepLine = True
                        keepCount += 1
                    else:
                        logging.info(str(samline[SAM_QNAME]) + " overhang on both ends of " + str(samline[SAM_RNAME]))
                        bothCount += 1
            # Optional telomeric motif filtering
            if motifList and keepLine and matchAnywhere:
                if any(s in samline[SAM_SEQ] for s in motifList):
                    sys.stdout.write(line)
                    motifCount += 1
                else:
                    removeCount += 1
            elif motifList and keepLine:
                if isMotifInClip(samline, motifList, leftClip, rightClip, leftClipLen, rightClipLen):
                    sys.stdout.write(line)
                    motifCount += 1
                else:
                    removeCount += 1
            elif keepLine:
                sys.stdout.write(line)
            else:
                removeCount += 1
    if motifList:
        logging.info(
            f"Processed {samlineCount} SAM records.\n"
            f"Found {keepCount} alignments soft-clipped at contig ends.\n"
            f"Output {motifCount} alignments containing motif matches.\n"
            f"Discarded {removeCount} terminal alignments after filtering."
        )
    else:
        logging.info(
            f"Processed {samlineCount} SAM records.\n"
            f"Found {keepCount} alignments soft-clipped at contig ends.\n"
            f"Found {bothCount} alignments spanning entire contigs.\n"
            f"Discarded {removeCount} terminal alignments after filtering."
        )

def splitCIGAR(SAM_CIGAR):
    """
    Split CIGAR string into list of tuples with format (len, operator)
    """
    CIGARlist = []
    for x in re.findall("[0-9]*[A-Z|=]", SAM_CIGAR):
        CIGARlist.append((int(re.findall("[0-9]*", x)[0]), re.findall("[A-Z]|=", x)[0]))
    return CIGARlist

def checkClips(SAM_CIGAR):
    """
    Get lengths of soft-clipped blocks from either end of an alignment given a CIGAR string.
    """
    leftClipLen = None
    rightClipLen = None
    CIGARlist = splitCIGAR(SAM_CIGAR)
    if CIGARlist[0][1] == "S":
        leftClipLen = int(CIGARlist[0][0])
    if CIGARlist[-1][1] == "S":
        rightClipLen = int(CIGARlist[-1][0])
    return (leftClipLen, rightClipLen)

def lenCIGAR(SAM_CIGAR):
    """
    Calculate alignment length on reference as sum of M, D, N, X, = operators.
    """
    alnLen = 0
    CIGARlist = splitCIGAR(SAM_CIGAR)
    for x in CIGARlist:
        if x[1] in {"D", "M", "N", "X", "="}:
            alnLen += x[0]
    return alnLen

def StreamingSamFilter(samfile=None, contigs=None, maxBreak=50, minClip=1):
    """Rewrite loadSam() as generator."""
    SAM_QNAME = 0
    SAM_RNAME = 2
    SAM_POS = 3
    SAM_CIGAR = 5
    SAM_SEQ = 9
    for line in samfile:
        if line[0][0] == "@":
            continue
        samline = line.split("\t")
        if "S" in samline[SAM_CIGAR] and not "H" in samline[SAM_CIGAR]:
            leftClipLen, rightClipLen = checkClips(samline[SAM_CIGAR])
            alnLen = lenCIGAR(samline[SAM_CIGAR])
            if leftClipLen:
                if (int(samline[SAM_POS]) <= maxBreak) and (leftClipLen >= (int(samline[SAM_POS]) + minClip)):
                    alnEnd = int(samline[SAM_POS]) + alnLen
                    try:
                        yield (samline[SAM_POS], alnEnd, leftClipLen, samline[SAM_SEQ], samline[SAM_QNAME], samline[SAM_RNAME], "L")
                    except:
                        logging.warning("Reference sequence not found in FAI file: " + str(samline[SAM_RNAME]))
            if rightClipLen:
                try:
                    ContigLen = contigs[str(samline[SAM_RNAME])]
                except:
                    logging.warning("Reference sequence not found in FAI file: " + str(samline[SAM_RNAME]))
                alnEnd = int(samline[SAM_POS]) + alnLen
                if ((ContigLen - alnEnd) <= maxBreak) and (alnEnd + rightClipLen >= ContigLen + 1):
                    yield (samline[SAM_POS], alnEnd, rightClipLen, samline[SAM_SEQ], samline[SAM_QNAME], samline[SAM_RNAME], "R")

def isMotifInClip(samline, motifList, leftClip, rightClip, leftClipLen, rightClipLen):
    SAM_SEQ = 9
    leftcheck = False
    rightcheck = False
    if leftClip:
        leftcheck = check_sequence_for_patterns(samline[SAM_SEQ][0:leftClipLen], motifList)
    if rightClip:
        rightcheck = check_sequence_for_patterns(samline[SAM_SEQ][-rightClipLen:], motifList)
    return any([leftcheck, rightcheck])

# ------------------ Part 2: Auxiliary tool functions ------------------

def read_fai(fai):
    """
    Read the FAI index file and return a dictionary with the key being the reference sequence name and the value     being the sequence length.
    """
    ContigDict = dict()
    with open(fai, "r") as f:
        for line in f.readlines():
            li = line.strip().split()
            ContigDict[li[0]] = int(li[1])
    return ContigDict

def crunchHomopolymers(motifList):
    crunchList = list()
    for motif in motifList:
        noReps = list()
        for base in motif:
            if not noReps:
                noReps.append(base)
            elif base != noReps[-1]:
                noReps.append(base)
        crunchList.append("".join(noReps))
    return list(set(crunchList))

def check_sequence_for_patterns(dna_sequence, regex_patterns):
    return any(re.search(pattern, dna_sequence) for pattern in regex_patterns)

# ------------------ Main process: parallel running of ONT/HiFi processes ------------------

def run_pipeline(read_label, genome, read_file, preset, threads, contig_dict, motifs, max_break, min_clip, output_bam):
    """
	For single read data:
	1. Use minimap2 to align and output SAM file;
	2. Use teloclip to process SAM file and output filtered SAM file;
	3. Convert filtered SAM to BAM file.
    """
    # Create a temporary SAM file output by minimap2
    with tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".sam") as tmp_minimap2:
        minimap2_sam = tmp_minimap2.name
    minimap2_cmd = ["minimap2", "-t", str(threads), "-ax", preset, genome, read_file]
    logging.info(f"Running minimap2 for {read_label}: " + " ".join(minimap2_cmd))
    with open(minimap2_sam, "w") as sam_out:
        subprocess.run(minimap2_cmd, stdout=sam_out, check=True)
    logging.info(f"{read_label} alignment completed. SAM output saved to {minimap2_sam}")

    # Filter: Call teloComp to process the SAM file and output the filtered SAM file
    with tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".sam") as tmp_filtered:
        filtered_sam = tmp_filtered.name
    with open(minimap2_sam, "r") as in_sam, open(filtered_sam, "w") as out_filtered:
        original_stdout = sys.stdout
        sys.stdout = out_filtered
        processSamlines(in_sam, contig_dict, motifList=motifs, maxBreak=max_break, minClip=min_clip)
        sys.stdout = original_stdout
    logging.info(f"Teloclip filtering for {read_label} completed. Filtered SAM saved to {filtered_sam}")

    # SAM to BAM
    logging.info(f"Converting filtered SAM to BAM for {read_label}.")
    with pysam.AlignmentFile(filtered_sam, "r") as samfile, \
         pysam.AlignmentFile(output_bam, "wb", template=samfile) as bamfile:
        for read in samfile:
            bamfile.write(read)
    logging.info(f"BAM file for {read_label} written to {output_bam}")

    # Clean up temporary files
    os.remove(minimap2_sam)
    os.remove(filtered_sam)

def main():
    # 自定义 Formatter，不显示 metavar
    class NoMetavarFormatter(argparse.HelpFormatter):
        def _format_action_invocation(self, action):
            if not action.option_strings:
                return super()._format_action_invocation(action)

            # 对于带参数的选项，只显示选项名，不显示 metavar
            return ', '.join(action.option_strings)

    parser = argparse.ArgumentParser(
        description="Pipeline: Align the genome FASTA and its FAI index with ONT/HiFi data, process the SAM file using teloclip, and output BAM files separately.",
        formatter_class=NoMetavarFormatter
    )
    parser.add_argument("--genome", required=True, help="Input genome FASTA file.")
    parser.add_argument("--fai", required=True, help="Input genome index (FAI) file.")
    parser.add_argument("--ont", help="Input ONT data file (optional).")
    parser.add_argument("--hifi", help="Input HiFi data file (optional).")
    parser.add_argument("--threads", type=int, default=5, help="Number of threads to use with minimap2.")
    parser.add_argument("--motifs", nargs="*", default=[], help="A list of telomeric repeat motifs to use for filtering (optional).")
    parser.add_argument("--max_break", type=int, default=50, help="Maximum tolerable fracture length for soft shear.")
    parser.add_argument("--min_clip", type=int, default=1, help="Minimum cutting length.")
    # Specify the output BAM files for ONT and HiFi respectively
    parser.add_argument("--Ob", help="BAM output path after ONT filtering.")
    parser.add_argument("--Hb", help="HiFi filtered BAM output path.")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    # Read the FAI file to obtain the reference sequence length information
    contig_dict = read_fai(args.fai)

    tasks = []
    with ThreadPoolExecutor(max_workers=2) as executor:
        if args.ont:
            output_bam_ont = args.Ob if args.Ob else "ont_output.bam"
            tasks.append(executor.submit(run_pipeline, "ONT", args.genome, args.ont, "map-ont", args.threads, contig_dict, args.motifs, args.max_break, args.min_clip, output_bam_ont))
        if args.hifi:
            output_bam_hifi = args.Hb if args.Hb else "hifi_output.bam"
            tasks.append(executor.submit(run_pipeline, "HiFi", args.genome, args.hifi, "map-hifi", args.threads, contig_dict, args.motifs, args.max_break, args.min_clip, output_bam_hifi))
        for future in as_completed(tasks):
            future.result()

if __name__ == "__main__":
    main()
