#!/usr/bin/env python


import multiprocessing
from multiprocessing import Process
from Bio import SeqIO
import argparse
import subprocess
import math
import os


#Step 1: Extract end alignment read and rename output data 
def rename_ONT_in_directory(directory, data_type):
    for filename in os.listdir(directory):
        parts = filename.split('_')
        if len(parts) >= 2:
            new_filename = f"{parts[0]}_{data_type}_{parts[1]}"
            os.rename(os.path.join(directory, filename), os.path.join(directory, new_filename))


def rename_HiFi_in_directory(directory, data_type):
    for filename in os.listdir(directory):
        parts = filename.split('_')
        if len(parts) >= 2:
            new_filename = f"{parts[0]}_{data_type}_{parts[1]}"
            os.rename(os.path.join(directory, filename), os.path.join(directory, new_filename))


def run_ont_alignments_extraction(genome, threads, fastq_ont, out_ONT_bam, data_type, ref_genome, motif,
                                  algn_output_ont):
    combined_command = [
        f"minimap2 -t {threads} -ax map-ont {genome} {fastq_ont} | teloclip --ref {ref_genome} | samtools sort > {out_ONT_bam}",
        f"samtools view -h {out_ONT_bam} | teloclip --ref {ref_genome} --motifs {motif} | teloclip-extract --refIdx {ref_genome} --extractReads --extractDir {algn_output_ont}"
    ]

    for cmd in combined_command:
        subprocess.run(cmd, shell=True, check=True)
    rename_ONT_in_directory(algn_output_ont, data_type)
    return algn_output_ont


def run_hifi_alignments_extraction(genome, threads, fastq_hifi, out_HiFi_bam, data_type, ref_genome, motif,
                                   algn_output_hifi):
    combined_command = [
        f"minimap2 -t {threads} -ax map-ont {genome} {fastq_hifi} | teloclip --ref {ref_genome} | samtools sort > {out_HiFi_bam}",
        f"samtools view -h {out_HiFi_bam} | teloclip --ref {ref_genome} --motifs {motif} | teloclip-extract --refIdx {ref_genome} --extractReads --extractDir {algn_output_hifi}"
    ]
    for cmd in combined_command:
        subprocess.run(cmd, shell=True, check=True)
    rename_HiFi_in_directory(algn_output_hifi, data_type)
    return algn_output_hifi


# Step 1': Extract end alignment reads for BAM
def run_BAM_ont_alignments_extraction(out_ONT_bam, data_type, ref_genome, motif, algn_output_ont):
    command1 = f"samtools view -h {out_ONT_bam} | teloclip --ref {ref_genome} --motifs {motif} | teloclip-extract --refIdx {ref_genome} --extractReads --extractDir {algn_output_ont}"
    subprocess.run(command1, shell=True, check=True)
    rename_ONT_in_directory(algn_output_ont, data_type)
    return algn_output_ont


def run_BAM_hifi_alignments_extraction(out_HiFi_bam, data_type, ref_genome, motif, algn_output_hifi):
    command2 = f"samtools view -h {out_HiFi_bam} | teloclip --ref {ref_genome} --motifs {motif} | teloclip-extract --refIdx {ref_genome} --extractReads --extractDir {algn_output_hifi}"
    subprocess.run(command2, shell=True, check=True)
    rename_HiFi_in_directory(algn_output_hifi, data_type)
    return algn_output_hifi


### Step 2: Merge and Trim Sequences
def merge_sequences(ont_file, hifi_file, merged_file):
    with open(merged_file, 'w') as out_handle:
        ont_records = {record.id: record for record in SeqIO.parse(ont_file, 'fasta')}
        hifi_records = {record.id: record for record in SeqIO.parse(hifi_file, 'fasta')}
        for record_id, ont_record in ont_records.items():
            out_handle.write(f'>{record_id}_ONT\n{ont_record.seq}\n')
            if record_id in hifi_records:
                hifi_record = hifi_records[record_id]
                out_handle.write(f'>{record_id}_HiFi\n{hifi_record.seq}\n')
        for record_id, hifi_record in hifi_records.items():
            if record_id not in ont_records:
                out_handle.write(f'>{record_id}_HiFi\n{hifi_record.seq}\n')


def calculate_lengths(input_file):
    lengths = {'lowercase': {}}
    for record in SeqIO.parse(input_file, 'fasta'):
        sequence = str(record.seq)
        lowercase_sequence = ''.join(char for char in sequence if char.islower())
        lengths['lowercase'][record.id] = len(lowercase_sequence)
    return lengths

# The length of coverage is selected
def calculate_lowercase_length(sequence):
    """Count the length of the lowercase portion of the sequence."""
    return sum(1 for c in sequence if c.islower())

def get_read_length_by_coverage(input_fasta, coverage):
    """Get the length of reads corresponding to the specified coverage."""
    # Read all reads in the fasta file
    reads = list(SeqIO.parse(input_fasta, "fasta"))

    #Calculate the length of the lowercase part of each read
    reads_with_length = [(read, calculate_lowercase_length(read.seq)) for read in reads]

    # Sort by length of lowercase part from longest to shortest
    sorted_reads = sorted(reads_with_length, key=lambda x: x[1], reverse=True)

    # Calculate the number of reads required
    total_reads = len(reads)
    num_reads_to_select = math.ceil(total_reads * (coverage / 100.0))

    # Get the length of reads to which the coverage falls
    if num_reads_to_select > total_reads:
        num_reads_to_select = total_reads

    selected_length = sorted_reads[num_reads_to_select - 1][1]

    return selected_length, sorted_reads[:num_reads_to_select]

def coverage_read_selection(input_file, chromosome, direction, seq_type, coverage, temp_files):
    """Select reads based on coverage and save them into new FASTA files."""
    selected_length, selected_reads = get_read_length_by_coverage(input_file, coverage)

    output_fasta = f"{chromosome}_{coverage}_{direction}.fasta"
    with open(output_fasta, "w") as output_handle:
        for read, _ in selected_reads:
            SeqIO.write(read, output_handle, "fasta")

    print(f"Selected reads saved to {output_fasta}")

    if output_fasta:
        # If the selected length corresponds to the first read, execute the function to extract the lowercase sequence
        if selected_length == selected_reads[0][1]:
            output_dir = f"trim_{direction}"
            os.makedirs(output_dir, exist_ok=True)
            output_file = f"{chromosome}_one_{seq_type}_{direction}.fasta"
            with open(output_file, 'w') as out_handle:
                for record in SeqIO.parse(output_fasta, 'fasta'):
                    # Extract the lowercase sequence
                    lowercase_seq = ''.join([char for char in str(record.seq) if char.islower()])
                    out_handle.write(f'>{record.id}\n{lowercase_seq}\n')

            # Move the output file to the specified directory
            subprocess.run(f"mv {output_file} {output_dir}", shell=True, check=True)

        else:
            # If the first read is not selected, execute the trim function
            lengths = calculate_lengths(output_fasta)
            min_record_id = min(lengths['lowercase'], key=lengths['lowercase'].get)
            min_length = min(lengths['lowercase'].values()) if lengths['lowercase'] else 0
            print(f"Min {chromosome}_{direction} lowercase sequence length on {direction}: {min_length}")
            trim_sequences(output_fasta, direction, min_length, min_record_id, seq_type, temp_files)

        # Delete output_fasta file
        if os.path.exists(output_fasta):
            os.remove(output_fasta)
            print(f"Deleted temporary file {output_fasta}")


def trim_sequences(input_file, direction, min_Length, min_record_id, seq_type, temp_files):
    base_name = os.path.basename(input_file).rsplit('_', 2)[0]
    input_FILE = f"{base_name}_min_{seq_type}_{direction}.fasta"

    with open(input_FILE, 'w') as out_handle:
        for record in SeqIO.parse(input_file, 'fasta'):
            if record.id == min_record_id:
                record.id += '_min'
                record.description = ""
            SeqIO.write(record, out_handle, 'fasta')
    temp_files.append(input_FILE)

    output_dir = f"trim_{direction}"
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{base_name}_trimmed_{seq_type}_{direction}.fasta")

    with open(output_file, 'w') as out_handle:
        for record in SeqIO.parse(input_FILE, 'fasta'):
            sequence = str(record.seq)
            uppercase_sequence = ''.join(char for char in sequence if char.isupper())
            lowercase_sequence = ''.join(char for char in sequence if char.islower())

            if direction == 'R':
                trimmed_lowsequence = lowercase_sequence[:min_Length] if len(
                    lowercase_sequence) > min_Length else lowercase_sequence
                out_handle.write(f'>{record.id}\n{uppercase_sequence}{trimmed_lowsequence}\n')
            elif direction == 'L':
                trimmed_lowsequence = lowercase_sequence[-min_Length:] if len(
                    lowercase_sequence) > min_Length else lowercase_sequence
                out_handle.write(f'>{record.id}\n{trimmed_lowsequence}{uppercase_sequence}\n')


def process_chromosome(chromosome, trim_merged_file, algn_output_ont, algn_output_hifi, direction, coverage, temp_files):
    ont_filename = f"{chromosome}_ONT_{direction}.fasta"
    hifi_filename = f"{chromosome}_HiFi_{direction}.fasta"
    input_file = None
    seq_type = None

    if ont_filename in trim_merged_file and hifi_filename in trim_merged_file:
        seq_type = "merged"
        ont_file = os.path.join(algn_output_ont, ont_filename)
        hifi_file = os.path.join(algn_output_hifi, hifi_filename)
        merged_file = f"{chromosome}_merged_{direction}.fasta"
        merge_sequences(ont_file, hifi_file, merged_file)
        input_file = merged_file
        temp_files.append(input_file)
    elif ont_filename in trim_merged_file:
        seq_type = "ONT"
        ont_file = os.path.join(algn_output_ont, ont_filename)
        input_file = f"{chromosome}_ONT_{direction}.fasta"
        with open(input_file, 'w') as out_handle:
            for record in SeqIO.parse(ont_file, 'fasta'):
                out_handle.write(f'>{record.id}_ONT\n{record.seq}\n')
        temp_files.append(input_file)
    elif hifi_filename in trim_merged_file:
        seq_type = "HiFi"
        hifi_file = os.path.join(algn_output_hifi, hifi_filename)
        input_file = f"{chromosome}_HiFi_{direction}.fasta"
        with open(input_file, 'w') as out_handle:
            for record in SeqIO.parse(hifi_file, 'fasta'):
                out_handle.write(f'>{record.id}_HiFi\n{record.seq}\n')
        temp_files.append(input_file)

    if input_file:
        coverage_read_selection(input_file, chromosome, direction, seq_type, coverage, temp_files)


def main():
    parser = argparse.ArgumentParser(
        prog='TeloComp',
        usage='%(prog)s [options]',
        exit_on_error=False,
        description='A tool for telomere extraction and genome polishing.',
        epilog='Text at the bottom of help')

    # Step 1: Telomere extraction arguments 
    parser.add_argument('-G', '--genome', metavar='', help='Input genome file (FASTA format)')
    parser.add_argument('-O', '--ONT', metavar='', help='Input ONT data')
    parser.add_argument('-H', '--HiFi', metavar='', help='Input HiFi data')
    parser.add_argument('-B', '--ONTbam', metavar='', help='The output bam file. If the [--BamExtr] parameter is selected, it is the input bam file.')
    parser.add_argument('-b', '--HiFibam', metavar='', help='The output bam file. If the [--BamExtr] parameter is selected, it is the input bam file.')
    parser.add_argument('--BamExtr', action='store_true', help='Selecting this parameter does not perform genome alignment to obtain bam, but directly inputs the sorted bam file, and then screens the qualified reads.')
    parser.add_argument('-r', '--ref', metavar='', required=True, help='Index file of reference genome')
    parser.add_argument('-c', '--coverage', type=float, metavar='', required=True, help='Choose the coverage you think is appropriate, that is, the minimum length of reads')
    parser.add_argument('-m', '--motif', metavar='', required=True,
                        help='Telomeric repeats sequences, e.g., plant: CCCTAAA(TTTAGGG), animal: TTAGGG(CCCTAA), etc.')
    parser.add_argument('-t', '--threads', type=int, default=20, metavar='',
                        help='Number of threads to use (default: 20)')

    args = parser.parse_args()

    # Create directories for alignment outputs if they do not exist
    algn_output_ont = 'algn_output_ont'
    algn_output_hifi = 'algn_output_hifi'
    if not os.path.exists(algn_output_ont):
        os.makedirs(algn_output_ont)
    if not os.path.exists(algn_output_hifi):
        os.makedirs(algn_output_hifi)

    # Step 1: Run ONT or HiFi alignment extraction, or process existing BAM files
    processes=[]
    if args.BamExtr:
        if args.ONTbam:
            run_BAM_ont_alignments_extraction(args.ONTbam, 'ONT', args.ref, args.motif, algn_output_ont)

        if args.HiFibam:
            run_BAM_hifi_alignments_extraction(args.HiFibam, 'HiFi', args.ref, args.motif, algn_output_hifi)


    else:
        if args.ONT:
            run_ont_alignments_extraction(args.genome, args.threads, args.ONT, args.ONTbam, 'ONT', args.ref, args.motif, args.algn_output_ont)

        if args.HiFi:
            run_hifi_alignments_extraction(args.genome, args.threads, args.HiFi, args.HiFibam, 'HiFi', args.ref, args.motif, algn_output_hifi)  


    # Step 2: Trim alignments sequence by highest coverage
    trim_ont_files = sorted(os.listdir(algn_output_ont))
    trim_hifi_files = sorted(os.listdir(algn_output_hifi))
    trim_merged_file = trim_ont_files + trim_hifi_files

    processed_chromosomes = set()
    temp_files = []

    for chromosome_file in trim_merged_file:
        if not chromosome_file.endswith('.fasta'):
            continue
        chromosome = chromosome_file.split('_')[0]
        if chromosome in processed_chromosomes:
            continue
        process_chromosome(chromosome, trim_merged_file, algn_output_ont, algn_output_hifi, 'L', args.coverage, temp_files)
        process_chromosome(chromosome, trim_merged_file, algn_output_ont, algn_output_hifi, 'R', args.coverage, temp_files)
        processed_chromosomes.add(chromosome)

    # Remove temporary files after processing
    for temp_file in temp_files:
        os.remove(temp_file)


if __name__ == "__main__":
    main()