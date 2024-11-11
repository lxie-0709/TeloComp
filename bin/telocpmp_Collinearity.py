#!/usr/bin/env python


import sys
import os
import re
import argparse
import shutil
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def calculate_sequence_frequencies(sequence, bin, step_size):
    """
    Calculate the frequency of specific sequences in a given DNA sequence within non-overlapping windows.
    :param sequence: The DNA sequence to analyze.
    :param window_size: The size of the sliding window in base pairs.
    :param step_size: The step size for the sliding window in base pairs.
    :return: A list of tuples (start, end, frequency) for each window.
    """
    seq_length = len(sequence)
    frequencies = []

    for start in range(0, seq_length, bin):
        end = min(start + bin, seq_length)
        window_seq = sequence[start:end]

        # Calculate motif counts within the window
        count_ccctaaa = window_seq.count('CCCTAAA')
        count_tttaggg = window_seq.count('TTTAGGG')
        total_count = count_ccctaaa + count_tttaggg

        # Calculate frequency based on the window size
        frequency = total_count

        frequencies.append((start, end, frequency))  # +1 to make start 1-based

    return frequencies


def write_bed_file(filename, sequences, chr_name, direction):
    """
    Append the frequency data to a BED file.
    :param filename: The name of the BED file to write.
    :param sequences: List of (start, end, frequency) tuples.
    :param chr_name: The chromosome name.
    :param direction: The direction (L or R).
    """
    with open(filename, 'a') as bed_file:  # Append mode ('a')
        for start, end, frequency in sequences:
            bed_file.write(f"{chr_name}_{direction}\t{start}\t{end}\t{frequency:.2f}\n")


def extract_chr_and_direction(header):
    """
    Extract chromosome name and direction from the header.
    :param header: The header string, e.g., 'chr1_L'.
    :return: A tuple (chromosome name, direction) or (None, None) if format is invalid.
    """
    parts = header.rsplit('_', 1)
    if len(parts) == 2:
        return parts[0], parts[1]
    else:
        return None, None


# Calculate the total number of telomeres for the 100kb intercepted from both ends of the new and old genomes,
# and input them into the information files telomere1.num.info and telomere2.num.info in the form of chromosome, chr_length, and number respectively
def extract_telomere1_sequences(genome1, genome2, positions_txt, length, motif, bin, log_file):
    # Read the positions file
    telomere_positions = []
    specified_chromosomes = set()

    try:
        with open(positions_txt, 'r') as pos_file:
            header_skipped = False
            for line in pos_file:
                if line.strip():
                    if not header_skipped:
                        header_skipped = True
                        continue  # Skip the header line
                    parts = line.strip().split()
                    if len(parts) == 4:
                        chromosome = parts[0]
                        try:
                            start = int(parts[1])
                            end = int(parts[2])
                            length_value = int(parts[3])
                            direction = chromosome.split('_')[-1]
                            telomere_positions.append((chromosome, start, end, length_value, direction))
                            chr_name = chromosome.split('_')[0]
                            specified_chromosomes.add((chr_name, direction))
                        except ValueError:
                            print(f"Warning: Invalid line format or non-integer values in positions file: {line.strip()}")
                            continue
                    else:
                        print(f"Warning: Invalid line format in positions file: {line.strip()}")
    except FileNotFoundError:
        print(f"Error: The positions file {positions_txt} was not found.")
        sys.exit(1)

    # Read the genome sequences
    genome_sequences = {}
    try:
        for record in SeqIO.parse(genome2, "fasta"):
            genome_sequences[record.id] = record.seq
    except FileNotFoundError:
        print(f"Error: The genome FASTA file {genome2} was not found.")
        sys.exit(1)

    # Sort function for chromosome names
    def natural_sort_key(s):
        return [int(text) if text.isdigit() else text for text in re.split(r'(\d+)', s)]

    with open(log_file, 'w') as log_handle, open("telomere.new.num.info", 'w') as tel1_info:
        log_handle.write("Processed Chromosomes:\n")
        sequences = []

        # Process specified telomere positions
        telomere_positions.sort(key=lambda x: natural_sort_key(x[0]))

        for chromosome, start, end, length_value, direction in telomere_positions:
            chr_name = chromosome.split('_')[0]
            if chr_name not in genome_sequences:
                print(f"Warning: Chromosome {chr_name} not found in genome FASTA.")
                continue

            sequence = genome_sequences[chr_name]

            if direction == 'L':
                left_telomere = sequence[:end]
                right_end = min(len(sequence), end + length)
                right_telomere = sequence[end:right_end]
                combined_sequence = left_telomere + right_telomere
                sequences.append((f"{chromosome}", combined_sequence))
                print(f"Prepared {chromosome} for output")

                combined_sequences = calculate_sequence_frequencies(combined_sequence, bin, step_size=10)
                write_bed_file("telo_1.bed", combined_sequences, chr_name, direction)

                # Calculate the number of telomeres in a 100kb segment and write the information to a file
                total_telomere_count = count_telomere_motifs(left_telomere, motif) + count_telomere_motifs(right_telomere, motif)
                tel1_info.write(f"{chr_name}_L\t{len(sequence)}\t{total_telomere_count}\n")

            elif direction == 'R':
                left_start = max(0, start - length)
                left_telomere = sequence[left_start:start]
                right_telomere = sequence[start:end]
                combined_sequence = left_telomere + right_telomere
                sequences.append((f"{chromosome}", combined_sequence))
                print(f"Prepared {chromosome} for output")

                combined_sequences = calculate_sequence_frequencies(combined_sequence, bin, step_size=10)
                write_bed_file("telo_1.bed", combined_sequences, chr_name, direction)

                # Calculate the number of telomeres in a 100kb segment and write the information to a file
                total_telomere_count = count_telomere_motifs(left_telomere, motif) + count_telomere_motifs(right_telomere, motif)
                tel1_info.write(f"{chr_name}_R\t{len(sequence)}\t{total_telomere_count}\n")

            log_handle.write(f"{chr_name}_{direction}\n")

        # Ensure all entries have a valid format
        filtered_sequences = []
        for header, seq in sequences:
            chr_name, direction = extract_chr_and_direction(header)
            if chr_name and direction:
                filtered_sequences.append((header, seq, chr_name, direction))

        # Sort filtered sequences by chromosome name and direction
        filtered_sequences.sort(key=lambda x: (natural_sort_key(x[2]), x[3]))

        output1_fasta = f"new_genome_{length}.fasta"
        # Write sorted sequences to output FASTA file
        with open(output1_fasta, 'w') as fasta_out:
            for header, seq, _, _ in filtered_sequences:
                fasta_out.write(f">{header}\n{seq}\n")
                
        # Sort again and read all sequences into a list
        records = list(SeqIO.parse(output1_fasta, "fasta"))

        # Define sorting rules, sort by chromosome number (numeric part) first, then by _L and _R
        def sort_key(record):
            # Assuming the format is chr1_L, chrA01_R or Vu01_L, etc.
            name = record.id

            # Use regular expressions to separate the numeric part and the arm part of the chromosome name
            match = re.match(r"[a-zA-Z]+(\d+)(?:_(L|R))?", name)
            if match:
                chrom_num, arm = match.groups()
                # Convert the chromosome number part to an integer for sorting
                chrom_num = int(chrom_num)
                # If the end is _L, set it to 0, and if it is _R, set it to 1 (L is sorted before R)
                arm = 0 if arm == 'L' else 1
                return (chrom_num, arm)
            else:
                # If the regular expression does not match, the original name is returned by default to avoid abnormal situations
                return (name,)

        # Sort the records
        sorted_records = sorted(records, key=sort_key)

        # Write the sorted sequence to the output file
        output_file = f"new_genome_{length}.fasta"
        with open(output_file, 'w') as out_handle:
            SeqIO.write(sorted_records, out_handle, "fasta")

        output2_fasta = f"old_genome_{length}.fasta"
        extract_telomere2_sequences(genome1, output2_fasta, length, specified_chromosomes, motif, bin)
               

# Function to extract telomere sequences
# Here, the extracted original genome sequences are sorted by chromosome and direction
def extract_telomere2_sequences(genome1, output2_fasta, length, specified_chromosomes, motif, bin):
    telomere_data = []  # Used to save the extracted sequence information (chromosome name, direction, sequence)

    with open("telomere.old.num.info", 'w') as tel2_info:
        for record in SeqIO.parse(genome1, "fasta"):
            chromosome = record.id
            sequence = record.seq

            # Process only the specified chromosomes and directions
            if (chromosome, "L") in specified_chromosomes:
                # Extract the left end sequence of the chromosome
                left_telomere = sequence[:length]
                left_telomeres = calculate_sequence_frequencies(left_telomere, bin, step_size=10)
                write_bed_file("telo_2.bed", left_telomeres, chromosome, "L")
                telomere_data.append((chromosome, "L", left_telomere))  # Save sequence information

                # Calculate the number of telomeres and write information to the file
                total_telomere_count = count_telomere_motifs(left_telomere, motif)
                tel2_info.write(f"{chromosome}_L\t{len(sequence)}\t{total_telomere_count}\n")

            if (chromosome, "R") in specified_chromosomes:
                # Extract the right end sequence of the chromosome
                right_telomere = sequence[-length:]
                right_telomeres = calculate_sequence_frequencies(right_telomere, bin, step_size=10)
                write_bed_file("telo_2.bed", right_telomeres, chromosome, "R")
                telomere_data.append((chromosome, "R", right_telomere))  # Save sequence information

                # Calculate the number of telomeres and write information to the file
                total_telomere_count = count_telomere_motifs(right_telomere, motif)
                tel2_info.write(f"{chromosome}_R\t{len(sequence)}\t{total_telomere_count}\n")

    # Write the extracted sequences to output2_fasta (unsorted)
    with open(output2_fasta, 'w') as out_handle:
        for chromosome, direction, telomere_seq in telomere_data:
            out_handle.write(f">{chromosome}_{direction}\n{telomere_seq}\n")

    # Step 2: Sort the extracted telomere sequences
    sort_fasta_by_chromosome(output2_fasta, length)


# Sorting that applies to all chromosomes
def sort_fasta_by_chromosome(fasta_file, length):
    output_file = f"old_genome_{length}.fasta"
    # Read all sequences into a list
    records = list(SeqIO.parse(fasta_file, "fasta"))

    # Define sorting rules, sort by chromosome number (numeric part) first, then by _L and _R
    def sort_key(record):
        # Assuming the format is chr1_L, chrA01_R or Vu01_L, etc.
        name = record.id

        # Use regular expressions to separate the numeric part and the arm part of the chromosome name
        match = re.match(r"[a-zA-Z]+(\d+)(?:_(L|R))?", name)
        if match:
            chrom_num, arm = match.groups()
            # Convert the chromosome number part to an integer for sorting
            chrom_num = int(chrom_num)
            # If the end is _L, set it to 0, and if it is _R, set it to 1 (L is sorted before R)
            arm = 0 if arm == 'L' else 1
            return (chrom_num, arm)
        else:
            # If the regular expression does not match, the original name is returned by default to avoid abnormal situations
            return (name,)

    # Sort the records
    sorted_records = sorted(records, key=sort_key)

    # Write the sorted sequence to the output file
    with open(output_file, 'w') as out_handle:
        SeqIO.write(sorted_records, out_handle, "fasta")
        

def count_telomere_motifs(sequence, motif):
    """
    Calculate the number of telomeric motifs in a sequence, supporting motifs and their reverse complementary sequences
    :param sequence: sequence
    :param motif: telomeric motif
    :return: number of motifs
    """
    motif_rc = str(Seq(motif).reverse_complement())
    motif_count = sequence.count(motif) + sequence.count(motif_rc)
    return motif_count


def natural_sort_key(s):
    # Helper function to sort chromosome names with numbers naturally
    return [int(text) if text.isdigit() else text for text in re.split(r'(\d+)', s)]


def sort_telomere_info_file(input_file, output_file):
    # Sort the telomere info file based on the chromosome name and direction (L to R)
    with open(input_file, 'r') as infile:
        lines = infile.readlines()

    sorted_data = []
    for line in lines:
        if line.strip():
            parts = line.strip().split()
            if len(parts) == 3:
                chromosome = parts[0]
                chr_length = int(parts[1])
                number = int(parts[2])
                direction = chromosome.split('_')[-1]
                chr_name = chromosome.split('_')[0]
                sorted_data.append((chr_name, direction, chr_length, number))
            else:
                print(f"Warning: Invalid line format in {input_file}: {line.strip()}")

    # Sort by chromosome name first (natural sort), then by direction ('L' comes before 'R')
    sorted_data.sort(key=lambda x: (natural_sort_key(x[0]), x[1]))

    with open(output_file, 'w') as outfile:
        # Write the header line
        outfile.write("chromosome\tchr_length\tnumber\n")

        # Write the sorted data
        for chr_name, direction, chr_length, number in sorted_data:
            outfile.write(f"{chr_name}_{direction}\t{chr_length}\t{number}\n")

    print(f"Sorted data written to {output_file}")
    

# Process the bed file and sort it by the first column, with values ​​from small to large and direction from L to R
def sort_and_modify_bed_file(input_bed, output_bed):
    with open(input_bed, 'r') as bed_file:
        lines = bed_file.readlines()

    # Define a custom sorting key
    def sorting_key(line):
        chr_name, direction = line.split('\t')[0].rsplit('_', 1)
        chr_number = int(re.sub(r'\D', '', chr_name))  # Extract the numeric part of the chromosome
        return chr_number, direction

    # Sort the lines based on the chromosome number and direction (L, R)
    sorted_lines = sorted(lines, key=sorting_key)

    # Add a header to the output file
    header = "seqnames\tstart\tend\tscore\ttype\n"

    with open(output_bed, 'w') as sorted_bed_file:
        sorted_bed_file.write(header)  # Write the header
        for line in sorted_lines:
            sorted_bed_file.write(line.strip() + "\ttelomere_repeat\n")  # Add the "telomere_repeat" in the type column


# Use GenomeSyn to draw a collinearity diagram
def run_GenomeSyn1(old_extr_genome, new_extr_genome, telo_1_bed_sorted, telo_2_bed_sorted, tel_max, genome_color1, genome_color2, synteny_color, telo_color, length):
    GenomeSn1_command = [f"GenomeSyn -g1 {old_extr_genome} -g2 {new_extr_genome} -snp1 {telo_2_bed_sorted} -snp2 {telo_1_bed_sorted} -snp_max {tel_max} -c1 \"{genome_color1}\" -c2 \"{genome_color2}\" -c4 \"{synteny_color}\" -c11 \"{telo_color}\""]
    try:
        print(f"Running command: {GenomeSn1_command}")
        result = subprocess.run(GenomeSn1_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Define source and target directories
        source_dir = './'
        output_dir = 'genomeSyn_result'
        os.makedirs(output_dir, exist_ok=True)

        # List of file extensions to move
        file_extensions = ['.bed', 'delta', '.delta.filter', '.delta.filter.coords', '.coords', f'_{length}.fasta', '.svg', '.pdf', '.log']

        # Move files that match the extensions
        for filename in os.listdir(source_dir):
            if any(filename.endswith(ext) for ext in file_extensions):
                shutil.move(os.path.join(source_dir, filename), os.path.join(output_dir, filename))

        print(f"All result files have been moved to {output_dir}")
        
        if result.returncode != 0:
            raise subprocess.CalledProcessError(result.returncode, GenomeSn1_command, result.stdout, result.stderr)
        else:
            print(f"GenomeSn1_command executed successfully!")
            print(result.stdout.decode())
    except subprocess.CalledProcessError as e:
        print(f"Error running {GenomeSn1_command}: {e}")
        print(f"Stdout: {e.stdout}")
        print(f"Stderr: {e.stderr}")


def run_GenomeSyn2(old_extr_genome, new_extr_genome, telo_1_bed_sorted, telo_2_bed_sorted, tel_max, genome_color1, genome_color2, synteny_color, telo_color, length):
    GenomeSn2_command = f"GenomeSyn -g1 {old_extr_genome} -g2 {new_extr_genome} -snp1 {telo_2_bed_sorted} -snp2 {telo_1_bed_sorted} -snp_max {tel_max} -c1 \"{genome_color1}\" -c2 \"{genome_color2}\" -c4 \"{synteny_color}\" -c11 \"{telo_color}\""
    try:
        print(f"Running command: {GenomeSn2_command}")
        result = subprocess.run(GenomeSn2_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Define source and target directories
        source_dir = './'
        output_dir = 'genomeSyn_result'
        os.makedirs(output_dir, exist_ok=True)

        # List of file extensions to move
        file_extensions = ['.bed', 'delta', '.delta.filter', '.delta.filter.coords', '.coords', f'_{length}.fasta', '.svg', '.pdf', '.log']

        # Move files that match the extensions
        for filename in os.listdir(source_dir):
            if any(filename.endswith(ext) for ext in file_extensions):
                shutil.move(os.path.join(source_dir, filename), os.path.join(output_dir, filename))

        print(f"All result files have been moved to {output_dir}")

        if result.returncode != 0:
            raise subprocess.CalledProcessError(result.returncode, GenomeSn1_command, result.stdout, result.stderr)
        else:
            print(f"GenomeSn1_command executed successfully!")
            print(result.stdout.decode())
    except subprocess.CalledProcessError as e:
        print(f"Error running {GenomeSn2_command}: {e}")
        print(f"Stdout: {e.stdout.decode()}")
        print(f"Stderr: {e.stderr.decode()}")


# Define the required parameters for input
def main():
    parser = argparse.ArgumentParser(
        prog='TeloComp',
        usage='%(prog)s [options]',
        exit_on_error=False,
        description='A tool for telomere extraction and genome polishing.',
        epilog='Text at the bottom of help'
    )

    parser.add_argument('-G1', '--genome1', metavar='', required=True, help='Input the original genome file (FASTA format)')
    parser.add_argument('-G2', '--genome2', metavar='', required=True, help='Input new genome file (FASTA format)')
    parser.add_argument('-pos', metavar='', required=True, help='The position information of telomere installed on the original genome')
    parser.add_argument('-m', '--motif', metavar='', required=True,
                        help='Telomeric repeats sequences, e.g., plant: CCCTAAA(TTTAGGG), animal: TTAGGG(CCCTAA), etc.')
    parser.add_argument('--bin', metavar='', type=int, default=100,
                        help='The bins for calculating the number of telomeres are divided into. The bin size can be determined by the user, and the default value is 100.')
    parser.add_argument('--genomeSyn1', action='store_true', help='GenomeSyn1 is in the order of original genome first and new genome second')
    parser.add_argument('--genomeSyn2', action='store_true', help='genomeSyn2 is in the order of original genome second and new genome first')
    parser.add_argument('--length', metavar='', help='Extract genome length')
    parser.add_argument('-tel_max', metavar='', help='The maximum value of telomere frequency within the calculated length')
    parser.add_argument("-C1","--genome_color1",  metavar='', type=str, default="#04686b",
                        help='Set the color parameters of genome 1, the default is #04686b, or select other hexadecimal or RGB color codes, such as "#04686b"/"rgb(4,104,107)"')
    parser.add_argument("-C2","--genome_color2",  metavar='', type=str, default="#87c9c3",
                        help='Set the color parameters of genome 2, the default is #87c9c3, or select other hexadecimal or RGB color codes, such as "#87c9c3"/"rgb(135,201,195)"')
    parser.add_argument("-C3", "--synteny_color",  metavar='', type=str, default="#DFDFE1",
                        help='Set the color parameters of the synteny blocks, the default is #DFDFE1, or select other hexadecimal or RGB color codes, such as "#DFDFE1"/"rgb(223,223,225)"')
    parser.add_argument("-C4","--telo_color",  metavar='', type=str, default="#fcaf7c",
                        help='Set the telomere color parameters for displaying genomes, the default is #fcaf7c, or select other hexadecimal or RGB color codes, such as "#fcaf7c"/"rgb(252,175,124)"')
    parser.add_argument('-t', '--threads', type=int, default=20, metavar='',
                        help='Number of threads to use (default: 20)')
    args = parser.parse_args()

    # Convert all paths to absolute paths
    args.genome1 = os.path.abspath(args.genome1)
    args.genome2 = os.path.abspath(args.genome2)
    args.pos = os.path.abspath(args.pos)
    args.length = int(args.length)  # Convert length to integer
    args.tel_max = int(args.tel_max)  # Convert tel_max to integer

    # Process the files needed for drawing collinearity
    log_file = f"genomeSyn.log"
    extract_telomere1_sequences(args.genome1, args.genome2, args.pos, args.length, args.motif, args.bin, log_file)

    # Sort and modify the BED files before running GenomeSyn
    telo_1_bed_sorted = "telo_1_sorted.bed"
    telo_2_bed_sorted = "telo_2_sorted.bed"
    telo_1_bed = "telo_1.bed"
    telo_2_bed = "telo_2.bed"


    sort_and_modify_bed_file(telo_1_bed, telo_1_bed_sorted)
    sort_and_modify_bed_file(telo_2_bed, telo_2_bed_sorted)

    # Sort the telomere1.num.info file with header
    sort_telomere_info_file("telomere.new.num.info", "telomere.complement.num.info")
    # Sort the telomere2.num.info file with header
    sort_telomere_info_file("telomere.old.num.info", "telomere.original.num.info")
    # Remove the file
    os.remove("telomere.old.num.info")
    os.remove("telomere.new.num.info")

    # Draw collinearity plots
    if args.genomeSyn1:
        old_extr_genome = f"old_genome_{args.length}.fasta"
        new_extr_genome = f"new_genome_{args.length}.fasta"
        run_GenomeSyn1(old_extr_genome, new_extr_genome, telo_1_bed_sorted, telo_2_bed_sorted, args.tel_max, args.genome_color1, args.genome_color2, args.synteny_color, args.telo_color, args.length)

    elif args.genomeSyn2:
        old_extr_genome = f"old_genome_{args.length}.fasta"
        new_extr_genome = f"new_genome_{args.length}.fasta"
        run_GenomeSyn2(old_extr_genome, new_extr_genome, telo_1_bed_sorted, telo_2_bed_sorted, args.tel_max, args.genome_color1, args.genome_color2, args.synteny_color, args.telo_color, args.length)




if __name__ == "__main__":
    main()
    