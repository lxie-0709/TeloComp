#!/usr/bin/env python


import os
import sys
import re
import shutil
import subprocess
import argparse
from pathlib import Path
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from threading import Lock
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio.SeqUtils import nt_search

print_lock = Lock()

def align_fasta_to_genome(tmp_sam_dir, genome_path, trim_L, trim_R, tmp_getpos_dir, lgsreads, wgs1, wgs2, threads, NextPolish):
    # Create an output directory for the next step of alignment. 
    # By comparing, we know the position that needs to be cut, so as to obtain the sequence of the telomere part
    tmp_pysam_dir = "tmp_pysam" 
    os.makedirs(tmp_pysam_dir, exist_ok=True)
    
    # The directory where the assembled and polished contig files are located
    contig_files = sorted(os.listdir(tmp_sam_dir)) 
    for contig_file in contig_files:
        if not contig_file.endswith('.fasta'):
            continue
        chromosome = contig_file.split('_')[0]
        direction_sam = contig_file.split('_')[-1].split('.')[0]
        input_file = os.path.join(tmp_sam_dir, contig_file)

        # Define the output SAM file path
        sam_output = Path(tmp_pysam_dir) / f"{chromosome}_asm_{direction_sam}.sam"

        # Run minimap2 for comparison
        command_sam = f"minimap2 -ax map-pb -t {threads} {genome_path} {input_file} > {sam_output}"
        subprocess.run(command_sam, shell=True, check=True)

        with print_lock:
            print(f"Aligned {input_file} to {genome_path} and saved to {sam_output}")


    pysam_get_tel_position(tmp_pysam_dir, tmp_getpos_dir, genome_path, trim_L, trim_R, lgsreads, wgs1, wgs2, threads, NextPolish)   


def pysam_get_tel_position(tmp_pysam_dir, tmp_getpos_dir, genome_path, trim_L, trim_R, lgsreads, wgs1, wgs2, threads, NextPolish):
    # Get the chromosome length and sequence of the reference genome
    chr_sequences = {}
    chr_lengths = {}
    for record in SeqIO.parse(genome_path, "fasta"):
        chr_sequences[record.id] = record.seq
        chr_lengths[record.id] = len(record.seq)

    #Store the chromosome numbers that do not meet the conditions
    invalid_contigs = []

    # Compare the sam file to get the location information
    samfiles = os.listdir(tmp_pysam_dir)
    for samfile in samfiles:
        chromosome = os.path.basename(samfile).split('_')[0]
        direction_sam = os.path.basename(samfile).split('_')[-1].split('.')[0]
        samfile_path = os.path.join(tmp_pysam_dir, samfile)
        valid_contig_found = False  # Track whether a contig that meets the criteria is found

        with pysam.AlignmentFile(samfile_path, "r") as samfile_obj:
            # Traverse each read in the SAM file
            for read in samfile_obj:
                if not read.is_unmapped and read.reference_name == chromosome:
                    contig_start = read.reference_start  # Starting position on the reference genome
                    contig_end = read.reference_end  # End position on the reference genome
                    contig_name = read.query_name  # Name of the contig
                    contig_seq_start = read.query_alignment_start  # Starting position on the contig
                    contig_seq_end = read.query_alignment_end  # End position on the contig
                    reference_length = chr_lengths[chromosome]

                    print(f"Processing {contig_name} on {chromosome} direction: {direction_sam}")
                    print(
                            f"contig_start: {contig_start}, contig_end: {contig_end}, reference_length: {reference_length}")

                    # Confirm the reference genome endpoint corresponding to the contig
                    if (direction_sam == 'L' and contig_start == 0) or (
                            direction_sam == 'R' and contig_end == reference_length):
                        # Find contigs that meet the criteria
                        valid_contig_found = True
                        telomere_seq = None

                        # Check if read.query_sequence exists
                        if not read.query_sequence:
                            print(f"Read query sequence is missing for contig {contig_name}.")
                            invalid_contigs.append((chromosome, direction_sam))
                            continue

                        if direction_sam == 'L':
                            if contig_seq_start == 0:
                                # Align to the left end of contig
                                telomere_seq = Seq(read.query_sequence[contig_seq_end:]).reverse_complement()
                            else:
                                # Align to the right end of contig
                                telomere_seq = Seq(read.query_sequence[:contig_seq_start])
                        elif direction_sam == 'R':
                            if contig_seq_end != len(read.query_sequence):
                                # Align to the left end of contig
                                telomere_seq = Seq(read.query_sequence[contig_seq_end:])
                            else:
                                # Align to the right end of contig
                                telomere_seq = Seq(read.query_sequence[:contig_seq_start]).reverse_complement()

                        contig_position = (f"Contig {contig_name}: {contig_seq_start}-{contig_seq_end} aligns to "
                                           f"Reference {chromosome}: {contig_start}-{contig_end}")
                        print(contig_position)

                        # Check if telomere_seq exists and is not empty
                        if telomere_seq and len(telomere_seq) > 0:
                            print(f"Telomere sequence: {telomere_seq}")
                            # Output telomere_seq to fasta file
                            output_file = os.path.join(tmp_getpos_dir, f"{chromosome}_pysam_np_{direction_sam}.fasta")
                            seq_record = SeqRecord(telomere_seq, id=f"{contig_name}_telomere", description="")
                            SeqIO.write(seq_record, output_file, "fasta")
                        else:
                            print(f"Telomere sequence is missing or empty for contig {contig_name}.")
                            invalid_contigs.append((chromosome, direction_sam))
                            continue  # Skip the current read and continue processing the next read

                        # End processing of the current sam file
                        break

            if not valid_contig_found:
                print(
                    f"Contig does not meet the criteria for {chromosome}_{direction_sam} direction. Run to extract the shortest sequence")
                invalid_contigs.append((chromosome, direction_sam))

    # Process chromosome numbers that do not meet the conditions
    if invalid_contigs:
        # Process contigs that do not meet the conditions (subsequent processing functions can be called here)
        get_shortest_DIR = ['trim_L', 'trim_R']
        get_shortest_reads(invalid_contigs, get_shortest_DIR, lgsreads, wgs1, wgs2, threads, NextPolish)


def get_shortest_reads(invalid_contigs, get_shortest_DIR, lgsreads, wgs1, wgs2, threads, NextPolish):
    output_dir = "tmp_getpos"
    os.makedirs(output_dir, exist_ok=True)

    # Make sure all directories exist
    for directory in get_shortest_DIR:
        if not os.path.exists(directory):
            print(f"Warning: Directory {directory} does not exist. Skipping...")
            continue

    if len(invalid_contigs) == 1:
        chromosome, direction_sam = invalid_contigs[0]
        output_file = os.path.join(output_dir, f"{chromosome}_shortest_np_{direction_sam}.fasta")
        for directory in get_shortest_DIR:
            fasta_files = [f for f in os.listdir(directory) if f.endswith('.fasta')]
            for fasta_file in fasta_files:
                chr_name = os.path.basename(fasta_file).split('_')[0]
                direction_trim = os.path.basename(fasta_file).split('_')[-1].split('.')[0]
                if chr_name == chromosome and direction_sam == direction_trim:
                    file_path = os.path.join(directory, fasta_file)
                    with open(file_path, 'r') as infile:
                        lengths = {'lowercase': {}}
                        for record in SeqIO.parse(infile, 'fasta'):
                            sequence = str(record.seq)
                            lowercase_sequence = ''.join(char for char in sequence if char.islower())
                            lengths['lowercase'][record.id] = len(lowercase_sequence)

                        min_record_id = min(lengths['lowercase'], key=lengths['lowercase'].get)
                        with open(output_file, 'w') as out_handle:
                            for record in SeqIO.parse(file_path, 'fasta'):
                                if record.id == min_record_id:
                                    # Extract lowercase portions of the sequence
                                    lowercase_seq = ''.join([char for char in str(record.seq) if char.islower()])
                                    out_handle.write(f'>{record.id}\n{lowercase_seq}\n')

        asm_contig_racon(output_file, lgsreads, wgs1, wgs2, threads, NextPolish, split_fasta=False)

    else:
        output_merged_file = os.path.join(output_dir, "merged_chr_shortest.fasta")
        unique_sequence = set()
        with open(output_merged_file, 'w') as outfile:
            for chromosome, direction_sam in invalid_contigs:
                for directory in get_shortest_DIR:
                    fasta_files = [f for f in os.listdir(directory) if f.endswith('.fasta')]
                    for fasta_file in fasta_files:
                        chr_name = os.path.basename(fasta_file).split('_')[0]
                        direction_trim = os.path.basename(fasta_file).split('_')[-1].split('.')[0]
                        if chr_name == chromosome and direction_sam == direction_trim:
                            file_path = os.path.join(directory, fasta_file)
                            with open(file_path, 'r') as infile:
                                for record in SeqIO.parse(infile, 'fasta'):
                                    seq_id = os.path.splitext(os.path.basename(file_path))[0]
                                    lowercase_seq = ''.join(char for char in str(record.seq) if char.islower())
                                    if lowercase_seq and seq_id not in unique_sequence:
                                        unique_sequence.add(seq_id)
                                        outfile.write(f'>{seq_id}\n{lowercase_seq}\n')

        asm_contig_racon(output_merged_file, lgsreads, wgs1, wgs2, threads, NextPolish, split_fasta=True)
        # Remove the file
        os.remove(output_merged_file)


# # racon preprocessing
def asm_contig_racon(input_file, lgsreads, wgs1, wgs2, threads, NextPolish, split_fasta):
    tmpdir = "tmp_dir"
    os.makedirs(tmpdir, exist_ok=True)

    while True:  # Use loops instead of recursion
        racon_commands = [
            f"minimap2 -ax map-pb -t 50 {input_file} {lgsreads} | samtools sort - -m 2g --threads 20 -o {tmpdir}/genome.lgs.bam",
            f"samtools index {tmpdir}/genome.lgs.bam",
            f"ls {os.path.abspath(tmpdir)}/genome.lgs.bam > {tmpdir}/pb.map.bam.fofn",
            f"python {NextPolish}/lib/nextpolish2.py -g {input_file} -l {tmpdir}/pb.map.bam.fofn -r clr -p 50 -a -s -o {tmpdir}/merged_genome.lgspolish.fasta"
        ]

        try:
            for command in racon_commands:
                print(f"Running command: {command}")
                result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                if result.returncode != 0:
                    raise subprocess.CalledProcessError(result.returncode, command, result.stdout, result.stderr)
                else:
                    print(f"Command executed successfully: {command}")
                    print(result.stdout.decode())

            break  # If all commands are executed successfully, jump out of the loop

        except subprocess.CalledProcessError as e:
            print(f"Error running racon_commands: {e.stderr.decode()}")
            error_message = e.stderr.decode()

            if "Failed to correct sequence" in error_message:
                # Extract and clean strings using regular expressions
                match = re.search(r'Failed to correct sequence:\s+(\S+)', error_message)
                if match:
                    failed_sequence = match.group(1).strip()
                    # Further clean non-printing characters in the string
                    clean_sequence_id = re.sub(r'[^\w_\-]', '', failed_sequence)

                    print(f"Failed sequence ID: {clean_sequence_id}")

                    # Delete the corresponding sequence from the FASTA file
                    new_input_file = remove_sequence_from_fasta(input_file, clean_sequence_id)
                    if new_input_file == input_file:
                        print(f"Failed sequence {clean_sequence_id} was not found. Exiting.")
                        return
                    else:
                        input_file = new_input_file
                        print(f"Re-running racon with updated FASTA file: {input_file}")
                else:
                    print("Failed to parse the failed sequence ID. Exiting.")
                    return

            else:
                print("Unexpected error occurred. Exiting.")
                return
                

    # split_fasta 
    if split_fasta:
        racon_genome = f"{tmpdir}/merged_genome.lgspolish.fasta"
        output_dir = "tmp_to_NP"
        os.makedirs(output_dir, exist_ok=True)

        for record in SeqIO.parse(racon_genome, "fasta"):
            base_name = record.id
            chr_name = base_name.split('_')[0]
            extr_type = base_name.split('_')[1]
            direction = base_name.split('_')[3]
            output_file = os.path.join(output_dir, f"{chr_name}_{extr_type}_np_{direction}.fasta")
            with open(output_file, 'w') as out_handle:
                SeqIO.write(record, out_handle, "fasta")

        # Remove the file(Delete the newly generated merged file)
        tmpgetpos = "tmp_getpos"
        if not os.path.exists(tmpgetpos):
            os.makedirs(tmpgetpos, exist_ok=True)
        new_input_file = f"{tmpgetpos}/new_merged_chr_shortest.fasta"
        # Check if the file exists
        if os.path.exists(new_input_file):
            # If it exists, delete the file
            os.remove(new_input_file)
        else:
            print("Execute the nextPolish polishing command...")

        asm_contig_polish(output_dir, lgsreads, wgs1, wgs2, threads, NextPolish)

        # Delete temporary directories and files
        shutil.rmtree(tmpdir)
        shutil.rmtree(output_dir)

    else:
        asm_contig_polish(tmpdir, lgsreads, wgs1, wgs2, threads, NextPolish)
        shutil.rmtree(tmpdir)
        shutil.rmtree(output_dir)


def remove_sequence_from_fasta(fasta_file, sequence_id):
    output_Dir = "tmp_getpos"
    os.makedirs(output_Dir, exist_ok=True)

    output_new_merged_file = os.path.join(output_Dir, "new_merged_chr_shortest.fasta")
    sequence_removed = False

    # Use regular expressions to remove control characters in sequence IDs
    clean_sequence_id = re.sub(r'[^\w_\-]', '', sequence_id.strip())

    print(f"Processed sequence ID: {clean_sequence_id}")

    with open(fasta_file, "r") as input_handle, open(output_new_merged_file, "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            record_id = re.sub(r'[^\w_\-]', '', record.id.strip())

            # Print all recorded IDs
            print(f"Record ID in FASTA file: '{record_id}'")

            # Check for partial matches (alternatively you can use startswith or in methods
            if clean_sequence_id in record_id or record_id in clean_sequence_id:
                # Save failed sequences to a separate file
                failed_seq_file = os.path.join(output_Dir, f"{record_id}.fasta")
                with open(failed_seq_file, "w") as failed_handle:
                    SeqIO.write(record, failed_handle, "fasta")
                print(f"Saved failed sequence {record_id} to {failed_seq_file}")
                sequence_removed = True
            else:
                # Write the remaining sequences into a new FASTA file
                SeqIO.write(record, output_handle, "fasta")

    if sequence_removed:
        print(f"A new merged FASTA file has been created: {output_new_merged_file}")
        return output_new_merged_file
    else:
        print(f"Sequence ID '{clean_sequence_id}' not found in {fasta_file}. Returning original file.")
        return fasta_file


def asm_contig_polish(output_dir, lgsreads, wgs1, wgs2, threads, NextPolish):
    fasta_files = [f for f in os.listdir(output_dir) if f.endswith('.fasta')]
    for fasta_file in fasta_files:
        reads_name = os.path.splitext(fasta_file)[0]
        file_path = os.path.join(output_dir, fasta_file)
        genome = file_path
        FoFn_dir = os.path.join("Dir_fofn", reads_name)
        os.makedirs(FoFn_dir, exist_ok=True)

        # Only use the filename of the input file
        genome_final = os.path.basename(genome)
        shutil.copyfile(genome, os.path.join(FoFn_dir, genome_final))

        sgs_fofn = os.path.join(FoFn_dir, "sgs.fofn")
        lgs_fofn = os.path.join(FoFn_dir, "lgs.fofn")
        nextPolish = os.path.join(NextPolish, "nextPolish")
        t = threads
        outputname = reads_name

        # Get the current working directory
        current_dir = os.getcwd()
        work_dir = os.path.join(current_dir, f"Dir_{outputname}")
        os.makedirs(work_dir, exist_ok=True)

        final_dir = os.path.join(current_dir, "tmp_getpos")
        os.makedirs(final_dir, exist_ok=True)

        with open(sgs_fofn, 'w') as sgs_file:
            sgs_file.write(f"{wgs1}\n")
            sgs_file.write(f"{wgs2}\n")

        with open(lgs_fofn, 'w') as lgs_file:
            lgs_file.write(f"{lgsreads}\n")

        run_cgf_content = f"""
        [General]
        job_type = local
        job_prefix = nextPolish
        task = best
        rewrite = no
        rerun = 3
        parallel_jobs = 1
        multithread_jobs = {t}
        genome = {genome_final}
        genome_size = auto
        workdir = {work_dir}
        polish_options = -p {t}

        [sgs_option]
        sgs_fofn = sgs.fofn
        sgs_options = -max_depth 100 -bwa

        [lgs_option]
        lgs_fofn = lgs.fofn
        lgs_options = -min_read_len 5k -max_depth 100
        lgs_minimap2_options = -x map-pb
        """

        config_file_path = os.path.join(FoFn_dir, "run.cgf")
        with open(config_file_path, 'w') as cfg_file:
            cfg_file.write(run_cgf_content)

        print(f"Running nextpolish to polish assembled and extracted contigs for {reads_name}")

        polish_command = f"{nextPolish} {config_file_path}"
        try:
            subprocess.run(polish_command, shell=True, check=True, capture_output=True, text=True)
            np_fastas = [f for f in os.listdir(work_dir) if f.endswith('.fasta')]
            for fasta in np_fastas:
                # Rename to {reads_name}.fasta
                new_name = f"{reads_name}.fasta"
                new_path = os.path.join(final_dir, new_name)
                shutil.move(os.path.join(work_dir, fasta), new_path)
                #shutil.move(os.path.join(work_dir, fasta), final_dir)
        except subprocess.CalledProcessError as e:
            print(f"Error running nextPolish command for {reads_name}: {e}")
            print(f"Stdout: {e.stdout}")
            print(f"Stderr: {e.stderr}")
            print(
            f"{reads_name} quality is too low, polishing is unsuccessful, and the sequence is directly output to tmp_getpos")
            # If the command fails, copy the file to the files_NP directory
            shutil.copy2(file_path, final_dir)
        finally:
            # Delete work_dir and any backup files
            shutil.rmtree(work_dir)
            backup_files = [f for f in os.listdir(current_dir) if f.startswith(f"Dir_{outputname}.backup")]
            for backup in backup_files:
                backup_path = os.path.join(current_dir, backup)
                shutil.rmtree(backup_path)
            print(f"Deleted work directory {work_dir} and backup files")


# Filter the fasta files of motifs without corresponding motifs, 
# and set the threshold to load into the genome only when the total length of telomeres (number of telomeres * 7, 6) accounts for more than 10% of the total length
def post_tel_to_genome(tmp_getpos_dir, genome_path, motif):

    output_genome_path = "new_genome.fasta"  # New genome sequence output path
    output_positions_path = "telomere_positions.txt"  # Telomere location information output path

    # Calculate the reverse complementary sequence of the motif
    motif_rc = str(Seq(motif).reverse_complement())

    # Read the original genome sequence
    genome_sequences = {record.id: record for record in SeqIO.parse(genome_path, "fasta")}

    # Used to record telomere sequence information
    telomere_positions = ["chromosome\tstart\tend\tlength"]

    # Traverse the telomere sequence files in tmp_getpos_dir
    for telomere_file in os.listdir(tmp_getpos_dir):
        if telomere_file.endswith(".fasta"):
            telomere_path = os.path.join(tmp_getpos_dir, telomere_file)
            telomere_records = list(SeqIO.parse(telomere_path, "fasta"))

            if not telomere_records:
                continue

            # Get chromosome name and direction information
            chrom_name = telomere_file.split('_')[0]
            direction = os.path.basename(telomere_file).split('_')[-1].split('.')[0]

            # Get the telomere sequence and convert it to uppercase
            telomere_seq = telomere_records[0].seq.upper()

            # Count the number of times the motif appears
            motif_count = telomere_seq.count(motif) + telomere_seq.count(motif_rc)

            # Set the total length of motif to motif_count * len(motif) and calculate the proportion
            motif_total_length = motif_count * len(motif)
            telomere_total_length = len(telomere_seq)
            motif_percentage = (motif_total_length / telomere_total_length) * 100

            # If the total length of the motif is less than 10%, skip the sequence
            if motif_percentage < 10:
                continue

            # Determine the corresponding chromosome in the reference genome
            if chrom_name in genome_sequences:
                chrom_seq = genome_sequences[chrom_name].seq
                telomere_seq = telomere_records[0].seq.upper()  # Convert telomere sequences to uppercase

                # Initialize the new_seq variable
                new_seq = chrom_seq
                start_pos = None
                end_pos = None

                # Sequence integration according to direction
                if direction == "L":
                    new_seq = telomere_seq + chrom_seq
                    start_pos = 0
                    end_pos = len(telomere_seq) - 1
                elif direction == "R":
                    new_seq = chrom_seq + telomere_seq
                    start_pos = len(chrom_seq)
                    end_pos = start_pos + len(telomere_seq) - 1

                # Update chromosome sequence
                genome_sequences[chrom_name].seq = new_seq

                # Record telomere information
                telomere_positions.append(f"{chrom_name}_{direction}\t{start_pos}\t{end_pos}\t{len(telomere_seq)}")

    # Output new genome sequence
    with open(output_genome_path, "w") as output_handle:
        SeqIO.write(genome_sequences.values(), output_handle, "fasta")

    # Output telomere position information
    with open(output_positions_path, "w") as pos_handle:
        pos_handle.write("\n".join(telomere_positions))

    print(f"The new genome has been saved to {output_genome_path}")
    print(f"Telomere position information has been saved to {output_positions_path}")


def determine_window_size(length):

    # Determine the extraction window size based on the given length
    if length <= 10000:
        return 10000
    elif length <= 20000:
        return 20000
    else:
        return min(200000, 10000 * ((length // 10000) + 1))


# Extract the corresponding chromosome {window_size} part of the new genome
def extract_sequences(output_genome_path, output_positions_path):
    output_dir = "extract_teloSeq"
    os.makedirs(output_dir, exist_ok=True)

    genome_records = {record.id: record for record in SeqIO.parse(output_genome_path, "fasta")}

    with open(output_positions_path, "r") as pos_file:
        lines = pos_file.readlines()

    header = lines[0].strip().split("\t") # Read the header row
    data_lines = lines[1:] # Skip the header row

    for line in data_lines: # Use data_lines to iterate over the data lines
        chromosome, start_pos, end_pos, telomere_length = line.strip().split('\t')

        try:
            start_pos = int(start_pos)
            end_pos = int(end_pos)
            telomere_length = int(telomere_length)
        except ValueError as e:
            print(f"Skipping this line due to value error: {line}")
            continue

        chrom_name, direction = chromosome.rsplit('_', 1)
        window_size = determine_window_size(telomere_length)

        if chrom_name in genome_records:
            chrom_seq = genome_records[chrom_name].seq

            if direction == "L":
                extract_seq = chrom_seq[:window_size]
            elif direction == "R":
                extract_seq = chrom_seq[-window_size:]

            out_file = os.path.join(output_dir, f"{chrom_name}_{direction}_{window_size}BP.fasta")
            with open(out_file, "w") as out_handle:
                record = SeqIO.SeqRecord(extract_seq, id=f"{chrom_name}_{direction}_{window_size}BP", description="")
                SeqIO.write(record, out_handle, "fasta")

            print(f"Extract and save {window_size}BP sequences of {chrom_name}_{direction} to {out_file}")


# Calculate the type and sequence of telomere motif 7 in the telomere part
def window7_cut(sequence, windows_size, step):
    seq_length = len(sequence)
    cut_num = int((seq_length - windows_size) / step)
    rep = []
    for i in range(cut_num + 1):
        start_index = i * step
        end_index = start_index + windows_size
        rep_seq = sequence[start_index:end_index]
        rep.append(rep_seq)
    return rep


def top_rep7_seq(sequence, windows_size , step, top):
    rep_seqs = window7_cut(sequence, windows_size, step)
    rep_counts = {}
    for j in rep_seqs:
        if j in rep_counts:
            rep_counts[j] += 1
        else:
            rep_counts[j] = 1
    top_rep_seq = sorted(rep_counts.items(), key=lambda item: item[1], reverse=True)[:top]
    return top_rep_seq


# Calculate the type and sequence of telomere motif 6 in the telomere part
def window6_cut(sequence, windows_size, step):
    seq_length = len(sequence)
    cut_num = int((seq_length - windows_size) / step)
    rep = []
    for i in range(cut_num + 1):
        start_index = i * step
        end_index = start_index + windows_size
        rep_seq = sequence[start_index:end_index]
        rep.append(rep_seq)
    return rep


def top_rep6_seq(sequence, windows_size , step, top):
    rep_seqs = window7_cut(sequence, windows_size, step)
    rep_counts = {}
    for j in rep_seqs:
        if j in rep_counts:
            rep_counts[j] += 1
        else:
            rep_counts[j] = 1
    top_rep_seq = sorted(rep_counts.items(), key=lambda item: item[1], reverse=True)[:top]
    return top_rep_seq


# Calculation type previous step processing
# Calculate the 7-base sequence
def process_telomere7_sequences(output_dir, output_txt_path, window_size=7, step=1, top=100):
    with open(output_txt_path, 'w') as txt_handle:
        txt_handle.write("chr_name\tsequence\ttelomere_count\tstart_pos\tend_pos\n")  # Write the header row
        for fasta_file in os.listdir(output_dir):
            if fasta_file.endswith(".fasta"):
                chromosome = fasta_file.split('_')[0]
                direction = fasta_file.split('_')[-1].split('.')[0]
                fasta_path = os.path.join(output_dir, fasta_file)
                for record in SeqIO.parse(fasta_path, "fasta"):
                    sequence = str(record.seq).upper()  # Convert the sequence to uppercase
                    top_rep7_seqs = top_rep7_seq(sequence, window_size, step, top)
                    for seq, count in top_rep7_seqs:
                        start_pos = sequence.find(seq)
                        end_pos = start_pos + len(seq) - 1
                        if start_pos != -1:  # Make sure the sequence is found
                            txt_handle.write(f"{chromosome}_{direction}\t{seq}\t{count}\t{start_pos}\t{end_pos}\n")


# Calculate the 6-base motif
def process_telomere6_sequences(output_dir, output_txt_path, window_size=6, step=1, top=100):
    with open(output_txt_path, 'w') as txt_handle:
        txt_handle.write("chr_name\tsequence\ttelomere_count\tstart_pos\tend_pos\n")  # Write the header row
        for fasta_file in os.listdir(output_dir):
            if fasta_file.endswith(".fasta"):
                chromosome = fasta_file.split('_')[0]
                direction = fasta_file.split('_')[-1].split('.')[0]
                fasta_path = os.path.join(output_dir, fasta_file)
                for record in SeqIO.parse(fasta_path, "fasta"):
                    sequence = str(record.seq).upper()  # Convert the sequence to uppercase
                    top_rep6_seqs = top_rep6_seq(sequence, window_size, step, top)
                    for seq, count in top_rep6_seqs:
                        start_pos = sequence.find(seq)
                        end_pos = start_pos + len(seq) - 1
                        if start_pos != -1:  # Make sure the sequence is found
                            txt_handle.write(f"{chromosome}_{direction}\t{seq}\t{count}\t{start_pos}\t{end_pos}\n")


# Use python to draw
def count_pattern(sequence, pattern):
    # Count the number of times a pattern appears in a sequence
    return sum(1 for _ in nt_search(sequence, pattern)[1:])


# Processing 7 bases drawing
def plot_telomere7_data(KB_dir, output_positions_path, output_repeats_path, output_plot_dir):
    if not os.path.exists(output_plot_dir):
        os.makedirs(output_plot_dir)

    # Read telomere location information file
    telomere_positions = {}
    telomere_lengths = {}
    try:
        with open(output_positions_path, "r") as pos_file:
            lines = pos_file.readlines()
            header = lines[0].strip().split('\t')
            data_lines = lines[1:]  # Skip the header row
            for line in data_lines:
                if line.strip():
                    chromosome, start_pos, end_pos, telomere_length = line.strip().split('\t')
                    telomere_positions[chromosome] = int(end_pos)
                    telomere_lengths[chromosome] = int(telomere_length)
                    print(f"Loaded telomere position: {chromosome} -> {end_pos}")  # Debug information调试信息
    except FileNotFoundError:
        print(f"File not found: {output_positions_path}")
        return
    except ValueError as e:
        print(f"Error reading {output_positions_path}: {e}")
        return

    # Read telomere repeat sequence information file
    telomere_repeats = {}
    try:
        with open(output_repeats_path, "r") as repeat_file:
            lines = repeat_file.readlines()
            header = lines[0].strip().split('\t')
            data_lines = lines[1:]  # Skip the header row
            for line in data_lines:
                if line.strip():
                    chrom_name, seq, count, start_pos, end_pos = line.strip().split('\t')
                    key = f"{chrom_name}"
                    if key not in telomere_repeats:
                        telomere_repeats[key] = []
                    telomere_repeats[key].append((seq, int(count)))
                    print(f"Loaded repeat: {chrom_name} -> {seq}, {count}")  # Debug information
    except FileNotFoundError:
        print(f"File not found: {output_repeats_path}")
        return
    except ValueError as e:
        print(f"Error reading {output_repeats_path}: {e}")
        return

    blacklist = ["CCCTAAA", "TTTAGGG", "CCTAAA", "TTTAGG"]  # Filter out only motifs in the blacklist

    for fasta_file in os.listdir(KB_dir):
        if fasta_file.endswith(".fasta"):
            fasta_path = os.path.join(KB_dir, fasta_file)
            file_name = os.path.splitext(fasta_file)[0]
            fasta_name = '_'.join(file_name.split('_')[:-1])
            key = f"{fasta_name}"

            if key in telomere_positions and key in telomere_repeats:
                end_pos = telomere_positions[key]
                telomere_length = telomere_lengths[key]

                # Use CCCTAAA and TTTAGGG as parameters, except for the motifs in the blacklist, there is only one motif
                non_blacklist_repeats = [repeat for repeat in telomere_repeats[key] if repeat[0] not in blacklist]

                if len(non_blacklist_repeats) == 1:
                    filtered_repeats = [('CCCTAAA', 1), ('TTTAGGG', 1)]
                else:
                    filtered_repeats = [repeat for repeat in telomere_repeats[key] if repeat[0] not in blacklist]

                if len(filtered_repeats) < 2:
                    print(f"Not enough valid repeats for {key}, skipping...")
                    continue

                top_repeats = filtered_repeats[:2]
                param1, param2 = top_repeats[0][0], top_repeats[1][0]
                reverse_param1 = str(Seq(param1).reverse_complement())
                reverse_param2 = str(Seq(param2).reverse_complement())

                print(
                    f"Parameters for plotting:\nparam1: {param1}\nparam2: {param2}\nreverse_param1: {reverse_param1}\nreverse_param2: {reverse_param2}")

                # Extract and calculate chrom_segment_size
                try:
                    size_str = file_name.split('_')[-1]
                    if size_str.endswith('BP'):
                        chrom_segment_size = int(size_str[:-2])
                        scale_breaks = max(chrom_segment_size // 10, 1)
                    else:
                        raise ValueError("Size does not end with 'BP'")
                except ValueError as e:
                    print(f"Error calculating chrom_segment_size or scale_breaks: {e}")
                    continue

                direction = file_name.split('_')[1]
                if direction == "L":
                    plot_end_pos = end_pos
                elif direction == "R":
                    plot_end_pos = chrom_segment_size - telomere_length
                else:
                    print(f"Invalid direction {direction} in file {file_name}, skipping...")
                    continue

                with open(fasta_path, "r") as file:
                    record = SeqIO.read(file, "fasta")
                    chrom_seq = str(record.seq)

                # Define telomere sequence and its reverse complementary sequence
                telomere_seq = ["CCCTAAA", param1, param2, "TTTAGGG"]
                telomere_rc = [str(Seq(seq).reverse_complement()) for seq in telomere_seq]

                # Define sliding window size and step size
                window_size = 100
                step_size = 10
                positions = np.arange(0, len(chrom_seq) - window_size + 1, step_size)

                telomere_frequencies = {seq: [] for seq in telomere_seq}
                rc_frequencies = {seq: [] for seq in telomere_rc}

                for pos in positions:
                    window_seq = chrom_seq[pos:pos + window_size]
                    for tel_seq in telomere_seq:
                        telomere_frequencies[tel_seq].append(count_pattern(window_seq, tel_seq))
                    for rc_seq in telomere_rc:
                        rc_frequencies[rc_seq].append(count_pattern(window_seq, rc_seq))

                # If lengths don't match, truncate or pad telomere_frequencies to match positions
                for seq in telomere_seq:
                    if len(telomere_frequencies[seq]) != len(positions):
                        print(f"Adjusting length for {seq} to match positions length")
                        if len(telomere_frequencies[seq]) > len(positions):
                            telomere_frequencies[seq] = telomere_frequencies[seq][:len(positions)]
                        else:
                            telomere_frequencies[seq] += [0] * (len(positions) - len(telomere_frequencies[seq]))

                telomere_data = pd.DataFrame({'position': positions})

                # Check data length consistency
                for seq in telomere_seq:
                    if len(telomere_frequencies[seq]) != len(positions):
                        print(f"Length mismatch after adjustment for {seq}, expected {len(positions)} but got {len(telomere_frequencies[seq])}")

                for seq in telomere_seq:
                    telomere_data[f'{seq}_frequency'] = telomere_frequencies[seq]

                print(telomere_data.head())

                # Create a graph
                output_file = os.path.join(output_plot_dir, f"{file_name}_plot.pdf")

                plt.figure(figsize=(8, 6))

                # Fixed color settings
                colors = {
                    "CCCTAAA": "#F6CAE5",
                    "TTTAGGG": "#C4A5DE",
                    param1: "#CFEAF1",
                    param2: "#96CCCB"
                }

                for seq in telomere_seq:
                    if not telomere_data[f'{seq}_frequency'].empty:
                        plt.plot(telomere_data['position'], telomere_data[f'{seq}_frequency'], label=seq,
                                 color=colors.get(seq, '#000000'), linewidth=2)

                plt.axvline(x=plot_end_pos, linestyle='--', color='#B3B3B3')
                plt.title('Telomere Sequence Frequency along Chromosome')
                plt.xlabel('Position (bp)')
                plt.ylabel('Telomere Frequency')
                plt.legend(loc='upper left', bbox_to_anchor=(1, 1), ncol=1)
                plt.grid(True)

                plt.gca().spines['top'].set_visible(False)
                plt.gca().spines['right'].set_visible(False)
                plt.tick_params(axis='both', which='both', length=0)
                plt.xticks(np.arange(0, chrom_segment_size, scale_breaks))
                # Set the y-axis scale to integers
                max_y_value = max(telomere_data.iloc[:, 1:].max())
                y_ticks = np.arange(0, max_y_value + 1, 1)  # Step size is 1, adjust as needed
                plt.yticks(y_ticks)

                plt.savefig(output_file, bbox_inches='tight')
                plt.close()

                print(f"Successfully plotted the telomere map of {file_name}, output file: {output_file}")


# Processing 6 bases drawing
def plot_telomere6_data(KB_dir, output_positions_path, output_repeats_path, output_plot_dir):
    if not os.path.exists(output_plot_dir):
        os.makedirs(output_plot_dir)

    # Read telomere location information file
    telomere_positions = {}
    telomere_lengths = {}
    try:
        with open(output_positions_path, "r") as pos_file:
            lines = pos_file.readlines()
            header = lines[0].strip().split('\t')
            data_lines = lines[1:]  # Skip the header row
            for line in data_lines:
                if line.strip():
                    chromosome, start_pos, end_pos, telomere_length = line.strip().split('\t')
                    telomere_positions[chromosome] = int(end_pos)
                    telomere_lengths[chromosome] = int(telomere_length)
                    print(f"Loaded telomere position: {chromosome} -> {end_pos}")  # Debug information
    except FileNotFoundError:
        print(f"File not found: {output_positions_path}")
        return
    except ValueError as e:
        print(f"Error reading {output_positions_path}: {e}")
        return

    # Read telomere repeat sequence information file
    telomere_repeats = {}
    try:
        with open(output_repeats_path, "r") as repeat_file:
            lines = repeat_file.readlines()
            header = lines[0].strip().split('\t')
            data_lines = lines[1:]  # Skip the header row
            for line in data_lines:
                if line.strip():
                    chrom_name, seq, count, start_pos, end_pos = line.strip().split('\t')
                    key = f"{chrom_name}"
                    if key not in telomere_repeats:
                        telomere_repeats[key] = []
                    telomere_repeats[key].append((seq, int(count)))
                    print(f"Loaded repeat: {chrom_name} -> {seq}, {count}")  # Debug information
    except FileNotFoundError:
        print(f"File not found: {output_repeats_path}")
        return
    except ValueError as e:
        print(f"Error reading {output_repeats_path}: {e}")
        return

    blacklist = ["CCCTAAA", "TTTAGGG", "CCTAAA", "TTTAGG"]  # Filter out only motifs in the blacklist

    for fasta_file in os.listdir(KB_dir):
        if fasta_file.endswith(".fasta"):
            fasta_path = os.path.join(KB_dir, fasta_file)
            file_name = os.path.splitext(fasta_file)[0]
            fasta_name = '_'.join(file_name.split('_')[:-1])
            key = f"{fasta_name}"

            if key in telomere_positions and key in telomere_repeats:
                end_pos = telomere_positions[key]
                telomere_length = telomere_lengths[key]

                # Use CCCTAAA and TTTAGGG as parameters to handle the case where there is only one motif other than the motifs in the blacklist
                non_blacklist_repeats = [repeat for repeat in telomere_repeats[key] if repeat[0] not in blacklist]

                if len(non_blacklist_repeats) == 1:
                    filtered_repeats = [('CCTAAA', 1), ('TTTAGG', 1)]
                else:
                    filtered_repeats = [repeat for repeat in telomere_repeats[key] if repeat[0] not in blacklist]

                if len(filtered_repeats) < 2:
                    print(f"Not enough valid repeats for {key}, skipping...")
                    continue

                top_repeats = filtered_repeats[:2]
                param1, param2 = top_repeats[0][0], top_repeats[1][0]
                reverse_param1 = str(Seq(param1).reverse_complement())
                reverse_param2 = str(Seq(param2).reverse_complement())

                print(
                    f"Parameters for plotting:\nparam1: {param1}\nparam2: {param2}\nreverse_param1: {reverse_param1}\nreverse_param2: {reverse_param2}")

                # Extract and calculate chrom_segment_size
                try:
                    size_str = file_name.split('_')[-1]
                    if size_str.endswith('BP'):
                        chrom_segment_size = int(size_str[:-2])
                        scale_breaks = max(chrom_segment_size // 10, 1)
                    else:
                        raise ValueError("Size does not end with 'BP'")
                except ValueError as e:
                    print(f"Error calculating chrom_segment_size or scale_breaks: {e}")
                    continue

                direction = file_name.split('_')[1]
                if direction == "L":
                    plot_end_pos = end_pos
                elif direction == "R":
                    plot_end_pos = chrom_segment_size - telomere_length
                else:
                    print(f"Invalid direction {direction} in file {file_name}, skipping...")
                    continue

                with open(fasta_path, "r") as file:
                    record = SeqIO.read(file, "fasta")
                    chrom_seq = str(record.seq)

                # Define telomere sequence and its reverse complementary sequence
                telomere_seq = ["CCTAAA", param1, param2, "TTTAGG"]
                telomere_rc = [str(Seq(seq).reverse_complement()) for seq in telomere_seq]

                # Define sliding window size and step size
                window_size = 100
                step_size = 10
                positions = np.arange(0, len(chrom_seq) - window_size + 1, step_size)

                telomere_frequencies = {seq: [] for seq in telomere_seq}
                rc_frequencies = {seq: [] for seq in telomere_rc}

                for pos in positions:
                    window_seq = chrom_seq[pos:pos + window_size]
                    for tel_seq in telomere_seq:
                        telomere_frequencies[tel_seq].append(count_pattern(window_seq, tel_seq))
                    for rc_seq in telomere_rc:
                        rc_frequencies[rc_seq].append(count_pattern(window_seq, rc_seq))

                # If lengths don't match, truncate or pad telomere_frequencies to match positions
                for seq in telomere_seq:
                    if len(telomere_frequencies[seq]) != len(positions):
                        print(f"Adjusting length for {seq} to match positions length")
                        if len(telomere_frequencies[seq]) > len(positions):
                            telomere_frequencies[seq] = telomere_frequencies[seq][:len(positions)]
                        else:
                            telomere_frequencies[seq] += [0] * (len(positions) - len(telomere_frequencies[seq]))

                telomere_data = pd.DataFrame({'position': positions})

                # Check data length consistency
                for seq in telomere_seq:
                    if len(telomere_frequencies[seq]) != len(positions):
                        print(f"Length mismatch after adjustment for {seq}, expected {len(positions)} but got {len(telomere_frequencies[seq])}")

                for seq in telomere_seq:
                    telomere_data[f'{seq}_frequency'] = telomere_frequencies[seq]

                print(telomere_data.head())

                # Create a graph
                output_file = os.path.join(output_plot_dir, f"{file_name}_plot.pdf")

                plt.figure(figsize=(8, 6))

                # Fixed color settings
                colors = {
                    "CCTAAA": "#F6CAE5",
                    "TTTAGG": "#C4A5DE",
                    param1: "#CFEAF1",
                    param2: "#96CCCB"
                }

                for seq in telomere_seq:
                    if not telomere_data[f'{seq}_frequency'].empty:
                        plt.plot(telomere_data['position'], telomere_data[f'{seq}_frequency'], label=seq,
                                 color=colors.get(seq, '#000000'), linewidth=2)

                plt.axvline(x=plot_end_pos, linestyle='--', color='#B3B3B3')
                plt.title('Telomere Sequence Frequency along Chromosome')
                plt.xlabel('Position (bp)')
                plt.ylabel('Telomere Frequency')
                plt.legend(loc='upper left', bbox_to_anchor=(1, 1), ncol=1)
                plt.grid(True)

                plt.gca().spines['top'].set_visible(False)
                plt.gca().spines['right'].set_visible(False)
                plt.tick_params(axis='both', which='both', length=0)
                plt.xticks(np.arange(0, chrom_segment_size, scale_breaks))
                # Set the y-axis scale to integers轴刻度为整数
                max_y_value = max(telomere_data.iloc[:, 1:].max())
                y_ticks = np.arange(0, max_y_value + 1, 1)  # Step size is 1, adjust as needed
                plt.yticks(y_ticks)

                plt.savefig(output_file, bbox_inches='tight')
                plt.close()

                print(f"Successfully plotted the telomere map of {file_name}, output file: {output_file}")


# Define the required parameters for input
def main():
    parser = argparse.ArgumentParser(
        prog='TeloComp',
        usage='%(prog)s [options]',
        exit_on_error=False,
        description='A tool for telomere extraction and genome polishing.',
        epilog='Text at the bottom of help'
    )

    parser.add_argument('-G', '--genome', metavar='', required=True, help='Input genome file (FASTA format)')
    parser.add_argument('--dir_contigs', metavar='', help='Input polished contigs')
    parser.add_argument('--dir_trim_L', metavar='', help='Input the trimmed reads. If the conditions are not met, extract the shortest reads.')
    parser.add_argument('--dir_trim_R', metavar='', help='Input the trimmed reads. If the conditions are not met, extract the shortest reads.')
    parser.add_argument('-L', '--lgsreads', metavar='', help='Long-read sequencing data')
    parser.add_argument('-W', '--wgs1', metavar='', help='Path to WGS reads (read 1)')
    parser.add_argument('-w', '--wgs2', metavar='', help='Path to WGS reads (read 2)')
    parser.add_argument('-N', '--NextPolish', metavar='', help='Path to NextPolish tool')
    parser.add_argument('-m', '--motif', metavar='', required=True,
                        help='Telomeric repeats sequences, e.g., plant: CCCTAAA(TTTAGGG), animal: TTAGGG(CCCTAA), etc.')
    parser.add_argument('-M', '--motif_num', metavar='', help='Input the number of bases of the telomere motif')
    parser.add_argument('--Normal', action='store_true', help='Execute the command according to the general process')
    parser.add_argument('--dir_Max', action='store_true', help='Select the telomere reads obtained by polishing the longest reads to add to the genome')
    parser.add_argument('--dir_Min', action='store_true', help='Select the telomere reads obtained by polishing the shortest reads to add to the genome')
    parser.add_argument('-t', '--threads', type=int, default=20, metavar='',
                        help='Number of threads to use (default: 20)')
    args = parser.parse_args()

    # Parameter selection section
    if args.Normal:
        # Check if directories exist
        for dir_path in [args.dir_trim_L, args.dir_trim_R]:
            if not os.path.exists(dir_path):
                print(f"Error: Directory {dir_path} does not exist.")
                return

        # Convert all paths to absolute paths
        args.genome = os.path.abspath(args.genome)
        args.dir_contigs = os.path.abspath(args.dir_contigs)
        args.dir_trim_L = os.path.abspath(args.dir_trim_L)
        args.dir_trim_R = os.path.abspath(args.dir_trim_R)
        args.lgsreads = os.path.abspath(args.lgsreads)
        args.wgs1 = os.path.abspath(args.wgs1)
        args.wgs2 = os.path.abspath(args.wgs2)
        args.NextPolish = os.path.abspath(args.NextPolish)

        # step1:Obtain sam files
        tmp_sam_dir = "tmp_sam"  # Used to compare back to the original genome and obtain the sam file
        tmp_getpos_dir = "tmp_getpos"  # Extract the shortest sequence and polish it.
        os.makedirs(tmp_sam_dir, exist_ok=True)
        os.makedirs(tmp_getpos_dir, exist_ok=True)

        contigs = args.dir_contigs  # The directory where the assembled and polished contig files are located
        contig_files = os.listdir(contigs)  # Get contig file list

        # Check if there is a file containing 'asm'
        has_asm_file = any('asm' in file for file in contig_files)

        if has_asm_file:
            # Process contig files
            for contig_file in contig_files:
                contig_path = Path(contigs) / contig_file
                if 'asm' in contig_file:
                    # Target path
                    sam_output_path = Path(tmp_sam_dir) / contig_file
                    # Copy the files to the tmp_sam directory
                    shutil.copy2(contig_path, sam_output_path)
                else:
                    # Target path
                    destination_path = Path(tmp_getpos_dir) / contig_file
                    # Create the target directory (if it does not exist)
                    os.makedirs(destination_path.parent, exist_ok=True)
                    # Copy the file to the tmp_getpos directory
                    shutil.copy2(contig_path, destination_path)
            # Call align_fasta_to_genome function
            align_fasta_to_genome(tmp_sam_dir, args.genome, args.dir_trim_L, args.dir_trim_R, tmp_getpos_dir, args.lgsreads, args.wgs1, args.wgs2, args.threads, args.NextPolish)
        else:
            # If there is no 'asm' file, set the tmp_getpos_dir directory and save the file there
            for contig_file in contig_files:
                contig_path = Path(contigs) / contig_file
                # Target path
                destination_path = Path(tmp_getpos_dir) / contig_file
                # Create the target directory (if it does not exist)
                os.makedirs(destination_path.parent, exist_ok=True)
                # Copy the file to the tmp_getpos directory
                shutil.copy2(contig_path, destination_path)
                
        tmp_getpos_dir = "tmp_getpos"

    # Step 2: Complement, extract and calculate telomere types
    # Set the tmp_getpos_dir directory according to the selected parameters
    elif args.dir_Max:
        tmp_getpos_dir = "tmp_Max_getpos"
    elif args.dir_Min:
        tmp_getpos_dir = "tmp_Min_getpos"
    genome_path = args.genome
    post_tel_to_genome(tmp_getpos_dir, genome_path, args.motif)

    # Extract the 1KB sequence of the chromosome corresponding to the new genome
    output_genome_path = "new_genome.fasta"
    output_positions_path = "telomere_positions.txt"
    extract_sequences(output_genome_path, output_positions_path)

    # Process and output telomere motif information
    if args.dir_Max:
        output_dir = "tmp_Max_getpos"
    elif args.dir_Min:
        output_dir = "tmp_Min_getpos"
    else:
        output_dir = "tmp_getpos"
    output_txt_path = "telomere_repeat_info.txt"
    if args.motif_num == "7":
        process_telomere7_sequences(output_dir, output_txt_path)
    elif args.motif_num == "6":
        process_telomere6_sequences(output_dir, output_txt_path)
    print(f"Telomere motif information has been saved to {output_txt_path}")

    # Draw the telomere distribution diagram at both ends
    if args.motif_num == "7":
        output_positions_path = "telomere_positions.txt"
        output_repeats_path = "telomere_repeat_info.txt"
        output_plot_dir = "telomere_plots18"
        KB_dir = "extract_teloSeq"
        plot_telomere7_data(KB_dir, output_positions_path, output_repeats_path, output_plot_dir)
    elif args.motif_num == "6":
        output_positions_path = "telomere_positions.txt"
        output_repeats_path = "telomere_repeat_info.txt"
        output_plot_dir = "telomere_plots18"
        KB_dir = "extract_teloSeq"
        plot_telomere6_data(KB_dir, output_positions_path, output_repeats_path, output_plot_dir)


    # Delete intermediate files
    tmp_getpos_dir = "tmp_getpos"
    tmp_sam_dir = "tmp_sam"
    shutil.rmtree(tmp_getpos_dir)
    shutil.rmtree(tmp_sam_dir)


if __name__ == "__main__":
    main()