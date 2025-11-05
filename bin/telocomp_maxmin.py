#!/usr/bin/env python

import os
import shutil
import argparse
import subprocess
import re
from Bio import SeqIO


# Merge the L and R FASTA files that are present in both ONT and HiFi datasets
def merge_sequences(ont_filename, hifi_filename, merged_file):
    with open(merged_file, 'w') as out_handle:
        ont_records = {record.id: record for record in SeqIO.parse(ont_filename, 'fasta')}
        hifi_records = {record.id: record for record in SeqIO.parse(hifi_filename, 'fasta')}
        for record_id, ont_record in ont_records.items():
            out_handle.write(f'>{record_id}_ONT\n{ont_record.seq}\n')
            if record_id in hifi_records:
                hifi_record = hifi_records[record_id]
                out_handle.write(f'>{record_id}_HiFi\n{hifi_record.seq}\n')
        for record_id, hifi_record in hifi_records.items():
            if record_id not in ont_records:
                out_handle.write(f'>{record_id}_HiFi\n{hifi_record.seq}\n')


# Extract the longest lowercase sequence as the contig
def max_length_extract(input_file, chromosome, direction, seq_type):
    output_dir = f"MaxLength_{direction}"
    os.makedirs(output_dir, exist_ok=True)

    output_file = f"{chromosome}_max_{seq_type}_{direction}.fasta"

    lengths = {'lowercase': {}}
    for record in SeqIO.parse(input_file, 'fasta'):
        sequence = str(record.seq)
        lowercase_sequence = ''.join(char for char in sequence if char.islower())
        lengths['lowercase'][record.id] = len(lowercase_sequence)

    max_record_id = max(lengths['lowercase'], key=lengths['lowercase'].get)
    with open(output_file, 'w') as out_handle:
        for record in SeqIO.parse(input_file, 'fasta'):
            if record.id == max_record_id:
                # Extract lowercase portions of the sequence
                lowercase_seq = ''.join([char for char in str(record.seq) if char.islower()])
                out_handle.write(f'>{record.id}\n{lowercase_seq}\n')

    subprocess.run(f"mv {output_file} {output_dir}", shell=True, check=True)


# Extract the shortest lowercase sequence as the contig
def min_length_extract(input_file, chromosome, direction, seq_type):
    output_dir = f"MinLength_{direction}"
    os.makedirs(output_dir, exist_ok=True)

    output_file = f"{chromosome}_min_{seq_type}_{direction}.fasta"

    lengths = {'lowercase': {}}
    for record in SeqIO.parse(input_file, 'fasta'):
        sequence = str(record.seq)
        lowercase_sequence = ''.join(char for char in sequence if char.islower())
        lengths['lowercase'][record.id] = len(lowercase_sequence)

    min_record_id = min(lengths['lowercase'], key=lengths['lowercase'].get)
    with open(output_file, 'w') as out_handle:
        for record in SeqIO.parse(input_file, 'fasta'):
            if record.id == min_record_id:
                # Extract lowercase portions of the sequence
                lowercase_seq = ''.join([char for char in str(record.seq) if char.islower()])
                out_handle.write(f'>{record.id}\n{lowercase_seq}\n')

    subprocess.run(f"mv {output_file} {output_dir}", shell=True, check=True)


# Preprocess with Racon
def Max_contig_racon(merged_file, lgsreads, wgs1, wgs2, threads, NextPolish, polish):
    tmpdir = "tmp_dir"
    os.makedirs(tmpdir, exist_ok=True)

    try:
        if polish:
            while True:  # Use iteration instead of recursion
                racon_commands = [
                    f"minimap2 -ax map-pb -t 50 {merged_file} {lgsreads} | samtools sort - -m 2g --threads 20 -o {tmpdir}/genome.lgs.bam",
                    f"samtools index {tmpdir}/genome.lgs.bam",
                    f"ls {os.path.abspath(tmpdir)}/genome.lgs.bam > {tmpdir}/pb.map.bam.fofn",
                    f"python {NextPolish}/lib/nextpolish2.py -g {merged_file} -l {tmpdir}/pb.map.bam.fofn -r clr -p 50 -a -s -o {tmpdir}/merged_genome.lgspolish.fasta"
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

                    break  # Break the loop if all commands are successfully executed

                except subprocess.CalledProcessError as e:
                    print(f"Error running racon_commands: {e.stderr.decode()}")
                    error_message = e.stderr.decode()

                    if "Failed to correct sequence" in error_message:
                        # Use regular expressions to extract and clean the string
                        match = re.search(r'Failed to correct sequence:\s+(\S+)', error_message)
                        if match:
                            failed_sequence = match.group(1).strip()
                            # Further clean non-printable characters from the string
                            clean_sequence_id = re.sub(r'[^\w_\-]', '', failed_sequence)

                            print(f"Failed sequence ID: {clean_sequence_id}")

                            # Remove the corresponding sequences from the FASTA file
                            new_input_file = Max_remove_sequence_from_fasta(input_file, clean_sequence_id)
                            if new_input_file == merged_file:
                                print(f"Failed sequence {clean_sequence_id} was not found. Exiting.")
                                return
                            else:
                                merged_file = new_input_file
                                print(f"Re-running racon with updated FASTA file: {input_file}")
                        else:
                            print("Failed to parse the failed sequence ID. Exiting.")
                            return

                    else:
                        print("Unexpected error occurred. Exiting.")
                        return

            # Split the combined FASTA file into individual FASTA files
            racon_genome = f"{tmpdir}/merged_genome.lgspolish.fasta"
            output_dir = "MaxLength_NP"
            if not os.path.exists(output_dir):
                os.makedirs(output_dir, exist_ok=True)
            for record in SeqIO.parse(racon_genome, "fasta"):
                base_name = record.id
                output_file = os.path.join(output_dir, f"{base_name}.fasta")
                with open(output_file, 'w') as out_handle:
                    SeqIO.write(record, out_handle, "fasta")

            Max_contig_polish(output_dir, lgsreads, wgs1, wgs2, threads, NextPolish, polish)

            # Delete the temporary directory and its copies
            shutil.rmtree(tmpdir)
        else:
            # If not polishing, just copy the merged file to output directory
            output_dir = "MaxLength_NP"
            os.makedirs(output_dir, exist_ok=True)
            for record in SeqIO.parse(merged_file, "fasta"):
                output_file = os.path.join(output_dir, f"{record.id}.fasta")
                with open(output_file, 'w') as out_handle:
                    SeqIO.write(record, out_handle, "fasta")
            shutil.rmtree(tmpdir)
    finally:
        # 无论如何都会执行
        if os.path.isdir(tmpdir):
            shutil.rmtree(tmpdir)
            print(f"Temporary directory '{tmpdir}' deleted.")
        # Delete any _merged_*.fasta temporary files
        merged_files = [f for f in os.listdir('.') if re.match(r'.*_merged_.*\.fasta$', f)]
        for f in merged_files:
            try:
                os.remove(f)
                print(f"Deleted temporary merged file: {f}")
            except Exception as e:
                print(f"Warning: failed to delete {f}: {e}")


def Max_remove_sequence_from_fasta(fasta_file, sequence_id):
    output_Dir = "tmp_getpos"
    os.makedirs(output_Dir, exist_ok=True)

    output_new_merged_file = os.path.join(output_Dir, "new_merged_chr_shortest.fasta")
    sequence_removed = False

    # Use regular expressions to remove control characters from the sequence IDs
    clean_sequence_id = re.sub(r'[^\w_\-]', '', sequence_id.strip())

    print(f"Processed sequence ID: {clean_sequence_id}")

    with open(fasta_file, "r") as input_handle, open(output_new_merged_file, "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            record_id = re.sub(r'[^\w_\-]', '', record.id.strip())

            # Print all recorded IDs
            print(f"Record ID in FASTA file: '{record_id}'")

            # Check partial matches
            if clean_sequence_id in record_id or record_id in clean_sequence_id:
                # Save the failed sequences to a separate file
                failed_seq_file = os.path.join(output_Dir, f"{record_id}.fasta")
                with open(failed_seq_file, "w") as failed_handle:
                    SeqIO.write(record, failed_handle, "fasta")
                print(f"Saved failed sequence {record_id} to {failed_seq_file}")
                sequence_removed = True
            else:
                # Write the remaining sequences to a new FASTA file
                SeqIO.write(record, output_handle, "fasta")

    if sequence_removed:
        print(f"New merged FASTA file created: {output_new_merged_file}")
        return output_new_merged_file
    else:
        print(f"Sequence ID '{clean_sequence_id}' not found in {fasta_file}. Return the original file.")
        return fasta_file


def Max_contig_polish(output_dir, lgsreads, wgs1, wgs2, threads, NextPolish, polish):
    fasta_files = [f for f in os.listdir(output_dir) if f.endswith('.fasta')]
    for fasta_file in fasta_files:
        reads_name = os.path.splitext(fasta_file)[0]
        file_path = os.path.join(output_dir, fasta_file)
        genome = file_path
        FoFn_dir = os.path.join("Dir_Max_fofn", reads_name)
        os.makedirs(FoFn_dir, exist_ok=True)

        # Use only the filename of the input file
        genome_final = os.path.basename(genome)
        shutil.copyfile(genome, os.path.join(FoFn_dir, genome_final))

        if polish:
            sgs_fofn = os.path.join(FoFn_dir, "sgs.fofn")
            lgs_fofn = os.path.join(FoFn_dir, "lgs.fofn")
            nextPolish = os.path.join(NextPolish, "nextPolish")
            t = threads
            outputname = reads_name

            # Get the current working directory
            current_dir = os.getcwd()
            work_dir = os.path.join(current_dir, f"Dir_{outputname}")
            os.makedirs(work_dir, exist_ok=True)

            final_dir = os.path.join(current_dir, "MaxLength_NP")
            os.makedirs(final_dir, exist_ok=True)

            with open(sgs_fofn, 'w') as sgs_file:
                sgs_file.write(f"{wgs1}\n")
                sgs_file.write(f"{wgs2}\n")

            with open(lgs_fofn, 'w') as lgs_file:
                lgs_file.write(f"{lgsreads}\n")

            # Generate configuration file
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

            Max_polish_command = f"{nextPolish} {config_file_path}"
            try:
                subprocess.run(Max_polish_command, shell=True, check=True, capture_output=True, text=True)
                np_fastas = [f for f in os.listdir(work_dir) if f.endswith('.fasta')]
                for fasta in np_fastas:
                    # Rename to {reads_name}.fasta
                    new_name = f"{reads_name}.fasta"
                    new_path = os.path.join(final_dir, new_name)
                    shutil.move(os.path.join(work_dir, fasta), new_path)
            except subprocess.CalledProcessError as e:
                print(f"Error running nextPolish command for {reads_name}: {e}")
                print(f"Stdout: {e.stdout}")
                print(f"Stderr: {e.stderr}")
                print(
                f"{reads_name}.fasta quality is too low, polishing is unsuccessful, and the sequence is directly output to files_NP")
                # If the command fails, copy the file to the files_NP directory
                shutil.copy2(file_path, final_dir)
            finally:
                # Delete the work_dir and any backup files
                if os.path.isdir(work_dir):
                    shutil.rmtree(work_dir, ignore_errors=True)
                backup_files = [f for f in os.listdir(current_dir) if f.startswith(f"Dir_{outputname}.backup")]
                for backup in backup_files:
                    backup_path = os.path.join(current_dir, backup)
                    shutil.rmtree(backup_path, ignore_errors=True)
                print(f"Deleted work directory {work_dir} and backup files")
                
                # Delete the Dir_Max_fofn folder
                if os.path.isdir("Dir_Max_fofn"):
                    shutil.rmtree("Dir_Max_fofn", ignore_errors=True)
                    print("Deleted Dir_Max_fofn directory")

                # Delete all .log.info files
                for f in os.listdir(current_dir):
                    if f.endswith(".log.info"):
                        try:
                            os.remove(os.path.join(current_dir, f))
                            print(f"Deleted log file: {f}")
                        except Exception as e:
                            print(f"Failed to delete {f}: {e}")
        else:
            # If not polishing, just copy the file to final directory
            final_dir = os.path.join(os.getcwd(), "MaxLength_NP")
            os.makedirs(final_dir, exist_ok=True)
            shutil.copy2(file_path, os.path.join(final_dir, f"{reads_name}.fasta"))


# Preprocess with Racon
def Min_contig_racon(merged_file, lgsreads, wgs1, wgs2, threads, NextPolish, polish):
    tmpdir = "tmp_dir"
    os.makedirs(tmpdir, exist_ok=True)

    try:
        if polish:
            while True:  # Use iteration instead of recursion
                racon_commands = [
                    f"minimap2 -ax map-pb -t 50 {merged_file} {lgsreads} | samtools sort - -m 2g --threads 20 -o {tmpdir}/genome.lgs.bam",
                    f"samtools index {tmpdir}/genome.lgs.bam",
                    f"ls {os.path.abspath(tmpdir)}/genome.lgs.bam > {tmpdir}/pb.map.bam.fofn",
                    f"python {NextPolish}/lib/nextpolish2.py -g {merged_file} -l {tmpdir}/pb.map.bam.fofn -r clr -p 50 -a -s -o {tmpdir}/merged_genome.lgspolish.fasta"
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

                    break  # Break the loop if all commands are successfully executed

                except subprocess.CalledProcessError as e:
                    print(f"Error running racon_commands: {e.stderr.decode()}")
                    error_message = e.stderr.decode()

                    if "Failed to correct sequence" in error_message:
                        # Use regular expressions to extract and clean the string
                        match = re.search(r'Failed to correct sequence:\s+(\S+)', error_message)
                        if match:
                            failed_sequence = match.group(1).strip()
                            # Further clean non-printable characters from the string
                            clean_sequence_id = re.sub(r'[^\w_\-]', '', failed_sequence)

                            print(f"Failed sequence ID: {clean_sequence_id}")

                            # Remove the corresponding sequences from the FASTA file
                            new_input_file = Min_remove_sequence_from_fasta(input_file, clean_sequence_id)
                            if new_input_file == merged_file:
                                print(f"Failed sequence {clean_sequence_id} was not found. Exiting.")
                                return
                            else:
                                merged_file = new_input_file
                                print(f"Re-running racon with updated FASTA file: {input_file}")
                        else:
                            print("Failed to parse the failed sequence ID. Exiting.")
                            return

                    else:
                        print("Unexpected error occurred. Exiting.")
                        return

            # Split the combined FASTA file into individual FASTA files
            racon_genome = f"{tmpdir}/merged_genome.lgspolish.fasta"
            output_dir = "MinLength_NP"
            if not os.path.exists(output_dir):
                os.makedirs(output_dir, exist_ok=True)
            for record in SeqIO.parse(racon_genome, "fasta"):
                base_name = record.id
                output_file = os.path.join(output_dir, f"{base_name}.fasta")
                with open(output_file, 'w') as out_handle:
                    SeqIO.write(record, out_handle, "fasta")

            Min_contig_polish(output_dir, lgsreads, wgs1, wgs2, threads, NextPolish, polish)

            # Delete the temporary directory and its copies
            shutil.rmtree(tmpdir)
        else:
            # If not polishing, just copy the merged file to output directory
            output_dir = "MinLength_NP"
            os.makedirs(output_dir, exist_ok=True)
            for record in SeqIO.parse(merged_file, "fasta"):
                output_file = os.path.join(output_dir, f"{record.id}.fasta")
                with open(output_file, 'w') as out_handle:
                    SeqIO.write(record, out_handle, "fasta")
            shutil.rmtree(tmpdir)
    finally:
        # 无论如何都会执行
        if os.path.isdir(tmpdir):
            shutil.rmtree(tmpdir)
            print(f"Temporary directory '{tmpdir}' deleted.")

        # Delete any _merged_*.fasta temporary files
        merged_files = [f for f in os.listdir('.') if re.match(r'.*_merged_.*\.fasta$', f)]
        for f in merged_files:
            try:
                os.remove(f)
                print(f"Deleted temporary merged file: {f}")
            except Exception as e:
                print(f"Warning: failed to delete {f}: {e}")
                

def Min_remove_sequence_from_fasta(fasta_file, sequence_id):
    output_Dir = "tmp_getpos"
    os.makedirs(output_Dir, exist_ok=True)

    output_new_merged_file = os.path.join(output_Dir, "new_merged_chr_shortest.fasta")
    sequence_removed = False

    # Use regular expressions to remove control characters from the sequence IDs
    clean_sequence_id = re.sub(r'[^\w_\-]', '', sequence_id.strip())

    print(f"Processed sequence ID: {clean_sequence_id}")

    with open(fasta_file, "r") as input_handle, open(output_new_merged_file, "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            record_id = re.sub(r'[^\w_\-]', '', record.id.strip())

            # Print all recorded IDs
            print(f"Record ID in FASTA file: '{record_id}'")

            # Check partial matches (alternatively, you can use the startswith or in methods)
            if clean_sequence_id in record_id or record_id in clean_sequence_id:
                # Save the failed sequences to a separate file
                failed_seq_file = os.path.join(output_Dir, f"{record_id}.fasta")
                with open(failed_seq_file, "w") as failed_handle:
                    SeqIO.write(record, failed_handle, "fasta")
                print(f"Saved failed sequence {record_id} to {failed_seq_file}")
                sequence_removed = True
            else:
                # Write the remaining sequences to a new FASTA file
                SeqIO.write(record, output_handle, "fasta")

    if sequence_removed:
        print(f"New merged FASTA file created: {output_new_merged_file}")
        return output_new_merged_file
    else:
        print(f"Sequence ID '{clean_sequence_id}' not found in {fasta_file}. Return the original file.")
        return fasta_file


def Min_contig_polish(output_dir, lgsreads, wgs1, wgs2, threads, NextPolish, polish):
    fasta_files = [f for f in os.listdir(output_dir) if f.endswith('.fasta')]
    for fasta_file in fasta_files:
        reads_name = os.path.splitext(fasta_file)[0]
        file_path = os.path.join(output_dir, fasta_file)
        genome = file_path
        FoFn_dir = os.path.join("Dir_Min_fofn", reads_name)
        os.makedirs(FoFn_dir, exist_ok=True)

        # Use only the filename of the input file
        genome_final = os.path.basename(genome)
        shutil.copyfile(genome, os.path.join(FoFn_dir, genome_final))

        if polish:
            sgs_fofn = os.path.join(FoFn_dir, "sgs.fofn")
            lgs_fofn = os.path.join(FoFn_dir, "lgs.fofn")
            nextPolish = os.path.join(NextPolish, "nextPolish")
            t = threads
            outputname = reads_name

            # Get the current working directory
            current_dir = os.getcwd()
            work_dir = os.path.join(current_dir, f"Dir_{outputname}")
            os.makedirs(work_dir, exist_ok=True)

            final_dir = os.path.join(current_dir, "MinLength_NP")
            os.makedirs(final_dir, exist_ok=True)

            with open(sgs_fofn, 'w') as sgs_file:
                sgs_file.write(f"{wgs1}\n")
                sgs_file.write(f"{wgs2}\n")

            with open(lgs_fofn, 'w') as lgs_file:
                lgs_file.write(f"{lgsreads}\n")

            # Generate configuration file
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

            Min_polish_command = f"{nextPolish} {config_file_path}"
            try:
                subprocess.run(Min_polish_command, shell=True, check=True, capture_output=True, text=True)
                np_fastas = [f for f in os.listdir(work_dir) if f.endswith('.fasta')]
                for fasta in np_fastas:
                    # Rename to {reads_name}.fasta
                    new_name = f"{reads_name}.fasta"
                    new_path = os.path.join(final_dir, new_name)
                    shutil.move(os.path.join(work_dir, fasta), new_path)
            except subprocess.CalledProcessError as e:
                print(f"Error running nextPolish command for {reads_name}: {e}")
                print(f"Stdout: {e.stdout}")
                print(f"Stderr: {e.stderr}")
                print(
                f"{reads_name}.fasta quality is too low, polishing is unsuccessful, and the sequence is directly output to files_NP")
                # If the command fails, copy the file to the files_NP directory
                shutil.copy2(file_path, final_dir)
            finally:
#                # Delete the work_dir and any backup files
#                shutil.rmtree(work_dir)
#                backup_files = [f for f in os.listdir(current_dir) if f.startswith(f"Dir_{outputname}.backup")]
#                for backup in backup_files:
#                    backup_path = os.path.join(current_dir, backup)
#                    shutil.rmtree(backup_path)
#                print(f"Deleted work directory {work_dir} and backup files")
#        else:
#            # If not polishing, just copy the file to final directory
#            final_dir = os.path.join(os.getcwd(), "MinLength_NP")
#            os.makedirs(final_dir, exist_ok=True)
#            shutil.copy2(file_path, os.path.join(final_dir, f"{reads_name}.fasta"))
                # Delete the work_dir and any backup files
                if os.path.isdir(work_dir):
                    shutil.rmtree(work_dir, ignore_errors=True)
                backup_files = [f for f in os.listdir(current_dir) if f.startswith(f"Dir_{outputname}.backup")]
                for backup in backup_files:
                    backup_path = os.path.join(current_dir, backup)
                    shutil.rmtree(backup_path, ignore_errors=True)
                print(f"Deleted work directory {work_dir} and backup files")
                
                # Delete the Dir_Min_fofn folder
                if os.path.isdir("Dir_Min_fofn"):
                    shutil.rmtree("Dir_Min_fofn", ignore_errors=True)
                    print("Deleted Dir_Min_fofn directory")

                # Delete all .log.info files
                for f in os.listdir(current_dir):
                    if f.endswith(".log.info"):
                        try:
                            os.remove(os.path.join(current_dir, f))
                            print(f"Deleted log file: {f}")
                        except Exception as e:
                            print(f"Failed to delete {f}: {e}")
        else:
            # If not polishing, just copy the file to final directory
            final_dir = os.path.join(os.getcwd(), "MinLength_NP")
            os.makedirs(final_dir, exist_ok=True)
            shutil.copy2(file_path, os.path.join(final_dir, f"{reads_name}.fasta"))


def process_extract_Max_Min(chromosome, extr_merged_files, input_dir_ont, input_dir_hifi, direction, max_length,
                            min_length, temp_files):
    # Modify filename matching logic
    ont_filename = f"{chromosome}_{direction}.fasta"  # Change to chr{num}_L.fasta format
    hifi_filename = f"{chromosome}_{direction}.fasta"  # Change to chr{num}_R.fasta format
    
    print(f"\nProcessing chromosome {chromosome} {direction}")
    print(f"Looking for files: {ont_filename} in ONT, {hifi_filename} in HiFi")
    
    input_file = None
    seq_type = None

    # Check if the file exists
    ont_file = os.path.join(input_dir_ont, ont_filename)
    hifi_file = os.path.join(input_dir_hifi, hifi_filename)
    
    ont_exists = os.path.exists(ont_file)
    hifi_exists = os.path.exists(hifi_file)
    
    print(f"ONT file exists: {ont_exists} ({ont_file})")
    print(f"HiFi file exists: {hifi_exists} ({hifi_file})")

    if ont_exists and hifi_exists:
        print("Found both ONT and HiFi files")
        seq_type = "merged"
        merged_file = f"{chromosome}_merged_{direction}.fasta"
        print(f"Merging files to: {merged_file}")
        merge_sequences(ont_file, hifi_file, merged_file)
        input_file = merged_file
        temp_files.append(input_file)
    elif ont_exists:
        print("Found only ONT file")
        seq_type = "ONT"
        input_file = ont_file
    elif hifi_exists:
        print("Found only HiFi file")
        seq_type = "HiFi"
        input_file = hifi_file
    else:
        print(f"Warning: No input files found for {chromosome} {direction}")
        return

    if input_file:
        print(f"Processing {seq_type} file: {input_file}")
        if max_length:
            print("Extracting max length...")
            max_length_extract(input_file, f"{chromosome}", direction, seq_type)
        if min_length:
            print("Extracting min length...")
            min_length_extract(input_file, f"{chromosome}", direction, seq_type)


# Step 2: Merge and Trim Sequences 
def process_Max_merge_fasta_files(directories, output_file, lgsreads, wgs1, wgs2, threads, NextPolish, polish):
    unique_sequence = set()
    with open(output_file, 'w') as outfile:
        for directory in directories:
            fasta_files = [f for f in os.listdir(directory) if f.endswith('.fasta')]
            for fasta_file in fasta_files:
                chr_name = os.path.basename(fasta_file).split('_')[0]
                directorie_max = os.path.basename(fasta_file).split('_')[-1].split('.')[0]
                file_path = os.path.join(directory, fasta_file)
                with open(file_path, 'r') as infile:
                    for record in SeqIO.parse(infile, 'fasta'):
                        seq_id = f"{chr_name}_longest_np_{directorie_max}"
                        if record.seq and seq_id not in unique_sequence:
                            unique_sequence.add(seq_id)
                            outfile.write(f'>{seq_id}\n{record.seq}\n')

    merged_file = output_file
    Max_contig_racon(merged_file, lgsreads, wgs1, wgs2, threads, NextPolish, polish)

    # Remove the file
    os.remove(merged_file)


def process_Min_merge_fasta_files(directories, output_file, lgsreads, wgs1, wgs2, threads, NextPolish, polish):
    unique_sequence = set()
    with open(output_file, 'w') as outfile:
        for directory in directories:
            fasta_files = [f for f in os.listdir(directory) if f.endswith('.fasta')]
            for fasta_file in fasta_files:
                chr_name = os.path.basename(fasta_file).split('_')[0]
                directorie_min = os.path.basename(fasta_file).split('_')[-1].split('.')[0]
                file_path = os.path.join(directory, fasta_file)
                with open(file_path, 'r') as infile:
                    for record in SeqIO.parse(infile, 'fasta'):
                        seq_id = f"{chr_name}_shortest_np_{directorie_min}"
                        if record.seq and seq_id not in unique_sequence:
                            unique_sequence.add(seq_id)
                            outfile.write(f'>{seq_id}\n{record.seq}\n')

    merged_file = output_file
    Min_contig_racon(merged_file, lgsreads, wgs1, wgs2, threads, NextPolish, polish)

    # Remove the file
    os.remove(merged_file)


def main():
    parser = argparse.ArgumentParser(
        prog='TeloComp',
        usage='%(prog)s [options]',
        exit_on_error=False,
        description='A tool for telomere extraction and genome polishing.',
        epilog='Text at the bottom of help')

    parser.add_argument('--Max_length', action='store_true', help='Extract longest reads')
    parser.add_argument('--Min_length', action='store_true', help='Extract shortest reads')
    parser.add_argument('--dir_ont', metavar='', required=True, help='Directory containing ONT files')
    parser.add_argument('--dir_hifi', metavar='', required=True, help='Directory containing HiFi files')
    parser.add_argument('-L', '--lgsreads', metavar='', help='Long-read sequencing data')
    parser.add_argument('-W', '--wgs1', metavar='', help='Path to WGS reads (read 1)')
    parser.add_argument('-w', '--wgs2', metavar='', help='Path to WGS reads (read 2)')
    parser.add_argument('-N', '--NextPolish', metavar='', help='Path to NextPolish tool')
    parser.add_argument('-t', '--threads', type=int, default=20, metavar='',
                        help='Number of threads to use (default: 20)')
    parser.add_argument('--polish', action='store_true', help='Perform polishing with NextPolish')

    args = parser.parse_args()

    # Check if NextPolish path is provided when polishing is requested
    if args.polish and not args.NextPolish:
        parser.error("--NextPolish is required when --polish is specified")

    # Step1: Extract longest or shortest sequence
    input_dir_ont = args.dir_ont
    input_dir_hifi = args.dir_hifi
    
	#Scanning fasta
    extr_ont_files = [f for f in sorted(os.listdir(input_dir_ont))
                      if f.endswith('.fasta') and ('_L.fasta' in f or '_R.fasta' in f)]
    extr_hifi_files = [f for f in sorted(os.listdir(input_dir_hifi))
                       if f.endswith('.fasta') and ('_L.fasta' in f or '_R.fasta' in f)]

    print("\nFound ONT files:", extr_ont_files)
    print("Found HiFi files:", extr_hifi_files)

	# Extract the "chromosome name" or "group name", removing the trailing _L/_R.
    chromosomes = set()
    for f in extr_ont_files + extr_hifi_files:
        chrom = "_".join(f.split('_')[:-1])
        chromosomes.add(chrom)

    print("Chromosomes to process:", sorted(chromosomes))

    
    processed_chromosomes = set()
    temp_files = []

    for chromosome in sorted(chromosomes):
        if chromosome in processed_chromosomes:
            continue

        print(f"\nProcessing chromosome: {chromosome}")
        process_extract_Max_Min(chromosome, extr_ont_files + extr_hifi_files, input_dir_ont, input_dir_hifi, 
                               'L', args.Max_length, args.Min_length, temp_files)
        process_extract_Max_Min(chromosome, extr_ont_files + extr_hifi_files, input_dir_ont, input_dir_hifi, 
                               'R', args.Max_length, args.Min_length, temp_files)
        processed_chromosomes.add(chromosome)

    # Step2: Polish the longest or shortest sequence

    if args.Max_length:
        # Extract the directory of the longest sequence output
        polish_Maxmerged_DIR = ['MaxLength_L', 'MaxLength_R']
        output_file = f"merged_Max_chr.fasta"

        process_Max_merge_fasta_files(polish_Maxmerged_DIR, output_file, args.lgsreads, args.wgs1, args.wgs2,
                                      args.threads, args.NextPolish, args.polish)

    elif args.Min_length:
        # Extract the directory of the shortest sequence output
        polish_Minmerged_DIR = ['MinLength_L', 'MinLength_R']
        output_file = f"merged_Min_chr.fasta"

        process_Min_merge_fasta_files(polish_Minmerged_DIR, output_file, args.lgsreads, args.wgs1, args.wgs2,
                                      args.threads, args.NextPolish, args.polish)

    # Remove temporary files after processing
    for temp_file in temp_files:
        if os.path.exists(temp_file):
            os.remove(temp_file)


if __name__ == "__main__":
    main()
