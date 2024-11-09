#!/usr/bin/env python


import os
import re
import shutil
import argparse
import subprocess
from Bio import SeqIO


# Extract the shortest reads as contig
def extract_shortest_reads(fasta_path, output_path, seq_type, base_name, direction):
    if seq_type not in ["ONT", "HiFi", "merged"]:
        raise ValueError("seq_type must be one of 'ONT', 'HiFi', or 'merged'")

    ont_output_path = os.path.join(output_path, f"{base_name}_shortest_ONT_{direction}.fasta")
    hifi_output_path = os.path.join(output_path, f"{base_name}_shortest_HiFi_{direction}.fasta")
    ont_handle = open(ont_output_path, 'w') if seq_type == "merged" else None
    hifi_handle = open(hifi_output_path, 'w') if seq_type == "merged" else None
    out_handle = None

    if seq_type != "merged":
        shortest_reads_output = os.path.join(output_path, f"{base_name}_shortest_{seq_type}_{direction}.fasta")
        out_handle = open(shortest_reads_output, 'w')

    has_ont = False
    has_hifi = False

    for record in SeqIO.parse(fasta_path, 'fasta'):
        if "min" in record.id:
            if seq_type == "merged":
                if "ONT" in record.id:
                    has_ont = True
                    ont_handle.write(f">{record.id}\n{record.seq}\n")
                elif "HiFi" in record.id:
                    has_hifi = True
                    hifi_handle.write(f">{record.id}\n{record.seq}\n")
            else:
                out_handle.write(f">{record.id}\n{record.seq}\n")

    if seq_type == "merged":
        if ont_handle:
            ont_handle.close()
            if not has_ont:
                os.remove(ont_output_path)
        if hifi_handle:
            hifi_handle.close()
            if not has_hifi:
                os.remove(hifi_output_path)
    else:
        if out_handle:
            out_handle.close()

# Assemble reads into contigs
def run_alignments_assembly_for_directory(input_dir_L, input_dir_R, threads):
    inputdir_L = os.listdir(input_dir_L)
    inputdir_R = os.listdir(input_dir_R)
    input_merged_dir = inputdir_L + inputdir_R

    # Run assembly for all FASTA files in a directory and organize outputs based on direction
    for fastafile in input_merged_dir:
        if fastafile.endswith(".fasta"):
            base_name = os.path.basename(fastafile).split('_')[0]  # Get base name like chr1
            seq_type = os.path.basename(fastafile).split('_')[2]

            if "L" in fastafile:
                direction = "L"
                asm_out_dir = "asm_L"
                fasta_path = os.path.join(input_dir_L, fastafile)
            elif "R" in fastafile:
                direction = "R"
                asm_out_dir = "asm_R"
                fasta_path = os.path.join(input_dir_R, fastafile)
            else:
                print(f"Error: Cannot determine direction for {fastafile}")
                continue

            os.makedirs(asm_out_dir, exist_ok=True)
            asm_out = os.path.join(asm_out_dir, f"{base_name}_asm_{direction}")

            combined_command = (
                f"flye --meta --no-alt-contigs --nano-raw {fasta_path} --out-dir {asm_out} --threads {threads}"
            )

            try:
                subprocess.run(combined_command, shell=True, check=True)

                asm_file = os.path.join(asm_out, "assembly.fasta")
                new_asm_file = os.path.join(asm_out_dir, f"{base_name}_asm_{seq_type}_{direction}.fasta")
                if os.path.exists(asm_file):
                    os.rename(asm_file, new_asm_file)
            except subprocess.CalledProcessError as e:
                print(f"Error: Command '{e.cmd}' returned non-zero exit status {e.returncode}.")
                # If the command fails, extract the shortest reads
                extract_shortest_reads(fasta_path, asm_out_dir, seq_type, base_name, direction)
            finally:
                # Delete the asm_out directory
                if os.path.exists(asm_out):
                    shutil.rmtree(asm_out)


# racon preprocessing
def asm_contig_racon(input_file, lgsreads, wgs1, wgs2, threads, NextPolish):
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

    # After racon, the total fasta file is re-divided into individual fasta files with id as the file name
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

    # Remove the file (delete the newly generated merged file)
    new_input_file = f"new_merged_chr_shortest.fasta"
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
    os.remove(input_file)


def remove_sequence_from_fasta(fasta_file, sequence_id):
    output_Dir = "tmp_to_NP"
    os.makedirs(output_Dir, exist_ok=True)

    output_new_merged_file = f"new_merged_chr_shortest.fasta"
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


def asm_contig_polish(output_dir, lgsreads, reads1, reads2, threads, NextPolish):
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

        current_dir = os.getcwd()
        work_dir = os.path.join(current_dir, f"Dir_{outputname}")
        os.makedirs(work_dir, exist_ok=True)

        final_dir = os.path.join(current_dir, "files_NP")
        os.makedirs(final_dir, exist_ok=True)

        with open(sgs_fofn, 'w') as sgs_file:
            sgs_file.write(f"{reads1}\n")
            sgs_file.write(f"{reads2}\n")

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
        except subprocess.CalledProcessError as e:
            print(f"Error running nextPolish command for {reads_name}: {e}")
            print(f"Stdout: {e.stdout}")
            print(f"Stderr: {e.stderr}")
            print(
            f"{reads_name}.fasta quality is too low, polishing is unsuccessful, and the sequence is directly output to files_NP")
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


# Step 2: Nextpolish previous processing operations
def process_asm_merge_fasta_files(directories, output_file, lgsreads, reads1, reads2, threads, NextPolish):
    unique_sequences = set()
    with open(output_file, 'w') as outfile:
        for directory in directories:
            fasta_files = [f for f in os.listdir(directory) if f.endswith('.fasta')]
            for fasta_file in fasta_files:
                file_path = os.path.join(directory, fasta_file)
                if 'asm' in fasta_file:
                    with open(file_path, 'r') as infile:
                        for record in SeqIO.parse(infile, 'fasta'):
                            seq_id = os.path.splitext(fasta_file)[0]
                            if seq_id not in unique_sequences:
                                unique_sequences.add(seq_id)
                                outfile.write(f'>{seq_id}\n{str(record.seq)}\n')
                else:
                    with open(file_path, 'r') as infile:
                        for record in SeqIO.parse(infile, 'fasta'):
                            seq_id = os.path.splitext(fasta_file)[0]
                            lowercase_sequence = ''.join(char for char in str(record.seq) if char.islower())
                            if lowercase_sequence and seq_id not in unique_sequences:
                                unique_sequences.add(seq_id)
                                outfile.write(f'>{seq_id}\n{lowercase_sequence}\n')

    merged_file = output_file
    asm_contig_racon(merged_file, lgsreads, reads1, reads2, threads, NextPolish)


# Perform different processing based on the file name in the directory
def run_select_execution_function(dir_IN_L, dir_IN_R, lgsreads, wgs1, wgs2, threads, NextPolish):
    """Perform different processing based on the file name in the directory"""
    one_files = []
    trimmed_files = []

    # Collect file names
    for filename in os.listdir(dir_IN_L):
        if 'one' in filename:
            one_files.append(filename)
        elif 'trimmed' in filename:
            trimmed_files.append(filename)

    for filename in os.listdir(dir_IN_R):
        if 'one' in filename:
            one_files.append(filename)
        elif 'trimmed' in filename:
            trimmed_files.append(filename)

    # If all are one files
    if one_files and not trimmed_files:
        asm_L_dir = "asm_L"
        asm_R_dir = "asm_R"
        os.makedirs(asm_L_dir, exist_ok=True)
        os.makedirs(asm_R_dir, exist_ok=True)
        # Extract the assembled directory
        polish_asmmerged_DIR = ['asm_L', 'asm_R']
        for filename in one_files:
            # Move into asm_L or asm_R depending on direction
            if 'L' in filename:
                shutil.move(os.path.join(dir_IN_L, filename), os.path.join(polish_asmmerged_DIR[0], filename))
            elif 'R' in filename:
                shutil.move(os.path.join(dir_IN_R, filename), os.path.join(polish_asmmerged_DIR[1], filename))
        output_file = "asm_merged_All_chr.fasta"
        process_asm_merge_fasta_files(polish_asmmerged_DIR, output_file, lgsreads, wgs1, wgs2, threads, NextPolish)

    # If there are both one and trimmed files
    elif one_files and trimmed_files:
        asm_L_dir = "asm_L"
        asm_R_dir = "asm_R"
        os.makedirs(asm_L_dir, exist_ok=True)
        os.makedirs(asm_R_dir, exist_ok=True)
        # Process one file
        polish_asmmerged_DIR = ['asm_L', 'asm_R']
        for filename in one_files:
            # Move into asm_L or asm_R depending on direction
            if 'L' in filename:
                shutil.move(os.path.join(dir_IN_L, filename), os.path.join(polish_asmmerged_DIR[0], filename))
            elif 'R' in filename:
                shutil.move(os.path.join(dir_IN_R, filename), os.path.join(polish_asmmerged_DIR[1], filename))

        # Processing trimmed files
        run_alignments_assembly_for_directory(dir_IN_L, dir_IN_R, threads)

    # If there is no one file, only trimmed files
    elif trimmed_files and not one_files:
        run_alignments_assembly_for_directory(dir_IN_L, dir_IN_R, threads)

    # If there is neither one file nor trimmed file
    else:
        print("No matching files found")


def main():
    parser = argparse.ArgumentParser(
        prog='TeloComp',
        usage='%(prog)s [options]',
        exit_on_error=False,
        description='A tool for telomere extraction and genome polishing.',
        epilog='Text at the bottom of help'
    )

    parser.add_argument('--dir_IN_L', metavar='', help='Directory containing left-aligned reads (FASTA format)')
    parser.add_argument('--dir_IN_R', metavar='', help='Directory containing right-aligned reads (FASTA format)')
    parser.add_argument('-L', '--lgsreads', metavar='', required=True, help='Long-read sequencing data')
    parser.add_argument('-W', '--wgs1', metavar='', required=True, help='Path to WGS reads (read 1)')
    parser.add_argument('-w', '--wgs2', metavar='', required=True, help='Path to WGS reads (read 2)')
    parser.add_argument('-N', '--NextPolish', metavar='', required=True, help='Path to NextPolish tool')
    parser.add_argument('-t', '--threads', type=int, default=20, metavar='',
                        help='Number of threads to use (default: 20)')
    args = parser.parse_args()

    #Step1: Assemble reads into contigs 
    if args.dir_IN_L and args.dir_IN_R:
        run_select_execution_function(args.dir_IN_L, args.dir_IN_R, args.lgsreads, args.wgs1, args.wgs2, args.threads, args.NextPolish)

    # Step2: Polish the longest or shortest sequence 
    # Extract the assembled directory
    polish_asmmerged_DIR = ['asm_L', 'asm_R']
    output_file = "asm_merged_All_chr.fasta"

    process_asm_merge_fasta_files(polish_asmmerged_DIR, output_file, args.lgsreads, args.wgs1, args.wgs2, args.threads, args.NextPolish)


if __name__ == "__main__":
    main()
