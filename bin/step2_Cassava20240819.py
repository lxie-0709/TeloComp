#组装成功的直接打磨，组装不成功直接拿最短的小写序列出来，先都合在一起racon，然后分开单个循环打磨成功
import os
import re
import shutil
import argparse
import subprocess
from Bio import SeqIO


### Step 1: Extract end alignment read and rename output data ###
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
    # subprocess.run(combined_command, shell=True, check=True)
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
    # subprocess.run(combined_command, shell=True, check=True)
    rename_HiFi_in_directory(algn_output_hifi, data_type)
    return algn_output_hifi


# 提取最短的reads作为contig
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
                # 如果命令失败，则提取最短reads
                extract_shortest_reads(fasta_path, asm_out_dir, seq_type, base_name, direction)
            finally:
                # 删除asm_out目录
                if os.path.exists(asm_out):
                    shutil.rmtree(asm_out)


# 20240906 修改racon逻辑
def asm_contig_racon(input_file, lgsreads, wgs1, wgs2, threads, NextPolish):
    tmpdir = "tmp_dir"
    os.makedirs(tmpdir, exist_ok=True)

    while True:  # 使用循环而非递归
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

            break  # 如果所有命令都成功执行，跳出循环

        except subprocess.CalledProcessError as e:
            print(f"Error running racon_commands: {e.stderr.decode()}")
            error_message = e.stderr.decode()

            if "Failed to correct sequence" in error_message:
                # 使用正则表达式提取并清理字符串
                match = re.search(r'Failed to correct sequence:\s+(\S+)', error_message)
                if match:
                    failed_sequence = match.group(1).strip()
                    # 进一步清理字符串中的非打印字符
                    clean_sequence_id = re.sub(r'[^\w_\-]', '', failed_sequence)

                    print(f"Failed sequence ID: {clean_sequence_id}")

                    # 从FASTA文件中删除相应的序列
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

    #将racon过后，总的fasta文件以id为文件名，重新分成单个fasta文件
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

    # 20240906新修改
    # Remove the file(删除新生成的merged文件）
    new_input_file = f"new_merged_chr_shortest.fasta"
    # 判断文件是否存在
    if os.path.exists(new_input_file):
        # 存在，则删除文件
        os.remove(new_input_file)
    else:
        print("Execute the nextPolish polishing command...")

    asm_contig_polish(output_dir, lgsreads, wgs1, wgs2, threads, NextPolish)

    # 删除临时目录和文件
    shutil.rmtree(tmpdir)
    shutil.rmtree(output_dir)


def remove_sequence_from_fasta(fasta_file, sequence_id):
    output_Dir = "tmp_to_NP"
    os.makedirs(output_Dir, exist_ok=True)

    output_new_merged_file = f"new_merged_chr_shortest.fasta"
    sequence_removed = False

    # 使用正则表达式去除序列ID中的控制字符
    clean_sequence_id = re.sub(r'[^\w_\-]', '', sequence_id.strip())

    print(f"处理后的序列ID: {clean_sequence_id}")

    with open(fasta_file, "r") as input_handle, open(output_new_merged_file, "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            record_id = re.sub(r'[^\w_\-]', '', record.id.strip())

            # 调试输出：打印所有记录的ID
            print(f"FASTA文件中的记录ID: '{record_id}'")

            # 检查部分匹配（或者可以使用 startswith 或 in 方法）
            if clean_sequence_id in record_id or record_id in clean_sequence_id:
                # 保存失败的序列到一个单独的文件
                failed_seq_file = os.path.join(output_Dir, f"{record_id}.fasta")
                with open(failed_seq_file, "w") as failed_handle:
                    SeqIO.write(record, failed_handle, "fasta")
                print(f"已保存失败的序列 {record_id} 到 {failed_seq_file}")
                sequence_removed = True
            else:
                # 将剩余的序列写入新的FASTA文件
                SeqIO.write(record, output_handle, "fasta")

    if sequence_removed:
        print(f"新的合并FASTA文件已创建: {output_new_merged_file}")
        return output_new_merged_file
    else:
        print(f"在 {fasta_file} 中未找到序列ID '{clean_sequence_id}'。返回原始文件。")
        return fasta_file


def asm_contig_polish(output_dir, lgsreads, reads1, reads2, threads, NextPolish):
    fasta_files = [f for f in os.listdir(output_dir) if f.endswith('.fasta')]
    for fasta_file in fasta_files:
        reads_name = os.path.splitext(fasta_file)[0]
        file_path = os.path.join(output_dir, fasta_file)
        genome = file_path
        FoFn_dir = os.path.join("Dir_fofn", reads_name)
        os.makedirs(FoFn_dir, exist_ok=True)

        # 只使用输入文件的文件名
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
                # 重命名为 {reads_name}.fasta
                new_name = f"{reads_name}.fasta"
                new_path = os.path.join(final_dir, new_name)
                shutil.move(os.path.join(work_dir, fasta), new_path)
                #shutil.move(os.path.join(work_dir, fasta), final_dir)
        except subprocess.CalledProcessError as e:
            print(f"Error running nextPolish command for {reads_name}: {e}")
            print(f"Stdout: {e.stdout}")
            print(f"Stderr: {e.stderr}")
            print(
            f"{reads_name}.fasta quality is too low, polishing is unsuccessful, and the sequence is directly output to files_NP")
            # 如果命令失败，则将文件复制到 files_NP 目录
            shutil.copy2(file_path, final_dir)
        finally:
            # 删除work_dir和任何备份文件
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


#根据目录中的文件名执行不同的处理
def run_select_execution_function(dir_IN_L, dir_IN_R, lgsreads, wgs1, wgs2, threads, NextPolish):
    """根据目录中的文件名执行不同的处理"""
    one_files = []
    trimmed_files = []

    # 收集文件名
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

    # 如果全是one文件
    if one_files and not trimmed_files:
        asm_L_dir = "asm_L"
        asm_R_dir = "asm_R"
        os.makedirs(asm_L_dir, exist_ok=True)
        os.makedirs(asm_R_dir, exist_ok=True)
        # Extract the assembled directory
        polish_asmmerged_DIR = ['asm_L', 'asm_R']
        for filename in one_files:
            # 根据方向移入asm_L或asm_R
            if 'L' in filename:
                shutil.move(os.path.join(dir_IN_L, filename), os.path.join(polish_asmmerged_DIR[0], filename))
            elif 'R' in filename:
                shutil.move(os.path.join(dir_IN_R, filename), os.path.join(polish_asmmerged_DIR[1], filename))
        output_file = "asm_merged_All_chr.fasta"
        process_asm_merge_fasta_files(polish_asmmerged_DIR, output_file, lgsreads, wgs1, wgs2, threads, NextPolish)

    # 如果同时有one和trimmed文件
    elif one_files and trimmed_files:
        asm_L_dir = "asm_L"
        asm_R_dir = "asm_R"
        os.makedirs(asm_L_dir, exist_ok=True)
        os.makedirs(asm_R_dir, exist_ok=True)
        # 处理one文件
        polish_asmmerged_DIR = ['asm_L', 'asm_R']
        for filename in one_files:
            # 根据方向移入asm_L或asm_R
            if 'L' in filename:
                shutil.move(os.path.join(dir_IN_L, filename), os.path.join(polish_asmmerged_DIR[0], filename))
            elif 'R' in filename:
                shutil.move(os.path.join(dir_IN_R, filename), os.path.join(polish_asmmerged_DIR[1], filename))

        # 处理trimmed文件
        run_alignments_assembly_for_directory(dir_IN_L, dir_IN_R, threads)

    # 如果没有one文件，只有trimmed文件
    elif trimmed_files and not one_files:
        run_alignments_assembly_for_directory(dir_IN_L, dir_IN_R, threads)

    # 如果既没有one文件，也没有trimmed文件
    else:
        print("没有找到符合条件的文件")


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

    ### Step1: Assemble reads into contigs ###
    if args.dir_IN_L and args.dir_IN_R:
        run_select_execution_function(args.dir_IN_L, args.dir_IN_R, args.lgsreads, args.wgs1, args.wgs2, args.threads, args.NextPolish)

    ### Step2: Polish the longest or shortest sequence ###
    # Extract the assembled directory
    polish_asmmerged_DIR = ['asm_L', 'asm_R']
    output_file = "asm_merged_All_chr.fasta"

    process_asm_merge_fasta_files(polish_asmmerged_DIR, output_file, args.lgsreads, args.wgs1, args.wgs2, args.threads, args.NextPolish)

if __name__ == "__main__":
    main()
