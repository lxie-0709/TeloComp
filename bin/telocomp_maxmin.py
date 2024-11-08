#20240807 从最初的数据中提取最长或者最短的reads，只要小写部分，然后先合并racon，最后分开单独打磨
import os
import shutil
import argparse
import subprocess
from Bio import SeqIO


# import concurrent.futures


# 合并在ONT或HiFi中都存在的L、R的fasta文件
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


# 提取最长的小写序列作为contig
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


# 提取最短的小写序列作为contig
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


#def Max_contig_racon(merged_file, lgsreads, wgs1, wgs2, threads, NextPolish):
#    tmpdir = "tmp_dir"
#    os.makedirs(tmpdir, exist_ok=True)
#
#    racon_command = (
#        f"minimap2 -ax map-pb -t {threads} {merged_file} {lgsreads} | samtools sort - -m 2g --threads 20 -o {tmpdir}/genome.lgs.bam",
#        f"samtools index {tmpdir}/genome.lgs.bam",
#        f"ls `pwd`/{tmpdir}/genome.lgs.bam > {tmpdir}/pb.map.bam.fofn",
#        f"python {NextPolish}/lib/nextpolish2.py -g {merged_file} -l {tmpdir}/pb.map.bam.fofn -r clr -p {threads} -a -s -o {tmpdir}/merged_genome.lgspolish.fasta"
#    )
#
#    print(f"Processing operations before nextpolish")
#
#    try:
#        for command in racon_command:
#            subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
#    except subprocess.CalledProcessError as e:
#        print(f"Error running Max_polish_command: {e}")
#        print(f"Stdout: {e.stdout}")
#        print(f"Stderr: {e.stderr}")
#        return None
#20240822 修改racon逻辑
def Max_contig_racon(merged_file, lgsreads, wgs1, wgs2, threads, NextPolish):
    tmpdir = "tmp_dir"
    os.makedirs(tmpdir, exist_ok=True)

    while True:  # 使用循环而非递归
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

    # 将总的 fasta 文件分成单个 fasta 文件
    racon_genome = f"{tmpdir}/merged_genome.lgspolish.fasta"
    output_dir = "MaxLength_NP"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    for record in SeqIO.parse(racon_genome, "fasta"):
        base_name = record.id
        # chr_name = base_name.split('_')[0]
        # extr_type = base_name.split('_')[1]
        # direction = base_name.split('_')[3]
        # output_file = os.path.join(output_dir, f"{chr_name}_{extr_type}_np_{direction}.fasta")
        output_file = os.path.join(output_dir, f"{base_name}.fasta")
        with open(output_file, 'w') as out_handle:
            SeqIO.write(record, out_handle, "fasta")

    Max_contig_polish(output_dir, lgsreads, wgs1, wgs2, threads, NextPolish)

    # 删除临时目录和副本
    shutil.rmtree(tmpdir)
    shutil.rmtree(output_dir)


def Max_remove_sequence_from_fasta(fasta_file, sequence_id):
    output_Dir = "tmp_getpos"
    os.makedirs(output_Dir, exist_ok=True)

    output_new_merged_file = os.path.join(output_Dir, "new_merged_chr_shortest.fasta")
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


def Max_contig_polish(output_dir, lgsreads, wgs1, wgs2, threads, NextPolish):
    fasta_files = [f for f in os.listdir(output_dir) if f.endswith('.fasta')]
    for fasta_file in fasta_files:
        reads_name = os.path.splitext(fasta_file)[0]
        file_path = os.path.join(output_dir, fasta_file)
        genome = file_path
        FoFn_dir = os.path.join("Dir_Max_fofn", reads_name)
        os.makedirs(FoFn_dir, exist_ok=True)

        # 只使用输入文件的文件名
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

        final_dir = os.path.join(current_dir, "tmp_Max_getpos")
        os.makedirs(final_dir, exist_ok=True)

        with open(sgs_fofn, 'w') as sgs_file:
            sgs_file.write(f"{wgs1}\n")
            sgs_file.write(f"{wgs2}\n")

        with open(lgs_fofn, 'w') as lgs_file:
            lgs_file.write(f"{lgsreads}\n")

        # 生成配置文件
        # {t}为每一个任务的线程，默认是5；{genome}为初始基因组文件；
        # {workdir}为工作的输出路径，一般与单独不同的一个；{sgs_fofn}为含有二代reads路径的文本，每一行文件；
        # {lgs_fofn}为含有长reads路径的文本文件
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
                # 重命名为 {reads_name}.fasta
                new_name = f"{reads_name}.fasta"
                new_path = os.path.join(final_dir, new_name)
                shutil.move(os.path.join(work_dir, fasta), new_path)
                # shutil.move(os.path.join(work_dir, fasta), final_dir)
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


#def Min_contig_racon(merged_file, lgsreads, wgs1, wgs2, threads, NextPolish):
#    tmpdir = "tmp_dir"
#    os.makedirs(tmpdir, exist_ok=True)
#
#    racon_command = (
#        f"minimap2 -ax map-pb -t {threads} {merged_file} {lgsreads} | samtools sort - -m 2g --threads 20 -o {tmpdir}/genome.lgs.bam",
#        f"samtools index {tmpdir}/genome.lgs.bam",
#        f"ls `pwd`/{tmpdir}/genome.lgs.bam > {tmpdir}/pb.map.bam.fofn",
#        f"python {NextPolish}/lib/nextpolish2.py -g {merged_file} -l {tmpdir}/pb.map.bam.fofn -r clr -p {threads} -a -s -o {tmpdir}/merged_genome.lgspolish.fasta"
#    )
#
#    print(f"Processing operations before nextpolish")
#
#    try:
#        for command in racon_command:
#            subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
#    except subprocess.CalledProcessError as e:
#        print(f"Error running Max_polish_command: {e}")
#        print(f"Stdout: {e.stdout}")
#        print(f"Stderr: {e.stderr}")
#        return None
#20240822 修改racon逻辑
def Min_contig_racon(merged_file, lgsreads, wgs1, wgs2, threads, NextPolish):
    tmpdir = "tmp_dir"
    os.makedirs(tmpdir, exist_ok=True)

    while True:  # 使用循环而非递归
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

    # 将总的 fasta 文件分成单个 fasta 文件
    racon_genome = f"{tmpdir}/merged_genome.lgspolish.fasta"
    output_dir = "MinLength_NP"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    for record in SeqIO.parse(racon_genome, "fasta"):
        base_name = record.id
        # chr_name = base_name.split('_')[0]
        # extr_type = base_name.split('_')[1]
        # direction = base_name.split('_')[3]
        # output_file = os.path.join(output_dir, f"{chr_name}_{extr_type}_np_{direction}.fasta")
        output_file = os.path.join(output_dir, f"{base_name}.fasta")
        with open(output_file, 'w') as out_handle:
            SeqIO.write(record, out_handle, "fasta")

    Min_contig_polish(output_dir, lgsreads, wgs1, wgs2, threads, NextPolish)

    # 删除临时目录和副本
    shutil.rmtree(tmpdir)
    shutil.rmtree(output_dir)


def Min_remove_sequence_from_fasta(fasta_file, sequence_id):
    output_Dir = "tmp_getpos"
    os.makedirs(output_Dir, exist_ok=True)

    output_new_merged_file = os.path.join(output_Dir, "new_merged_chr_shortest.fasta")
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


def Min_contig_polish(output_dir, lgsreads, wgs1, wgs2, threads, NextPolish):
    fasta_files = [f for f in os.listdir(output_dir) if f.endswith('.fasta')]
    for fasta_file in fasta_files:
        reads_name = os.path.splitext(fasta_file)[0]
        file_path = os.path.join(output_dir, fasta_file)
        genome = file_path
        FoFn_dir = os.path.join("Dir_Min_fofn", reads_name)
        os.makedirs(FoFn_dir, exist_ok=True)

        # 只使用输入文件的文件名
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

        final_dir = os.path.join(current_dir, "tmp_Min_getpos")
        os.makedirs(final_dir, exist_ok=True)

        with open(sgs_fofn, 'w') as sgs_file:
            sgs_file.write(f"{wgs1}\n")
            sgs_file.write(f"{wgs2}\n")

        with open(lgs_fofn, 'w') as lgs_file:
            lgs_file.write(f"{lgsreads}\n")

        # 生成配置文件
        # {t}为每一个任务的线程，默认是5；{genome}为初始基因组文件；
        # {workdir}为工作的输出路径，一般与单独不同的一个；{sgs_fofn}为含有二代reads路径的文本，每一行文件；
        # {lgs_fofn}为含有长reads路径的文本文件
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
                # 重命名为 {reads_name}.fasta
                new_name = f"{reads_name}.fasta"
                new_path = os.path.join(final_dir, new_name)
                shutil.move(os.path.join(work_dir, fasta), new_path)
                # shutil.move(os.path.join(work_dir, fasta), final_dir)
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


def process_extract_Max_Min(chromosome, extr_merged_files, input_dir_ont, input_dir_hifi, direction, max_length,
                            min_length, temp_files):
    ont_filename = f"{chromosome}_ONT_{direction}.fasta"
    hifi_filename = f"{chromosome}_HiFi_{direction}.fasta"
    input_file = None
    seq_type = None

    if ont_filename in extr_merged_files and hifi_filename in extr_merged_files:
        seq_type = "merged"
        ont_file = os.path.join(input_dir_ont, ont_filename)
        hifi_file = os.path.join(input_dir_hifi, hifi_filename)
        merged_file = f"{chromosome}_merged_{direction}.fasta"
        merge_sequences(ont_file, hifi_file, merged_file)
        input_file = merged_file
        temp_files.append(input_file)
    elif ont_filename in extr_merged_files:
        seq_type = "ONT"
        input_file = os.path.join(input_dir_ont, ont_filename)
        # temp_files.append(input_file)
    elif hifi_filename in extr_merged_files:
        seq_type = "HiFi"
        input_file = os.path.join(input_dir_hifi, hifi_filename)
        # temp_files.append(input_file)

    if input_file:
        if max_length:
            max_length_extract(input_file, chromosome, direction, seq_type)
        if min_length:
            min_length_extract(input_file, chromosome, direction, seq_type)


### Step 2: Merge and Trim Sequences ###
def process_Max_merge_fasta_files(directories, output_file, lgsreads, wgs1, wgs2, threads, NextPolish):
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
    Max_contig_racon(merged_file, lgsreads, wgs1, wgs2, threads, NextPolish)

    # Remove the file
    os.remove(merged_file)


def process_Min_merge_fasta_files(directories, output_file, lgsreads, wgs1, wgs2, threads, NextPolish):
    unique_sequence = set()
    with open(output_file, 'w') as outfile:
        for directory in directories:
            fasta_files = [f for f in os.listdir(directory) if f.endswith('.fasta')]
            for fasta_file in fasta_files:
                chr_name = os.path.basename(fasta_file).split('_')[0]
                directorie_min = os.path.basename(fasta_file).split('_')[-1].split('.')[0]
                file_path = os.path.join(directory, fasta_file)
                with open(file_path, 'r') as infile:
                    # for record in SeqIO.parse(infile, 'fasta'):
                    #     record.id = os.path.splitext(fasta_file)[0]
                    #     record.description = ''
                    #     SeqIO.write(record, outfile, 'fasta')
                    for record in SeqIO.parse(infile, 'fasta'):
                        seq_id = f"{chr_name}_shortest_np_{directorie_min}"
                        if record.seq and seq_id not in unique_sequence:
                            unique_sequence.add(seq_id)
                            outfile.write(f'>{seq_id}\n{record.seq}\n')

    merged_file = output_file
    Min_contig_racon(merged_file, lgsreads, wgs1, wgs2, threads, NextPolish)

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
    parser.add_argument('-L', '--lgsreads', metavar='', required=True, help='Long-read sequencing data')
    parser.add_argument('-W', '--wgs1', metavar='', required=True, help='Path to WGS reads (read 1)')
    parser.add_argument('-w', '--wgs2', metavar='', required=True, help='Path to WGS reads (read 2)')
    parser.add_argument('-N', '--NextPolish', metavar='', required=True, help='Path to NextPolish tool')
    parser.add_argument('-t', '--threads', type=int, default=20, metavar='',
                        help='Number of threads to use (default: 20)')

    args = parser.parse_args()

    # Step1: Extract longest or shortest sequence
    input_dir_ont = args.dir_ont
    input_dir_hifi = args.dir_hifi

    extr_ont_files = [f for f in sorted(os.listdir(input_dir_ont)) if f.endswith('.fasta')]
    extr_hifi_files = [f for f in sorted(os.listdir(input_dir_hifi)) if f.endswith('.fasta')]
    extr_merged_files = extr_ont_files + extr_hifi_files

    processed_chromosomes = set()
    temp_files = []

    for chromosome_file in extr_merged_files:
        if not chromosome_file.startswith('chr'):
            continue
        chromosome = chromosome_file.split('_')[0]
        if chromosome in processed_chromosomes:
            continue

        process_extract_Max_Min(chromosome, extr_merged_files, input_dir_ont, input_dir_hifi, 'L', args.Max_length,
                                args.Min_length, temp_files)
        process_extract_Max_Min(chromosome, extr_merged_files, input_dir_ont, input_dir_hifi, 'R', args.Max_length,
                                args.Min_length, temp_files)
        processed_chromosomes.add(chromosome)

    # Step2: Polish the longest or shortest sequence

    if args.Max_length:
        # Extract the directory of the longest sequence output
        polish_Maxmerged_DIR = ['MaxLength_L', 'MaxLength_R']
        output_file = f"merged_Max_chr.fasta"

        process_Max_merge_fasta_files(polish_Maxmerged_DIR, output_file, args.lgsreads, args.wgs1, args.wgs2,
                                      args.threads, args.NextPolish)

    elif args.Min_length:
        # Extract the directory of the shortest sequence output
        polish_Minmerged_DIR = ['MinLength_L', 'MinLength_R']
        output_file = f"merged_Min_chr.fasta"

        process_Min_merge_fasta_files(polish_Minmerged_DIR, output_file, args.lgsreads, args.wgs1, args.wgs2,
                                      args.threads, args.NextPolish)

    # Remove temporary files after processing
    for temp_file in temp_files:
        if os.path.exists(temp_file):
            os.remove(temp_file)


if __name__ == "__main__":
    main()
    