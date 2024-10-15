# TeloComp
TeloComp is an efficient integrated software package for telomere extraction and complementation. It outputs new genome and telomere complementation information, and visualizes the complemented telomere parts with line graphs and colinearity graphs. It is more friendly to scientific researchers and is committed to more complete T2T genome assembly.

![image](https://github.com/lxie-0709/TeloComp/blob/main/Synteny.png)

# Install
TeloComp is an executable file written in Python that users can run directly, but before using TeloComp, you need to install all the dependencies of the software. You can download the current source distribution here or visit Github.

    git clone git@github.com:lxie-0709/TeloComp.git

# Dependencies
Please note that you must install the following versions of dependent software or higher before running:

    * samtools-1.18
    * minimap2-2.27
    * bwa-0.7.12
    * Flye-2.9.3
    * racon-v1.4.3
    * NextPolish-v1.4.1
    
The above software can be installed using conda(install_conda.sh) or by downloading the software source code from github(install.sh)

3.GenomeSyn
GenomeSyn cannot be installed using conda, you need to refer to its [github](https://github.com/jmsong2/GenomeSyn) tutorial for installation.


#Usage
    
    screen -L -dmS step1_Cassava bash -c "/usr/bin/time -v python step1_Cassava20240809.py -G ../GWHDEDE00000000.genome.fasta -O ../CRR800583_ont.fq.gz -H ../CRR780166_hifi.fastq.gz -B out_ONT_Casssava.bam -b out_HiFi_Casssava.bam -r ../GWHDEDE00000000.genome.fasta.fai -c 100 -m CCCTAAA -t 50"

    nohup /usr/bin/time -v python step2_Cassava20240819.py --dir_IN_L trim_L --dir_IN_R trim_R -L /pub5/Hshoub/ZH13_T2T/CRR70524_merged_hifi.fq.gz -W /pub5/Hshoub/ZH13_T2T/WGS/CRR705250_f1.fq.gz -w /pub5/Hshoub/ZH13_T2T/WGS/CRR705250_r2.fq.gz -N /pub4/huangshoubian1/NextPolish -t 50 &

    nohup /usr/bin/time -v python step3_Vigna_20240822.py --Normal -G ../../GWHDEDE00000000.genome.fasta --dir_contigs ../files_NP --dir_trim_L trim_L --dir_trim_R trim_R -L /pub5/Hshoub/Cassava/CRR780166_hifi.fastq.gz -W /pub5/Hshoub/Cassava/WGS/CRR780168_f1.fq.gz -w /pub5/Hshoub/Cassava/WGS/CRR780168_r2.fq.gz -t 20 -N /pub4/huangshoubian1/NextPolish -m CCCTAAA -M 7 &

    python step4_genomeSyn_20240830_color.py -G1 GCA_028455895.1_ASM2845589v1_genomic_upseq.fna -G2 new_genome_upseq.fasta -pos ../telomere_positions.txt --length 20000 -tel_max 10 -m CCCTAAA --genomeSyn2

    








