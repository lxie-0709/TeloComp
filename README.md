# TeloComp
TeloComp is an efficient integrated software package for telomere extraction and complementation. It finalizes the output of the new genome and telomere complementation information and visualizes the complemented telomere portion through line graphs and covariance plots. It is more friendly to researchers and works towards a more complete T2T genome assembly.
<div align="center">
<img src="https://github.com/lxie-0709/TeloComp/blob/1.0.0/example/TeloComp.png" width="688px">
</div>

# Install
TeloComp is an executable program written in **Python (3.11.6)** that can be run directly by the user, but all dependencies of the program need to be installed before using TeloComp.
## Dependencies
Please note that you must install the following versions of dependent software or higher before running:

* [samtools-1.18](https://github.com/lxie-0709/TeloComp/blob/main/Dependencies/samtools-1.18.tar.gz)
* [minimap2-2.27](https://github.com/lxie-0709/TeloComp/blob/main/Dependencies/minimap2-2.27.tar.gz)
* [bwa-0.7.17](https://github.com/lxie-0709/TeloComp/blob/main/Dependencies/bwa-0.7.17.tar.bz2)
* [Flye-2.9.4](https://github.com/lxie-0709/TeloComp/blob/main/Dependencies/Flye-2.9.4.tar.gz)
* [pilon-1.24](https://github.com/lxie-0709/TeloComp/blob/main/Dependencies/pilon-1.24.tar.gz)
* [NextPolish-v1.4.1](https://github.com/lxie-0709/TeloComp/blob/main/Dependencies/NextPolish-v1.4.1.tar.gz)
* [GenomeSyn-1.2.7](https://github.com/JM-SONG/GenomeSyn/tree/main)
    
The above software can be installed using conda, downloaded and installed from its github, or by running TeloComp's `install.sh`.

## Building on Linux
To use the software, you need to follow the following steps to install it.

#### &emsp;1.First, get the source code.

    git clone git@github.com:lxie-0709/TeloComp.git
    
    cd TeloComp

#### &emsp;2.Next, execute `install.sh` and `setup.sh` in `Dependencies` and `bin` respectively “ to install the software dependencies and configure the software. 
    （1）Installing dependencies
        $ sh install.sh 

    （2）Configuring TeloComp
        $ sh setup.sh

#### &emsp;3.Download [GenomeSyn](https://github.com/JM-SONG/GenomeSyn/archive/refs/heads/main.zip), add the GenomeSyn-1.2.7 file to your root directory, and perform the following steps：
        $ chmod -R 777 GenomeSyn-1.2.7
        $ echo "export PATH=$PATH:/yourPATH/GenomeSyn-1.2.7/bin" >> ~/.bashrc

#### &emsp;3.Finally, activate the environment variable and then verify that it is properly installed and executable with the following command:

    # Activation environment
      source ~/.bashrc
    
      telocomp_Filter_1 -h

      telocomp_Filter_2 -h
    
      telocomp_Assembly -h
    
      telocomp_maxmin -h
    
      telocomp_Complement -h

      telocomp_Collinearity -h

# Usage
Note: TeloComp requires that you run the telomere complement command in the same directory from start to finish!

## Filter

#### Options_1:
      -h, --help   show this help message and exit
      --genome     Input genome FASTA file.
      --fai        Input genome index (FAI) file.
      --ont        Input ONT data file (optional).
      --hifi       Input HiFi data file (optional).
      --threads    Number of threads to use with minimap2.
      --motifs     A list of telomeric repeat motifs to use for filtering (optional).
      --max_break  Maximum tolerable fracture length for soft shear.
      --min_clip   Minimum cutting length.
      --Ob         BAM output path after ONT filtering.
      --Hb         HiFi filtered BAM output path.

#### Options_2:
      -h, --help       show this help message and exit
      --ont_bam        ONT BAM
      --hifi_bam       HiFi BAM
      -o, --out_dir    Directory where both step1_2 and step1_3 should write their outputs
      -c, --coverage   The coverage parameter ranges from 0 to 100 and is used to trim reads according to the selected                           coverage level.

      -p, --parallels  Parameter for parallel processing of reads, with a default value of 4.
      --min_ratio      The proportion of the original genome sequence to the length of the reads, default=0.2

#### Run:
（1）Get the bam file containing the end software cutting sequence

    telocomp_Filter_1 --genome /PATH/genome.fasta --fai /PATH/genome.fasta.fai --ont /PATH/ont.fastq.gz --hifi /PATH/hifi.fastq.gz --threads 50 --Ob /PATH/ont_out.bam --Hb /PATH/hifi_out.bam

（2）Detection, extraction, and processing of reads.Start by importing the bam file（Here, run the test using this procedure.）：
    
    telocomp_Filter_2 --ont_bam /PATH/ont_out.bam --hifi_bam /PATH/hifi_out.bam -o PATH/output_dir/ -c 100 -p 10 --min_ratio 0.2




First, this step mainly detects and filters out reads containing telomeres outside the ends of the genome, trims the reads according to the coverage, and finally outputs the final results to the `trim_L` and `trim_R` directories according to the direction.

## Assembly

#### Options:
      -h, --help          show this help message and exit
      --dir_IN_L          Directory containing left-aligned reads (FASTA format)
      --dir_IN_R          Directory containing right-aligned reads (FASTA format)
      --flye              Flye assembly module
      -a, --assemble      Alternative assemble module
      -t , --threads      Threads (default:20)
      --min_overlap       Min overlap (default:50)
      --error_rate        Error rate (default:0.15)
      --kmer_size         K-mer size (default:15)
      -L , --lgsreads     Long-read sequencing data
      -W , --wgs1         Path to WGS reads (read 1)
      -w , --wgs2         Path to WGS reads (read 2)
      -N , --NextPolish   Path to NextPolish tool
      --jar               Path to jar in Pilon



#### Run:
    （1）Flye assembly module (default assembly)
     telocomp_Assembly --dir_IN_L trim_L --dir_IN_R trim_R -L /PATH/test_HiFi.fq.gz -W /PATH/test_WGS_f1.fq.gz -w /PATH/test_WGS_r2.fq.gz -N /PATH/NextPolish -t 50 
    （2）assemble assembly module
     telocomp_Assembly --dir_IN_L trim_L --dir_IN_R trim_R -L /PATH/test_HiFi.fq.gz -W /PATH/test_WGS_f1.fq.gz -w /PATH/test_WGS_r2.fq.gz -N /PATH/NextPolish -t 50 --assemble
    
Next,the screened and processed reads are assembled and polished, and the final results are output to the directory `files_NP`.

## Extract the longest or shortest reads

If you choose to directly extract the longest or shortest reads, you can skip the assembly step and run this step directly.

#### Options:

      -h, --help          show this help message and exit
      --Max_length        Extract longest reads
      --Min_length        Extract shortest reads
      --dir_ont           Directory containing ONT files
      --dir_hifi          Directory containing HiFi files
      -L , --lgsreads     Long-read sequencing data
      -W , --wgs1         Path to WGS reads (read 1)
      -w , --wgs2         Path to WGS reads (read 2)
      -N , --NextPolish   Path to NextPolish tool
      -t , --threads      Number of threads to use (default: 20)
      --polish            Perform polishing with NextPolish

#### Run:
##### Extract reads
    （1）Extract the longest reads
    telocomp_maxmin --Max_length --dir_ont /PATH/algn_output_ont --dir_hifi /PATH/algn_output_hifi -L /PATH/test_HiFi.fq.gz -W /PATH/test_WGS_f1.fq.gz -w /PATH/test_WGS_r2.fq.gz -N /PATH/NextPolish -t 50 --polish


    （2）Extract the shortest reads
    telocomp_maxmin --Min_length --dir_ont /PATH/algn_output_ont --dir_hifi /PATH/algn_output_hifi -L /PATH/test_HiFi.fq.gz -W /PATH/test_WGS_f1.fq.gz -w /PATH/test_WGS_r2.fq.gz -N /PATH/NextPolish -t 50 --polish

    （3）Without polishing, run the following code
    telocomp_maxmin --Min_length --dir_ont /PATH/algn_output_ont --dir_hifi /PATH/algn_output_hifi
    telocomp_maxmin --Max_length --dir_ont /PATH/algn_output_ont --dir_hifi /PATH/algn_output_hifi

Here you need to enter the untrimmed end alignment reads in the Filter, `algn_output_ont` and `algn_output_hifi` respectively, and finally output the polished reads to the directory `tmp_Min_getpos` and `tmp_Max_getpos`.

##### （2）Telomere complement 
    telocomp_Complement --dir_Min -G /PATH/test_sequence.fasta -m CCCTAAA -M 7 

This is the same as the Telomere complement below, both of which complete the telomere part to the original genome, but the running command is different.


## Telomere complement

#### Option:
      -G , --genome       Input genome file (FASTA format)
      --dir_contigs       Input polished contigs
      --dir_trim_L        Input the trimmed reads. If the conditions are not
                          met, extract the shortest reads.
      --dir_trim_R        Input the trimmed reads. If the conditions are not
                          met, extract the shortest reads.
      -L , --lgsreads     Long-read sequencing data
      -W , --wgs1         Path to WGS reads (read 1)
      -w , --wgs2         Path to WGS reads (read 2)
      -N , --NextPolish   Path to NextPolish tool
      -m , --motif        Telomeric repeats sequences, e.g., plant:
                          CCCTAAA(TTTAGGG), animal: TTAGGG(CCCTAA), etc.
      -M , --motif_num    Input the number of bases of the telomere motif
      --Normal            Execute the command according to the general process
      --dir_Max           Select the telomere reads obtained by polishing the
                          longest reads to add to the genome
      --dir_Min           Select the telomere reads obtained by polishing the
                          shortest reads to add to the genome
      -t , --threads      Number of threads to use (default: 20)
      
#### Run:
    telocomp_Complement  --Normal -G /PATH/test_sequence.fasta --dir_contigs /PATH/files_NP --dir_trim_L trim_L --dir_trim_R trim_R -L /PATH/test_HiFi.fq.gz -W /PATH/test_WGS_f1.fq.gz -w /PATH/test_WGS_r2.fq.gz -t 50 -N /PATH/NextPolish -m CCCTAAA -M 7 

Finally，complement the assembled and polished reads to the original genome, then output the complemented genome(`new_genome.fasta`), and output the telomere position file(`telomere_positions.txt`) and telomere type file(`telomere_repeats_info.txt`), and the density distribution map of telomeres at each chromosome end(stored in the folder `telomere_plots`).

####  Example: the following pictures show the left and right ends of chromosome 2 of *Morus* *notabilis*.
Telomere density distribution diagram of chromosome ends with telomere complementation（The following pictures show the left and right ends of chromosome 2).
<div align="center">
    <img src="https://github.com/lxie-0709/TeloComp/blob/1.0.0/example/test_L_telomere.png" width="388px"/>
    <img src="https://github.com/lxie-0709/TeloComp/blob/1.0.0/example/test_R_telomere.png" width="388px"/>
</div>


## Collinearity analysis

#### Option:
      -G1 , --genome1       Input the original genome file (FASTA format)
      -G2 , --genome2       Input new genome file (FASTA format)
      -pos                  The position information of telomere installed on
                            the original genome
      -m , --motif          Telomeric repeats sequences, e.g., plant:
                            CCCTAAA(TTTAGGG), animal: TTAGGG(CCCTAA), etc.
      --bin                 The bins for calculating the number of telomeres are
                            divided into. The bin size can be determined by the
                            user, and the default value is 100.
      --genomeSyn1          GenomeSyn1 is in the order of original genome first
                            and new genome second
      --genomeSyn2          genomeSyn2 is in the order of original genome second
                            and new genome first
      --length              Extract genome length
      -tel_max              The maximum value of telomere frequency within the
                            calculated length
      -C1 , --genome_color1 
                            Set the color parameters of genome 1, the default is
                            #04686b, or select other hexadecimal or RGB color
                            codes, such as "#04686b"/"rgb(4,104,107)"
      -C2 , --genome_color2 
                            Set the color parameters of genome 2, the default is
                            #87c9c3, or select other hexadecimal or RGB color
                            codes, such as "#87c9c3"/"rgb(135,201,195)"
      -C3 , --synteny_color 
                            Set the color parameters of the synteny blocks, the
                            default is #DFDFE1, or select other hexadecimal or
                            RGB color codes, such as
                            "#DFDFE1"/"rgb(223,223,225)"
      -C4 , --telo_color    Set the telomere color parameters for displaying
                            genomes, the default is #fcaf7c, or select other
                            hexadecimal or RGB color codes, such as
                            "#fcaf7c"/"rgb(252,175,124)"
      -t , --threads        Number of threads to use (default: 20)

#### Run:
    telocomp_Collinearity -G1 /PATH/test_sequence.fasta -G2 new_genome.fasta -pos /PATH/telomere_positions.txt --length 20000 -tel_max 10 -m CCCTAAA --genomeSyn2

In addition, the left and right ends of the chromosomes of the original genome and the complemented genome(`new_genome.fasta`) are extracted for collinear comparison(`genomeSyn_result`), and the number of telomeres at the corresponding chromosome ends of the two genomes is output,including `telomere.original.num.info` and `telomere.complement.num.info`.

####  Collinearity result
Collinear alignment between the ***Morus*** ***notabilis*** original genome and the telomeric complemented genome (only chromosomes with complementary telomeres are shown).
<div align=center>
<img src="https://github.com/lxie-0709/TeloComp/blob/1.0.0/example/test_Collinearity.png" width="888px">
</div>








