# TeloComp
TeloComp is an efficient integrated software package for telomere extraction and complementation. It outputs new genome and telomere complementation information, and visualizes the complemented telomere parts with line graphs and colinearity graphs. It is more friendly to scientific researchers and is committed to more complete T2T genome assembly.
<div align="center">
<img src="https://github.com/lxie-0709/TeloComp/blob/1.0.0/example/TeloComp.png" width="688px">
</div>

# Install
TeloComp is an executable file written in `Python(3.11.6)` that users can run directly, but before using TeloComp, you need to install all the dependencies of the software. 
## Dependencies
Please note that you must install the following versions of dependent software or higher before running:

* [samtools-1.18](https://github.com/lxie-0709/TeloComp/blob/main/Dependencies/samtools-1.18.tar.gz)
* [minimap2-2.27](https://github.com/lxie-0709/TeloComp/blob/main/Dependencies/minimap2-2.27.tar.gz)
* [bwa-0.7.17](https://github.com/lxie-0709/TeloComp/blob/main/Dependencies/bwa-0.7.17.tar.bz2)
* [Flye-2.9.4](https://github.com/lxie-0709/TeloComp/blob/main/Dependencies/Flye-2.9.4.tar.gz)
* [racon-v1.4.3](https://github.com/isovic/racon/releases/download/1.4.3/racon-v1.4.3.tar.gz)
* [NextPolish-v1.4.1](https://github.com/lxie-0709/TeloComp/blob/main/Dependencies/NextPolish-v1.4.1.tar.gz)
* teloclip-0.0.4
* GenomeSyn-1.2.7
    
The above software can be installed using conda, downloaded and installed from its github, or installed by running TeloComp's `install.sh`.
In addition, GenomeSyn needs to check its [github](https://github.com/JM-SONG/GenomeSyn) and install it according to the software introduction.

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
    
      telocomp_Filter -h
    
      telocomp_Assembly -h
    
      telocomp_maxmin -h
    
      telocomp_Complement -h

      telocomp_Collinearity -h

# Usage
Note: TeloComp needs to be run in the same directory from beginning to end!

## Filter

#### Options:
      -h, --help        show this help message and exit
      -G , --genome     Input genome file (FASTA format)
      -O , --ONT        Input ONT data
      -H , --HiFi       Input HiFi data
      -B , --ONTbam     The output bam file. If the [--BamExtr] parameter is
                        selected, it is the input bam file.
      -b , --HiFibam    The output bam file. If the [--BamExtr] parameter is
                        selected, it is the input bam file.
      --BamExtr         Selecting this parameter does not perform genome
                        alignment to obtain bam, but directly inputs the sorted
                        bam file, and then screens the qualified reads.
      -r , --ref        Index file of reference genome
      -c , --coverage   Choose the coverage you think is appropriate, that is,
                        the minimum length of reads
      -m , --motif      Telomeric repeats sequences, e.g., plant:
                        CCCTAAA(TTTAGGG), animal: TTAGGG(CCCTAA), etc.
      -t , --threads    Number of threads to use (default: 20)

#### Run:
    Direct input of genomic data：
    telocomp_Filter -G /PATH/test_sequence.fasta.gz -O /PATH/test_ONT.fq.gz -H /PATH/test_HiFi.fq.gz -B out_ONT.bam -b out_HiFi.bam -r /PATH/test_sequence.fasta.fai -c 100 -m CCCTAAA -t 50
    
    Start by importing the bam file（Here, use this method to test and run）：
    telocomp_Filter --BamExtr -B /PATH/test_ONT.bam -b /PATH/test_HiFi.bam -r /PATH/test_sequence.fasta.fai -c 100 -m CCCTAAA


First,this step mainly screens out reads containing telomeres beyond the end of the genome, trims reads according to coverage, and outputs the final results to the directories `trim_L` and `trim_R` according to the direction.

## Assembly

#### Options:
      --dir_IN_L          Directory containing left-aligned reads (FASTA format)
      --dir_IN_R          Directory containing right-aligned reads (FASTA
                          format)
      -L , --lgsreads     Long-read sequencing data
      -W , --wgs1         Path to WGS reads (read 1)
      -w , --wgs2         Path to WGS reads (read 2)
      -N , --NextPolish   Path to NextPolish tool
      -t , --threads      Number of threads to use (default: 20)

#### Run:
     telocomp_Assembly --dir_IN_L trim_L --dir_IN_R trim_R -L /PATH/test_HiFi.fq.gz -W /PATH/test_WGS_f1.fq.gz -w /PATH/test_WGS_r2.fq.gz -N /PATH/NextPolish -t 50 

Next,the screened and processed reads are assembled and polished, and the final results are output to the directory `fils_NP`.

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

#### Run:
##### Extract reads
    （1）Extract the longest reads
    telocomp_maxmin --Max_length --dir_ont /PATH/algn_output_ont --dir_hifi /PATH/algn_output_hifi -L /PATH/test_HiFi.fq.gz -W /PATH/test_WGS_f1.fq.gz -w /PATH/test_WGS_r2.fq.gz -N /PATH/NextPolish -t 50 


    （2）Extract the shortest reads
    telocomp_maxmin --Min_length --dir_ont /PATH/algn_output_ont --dir_hifi /PATH/algn_output_hifi -L /PATH/test_HiFi.fq.gz -W /PATH/test_WGS_f1.fq.gz -w /PATH/test_WGS_r2.fq.gz -N /PATH/NextPolish -t 50

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

Finally，complement the assembled and polished reads to the original genome, then output the complemented genome(`new_genome.fasta`), and output the telomere position file(`telomere position`) and telomere type file(`telomere_repeats_info.txt`), and the density distribution map of telomeres at each chromosome end(stored in the folder `telomere_plots`).

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








