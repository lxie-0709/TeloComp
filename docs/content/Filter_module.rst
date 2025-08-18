Filter module
=============

The **Filter module** primarily performs the following tasks: extracting soft-clipped sequences, detecting telomere motifs, extracting reads containing the telomere motifs, and performing pre-assembly processing on the obtained reads. **TeloComp Filter_1** outputs `BAM` files containing soft-clipped sequences that extend beyond the chromosomal ends, for both ONT and HiFi reads. **TeloComp Filter_2** first identifies the main telomere sequence types, displaying the top 10 on the screen and saving the remaining types to a `TXT` file. After the user selects the desired telomere types, Filter2 extracts and outputs the corresponding reads in FASTA format, stored separately in the `ONT` and `HiFi` directories. Finally, the processed data are output to the `trim_L` and `trim_R` directories.

TeloComp Filter_1
-----------------

The first step of `Filter module` is intended to extract soft-clipped sequences located beyond the chromosomal ends of the genome.

.. code:: bash

    # optional arguments:
    #   -h, --help   show this help message and exit
    #   --genome     Input genome FASTA file.
    #   --fai        Input genome index (FAI) file.
    #   --ont        Input ONT data file (optional).
    #   --hifi       Input HiFi data file (optional).
    #   --threads    Number of threads to use with minimap2.
    #   --motifs     A list of telomeric repeat motifs to use for filtering (optional).
    #   --max_break  Maximum tolerable fracture length for soft shear.
    #   --min_clip   Minimum cutting length.
    #   --Ob         BAM output path after ONT filtering.
    #   --Hb         HiFi filtered BAM output path.

    $ telocomp_Filter_1 --genome genome.fasta \
                        --fai genome.fasta.fai \
                        --ont ont.fq.gz \
                        --hifi hifi.fastq.gz \
                        --threads 50 \
                        --Ob ont_out.bam --Hb hifi_out.bam 


TeloComp Filter_2
-----------------

The second step of the `Filter module` is designed to detect, extract, and process reads containing the predefined telomere motifs of interest, starting with the import of the BAM file.

.. code:: bash

    # optional arguments:
    #    -h, --help       show this help message and exit
    #    --ont_bam        Input ONT BAM
    #    --hifi_bam       Input HiFi BAM
    #    -o, --out_dir    Output directory
    #    -c, --coverage   The coverage parameter ranges from 0 to 100 and is used to trim reads
    #                     according to the selected coverage level
    #    -p, --parallels  Parameter for parallel processing of reads, with a default value of 5
    #    --min_ratio      The proportion of the original genome sequence to the length of the
    #                     reads, default=0.2

    $ telocomp_Filter_2 --ont_bam ont_out.bam \
                        --hifi_bam hifi_out.bam \
                        -o output_dir/ \
                        -c 100 -p 10 --min_ratio 0.2




