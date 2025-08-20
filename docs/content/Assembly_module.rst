Assembly Module  
===============

The **Assembly module** primarily assembles the reads processed by the **Filter module**, and outputs the final assembly results to the ``Files_NP`` directory.
Note: If reads at a chromosome end cannot be assembled, the shortest reads are selected to supplement and integrate into the reference genome.

.. code:: bash

    # optional arguments:
    # -h, --help          show this help message and exit
    # --dir_IN_L          Directory containing left-aligned reads (FASTA format)
    # --dir_IN_R          Directory containing right-aligned reads (FASTA format)
    # --flye              Flye assembly module (default: True)
    # --assemble, -a      Assemble using an alternative assembly module
    # -L , --lgsreads     Long-read sequencing data
    # -W , --wgs1         Path to WGS reads (read 1)
    # -w , --wgs2         Path to WGS reads (read 2)
    # -N , --NextPolish   Path to NextPolish tool
    # -t , --threads      Number of threads to use (default: 20)

    # The Assemble module can be adjusted using the following parameters
    # --min_overlap       Minimum overlap length (default: 50)
    # --error_rate        Error rate for assembly (default: 0.15)
    # --kmer_size         K-mer size (default: 15)


    # Flye assembly module (default assembly)
    $ telocomp_Assembly --dir_IN_L trim_L \
                        --dir_IN_R trim_R \
                        -L HiFi.fq.gz \
                        -W WGS_f1.fq.gz \
                        -w WGS_r2.fq.gz \
                        -N /PATH/NextPolish -t 50 

    # assemble assembly module
    $ telocomp_Assembly --dir_IN_L trim_L \
                        --dir_IN_R trim_R \
                        -L HiFi.fq.gz \
                        -W WGS_f1.fq.gz \
                        -w WGS_r2.fq.gz \
                        -N /PATH/NextPolish \
                        -t 50 --assemble

