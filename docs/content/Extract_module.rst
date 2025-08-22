Extract Module 
==============

Assembly-free Module for Extracting the Longest or Shortest Telomeric Reads
---------------------------------------------------------------------------

The **Extract module** is designed for use cases where direct extraction is required without assembly. It extracts telomeric sequences located beyond the chromosomal ends and containing telomere motifs, and directly integrates the extracted sequences into the original genome.

.. code:: bash

    # optional arguments:
    #  -h, --help          show this help message and exit
    #  --Max_length        Extract longest reads
    #  --Min_length        Extract shortest reads
    #  --dir_ont           Directory containing ONT files
    #  --dir_hifi          Directory containing HiFi files
    #  -L , --lgsreads     Long-read sequencing data
    #  -W , --wgs1         Path to WGS reads (read 1)
    #  -w , --wgs2         Path to WGS reads (read 2)
    #  -N , --NextPolish   Path to NextPolish tool
    #  -t , --threads      Number of threads to use (default: 20)
    # --polish            Perform polishing with NextPolish

    # Parameters of telocomp_Complement
    # --dir_Max           Select the telomere reads obtained by polishing the longest reads to
    #                     add to the genome
    # --dir_Min           Select the telomere reads obtained by polishing the shortest reads
    #                     to add to the genome
    # -m , --motif        Telomeric repeats sequences, e.g., plant: CCCTAAA(TTTAGGG), animal:
    #                     TTAGGG(CCCTAA), etc.
    # -M , --motif_num    Input the number of bases of the telomere motif


Extracting the longest telomeric reads
--------------------------------------

In this step, the longest reads are directly extracted and the results are saved to the ``MaxLength_L`` and ``MaxLength_R`` directories. The sequences used for genome completion are obtained by merging the ``FASTA`` files from ``MaxLength_L`` and ``MaxLength_R`` into ``MaxLength_NP``.

.. code:: bash

    # No polishing
    $ telocomp_maxmin --Max_length \
                      --dir_ont algn_output_ont \
                      --dir_hifi algn_output_hifi \


    # Polishing
    $ telocomp_maxmin --Max_length \
                      --dir_ont algn_output_ont \
                      --dir_hifi algn_output_hifi \
                      -L HiFi.fastq.gz \
                      -W WGS_f1.fq.gz \
                      -w WGS_r2.fq.gz \
                      --polish \
                      -N /PATH/NextPolish -t 50 


Extracting the shortest telomeric reads
---------------------------------------

In this step, the shortest reads are directly extracted and the results are saved to the ``MinLength_L`` and ``MinLength_R`` directories. The sequences used for genome completion are obtained by merging the ``FASTA`` files from ``MinLength_L`` and ``MinLength_R`` into ``MinLength_NP``.

.. code:: bash

    # No polishing
    $ telocomp_maxmin --Min_length \
                      --dir_ont algn_output_ont \
                      --dir_hifi algn_output_hifi \


    # Polishing
    $ telocomp_maxmin --Min_length \
                      --dir_ont algn_output_ont \
                      --dir_hifi algn_output_hifi \
                      -L HiFi.fastq.gz \
                      -W WGS_f1.fq.gz \
                      -w WGS_r2.fq.gz \
                      --polish \
                      -N /PATH/NextPolish -t 50 

Complement part
---------------

The **Complement part** operates on the longest and shortest reads extracted in the previous step, and integrates these reads into the corresponding positions of the genome.

.. code:: bash

    $ telocomp_Complement --dir_Max -G /PATH/test_sequence.fasta -m CCCTAAA -M 7

    $ telocomp_Complement --dir_Min -G /PATH/test_sequence.fasta -m CCCTAAA -M 7 




