Filter module
=============

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
                        --Ob ont_out.bam --Hb hifi_out.bam \


TeloComp Filter_1
-----------------

To use the software, you need to follow the following steps to install it.

1.Obtain software package from GitHub:
Open the software’s GitHub repository, e.g., https://github.com/lxie-0709/TeloComp.
Click the “Code” button and select “Download ZIP” to get the package, or copy the repository URL for git clone.
To clone via command line:

.. code:: bash

    $ git clone git@github.com:lxie-0709/TeloComp.git
    $ cd TeloComp

2.Install dependencies and configure the software.
Please install the required dependencies under the ``Dependencies/`` directory and configure the executable programs in the ``bin/``directory, respectively.
The Dependencies folder is intended for third-party dependency packages, whereas the bin directory contains or links to the actual tools to be executed.

(1)Installing dependencies

.. code:: bash

    $ sh install.sh 

(2)Configuring TeloComp

.. code:: bash

    $ sh setup.sh

3.Install GenomeSyn

Download GenomeSyn and place the uncompressed GenomeSyn-1.2.7 directory under your root path (/yourPATH/).
Set the execution permission and add the binaries to your system PATH:

.. code:: bash
  
    $ chmod -R 777 GenomeSyn-1.2.7
    $ echo "export PATH=\$PATH:/yourPATH/GenomeSyn-1.2.7/bin" >> ~/.bashrc
    $ source ~/.bashrc




