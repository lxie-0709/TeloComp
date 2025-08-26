Installation
============

This tool requires a number of third-party software packages.  
Before installation, please make sure you are working within a Linux environment and that
`git`, `conda` or `mamba` are already installed.

.. contents::
    :local:

1. Dependency installation
--------------------------

* python = 3.11
* pysam = 0.22.0
* biopython = 1.81
* samtools = 1.18
* minimap2 = 2.27
* bwa = 0.7.17
* Flye = 2.9.4
* pilon = 1.2.4
* NextPolish = 1.4.1
* GenomeSyn = 1.2.7

1.1 NextPolish (recommended: install from source)

.. code-block:: bash

    $ git clone https://github.com/shilong-liu/NextPolish.git
    $ cd NextPolish
    $ make

Add NextPolish to your ``PATH`` (recommended):

.. code-block:: bash

    $ echo "export PATH=\$PATH:/path/to/NextPolish" >> ~/.bashrc
    $ source ~/.bashrc

1.2 Other dependencies: They can be installed via conda, or downloaded and installed from their GitHub repositories, or installed automatically by running the ``install.sh`` script provided with TeloComp.

2. Linux Installation
---------------------

To use the software, you need to follow the following steps to install it.

2.1 Obtain software package from GitHub:
Open the software’s GitHub repository, e.g., https://github.com/lxie-0709/TeloComp.
Click the “Code” button and select “Download ZIP” to get the package, or copy the repository URL for git clone.
To clone via command line:

.. code:: bash

    $ git clone git@github.com:lxie-0709/TeloComp.git
    $ cd TeloComp

2.2 Install dependencies and configure the software.
Please install the required dependencies under the ``Dependencies/`` directory and configure the executable programs in the ``bin/``directory, respectively.
The Dependencies folder is intended for third-party dependency packages, whereas the bin directory contains or links to the actual tools to be executed.

(1)Installing dependencies

.. code:: bash

    $ sh install.sh 

(2)Configuring TeloComp

.. code:: bash

    $ sh setup.sh

3.Install GenomeSyn
-------------------

Download `GenomeSyn <https://github.com/JM-SONG/GenomeSyn/archive/refs/heads/main.zip>`_ and place the uncompressed GenomeSyn-1.2.7 directory under your root path (/yourPATH/).
Set the execution permission and add the binaries to your system PATH:

.. code:: bash
  
    $ chmod +x GenomeSyn-1.2.7
    $ echo "export PATH=\$PATH:/yourPATH/GenomeSyn-1.2.7/bin" >> ~/.bashrc
    $ source ~/.bashrc

4. Verify Installation
----------------------

Check that all required tools are available:

.. code-block:: bash

    # Activation environment
    $ source ~/.bashrc
    
    $ telocomp_Filter_1 -h

    $ telocomp_Filter_2 -h
    
    $ telocomp_Assembly -h
    
    $ telocomp_maxmin -h
    
    $ telocomp_Complement -h

    $ telocomp_Collinearity -h

If the help message of each command is printed, the installation has been completed successfully.


.. tip:: Please open an issue on `GitHub <https://github.com/lxie-0709/TeloComp>`_ for feature requests or bug reports.
