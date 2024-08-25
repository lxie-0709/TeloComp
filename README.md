# TeloComp
TeloComp is an efficient integrated software package for telomere extraction and complementation. It outputs new genome and telomere complementation information, and visualizes the complemented telomere parts with line graphs and colinearity graphs. It is more friendly to scientific researchers and is committed to more complete T2T genome assembly.

![image](https://github.com/lxie-0709/TeloComp/blob/main/Synteny.png)

# Install
TeloComp is an executable file written in Python that users can run directly, but before using TeloComp, you need to install all the dependencies of the software. You can download the current source distribution here or visit Github, or you can install it directly using [conda](https://anaconda.org/telocomp).

    git clone git@github.com:lxie-0709/TeloComp.git

# Dependencies

1.Flye

You can download it [here](https://github.com/fenderglass/Flye/archive/refs/tags/2.9.3.tar.gz). We use Flye version 2.9.3, you can use       this version or higher.

    $ tar zxvf Flye-2.9.3.tar.gz
    $ cd Flye-2.9.3/
    $ python setup.py install
    #测试
    $ flye -h
