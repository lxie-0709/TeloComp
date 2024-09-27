# TeloComp
TeloComp is an efficient integrated software package for telomere extraction and complementation. It outputs new genome and telomere complementation information, and visualizes the complemented telomere parts with line graphs and colinearity graphs. It is more friendly to scientific researchers and is committed to more complete T2T genome assembly.

![image](https://github.com/lxie-0709/TeloComp/blob/main/Synteny.png)

# Install
TeloComp is an executable file written in Python that users can run directly, but before using TeloComp, you need to install all the dependencies of the software. You can download the current source distribution here or visit Github, or you can install it directly using [conda](https://anaconda.org/telocomp).

    git clone git@github.com:lxie-0709/TeloComp.git

# Dependencies

1.Flye

You can download it [here](https://github.com/fenderglass/Flye/archive/refs/tags/2.9.3.tar.gz). We use Flye version 2.9.4, you can use this version or higher.

    $ tar zxvf Flye-2.9.3.tar.gz
    $ cd Flye-2.9.3/
    $ python setup.py install
    #测试
    $ flye -h

2.NextPolish and Racon

(1) NextPolish 

You can download it [here](https://github.com/Nextomics/NextPolish/releases/download/v1.4.1/NextPolish.tgz). We use NextPolish version 1.4.1, you can use this version or higher.

    $ tar -zxvf NextPolish.tgztar
    $ cd NextPolish &&make-j 10
    $export PATH=$PATH:/pub4/huangshoubian1/NextPolis

(2) Racon

You can download it [here](https://github.com/isovic/racon/archive/refs/tags/1.4.3.tar.gz). We use Racon version 1.4.3, you can use this version or higher.

    $ tar -zxvf racon-v1.4.3.tar.gz
    $ cd racon-v1.4.3/
    $ conda install cmake 
    $ cmake -DCMAKE_BUILD_TYPE=Release..
    $ make&&make install
    $ export PATH=$PATH:/pub4/huangshoushoubian1/racon/bin #写入环境变量
    #测试
    $ racon -h 
    $ racon --version 
