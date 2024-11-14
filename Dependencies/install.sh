#!bin/bash
##1.samtools
current_path=$(pwd)
tar jxf samtools-1.18.tar.bz2
cd samtools-1.18/
./configure --prefix=$current_path/samtools-1.18
make
make install
#Add samtools to your PATH
echo "export PATH=\$PATH:$current_path/samtools-1.18" >> ~/.bashrc
##2.minimap2
cd $current_path
tar -zxvf minimap2-2.27.tar.gz
cd minimap2-2.27/
make
#Add minimap2 to your PATH
echo "export PATH=\$PATH:$current_path/minimap2-2.27" >> ~/.bashrc
##3.bwa
cd $current_path
tar jxf bwa-0.7.17.tar.bz2
cd bwa-0.7.17/
make
#Add bwa to your PATH
echo "export PATH=\$PATH:$current_path/bwa-0.7.17" >> ~/.bashrc
##4.Flye
cd $current_path
tar zxvf Flye-2.9.4.tar.gz
cd Flye-2.9.4/
python setup.py install
##5.racon
cd $current_path
tar -zxvf racon-v1.4.3.tar.gz
cd racon-v1.4.3/
conda install cmake
cmake -DCMAKE_BUILD_TYPE=Release..
make&&make install
#Add racon to your PATH
echo "export PATH=\$PATH:$current_path/racon-v1.4.3/bin" >> ~/.bashrc
##6.NextPolish
cd $current_path
tar -zxvf NextPolish.tgz
cd NextPolish &&make-j 10
pip install paralleltask
#Add NextPolish to your PATH
echo "export PATH=\$PATH:$current_path/NextPolish" >> ~/.bashrc
##7.SVG of Perl module
cd $current_path
tar -zxvf SVG-2.85.tar.gz
cd SVG-2.85/
perl Makefile.PL
make
make test
