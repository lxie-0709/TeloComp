#!bin/bash
##1.samtools
current_path=`pwd`
tar jxvf samtools-1.18.tar.gz
cd samtools-1.18/
./configure --prefix=/home/vip47/biosoft/samtools-1.18
make
make install
#Add bwa to your PATH
bwa_path="export PATH=$current_path/samtools-1.18:\$PATH"
echo $bwa_path >> ~/.bashrc
##2.minimap2
tar -zxvf minimap2-2.27.tar.gz
cd minimap2-2.27/
make
##3.bwa
cd $current_path
tar jxf bwa-0.7.12.tar.bz2
cd bwa-0.7.12/
make
#Add bwa to your PATH
bwa_path="export PATH=$current_path/bwa-0.7.12:\$PATH"
echo $bwa_path >> ~/.bashrc
##4.Flye
cd $current_path
tar zxvf Flye-2.9.4.tar.gz
cd Flye-2.9.4/
setup.py install
##5.racon
cd $current_path
tar -zxvf racon-v1.4.3.tar.gz
cd racon-v1.4.3/
make&&make install
#Add bwa to your PATH
racon_path="export PATH=$PATH:$current_path/racon/bin:\$PATH"
echo $bwa_path >> ~/.bashrc
##6.NextPolish
cd $current_path
tar -zxvf NextPolish.tgztar
cd NextPolish &&make-j 10
#Add bwa to your PATH
nextPolish_path="export PATH=$PATH:$current_path/NextPolis:\$PATH"
echo $bwa_path >> ~/.bashrc

