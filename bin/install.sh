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
tar zxvf Flye-2.9.3.tar.gz
cd Flye-2.9.3/
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

###GenomeSyn
##1.MUMmer
current_path=`pwd`
tar -zxvf mummer-4.0.0beta2.tar.gz
cd mummer-4.0.0beta2
./configure --prefix=`pwd`
make
# Add MUMmer tools to your PATH
mummer_path="export PATH=$current_path/mummer-4.0.0beta2:\$PATH"
echo $mummer_path >> ~/.bashrc
##2.SVG of Perl module
cd $current_path
tar -zxvf SVG-2.85.tar.gz
cd SVG-2.85/
perl Makefile.PL
make
make test
##3.BioPerl of Perl module
cd $current_path
tar -zxvf BioPerl-1.7.8.tar.gz
cd BioPerl-1.7.8/
perl Makefile.PL
make
make test
# Add Perl module to your PATH
perl_path="export PERL5LIB=$current_path/SVG-2.85/lib:$current_path/BioPerl-1.7.8/lib:\$PERL5LIB"
echo $perl_path >> ~/.bashrc
##4.Python and python module
cd $current_path
tar -zxvf Python-3.9.4.tgz -C ./
cd ./Python-3.9.4
./configure --prefix=$current_path/Python-3.9.4/localpython
make
make install
#Add Python to your PATH
python_path="export PATH=$current_path/Python-3.9.4/python:\$PATH"
echo $python_path >> ~/.bashrc
##4.1 svglib of python module
cd $current_path
tar -zxvf svglib-1.1.0.tar.gz -C ./
cd svglib-1.1.0/
$current_path/Python-3.9.4/python setup.py install
#Add GenomeSyn to your PATH
genomesyn_path="export PATH=$current_path/bin:\$PATH"
echo $genomesyn_path >> ~/.bashrc
#Refers to the tool execution name
chmod +x telocomp_Filter.py telocomp_Assembly.py telocomp_maxmin.py telocomp_Complement.py telocomp_Collinearity.py
echo "alias telocomp_Filter='python \"$(pwd)/bin/telocomp_Filter.py\"'" >> ~/.bashrc
echo "alias telocomp_Assembly='python \"$(pwd)/bin/telocomp_Assembly.py\"'" >> ~/.bashrc
echo "alias telocomp_maxmin='python \"$(pwd)/bin/telocomp_maxmin.py\"'" >> ~/.bashrc
echo "alias telocomp_Complement='python \"$(pwd)/bin/telocomp_Complement.py\"'" >> ~/.bashrc
echo "alias telocomp_Collinearity='python \"$(pwd)/bin/telocomp_Collinearity.py\"'" >> ~/.bashrc
source ~/.bashrc








