#!bin/bash
#Refers to the tool execution name
chmod +x telocomp_Filter.py telocomp_Assembly.py telocomp_maxmin.py telocomp_Complement.py telocomp_Collinearity.py
echo "alias telocomp_Filter='python \"$(pwd)/bin/telocomp_Filter.py\"'" >> ~/.bashrc
echo "alias telocomp_Assembly='python \"$(pwd)/bin/telocomp_Assembly.py\"'" >> ~/.bashrc
echo "alias telocomp_maxmin='python \"$(pwd)/bin/telocomp_maxmin.py\"'" >> ~/.bashrc
echo "alias telocomp_Complement='python \"$(pwd)/bin/telocomp_Complement.py\"'" >> ~/.bashrc
echo "alias telocomp_Collinearity='python \"$(pwd)/bin/telocomp_Collinearity.py\"'" >> ~/.bashrc
source ~/.bashrc








