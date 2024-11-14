#!bin/bash
#Refers to the tool execution name
chmod +x telocomp_Filter.py telocomp_Assembly.py telocomp_maxmin.py telocomp_Complement.py telocomp_Collinearity.py
echo "alias telocomp_Filter='python \"$(pwd)/telocomp_Filter.py\"'" >> ~/.bashrc
echo "alias telocomp_Assembly='python \"$(pwd)/telocomp_Assembly.py\"'" >> ~/.bashrc
echo "alias telocomp_maxmin='python \"$(pwd)/telocomp_maxmin.py\"'" >> ~/.bashrc
echo "alias telocomp_Complement='python \"$(pwd)/telocomp_Complement.py\"'" >> ~/.bashrc
echo "alias telocomp_Collinearity='python \"$(pwd)/telocomp_Collinearity.py\"'" >> ~/.bashrc

#Installing python modules

packages=("numpy" "pandas" "pysam" "biopython" "matplotlib")
echo "Installing Python packages..."

for package in "${packages[@]}"; do
    echo "Attempting to install $package with pip..."
    if ! pip3 install "$package"; then
        echo "pip installation for $package failed. Attempting to install with conda..."
        conda install -y -c bioconda "$package"
    else
        echo "$package installed successfully with pip."
    fi
done

echo "Python packages installation complete."





