#!bin/bash
#Refers to the tool execution name
ln -s $(pwd)/telocomp_Filter.py ./telocomp_Filter
ln -s $(pwd)/telocomp_Assembly.py ./telocomp_Assembly
ln -s $(pwd)/telocomp_maxmin.py ./telocomp_maxmin
ln -s $(pwd)/telocomp_Complement.py ./telocomp_Complement
ln -s $(pwd)/telocomp_Collinearity.py ./telocomp_Collinearity
echo "export PATH=\$PATH:$(pwd)" >> ~/.bashrc

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





