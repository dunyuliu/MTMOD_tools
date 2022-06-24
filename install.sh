## Install Anaconda on Linux
#wget https://repo.anaconda.com/archive/Anaconda3-2022.05-Linux-x86_64.sh
#bash Anaconda3-2022.05-Linux-x86_64.sh

# Creating a directory mtmod to host all the tools.
mkdir mtmod
cd mtmod

# Creating a conda environment 'mtmod' using the MTMOD_env.yml.
conda env create -f MTMOD_env.yml
# Activate mtmod.
conda activate mtmod

# Downloading Dr. Ben Thompson's okada_wrapper. (It is already installed in the MTMOD_env.yml.)
git clone https://github.com/tbenthompson/okada_wrapper.git

https://github.com/tbenthompson/okada_wrapper.git
# Downloading Dr. Kathryn Materna's Elastic_stresses_py for Coulomb Stress.
git clone https://github.com/kmaterna/Elastic_stresses_py.git
cd Elastic_stresses_py
python setup.py install

