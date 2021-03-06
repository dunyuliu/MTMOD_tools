## First, you need to install Anaconda on Linux
# wget https://repo.anaconda.com/archive/Anaconda3-2022.05-Linux-x86_64.sh
# bash Anaconda3-2022.05-Linux-x86_64.sh

# Creating a conda environment 'mtmod' using the MTMOD_env.yml.
conda env create -f MTMOD_env.yml
# Activate mtmod.
conda activate mtmod
# 1 #
# Downloading Dr. Ben Thompson's okada_wrapper. (It is already installed in the MTMOD_env.yml.)
rm -rf okada_wrapper
echo 'GIT CLONENING Dr. Ben Thompson okada_wrapper.'
git clone https://github.com/tbenthompson/okada_wrapper.git
echo '------------------------'
# Run the test.
python okada_wrapper/test_okada.py && echo 'okada_wrapper installed and tested SUCESSFULLY' || echo FAIL
echo '------------------------'

# 2 #
echo 'GIT CLONENING Dr. Ben Thompson cutde.'
rm -rf cutde
git clone https://github.com/tbenthompson/cutde.git
# Run the test
cd test
cd cutde
python test_cutde.py && echo 'cutde installed and tested SUCCESSFULLY' || echo FAIL
cd ../..

# 3 #
# Downloading Dr. Kathryn Materna's Elastic_stresses_py for Coulomb Stress.
rm -rf Elastic_stresses_py
git clone https://github.com/kmaterna/Elastic_stresses_py.git
echo '------------------------'
echo 'INSTALLING Dr. Kathryn Materna Elastic_stresses_py.'
cd Elastic_stresses_py
python setup.py install && echo 'Elastic_stresses_py installed SUCCESSFULLY' || echo FAIL
elastic_stresses_driver.py examples/example_config.txt && echo 'Example run SUCCESS' ||echo FAIL
cd ..

# 4 #
# Downloading Dr. Camilla Cattania's crs to calculate rate-state seismicity evolution induced by time dependent, heterogeneous Coulomb stress changes.
git clone https://github.com/dunyuliu/crs_mtmod.git
cd crs_mtmod
echo '------------------------'
echo 'INSTALLING Dr. Camilla Cattania crs.'
source crs_install.sh && echo 'crs installed and tested SUCCESSFULLY' || echo FAIL
cd ..

# Exiting ... #
echo '------------------------'
echo 'MTMOD tools READY to use. Cheers!'
