#! /bin/bash

# 4 #
# Downloading Dr. Camilla Cattania's crs to calculate rate-state seismicity evolution induced by time dependent, heterogeneous Coulomb stress changes.
git clone https://github.com/dunyuliu/crs_mtmod.git
cd crs_mtmod
echo '------------------------'
echo 'INSTALLING Dr. Camilla Cattania crs.'
source crs_install.sh && echo 'crs installed and tested SUCCESSFULLY' || echo FAIL
cd ..
