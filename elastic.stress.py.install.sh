#! /bin/bash
rm -rf Elastic_stresses_py
git clone https://github.com/kmaterna/Elastic_stresses_py.git
echo '------------------------'
echo 'INSTALLING Dr. Kathryn Materna Elastic_stresses_py.'
cd Elastic_stresses_py
python setup.py install && echo 'Elastic_stresses_py installed SUCCESSFULLY' || echo FAIL
elastic_stresses_driver.py examples/example_config.txt && echo 'Example run SUCCESS' ||echo FAIL
cd ..
