# This is a recreation of MTMOD_env.yml to shell script. 
# 20230731
conda install -c conda-forge python>=3 matplotlib numpy scipy xarray
conda install -c conda-forge gmt pygmt 
conda install -c conda-forge gfortran gsl jupyter pip
pip install okada-wrapper tectonic-utils # used by elastic.stress.py
