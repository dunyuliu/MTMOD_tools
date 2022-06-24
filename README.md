# MTMOD_tools
This is the repository hosing information to install necessary tools for the MTMOD project.

1. Preparing suggested operating environments.
Anaconda (on Linux) is a good cross-platoform computing environment for the majority of tools used in MTMOD. Docker is needed for earthquake dynamic rupture software SeisSol. 

2. To install Anaconda on Linux (or Ubuntu subsystem on Windows),
```
wget https://repo.anaconda.com/archive/Anaconda3-2022.05-Linux-x86_64.sh
bash Anaconda3-2022.05-Linux-x86_64.sh
```
3. To install Docker on Windows, click this [link](https://desktop.docker.com/win/main/amd64/Docker%20Desktop%20Installer.exe?utm_source=docker&utm_medium=webreferral&utm_campaign=dd-smartbutton&utm_location=module) for the installer. 

4. Run the shell script install.sh
```
source sourceme_install.sh
```
to install the necessary python packages, okada_wrapper, and Elastic_stresses_py.

Misc.
useful conda commands:
```
conda update -n base -c defaults conda
conda env remove -n mtmod
```
