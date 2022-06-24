# MTMOD_tools
This is the repository hosing information to install necessary tools for the MTMOD project.

## 1. Preparing suggested operating environments.
Anaconda (on Linux) is a good cross-platoform computing environment for the majority of tools used in MTMOD. Docker is needed for earthquake dynamic rupture software SeisSol. 

Disk usage: Anaconda may need a minimum of 6GB (13 GB on Dunyu's Ubuntu Windows subsystem). the Anaconda environment 'mtmod' may need 2GB. The MTMOD_tools folder may need 200 MB. 

## 2. To install Anaconda on Linux (or Ubuntu subsystem on Windows),
```
wget https://repo.anaconda.com/archive/Anaconda3-2022.05-Linux-x86_64.sh
bash Anaconda3-2022.05-Linux-x86_64.sh
```
Right after the installation of Anaconda, you will be asked if to run 'conda init'. We suggested 'yes'.
Or, you need to use 
```
source path_to_anaconda/bin/activate
conda init
```
to initiate an Anaconda environment. 

## 3. To install Docker on Windows, click this [link](https://desktop.docker.com/win/main/amd64/Docker%20Desktop%20Installer.exe?utm_source=docker&utm_medium=webreferral&utm_campaign=dd-smartbutton&utm_location=module) for the installer. 

## 4. Source the sourceme_install.sh to install the mtmod conda env. 
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
