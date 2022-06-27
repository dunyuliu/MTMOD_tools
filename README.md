# MTMOD_tools
This is the repository hosing information to install necessary tools for the MTMOD project.

## 1. Summary of operating environments, disk usage, and others.
[Anaconda](https://www.anaconda.com/) is a good cross-platform computing environment to perform Python/R data science and machine learning applications with thousands of open-source packages and libraries. And Anaconda on Linux works well for the majority of tools used in the [MTMOD: Megathrust Modeling Framework](https://sites.utexas.edu/mtmod/) project's 2022 summer school. Docker is needed for earthquake dynamic rupture software [SeisSol](https://www.seissol.org/). (Not sure about Anaconda on Mac, update this later.)

### Disk usage: 
A. We need to install Anaconda, which may need a minimum of 6GB (13 GB on Dunyu's Ubuntu Windows subsystem). 
B. We need to clone this repositroy (named MTMOD_tools by default) to your local machine and it will finally use about 200 MB.
C. We will install the conda environment 'mtmod' with the instruction in the section 4. It may need 2GB.

### Installation time
The installation will take 20 minutes. So, it is suggested to get everything installed before attending the summer school.

## 2. To install Anaconda on Linux (or Ubuntu subsystem on Windows),
```
wget https://repo.anaconda.com/archive/Anaconda3-2022.05-Linux-x86_64.sh
bash Anaconda3-2022.05-Linux-x86_64.sh
```
Right after the installation of Anaconda, you will be asked if to run 'conda init'. We suggest to type in 'yes'.
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
