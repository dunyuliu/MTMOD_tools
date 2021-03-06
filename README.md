# MTMOD_tools
The repository contains the information on how to create the computing enviroment for the [MTMOD: Megathrust Modeling Framework](https://sites.utexas.edu/mtmod/) project's 2022 summer school.

## 1. Summary of operating environments, disk usage, and others.
[Anaconda](https://www.anaconda.com/) is a good cross-platform computing environment to perform Python/R data science and machine learning applications with thousands of open-source packages and libraries. And Anaconda on Linux works well for the majority of tools used in the MTMOD summer school. <br/>
[Docker](https://www.docker.com/) is needed for earthquake dynamic rupture software [SeisSol](https://www.seissol.org/). For details, please refer to section 3 in this doc. <br/> 
[MATLAB](https://www.mathworks.com/products/matlab.html) is needed for RateState and plotting crs results. <br/>

### Working systems
Ubuntu windows subsystem <br/>
(Currently testing on Mac, will update this later.)

### Tools currently included
  [okada_wrapper](https://github.com/tbenthompson/okada_wrapper.git) <br />
  [Elastic_stresses_py](https://github.com/kmaterna/Elastic_stresses_py.git) <br />
  [cutde](https://github.com/tbenthompson/cutde.git) <br />
  [crs](https://github.com/dunyuliu/crs_mtmod) <br />
  Earthquake-Python-Examples <br/>
  [SeisSol](https://www.seissol.org/) and [SeisSol Tutorial](https://github.com/SeisSol/Training) <br/>
  [RateState](https://drive.google.com/drive/folders/15nl880SFTFe61iJDIw38vunxTQdSBYZY?usp=sharing)
  
### Disk usage
A. The installation of Anaconda may need a minimum of 4.5GB. <br />
B. To clone this repositroy (named MTMOD_tools by default) to your local machine and it will finally use about 200 MB. <br />
C. To install the conda environment 'mtmod' with the instruction in the section 4. It may need 2GB. <br />

### Installation time
The installation of Anaconda may take 15 minutes. <br />
The creation of mtmod conda environment may take 20 minutes. <br />
So, it is suggested to get everything installed before attending the summer school.

## 2. To install Anaconda on Linux (or Ubuntu subsystem on Windows),
### Download and install
```
wget https://repo.anaconda.com/archive/Anaconda3-2022.05-Linux-x86_64.sh
bash Anaconda3-2022.05-Linux-x86_64.sh
```
### Initiate conda with different shells

Right after the installation of Anaconda, you will be asked to initiate conda with the command
```
conda init
```
Or, you can use 
```
conda init --all
```
to initiate conda for all the shells (bash, tcsh, fish, xonsh, zsh, powershell, etc.. NOTE, the command will modify .bashrc and/or .tcshrc). 

### Alternatives
You may want to initiate conda with the following lines
```
source $Anaconda_root_path"/etc/profile.d/conda.csh"
```
for tcsh or, 
```
source $Anaconda_root_path"etc/profile.d/conda.sh"
```
for bash shell.

## 3. To install and run [SeisSol tutorial](https://github.com/SeisSol/Training), please install Docker.
### To install Docker on Windows, click this [link](https://desktop.docker.com/win/main/amd64/Docker%20Desktop%20Installer.exe?utm_source=docker&utm_medium=webreferral&utm_campaign=dd-smartbutton&utm_location=module) for the installer. 
### On Windows, to run SeisSol tutorials,
Frist, Open a cmd terminal. <br/>
Then, run the following command
```
docker run -p 53155:53155 alicegabriel/seissol-training
```
Then, you will get a link after a few minutes downloading. <br/>
Finally, copy the link to a web browser and you can navigate to the jupyter notebook TPV13.ipynb for the benchmark problem TPV13 from the [SCEC dynamic rupture code verification project](https://strike.scec.org/cvws/tpv12_13docs.html).

## 4. Install the mtmod conda environment. 
First, clone this repository to your local machine with the command,
```
git clone https://github.com/dunyuliu/MTMOD_tools.git
cd MTMOD_tools
```
Then, use the following command line to create conda env mtmod,
```
conda env create -f MTMOD_env.yml
```
Then, use the following command line to download and install relevant tools
```
source sourceme_install.sh
```

## 5. Install RateState
RateState is currently available through this [google drive link](https://drive.google.com/drive/folders/15nl880SFTFe61iJDIw38vunxTQdSBYZY?usp=sharing). <br/>
It runs mainly in MATLAB but relis on hmmvp, which needs to be compiled with MATLAB mex C++ compiler. <br/> 

First, you need to install MinGW add-on for MATLAB, which can be found on MATLAB HOME/add-on. <br/>
Second, install [gnumex](https://sourceforge.net/projects/gnumex/files/latest/download). <br/>
```
cd gnumex
getenv MW_MINGW64_LOC
gnumex
```
Then, you will see a graphic interface.  <br/>
Copy and paste the path from getenv MW_MINGW64_LOC to MinGW root directory.  <br/>
Then, run the interface.  <br/>

Next, you need to install hmmvp. <br/>
```
cd RateState/hmmvp (or hmmvp_win for Windows users)
make
```
After hmmvp is successfully installed, you can run the test_run.m.  <br/>
For more details, please refer to Brief_Instruction.txt.

## Misc.
Useful conda commands,
```
conda update -n base -c defaults conda
conda env remove -n mtmod
du -hs $(conda info --base)/envs/mtmod
```

## On choosing a proper license for your open-source project
This [link](https://choosealicense.com/about/) provides very useful information on licenses.

