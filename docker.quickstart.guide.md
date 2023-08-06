## Docker Quick Start Guide

First, you need to install [Docker Desktop](https://www.docker.com/products/docker-desktop/). For Windows users, please download Docker Desktop from previous link and install it. For Mac users, please choose between Apple Chip or Intel Chip for the installer. <br/>

After Docker Desktop is installed, you need to run it with admin access. Then, you need to open a terminal for Mac users or a Powershell terminal with admin access for Windows users, and type the command 
```
docker pull dunyuliu/mtmod:2023
```
where the mtmod docker image that contains Anaconda ```mtmod``` env is available and will be pulled and run as a container. <br/>

For Windows OS with WSL available, Docker Desktop will be installed in the C:/ drive or system drive. <br/>
The external drives could be mounted to the docker container, which allows editing files and viewing results. <br/>

Currently, the image is ~10 GBs. <br>

Below are detailed commands to pull and run the mtmod docker image and how to mount external drives for Windows users as an example. <br/>

Open Powershell with admin access, then type the following command in the terminal.
```
docker run -it --name $mtmod -v $mtmod.external.drive:$proxy.path dunyuliu/mtmod:2023
```
where $mtmod is the name of the container shown in the Docker Desktop, for example mtmod, and <br/>
$mtmod.external.drive is the path of the external drive to be mounted to the container, say D:/mtmod.data, and <br/>
$proxy.path is the proxy path of the external drive in the container, say /mount  . <br/>

Here is an example:
```
docker run -it --name mtmod -v D:/mtmod.data:/mount dunyuliu/mtmod:2023
```

NOTE: If you don't want to mount external drives to the container, please don't add the flag -v and the following path. In this case, you will run everything inside the container, which is the C:/ drive and it may eat up disk space <br/>

AFter running the command, you will see files downloaded, if the image is not available locally, and a container created. <br/>

Then, you can open a terminal in the container from the Docker Desktop and run it as a normal Linux Anaconda system.

The tools are available under the path /home/MTMOD_tools/.
