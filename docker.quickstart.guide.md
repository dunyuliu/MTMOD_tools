## Docker Quick Start Guide

The mtmod docker image that contains Anaconda ```mtmod``` env is available via ```docker pull dunyuliu/mtmod```. <br/>
If [Docker Desktop](https://www.docker.com/products/docker-desktop/) is installed locally, users can pull the image and run it as a container. <br/>

For Windows OS with WSL available, Docker Desktop will be installed in the C:/ drive. <br/>
The external drives could be mounted to the docker container, which allows editing files and viewing results. <br/>

Currently, the image is ~10 GBs. <br>

Below are detailed commands to pull and run the mtmod docker image and how to mount external drives. <br/>

# Windows OS

Open Powershell with admin access, then type the following command in the terminal.
```
docker run -it --name $mtmod -v $mtmod.external.drive:$proxy.path dunyuliu/mtmod
```
where $mtmod is the name of the container shown in the Docker Desktop, for example mtmod, and <br/>
$mtmod.external.drive is the path of the external drive to be mounted to the container, say D:/mtmod.data, and <br/>
$proxy.path is the proxy path of the external drive in the container, say /mount  . <br/>

Here is an example:
```
docker run -it --name mtmod -v D:/mtmod.data:/mount dunyuliu/mtmod
```

NOTE: If you don't want to mount external drives to the container, please don't add the flag -v and the following path. In this case, you will run everything inside the container, which is the C:/ drive and it may eat up disk space <br/>

AFter running the command, you will see files downloaded, if the image is not available locally, and a container created. <br/>

Then, you can open a terminal in the container from the Docker Desktop and run it as a normal Linux Anaconda system.

The tools are available under the path /home/MTMOD_tools/.
