rm -rf FaultDynamics.tar FaultDynamics
mkdir FaultDynamics

wget --no-check-certificate "https://utexas-my.sharepoint.com/:u:/g/personal/dliu_ig_utexas_edu/ETNgb4u2othCqEzG9K8qsdcB7xMFbGRV5ByQofcrHXE9_A?e=v2Yu1c&download=1" -O FaultDynamics.tar
#wget --no-check-certificate "https://utexas-my.sharepoint.com/:u:/g/personal/dliu_ig_utexas_edu/EYvuFuoVZm5HmgAqMPNVk5kBDOjqYDR1HjKivchlHFKMqA?e=Mcmfh0&download=1" -O FaultDynamics.tar
# need to add $download=1 in the OneDrive link

tar -xvf FaultDynamics.tar -C FaultDynamics
cd FaultDynamics
python fault2D_SEM.py
cd ..
