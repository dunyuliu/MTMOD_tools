This README describes how to install Julia inside the Anaconda mtmod env and run stress.inversion.julia. <br/>
First, activate the conda mtmod env by 
```
conda activate mtmod
```
Then, you may want to find a location to download and install Julia. <br/> 
For example, you may create a folder called Julia/ inside MTMOD_tools/ by
```
cd Julia
wget https://julialang-s3.julialang.org/bin/linux/x64/1.8/julia-1.8.0-linux-x86_64.tar.gz
tar zxvf julia-1.8.0-linux-x86_64.tar.gz
```
Then, you need to add Julia's path to the environment variable PATH by adding
```
export PATH="$PATH:<Julia directory>/bin"
```
to .bashrc if you are using bash shell. <br/>

Then, restart a new terminal, and activate mtmod again. <br/>
Enter Julia by typing
```
julia
```
in the terminal. <br/>
Then, type 
```
]
```
Then, install and build packages by 
```
add IJulia
build IJulia
add StatsBase LinearAlgebra DelimitedFiles DataFrames LazyGrids CSV
```
After you cd into stress.inversion.julia/, you should be able to run the julia script for stress inversion by 
```
julia serial_uchide_inversion.jl
```
