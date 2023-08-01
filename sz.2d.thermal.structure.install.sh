#! /bin/bash
echo 'Cloning SZ_2D_thermal_structure ... ...'
echo 'It requires fenics-2019, petsc, gmsh, and ParametricModelUtils ... ...'
echo ' '
rm -rf SZ_2D_thermal_structure
git clone https://github.com/gabriellemhobson/SZ_2D_thermal_structure

echo 'Cloning ParametricModelUtils ... ...'
cd SZ_2D_thermal_structure
git clone https://github.com/hpc4geo/ParametricModelUtils.git
cd ..

echo 'Now back to SZ_2D_thermal_structure root directory ... ...'
echo 'Running a test example ... ...'





