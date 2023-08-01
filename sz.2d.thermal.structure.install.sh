#! /bin/bash
echo 'Cloning SZ_2D_thermal_structure ... ...'
echo 'It requires fenics-2019, petsc, gmsh, and ParametricModelUtils ... ...'
echo ' '
echo 'It also requires gmsh, which is installed via apt-get ... ...'
apt-get install gmsh
echo 'Check gmsh version ... ...'
gmsh --info 
echo ' '
  
rm -rf SZ_2D_thermal_structure
git clone https://github.com/gabriellemhobson/SZ_2D_thermal_structure

echo 'Cloning ParametricModelUtils ... ...'
cd SZ_2D_thermal_structure
git clone https://github.com/hpc4geo/ParametricModelUtils.git
cd ..

echo 'Now back to SZ_2D_thermal_structure root directory ... ...'
echo 'Running a test example ... ...'

echo 'Check the repo is up to date and paths are correct'

source setup.sh

echo 'Enter the generate_mesh subdirectory'

cd generate_meshes

echo 'Check that the slab profile info in the csv input file is correct.'
echo 'Generate the geometry and mesh files:'

python3 driver_generate_mesh_generic.py --profile_csv "cascadia_start_end_points.csv" --slab_name "cascadia" --corner_depth -35.0 --output_path "output"

echo 'Create the required mesh files for fenics usage and post-processing steps.'
cd ..

python3 convert_msh_to_fenics_files.py --mesh_dir 'generate_meshes/output/cascadia_profile_B' --profile_name 'cascadia_profile_B' 

echo 'Set the desired ranges of input parameters in input_param.csv.'
echo 'Set the forward model running.'

python3 schedule_script.py --profile_name "cascadia_profile_B" --mesh_dir "generate_meshes/output/cascadia_profile_B" --output_path "output/cascadia_profile_B/example" --sample_method "halton" --n1 1 --n2 1 --seed 92014 --jobs_csv "cascadia_profile_B_example_log.csv" --viscosity_type "isoviscous" 

echo 'Once the forward model is done, perform the post-processing steps to create plots and compute isotherm-slab interface intersection locations.'

python3 post_process.py --jobs_csv "cascadia_profile_B_example_log.csv" --mesh_path "generate_meshes/output/cascadia_profile_B" --profile_name "cascadia_profile_B" --include "halton"

echo 'Look at your plots and be proud that you've run this code!'



