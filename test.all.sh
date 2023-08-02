#! /bin/bash
echo 'Usage: source test.all.sh'
echo 'Script to test all the software packages ... ...'

echo 'Testing elastic.stress.py ... ...'
echo 'source elastic.stress.py.install.sh'
source elastic.stress.py.install.sh
echo ' '
echo ' '

echo 'Testing crs.mtmod ... ...'
echo 'source crs.mtmod.install.sh'
source crs.mtmod.install.sh
echo ' '
echo ' '

echo 'Testing sz.2d.thermal.structure ... ...'
echo 'source sz.2d.thermal.structure.install.sh'
source sz.2d.thermal.structure.install.sh
echo ' '
echo ' '
