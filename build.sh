#!/bin/bash
# EXECUTE AS:
# ./build.sh N= T= init='' Beta= J=
# Both versions (IsingModel_v1 - fixed boundaries creates image of final lattice; IsingModel_v2 - periodic boundaries, creates animation of lattice evolution) 
# are compiled and executed for the same arguments 
# Purpose of using both versions here is that despite the differences in the way the lattice is defined, 
# the write_netcdf module works well in both cases. 
# write_netcdf can also be used to animate the lattice evolution as demonstrated in IsingModel_v2.

#List here all the f90 files to compile, separated by spaces
#THERE MUST NOT be spaces around '=' in bash


myprogramfiles1="command_line.f90 write_netcdf.f90 IsingModel_v1.f90"
myprogramfiles2="command_line.f90 write_netcdf.f90 IsingModel_v2.f90"

#Name of compiled file
outfile1="ising_netcdf1"
outfile2="ising_netcdf2"

#Name of Python visualisation script
#There is only one visualisation script for both the versions
python_file='Visualise.py'

#Name of compiler
fc=gfortran
#Use nf-config to grab the compile and link flags. Backticks run command and grab output
fflags=`nf-config --fflags`
flibs=`nf-config --flibs`

#Actual compile line. Other flags etc can be added
$fc -g $fflags $myprogramfiles1 $flibs -o $outfile1

#Checks for compilation failure and exits
if [ $? -ne 0 ]; then
    echo "Error during compilation"
    exit 1
fi
echo 'Version 1 compiled.'

$fc -g $fflags $myprogramfiles2 $flibs -o $outfile2
if [ $? -ne 0 ]; then
    echo "Error during compilation"
    exit 1
fi
echo 'Version 2 compiled.'

#Run ising_netcdf with arguments 
./$outfile1 "$@"
./$outfile2 "$@"

if [ $? -ne 0 ]; then
    echo "Error executing FORTRAN"
    exit 1
fi

#Run python file
python3 $python_file

if [ $? -ne 0 ]; then
    echo "Error visualising"
    exit 1
fi
