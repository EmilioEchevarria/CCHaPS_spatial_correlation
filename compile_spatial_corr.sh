#!/bin/bash

module load intel-fc/2020.4.304
module load netcdf/4.8.1-intel20
module load openmpi/4.1.0-ofed51-intel20

source /apps/intel/fc/2020.4.304/bin/compilervars.sh intel64

export NC_INC="-I/apps/netcdf/4.8.1-intel20/include"
export NC_LIB="-L/apps/netcdf/4.8.1-intel20/lib64/ -lnetcdff -lnetcdf"

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/apps/netcdf/4.8.1-intel20/lib64/:/apps/openmpi/4.1.0-ofed51-intel20/lib64/

mpifort -c $NC_INC -O3 -ip -g -traceback -assume byterecl spatial_correlation.f90
mpifort -lcurl -o spatial_correlation spatial_correlation.o $NC_LIB

rm spatial_correlation.o