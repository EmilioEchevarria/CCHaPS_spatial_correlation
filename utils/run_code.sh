#!/bin/bash 
#SBATCH --ntasks=12
#SBATCH --job-name=corr_
#SBATCH --time=00:25:00
#SBATCH --mem=400GB
#SBATCH -A OD-215204

module load intel-fc/2020.4.304
module load netcdf/4.8.1-intel20
module load openmpi/4.1.0-ofed51-intel20

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/apps/netcdf/4.8.1-intel20/lib64/:/apps/openmpi/4.1.0-ofed51-intel20/lib64/
source /apps/intel/fc/2020.4.304/bin/compilervars.sh intel64

rm -f output_*.log
echo $SLURM_NTASKS

/apps/openmpi/4.1.0-ofed51-intel20/bin/mpirun -n $SLURM_NTASKS ./spatial_correlation
