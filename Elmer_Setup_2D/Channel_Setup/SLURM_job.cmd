#!/bin/bash
#SBATCH -o  /home/hpc/a2901/di49sux/elmer_simulations/IceShelfChannels/SIMULATIONS_IN_PAPER/AsymetricSMB_and_BMB/Channel_GrdFactorSurf0.0_GrdFactorBase-2.631/SLURM_job.%j.%N.out
#SBATCH -D /home/hpc/a2901/di49sux/elmer_simulations/IceShelfChannels/SIMULATIONS_IN_PAPER/AsymetricSMB_and_BMB/Channel_GrdFactorSurf0.0_GrdFactorBase-2.631/  
#SBATCH -J RdrewsTest
#SBATCH --get-user-env
#SBATCH --ntasks=8
#SBATCH --mail-type=none
#SBATCH --mail-user=rdrews@benicetoice.eu
#SBATCH --export=NONE
#SBATCH --time=02:00:00
#SBATCH --clusters=mpp2
source /etc/profile.d/modules.sh

module unload mpi.ibm
module unload mpi.ompi
module load mpi.intel
module load scalapack
module load mumps
module load metis
#
#
export ELMER_HOME="/home/hpc/a2901/di36hov/INSTALL/InstallElmerv83WithMMGRecomp/elmerice/Elmer_devel_12-04-18/"
export ELMER_SOLVER_HOME="/home/hpc/a2901/di36hov/INSTALL/InstallElmerv83WithMMGRecomp/elmerice/Elmer_devel_12-04-18/bin"
export PATH=/home/hpc/a2901/di36hov/INSTALL/InstallElmerv83WithMMGRecomp/elmerice/Elmer_devel_12-04-18/bin:$PATH
export LD_LIBRARY_PATH=/home/hpc/a2901/di36hov/INSTALL/InstallElmerv83WithMMGRecomp/elmerice/Elmer_devel_12-04-18/share/elmersolver/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/hpc/a2901/di36hov/INSTALL/InstallElmerv83WithMMGRecomp/elmerice/Elmer_devel_12-04-18/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/hpc/a2901/di36hov/INSTALL/InstallElmerv83WithMMGRecomp/elmerice/Elmer_devel_12-04-18/lib/elmersolver/:$LD_LIBRARY_PATH

#export PATH=/lrz/sys/io_tools/netcdf/4.3.3/serial_icc17/lib/lib/libnetcdf.so:$PATH
#export PATH=/lrz/sys/io_tools/netcdf/4.4.1/intel/impi_2017/lib/libnetcdf.so:$PATH
export LD_LIBRARY_PATH=/home/hpc/a2901/di36hov/INSTALL/mmg-5.3.10/MMG_devel/lib:$LD_LIBRARY_PATH
export PATH=/home/hpc/a2901/di36hov/INSTALL/mmg-5.3.10/MMG_devel/bin:$PATH
export LD_LIBRARY_PATH=/home/hpc/a2901/di36hov/INSTALL/mmg-5.3.10/MMG_devel/lib:$LD_LIBRARY_PATH
export PATH=/home/hpc/a2901/di36hov/INSTALLMPP3/gmsh-2.12.0-source/build:$PATH
export LD_LIBRARY_PATH=/home/hpc/a2901/di36hov/INSTALL/mumps_LRZ/MUMPS_5.0.2/lib:$LD_LIBRARY_PATH


make clean grid submit

### Check job status
### scontrol --clusters=mpp2 show jobid=109938
