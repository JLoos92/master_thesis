#!/bin/bash
#SBATCH -o  /home/hpc/a2901/di49sux/elmer_simulations/IceShelfChannels/AsymetricSMB_and_BMB/Channel_GrdFactorSurf20.0_GrdFactorBase0.0/SLURM_job.%j.%N.out
#SBATCH -D  /home/hpc/a2901/di49sux/elmer_simulations/IceShelfChannels/AsymetricSMB_and_BMB/Channel_GrdFactorSurf20.0_GrdFactorBase0.0/ 
#SBATCH -J RdrewsTest
#SBATCH --get-user-env
#SBATCH --ntasks=6
#SBATCH --mail-type=none
#SBATCH --mail-user=rdrews@benicetoice.eu
#SBATCH --export=NONE
#SBATCH --time=01:00:00
#SBATCH --clusters=mpp2
source /etc/profile.d/modules.sh

module unload mpi.ibm
module unload mpi.ompi
module load mpi.intel
module load scalapack
module load mumps
module load metis


export ELMER_HOME="/home/hpc/a2901/di49sux/soft_elmer/elmerice/Elmer_devel_installed_mumps/"
export ELMER_SOLVER_HOME="/home/hpc/a2901/di49sux/soft_elmer/elmerice/Elmer_devel_installed_mumps/share/elmersolver/"

make clean grid submit

### Check job status
### scontrol --clusters=mpp2 show jobid=109938
