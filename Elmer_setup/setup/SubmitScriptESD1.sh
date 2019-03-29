#!/bin/bash -l
###----------------------------------------------------------------------------#
### Run script for ESD (1 or 2) refined
###----------------------------------------------------------------------------#

#PBS -N Bump100_5001000_0
#PBS -m abe
#PBS -M julius-loos@gmx.de
#PBS -o ${PBS_JOBNAME}${PBS_JOBID}.out
#PBS -e ${PBS_JOBNAME}${PBS_JOBID}.err
#PBS -l walltime=999:99:00
#PBS -l nodes=1:ppn=20:esd2
#PBS -q esd2
##PBS -l pmem=150gb
#=================================================================================================================
module load chains/INTEL-17.0
module load compiler/intel/17.0
module load mpi/impi/2017
module load numlib/mkl/2017.6.256

export ELMER_HOME="/home-link/epioi01/elmerice/Elmer_devel_04-17-18"
export ELMER_SOLVER_HOME="$ELMER_HOME/bin"

export PATH=/home-link/epioi01/elmerice/Elmer_devel_04-17-18/bin:$PATH
export LD_LIBRARY_PATH=/home-link/epioi01/elmerice/Elmer_devel_04-17-18/share/elmersolver/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home-link/epioi01/elmerice/Elmer_devel_04-17-18/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home-link/epioi01/elmerice/Elmer_devel_04-17-18/lib/elmersolver/:$LD_LIBRARY_PATH

#export PATH=/home-link/epioi01/INSTALL/mmg-5.3.10/MMG_devel/bin:$PATH
#export LD_LIBRARY_PATH=/home-link/epioi01/INSTALL/mmg-5.3.10/MMG_devel/lib:$LD_LIBRARY_PATH

export LD_LIBRARY_PATH=/home-link/epioi01/INSTALL/MUMPS_5.1.2/lib:$LD_LIBRARY_PATH
export PATH=/home-link/epioi01/INSTALL/gmsh-3.0.6-source/build:$PATH
#export LD_LIBRARY_PATH=/beegfs/home/tu/epioi01/INSTALL/MUMPS_5.1.2/lib:$LD_LIBRARY_PATH

#=================================================================================================================
echo Here comes the partition the job runs in:
echo $PBS_QUEUE
cd $PBS_O_WORKDIR

#ELMEREXEC="$ELMER_HOME/bin/ElmerSolver_mpi"
#cat $PBS_NODEFILE | sort | uniq > nodelist.tmp

#while read line; do
	#echo $line
	#scp $ELMEREXEC $line:/scratch/${PBS_JOBID}/ElmerSolver_mpi
	##scp $ELMEREXEC $line:/scratch/Elmer
#done < nodelist.tmp

#sleep 3s

ncpus=`wc -l < ${PBS_NODEFILE}`
echo Here comes the Nodelist:${PBS_NODEFILE}}
echo Number of Nodes: ${ncpus}}
#echo $PBS_NODEFILE


cp $ELMER_HOME/share/elmersolver/lib/FreeSurfaceSolver.so src/MyFreeSurfaceSolver.so
make compile
make ini
make grid
#mpirun -np $ncpus -machinefile $PBS_NODEFILE ElmerSolver_mpi
mpirun -n $ncpus ElmerSolver_mpi
#make submit
#mpirun -n 60 /scratch/Elmer
#mpirun -n 60 /scratch/${PBS_JOBID}/ElmerSolver_mpi



