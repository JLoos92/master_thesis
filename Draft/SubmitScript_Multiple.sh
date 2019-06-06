#!/bin/bash -l

# define number of runs for node and cpu
NUMBERS=$(seq 1 10)
AMPLITUDES=$(seq 50 10 500)
for num in ${NUMBERS}
do
    for amp in ${AMPLITUDES}
    do
        #PBS -N Channel_2D_${amp}
        #PBS -m abe
        #PBS -M julius-loos@gmx.de
        #PBS -o ${PBS_JOBNAME}${PBS_JOBID}.out
        #PBS -e ${PBS_JOBNAME}${PBS_JOBID}.err
        #PBS -l walltime=999:99:00
        #PBS -l nodes=1:ppn=2:esd1
        #PBS -q esd1
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

        export PATH=/home-link/epioi01/INSTALL/mmg-5.3.10/MMG_devel/bin:$PATH
        export LD_LIBRARY_PATH=/home-link/epioi01/INSTALL/mmg-5.3.10/MMG_devel/lib:$LD_LIBRARY_PATH

        export LD_LIBRARY_PATH=/home-link/epioi01/INSTALL/MUMPS_5.1.2/lib:$LD_LIBRARY_PATH
        export PATH=/home-link/epioi01/INSTALL/gmsh-3.0.6-source/build:$PATH
        export LD_LIBRARY_PATH=/beegfs/home/tu/epioi01/INSTALL/MUMPS_5.1.2/lib:$LD_LIBRARY_PATH

        #=================================================================================================================
        echo Here comes the partition the job runs in:
        echo $PBS_QUEUE
        cd $PBS_O_WORKDIR


        ncpus=`wc -l < ${PBS_NODEFILE}`
        echo Here comes the Nodelist:${PBS_NODEFILE}}
        echo Number of Nodes: ${ncpus}}


        find . -type f -name "channel2dIni.sif" | xargs sed -i '' '26s/.*/  _lower_surface = -rhoi\/rhoo*('${amp}' - 100.0*exp(-(x*x)\/1000000.0))\\/g' channel2dIni.sif
        find . -type f -name "channel2dIni.sif" | xargs sed -i '' '30s/.*/  _upper_surface = (1.0-rhoi\/rhoo)*('${amp}' - 100.0*exp(-(x*x)\/1000000.0))\\/g' channel2dIni.sif


#cp $ELMER_HOME/share/elmersolver/lib/FreeSurfaceSolver.so src/MyFreeSurf:xaceSolver.so
        make compile
        make ini
        make grid
        #mpirun -np $ncpus -machinefile $PBS_NODEFILE ElmerSolver_mpi
        mpirun -n $ncpus ElmerSolver_mpi
        #make submit
        #mpirun -n 60 /scratch/Elmer
        #mpirun -n 60 /scratch/${PBS_JOBID}/ElmerSolver_mpi
        done
done

