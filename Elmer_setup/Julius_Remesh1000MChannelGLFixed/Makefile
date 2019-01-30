
# Makefile for Elmer
# ----------------------------------------
# Use external Geometry to create mesh
# Calculate Depth and Height for Paraview

EXECUTABLES = src/GroundedSolver src/Grid2DInterpolator src/USF_Contact src/USFs src/GroundedMaskSolver src/FlotingShelves src/CostSolver_DiscRD src/Optimize_m1qn3Parallel_Multi src/Optimize_m1qn3Parallel




NumProcs=20
NumProcx=4
NumProcy=5

#InputSif=Mismip3DInfluxRemesh.sif
InputSif=Mismip3DInfluxRestartBedrockBumpInElmer.sif



.SUFFIXES: .f90

all: clean ini grid submit

grid:
ifeq ($(NumProcs),1)
	gmsh Mesh.geo -1 -2
	ElmerGrid 14 2 Mesh.msh  -autoclean
else
	gmsh Mesh.geo -1 -2
	ElmerGrid 14 2 Mesh.msh  -autoclean
	ElmerGrid 2 2 Mesh -partition $(NumProcx) $(NumProcy) 1
endif


submit: ini

	mpirun -n $(NumProcs) ElmerSolver_mpi

drag:
	sed '1d' results/OutputInversionNorm1e51e5.csv > tmp
	matlab -nosplash -nodesktop -r InversionInterpolationScript
	rm tmp
#compile:  $(EXECUTABLES)
compile:
	#make --file MakefileJoh all
	#elmerf90 src/FlotingShelves.f90 -o src/FlotingShelves
	#elmerf90 src/MyFreeSurfaceSolver.f90 -o src/MyFreeSurfaceSolver
	elmerf90 src/DistanceSolveRD.f90 -o src/DistanceSolveRD
	elmerf90 src/USF_BMB.f90 -o src/USF_BMB
	elmerf90 src/BedrockBump.f90 -o src/BedrockBump
	elmerf90 src/SurfaceBump.f90 -o src/SurfaceBump

clean:
	rm -fr results/*

ini:
	echo $(InputSif) > ELMERSOLVER_STARTINFO

.f90:
	elmerf90  -o $@ $<
.c:
	gcc  -o $@ $< -lm

