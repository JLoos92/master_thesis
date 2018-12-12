#!/bin/zsh -l

###----------------------------------------------------------------------------#
### Copy, read and write for new gauss-bump topography 
###----------------------------------------------------------------------------#

#scp -r /Users/Schmulius/Desktop/UniMaster/Thesis/Clemens_Grid_Bump/BED_Bump400.xyz Cluster:/beegfs/work/zxmjf89/DEM/ 

#find ~ -name Mismip3DInflux.sif and type name of bump.xyz
#DIR = "$/Users/Schmulius/Desktop/UniMaster/Thesis/Clemens_Grid_Bump/"


old_IFS="$IFS"
IFS=\;


echo "----------------------------------------------------------------------------------
 Type in the amplitude (from 140 to 500), second the distribution in x direction, third the distribution in y direction and fourth the distance (m)\n
 towards icesheet from the groundingline (0 is peak at GL). Seperation by ; ... (e.g. 200;1000;2000;0)
 ----------------------------------------------------------------------------------"
 
read a x y dl
IFS=$old_IFS

echo "Amplitude = $a for the given BED topography. X and Y direction $x + $y for the distribution. Distance to GL = $dl"
Mismip3DSetUpSteadyState_UniformMesh

# Find and replace line in .sif file for given amplitude $x (changes BED)
find . -type f -name "Mismip3DSetUpSteadyState_UniformMesh/Mismip3DInflux.sif" | xargs sed -i '' '181s/.*/Variable 2 data file = File ".\/DEM\/BED_bump'$a''$x''$y'_'$dl'.xyz"/g' Mismip3DInflux.sif 
find . -type f -name "Mismip3DSetUpSteadyState_UniformMesh/Mismip3DInflux.sif" | xargs sed -i '' '198s/.*/Variable 2 data file = File ".\/DEM\/ZB_bump'$a''$x''$y'_'$dl'.xyz"/g' Mismip3DInflux.sif | scp -r Mismip3DSetUpSteadyState_UniformMesh/Mismip3DInflux.sif Cluster:/beegfs/work/zxmjf89/Mismip3DSetUpSteadyState_UniformMesh/
#find . -type d -name "Cluster:/beegfs/work/zxmjf89/Mismip3DInflux.sif" | xargs sed -i '' '181s/.*/Variable 1 data file = File ".\/DEM\/BED_bump'$x'.xyz"/g' Mismip3DInflux.sif
#ssh Cluster "find /beegfs/work/zxmjf89 -name Mismip3DInflux.sif" | xargs sed -i '' '181s/.*/Variable 2 data file = File ".\/DEM\/ZB_bump'$x'.xyz"/g' SubmitScriptESD1.sh

find . -type f -name "Mismip3DSetUpSteadyState_UniformMesh/SubmitScriptESD1.sh" |  xargs sed -i '' '7s/.*/#PBS -N Bump'$a''$x''$y'_msh/g' SubmitScriptESD1.sh | scp -r SubmitScriptESD1.sh Cluster:/beegfs/work/zxmjf89/


# Find and replace line in .mismip_geoscript.py and creates new .py and .xyz BED file
find . -type f -name "Mismip3DSetUpSteadyState_UniformMesh/Rearrange_bump.py" | xargs sed -i '' '81s/.*/maxAmplitude = '$a'/g' Rearrange_bump.py # > mismip_geoscript$x.py
find . -type f -name "Mismip3DSetUpSteadyState_UniformMesh/Rearrange_bump.py" | xargs sed -i '' '82s/.*/sigmax = '$x'/g' Rearrange_bump.py
find . -type f -name "Mismip3DSetUpSteadyState_UniformMesh/Rearrange_bump.py" | xargs sed -i '' '83s/.*/sigmay = '$y'/g' Rearrange_bump.py # > mismip_geoscript$x.py
find . -type f -name "Mismip3DSetUpSteadyState_UniformMesh/Rearrange_bump.py" | xargs sed -i '' '86s/.*/dl = '$dl'/g' Rearrange_bump.py # > mismip_geoscript$x.py

python Mismip3DSetUpSteadyState_UniformMesh/Rearrange_bump.py

rsync -ruvt Mismip3DSetUpSteadyState_UniformMesh/DEM/BED_bump$a$x$y_$dl.xyz Cluster:/beegfs/work/zxmjf89/Mismip3DSetUpSteadyState_UniformMesh/DEM/
rsync -ruvt Mismip3DSetUpSteadyState_UniformMesh/DEM/ZB_bump$a$x$y_$dl.xyz Cluster:/beegfs/work/zxmjf89/Mismip3DSetUpSteadyState_UniformMesh/DEM/
#scp -r Mismip3DInflux.sif Cluster:/beegfs/work/zxmjf89/


#scp -r Cluster:/beegfs/work/zxmjf89/Mesh /Users/Schmulius/Desktop




