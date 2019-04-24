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
Please type in the amplitude (from 140 to 500), second the distribution in x direction, third the distribution in y direction and fourth the distance (m)\n
 towards icesheet from the groundingline (0 is peak at GL). Seperation by ; ... (e.g. 200;1000;2000;0)
 ----------------------------------------------------------------------------------"
 
read a x y dl
IFS=$old_IFS

echo "Amplitude = $a for the given BED topography. X and Y direction $x + $y for the distribution. Distance to GL = $dl"


# Find and replace line in .sif file for given amplitude $x (changes BED)
find . -type f -name "Mismip3DInfluxRemesh.sif" | xargs sed -i '' '179s/.*/   Variable 1 data file = File ".\/DEM\/BED_bump'$a'_'$x''$y'_'$dl'.xyz"/g' Mismip3DInfluxRemesh.sif 
find . -type f -name "Mismip3DInfluxRemesh.sif" | xargs sed -i '' '195s/.*/   Variable 2 data file = File ".\/DEM\/ZB_bump'$a'_'$x''$y'_'$dl'.xyz"/g' Mismip3DInfluxRemesh.sif | rsync -ruvt 'Mismip3DInfluxRemesh.sif' 'Cluster:/beegfs/work/zxmjf89/Mismip3DSetUpSteadyState_Remesh/'
#find . -type d -name "Cluster:/beegfs/work/zxmjf89/Mismip3DInflux.sif" | xargs sed -i '' '181s/.*/Variable 1 data file = File ".\/DEM\/BED_bump'$x'.xyz"/g' Mismip3DInflux.sif
#ssh Cluster "find /beegfs/work/zxmjf89 -name Mismip3DInflux.sif" | xargs sed -i '' '181s/.*/Variable 2 data file = File ".\/DEM\/ZB_bump'$x'.xyz"/g' SubmitScriptESD1.sh

find . -type f -name "SubmitScriptESD1.sh" |  xargs sed -i '' '6s/.*/#PBS -N Bump'$a'_'$x''$y'_'$dl'/g' SubmitScriptESD1.sh | rsync -ruvt 'SubmitScriptESD1.sh' 'Cluster:/beegfs/work/zxmjf89/Mismip3DSetUpSteadyState_Remesh/'

# Find and replace line in .mismip_geoscript.py and creates new .py and .xyz BED file
find . -type f -name "Rearrange_bump.py" | xargs sed -i '' '81s/.*/maxAmplitude = '$a'/g' Rearrange_bump.py # > mismip_geoscript$x.py
find . -type f -name "Rearrange_bump.py" | xargs sed -i '' '82s/.*/sigmax = '$x'/g' Rearrange_bump.py
find . -type f -name "Rearrange_bump.py" | xargs sed -i '' '83s/.*/sigmay = '$y'/g' Rearrange_bump.py # > mismip_geoscript$x.py
find . -type f -name "Rearrange_bump.py" | xargs sed -i '' '86s/.*/dl = '$dl'/g' Rearrange_bump.py # > mismip_geoscript$x.py

python3 Rearrange_bump.py

rsync -ruvt 'BED_bump'$a'_'$x''$y'_'$dl'.xyz' 'Cluster:/beegfs/work/zxmjf89/Mismip3DSetUpSteadyState_Remesh/DEM/'
rsync -ruvt 'ZB_bump'$a'_'$x''$y'_'$dl'.xyz' 'Cluster:/beegfs/work/zxmjf89/Mismip3DSetUpSteadyState_Remesh/DEM/'
#scp -r Mismip3DInflux.sif Cluster:/beegfs/work/zxmjf89/


#scp -r Cluster:/beegfs/work/zxmjf89/Mesh /Users/Schmulius/Desktop




