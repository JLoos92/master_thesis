#!/bin/zsh -l


###----------------------------------------------------------------------------#
### Structure model runs of Elmer/Ice
###----------------------------------------------------------------------------#


# Open FileServer ESD1, make new dir and copy mesh
open 'smb://134.2.5.43/esd01'


amp="$AMP"
AMP=:
echo "\n Which amplitude and sigma did you use? Please define number for dirname (e.g. 200 for maxAmplitude 1000 for sigma = Mesh2001000)"
read x
AMP=$amp

PS3='Please enter your choice: '
options=("Option 1: Pick uniform-mesh folder" "Option 2: Pick refined-mesh folder" "Quit")
select opt in "${options[@]}"
do
    case $opt in
        "Option 1: Save in uniform-mesh folder")
            #echo "Folder for uniform-mesh"
			rsync -ruvt Cluster:/beegfs/work/zxmjf89/Mismip3DSetUpSteadyState_UniformMesh/Mesh_uni //volumes/esd01/docs/jloos/data_small/runs_elmerice_uniform/Mesh$x/ 
            ;;
        "Option 2: Save in refined-mesh folder")
            echo "Folder for refined-mesh"
			rsync -ruvt Cluster:/beegfs/work/zxmjf89/Mismip3DSetUpSteadyState_UniformMesh/Mesh_uni //volumes/esd01/docs/jloos/data_small/runs_elmerice_refined/Mesh$x/ 
            ;;
                "Quit")
            break
            ;;
        *) echo "invalid option $REPLY";;
    esac
done



mkdir -p //volumes/esd01/docs/jloos/data_small/runs_elmerice_uniform/Mesh$x 

echo "----------------------------------------------------------------------------------
	New folder Mesh$x can be found on ESD1 FileServer under data_small
------------------------------------------------------------------------------------"

#rsync -ruvt Cluster:/beegfs/work/zxmjf89/Mismip3DSetUpSteadyState_UniformMesh/Mesh_uni //volumes/esd01/docs/jloos/data_small/runs_elmerice_uniform/Mesh$x/ 
