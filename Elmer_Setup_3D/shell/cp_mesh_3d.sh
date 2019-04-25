#!/bin/zsh -l


###----------------------------------------------------------------------------#
###
###----------------------------------------------------------------------------#


# Open FileServer ESD1, make new dir and copy mesh
#open 'smb://134.2.5.43/esd01'


amp="$AMP"
AMP=" "
echo "Amplitude,sigma and GL position for modelrun Please define number for dirname (e.g. 250 for maxAmplitude,  for sigma 200 and 200 = Mesh250_500500_0)"
read a x y dl
AMP=$amp

PS3='Please choose directory:'
options=("fixed_GL folder" "uniform-mesh folder" "standard/test" "Quit")
select opt in "${options[@]}"
do
    case $opt in
        "fixed_GL folder")
			mkdir -p //volumes/esd01/docs/jloos/data_small/runs_elmerice_fixed/Mesh"$a"_"$x""$y"_"$dl" | rsync -ruvt Cluster:/beegfs/work/zxmjf89/Remesh1000mChannelGLFixed/Mesh //volumes/esd01/docs/jloos/data_small/runs_elmerice_fixed/Mesh"$a"_"$x""$y"_"$dl"/
            break
            ;;
        "refined folder")
			mkdir -p //volumes/esd01/docs/jloos/data_small/runs_elmerice_refined/Mesh"$a"_"$x""$y"_"$dl" | rsync -ruvt Cluster:/beegfs/work/zxmjf89/Remesh1000mChannelGLFixed/Mesh //volumes/esd01/docs/jloos/data_small/runs_elmerice_refined/Mesh"$a"_"$x""$y"_"$dl"/
            break
            ;;
        "standard/test")
            mkdir -p //volumes/esd01/docs/jloos/data_small/runs_elmerice/Mesh"$a"_"$x""$y"_"$dl" | rsync -ruvt Cluster:/beegfs/work/zxmjf89/Remesh1000mChannelGLFixed/Mesh //volumes/esd01/docs/jloos/data_small/runs_elmerice_uniform/Mesh"$a"_"$x""$y"_"$dl"/
            break
            ;;
                "Quit")
            break
            ;;
        *) echo "invalid option $REPLY";;
    esac
done



echo "----------------------------------------------------------------------------------
	New folder Mesh"$a"_"$x""$y"_"$dl" can be found on ESD01 FileServer under data_small
------------------------------------------------------------------------------------"

 
