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
mkdir -p //volumes/esd01/docs/jloos/data_small/runs_elmerice/Mesh$x 

echo "----------------------------------------------------------------------------------
	New folder Mesh$x can be found on ESD1 FileServer
------------------------------------------------------------------------------------"

scp -r Cluster:/beegfs/work/zxmjf89/Mesh //volumes/esd01/docs/jloos/data_small/runs_elmerice/Mesh$x/ 
