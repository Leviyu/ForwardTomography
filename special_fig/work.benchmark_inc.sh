#!/bin/tcsh

set PWD = `pwd`

set INFILE = $PWD/INFILE


foreach inc (0.005 0.01 0.02 0.05 0.1 0.2 0.3)
#foreach inc ( 0.005 )

# 0. change inc
sed -i "s/taup.path.maxPathInc.*/taup.path.maxPathInc=${inc}/g" ~/.taup

# 1. go to TAUP path dir and make taup path
cd /DATA1/ForwardTomography/TAUP_PATH_DIR
csh mother.sh 

sleep 5m

# 2. go to forward tomography and run
cd ~/ForwardTomography
set ID = INC4${inc}
csh mother.sh $ID &

end 


