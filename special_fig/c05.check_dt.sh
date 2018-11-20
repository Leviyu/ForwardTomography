#!/bin/tcsh
# plot the difference between taup+inc and travel time 

set PWD = `pwd`


#set DIR = /DATA1/ForwardTomography/WORKDIR/

foreach inc (0.005 0.01 0.02 0.05 0.1 0.2 0.3)
set ID = INC4${inc}

set log = ~/ForwardTomography/LOG/logfile.${ID}

#set tmp = $PWD/WORKDIR/tmp.${inc}

echo "---> inc $inc"
cat $log |grep EQ_STA_PHASE  |awk '{print $9}'
#cat $log |grep "max dl"




end 



