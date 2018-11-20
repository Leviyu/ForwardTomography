#!/bin/tcsh

# we use the prediction travel time made from a given to replace the dt of its input
set PWD = `pwd`

#set ID = CHECK10time
set ID = CHECK20time
#set ID = CHECK30time


# line 6 is time
set pred = /DATA1/ForwardTomography/WORKDIR/${ID}/prediction.data
# 19 is time
set input_file = /DATA1/ForwardTomography/WORKDIR/$ID/LSM_record_input


set file_line = 47

set outfile = $PWD/event.${ID}

set tmp1 = $PWD/.tmp1
cat $pred|awk '{print $6}' >! $tmp1

paste $input_file $tmp1 |awk '{$19=$48;print $0}'|awk '{$48="";print $0}' >! $outfile


