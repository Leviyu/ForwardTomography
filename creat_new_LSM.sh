#!/bin/tcsh

set PWD = `pwd`
set bigevent = $PWD/back/eventinfo.new

set event = $PWD/LSM_record_input
cat /dev/null >! $event


foreach PHASE (S SS SSS ScS ScSScS Sdiff)
	echo "---> On $PHASE"
	cat $bigevent |grep -w $PHASE |head -n 100 >> $event
end 

