#!/bin/tcsh


set in = ../LSM_record_input.1000
set out = ../LSM_record_input


cat /dev/null >! $out

foreach PHASE (S SS SSS ScSScS Sdiff ScS)

cat $in |grep -w $PHASE |head -n 20 >> $out
echo "$PHASE"
cat $in |grep -w $PHASE |head -n 20 |wc -l

end 
