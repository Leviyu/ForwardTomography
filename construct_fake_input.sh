#!/bin/tcsh



set PWD = `pwd`
set event = $PWD/LSM_record_input


set tmp = $PWD/.tmp.event
mv $event $tmp



set outfile = /DATA1/ForwardTomography/WORKDIR/T21/prediction.data

set NUM = 1
set NUMMAX = `cat $tmp|wc -l`

while($NUM <= $NUMMAX)
set model_dt = `cat $outfile |awk 'NR=='$NUM' {print $6}'`
	echo "$NUM / $NUMMAX $model_dt"
cat $tmp |awk -v dd=$model_dt 'NR=='$NUM' {$19=dd; print $0}' >> $event

@ NUM ++
end
