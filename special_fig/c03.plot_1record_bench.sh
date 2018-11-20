#!/bin/tcsh


set PWD = `pwd`

set event = ../LSM_record_input.all_phase

set NUM = 1
set NUMMAX = 120

set outfile = $PWD/WORKDIR/out.c03
cat /dev/null >! $outfile


while($NUM <= $NUMMAX )


	set dt1 = `cat $event|awk 'NR=='$NUM' {print $19}'`
	set dt1_cor = `cat ../LOG/logfile.P${NUM} |grep correction |awk 'NR==1 {print $4}'`
	set dt2 = `cat ../LOG/logfile.Q${NUM} |grep correction |awk 'NR==1 {print $4}'`


	set dif = `echo "$dt1 - $dt1_cor - $dt2"|bc -l`

	#echo $dt1 $dt1_cor $dt2 $dif

	echo $dt1 $dt2 >> $outfile
	

	#sleep 2s



@ NUM ++
end

