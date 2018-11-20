#!/bin/tcsh




set PWD = `pwd`
set big_file = $PWD/back/eventinfo.july18.csv


cd $PWD/tmp
cp $big_file eventinfo


divide_file eventinfo 7


set NUM = 0
set NUMMAX = 7

while( $NUM < $NUMMAX)
echo "-> on $NUM / $NUMMAX"
	cd $PWD
	cp $PWD/tmp/.eventinfo.${NUM} $PWD/LSM_record_input
	csh mother.sh Predict${NUM} &
	sleep 2s

@ NUM ++
end 






