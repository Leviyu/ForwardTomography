#!/bin/tcsh

set PWD = `pwd`

set INFILE = $PWD/INFILE

set event = $PWD/LSM_record_input.all_phase

set NUM = 1
while( $NUM <= 120)

	# put record inside 
	cat $event |awk 'NR=='$NUM' {print $0}' > ! $PWD/LSM_record_input
	sed -i 's/<MODEL_NAME>.*/<MODEL_NAME>     TEST/g' $PWD/INFILE
	csh $PWD/mother.sh P${NUM}

	mkdir -p ~/Tomography/P${NUM}
	set target = ~/Tomography/P${NUM}/P${NUM}.standard.model
	cp /DATA1/ForwardTomography/WORKDIR/P${NUM}/TEST.standard.model.2 $target

	cd $PWD
	sed -i 's/<MODEL_NAME>.*/<MODEL_NAME>     P'${NUM}'/g' $PWD/INFILE
	csh $PWD/mother.sh Q${NUM} &
	sleep 2s


@ NUM ++
end 


