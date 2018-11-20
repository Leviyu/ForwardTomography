#!/bin/tcsh


set PWD = `pwd`


#set MODELS = ( CHECK10time CHECK20time CHECK10time)
set MODELS = ( CHECK10time  CHECK10time)
set old_event = ../back/eventinfo.new

foreach MODEL (`echo $MODELS`)
	echo "--> on $MODEL"


	set prediction_file = /DATA1/ForwardTomography/WORKDIR/${MODEL}/prediction.data
	set p_size = `cat $prediction_file|wc -l`
	echo "size $p_size"
	set tmp_time = .tmp.tt
	cat $prediction_file|awk '{print $6}' >! $tmp_time

	set new_event = ./event.${MODEL}

	#cat $prediction_file
	paste $old_event $tmp_time |awk ''

	

end 

























set MODEL = S40RTS
set old_event = ../back/eventinfo.new
set s1 = `cat $old_event|wc -l`
echo "old size $s1"
