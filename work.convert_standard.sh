#!/bin/tcsh

set PWD = `pwd`

set INFILE = $PWD/INFILE





##foreach MODEL ( GYPSUM_S SAW24B16 HMSL_S06 S40RTS S362WANI )
#foreach MODEL (S40RTS  SAW24B16 HMSL_S06  S362WANI GYPSUM_S )
#foreach MODEL (  SAW24B16 )
foreach MODEL (  S40RTS )
	echo "--> convert format for $MODEL"
	set BIG_ID = Convert
	set ID = ${BIG_ID}_${MODEL}
	sed -i '/MODEL_NAME/c\<MODEL_NAME> 	'${MODEL}'' $PWD/INFILE
	csh $PWD/mother.sh $ID
	
	set target = ~/Tomography/${MODEL}/${MODEL}.standard.model
	/bin/rm -rf $target
	cp /DATA1/ForwardTomography/WORKDIR/${ID}/${MODEL}.standard.model $target

end 

