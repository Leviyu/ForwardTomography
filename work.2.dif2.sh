#!/bin/tcsh

set PWD = `pwd`



##foreach MODEL ( GYPSUM_S SAW24B16 HMSL_S06 S40RTS S362WANI )
foreach MODEL (  SAW24B16 HMSL_S06 S40RTS S362WANI GYPSUM_S )
##foreach MODEL (  SAW24B16 )


	cat $PWD/INFILE_try |awk -v dd=$MODEL '{if($1 == "<MODEL_NAME>") print $1, dd; else print $0}' > ! $PWD/INFILE

	set ID = Sep16_MODEL_${MODEL}

	csh $PWD/c00.forward_tomography.sh $ID > & logfile.${ID}

end 

