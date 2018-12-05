#!/bin/tcsh


set PWD = `pwd`

#########################################################
## Test Run with Data
#########################################################

set ID = P4
set NUM = 100
set starting_model = GYPSUM_P
cat $PWD/back/eventinfo.P.G123 |head -n ${NUM} > $PWD/LSM_record_input
sed -i "/MODEL_NAME/c\<MODEL_NAME> ${starting_model}" $PWD/INFILE
csh $PWD/mother.sh ${ID} &
echo "--------------> Working on ID $ID RecordNUM: $NUM starting model: $starting_model"
sleep 1s
exit 0


