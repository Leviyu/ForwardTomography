#!/bin/tcsh


set PWD = `pwd`

#########################################################
## Test Run with Data
#########################################################

set ID = T400
set NUM = 300000
set starting_model = S40RTS
cat $PWD/back/eventinfo.data |head -n ${NUM} > $PWD/LSM_record_input
sed -i "/MODEL_NAME/c\<MODEL_NAME> ${starting_model}" $PWD/INFILE
csh $PWD/mother.sh ${ID} &
echo "--------------> Working on ID $ID RecordNUM: $NUM starting model: $starting_model"
sleep 5s



set ID = T401
set starting_model = GYPSUM_S
cat $PWD/back/eventinfo.data |head -n ${NUM} > $PWD/LSM_record_input
sed -i "/MODEL_NAME/c\<MODEL_NAME> ${starting_model}" $PWD/INFILE
csh $PWD/mother.sh ${ID} &
echo "--------------> Working on ID $ID RecordNUM: $NUM starting model: $starting_model"
sleep 5s

set ID = T402
set starting_model = SEMUCB_WM1
cat $PWD/back/eventinfo.data |head -n ${NUM} > $PWD/LSM_record_input
sed -i "/MODEL_NAME/c\<MODEL_NAME> ${starting_model}" $PWD/INFILE
csh $PWD/mother.sh ${ID} &
echo "--------------> Working on ID $ID RecordNUM: $NUM starting model: $starting_model"
sleep 5s
exit 0

#########################################################
## Benchmark test
#########################################################

set ID = B25
set NUM = 1000
set starting_model = EMPTY
cat $PWD/back/event.CHECK20time |head -n ${NUM} > $PWD/LSM_record_input
sed -i "/MODEL_NAME/c\<MODEL_NAME> ${starting_model}" $PWD/INFILE
csh $PWD/mother.sh ${ID} &
echo "--------------> Working on ID $ID RecordNUM: $NUM starting model: $starting_model"
sleep 1s

######################################
set ID = K4
set NUM = 20000
set starting_model = S40RTS
cat $PWD/back/eventinfo.data |grep -w Sdiff |head -n ${NUM} > $PWD/LSM_record_input
sed -i "/MODEL_NAME/c\<MODEL_NAME> ${starting_model}" $PWD/INFILE
csh $PWD/mother.sh ${ID} &
echo "--------------> Working on ID $ID RecordNUM: $NUM starting model: $starting_model"
sleep 1s
exit 0
#######################################

exit 0

