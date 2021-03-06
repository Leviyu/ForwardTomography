#!/bin/tcsh


set PWD = `pwd`

# ==========================================================
set ID = H13
set NUM = 300000
set starting_model = S40RTS
cat $PWD/back/eventinfo.data |head -n ${NUM} > $PWD/LSM_record_input
sed -i "/MODEL_NAME/c\<MODEL_NAME> ${starting_model}" $PWD/INFILE
csh $PWD/mother.sh ${ID} &
echo "--------------> Working on ID $ID RecordNUM: $NUM starting model: $starting_model"
sleep 5s

set ID = H14
set starting_model = GYPSUM_S
cat $PWD/back/eventinfo.data |head -n ${NUM} > $PWD/LSM_record_input
sed -i "/MODEL_NAME/c\<MODEL_NAME> ${starting_model}" $PWD/INFILE
csh $PWD/mother.sh ${ID} &
echo "--------------> Working on ID $ID RecordNUM: $NUM starting model: $starting_model"
sleep 5s

set ID = H15
set starting_model = SEMUCB_WM1
cat $PWD/back/eventinfo.data |head -n ${NUM} > $PWD/LSM_record_input
sed -i "/MODEL_NAME/c\<MODEL_NAME> ${starting_model}" $PWD/INFILE
csh $PWD/mother.sh ${ID} &
echo "--------------> Working on ID $ID RecordNUM: $NUM starting model: $starting_model"
sleep 5s
exit 0
# ================================================================

# ==========================================================
set ID = H2
set NUM = 20000
set iteration = 3
set starting_model = S40RTS
##cat $PWD/back/eventinfo.data |grep -w Sdiff |head -n ${NUM} > $PWD/LSM_record_input
cat $PWD/back/eventinfo.data | head -n ${NUM} > $PWD/LSM_record_input
sed -i "/MODEL_NAME/c\<MODEL_NAME> ${starting_model}" $PWD/INFILE
sed -i "/Iteration_M/c\<Iteration_M> ${iteration}" $PWD/INFILE
csh $PWD/mother.sh ${ID} &
echo "--------------> Working on ID $ID RecordNUM: $NUM starting model: $starting_model"
sleep 1s
exit 0


#########################################################
## Test Run with Data
#########################################################
set ID = T809Sdiff
set NUM = 1000
set starting_model = S40RTS
cat $PWD/back/eventinfo.data |grep -w Sdiff |head -n ${NUM} > $PWD/LSM_record_input
sed -i "/MODEL_NAME/c\<MODEL_NAME> ${starting_model}" $PWD/INFILE
csh $PWD/mother.sh ${ID} &
echo "--------------> Working on ID $ID RecordNUM: $NUM starting model: $starting_model"
sleep 1s
exit 0

#########################################################
## Benchmark test
#########################################################

set ID = B40
set NUM = 1000
set starting_model = EMPTY
cat $PWD/back/event.CHECK20time |head -n ${NUM} > $PWD/LSM_record_input
sed -i "/MODEL_NAME/c\<MODEL_NAME> ${starting_model}" $PWD/INFILE
csh $PWD/mother.sh ${ID} &
echo "--------------> Working on ID $ID RecordNUM: $NUM starting model: $starting_model"
sleep 1s
exit 0




set ID = T301
set starting_model = GYPSUM_S
cat $PWD/back/eventinfo.data |head -n ${NUM} > $PWD/LSM_record_input
sed -i "/MODEL_NAME/c\<MODEL_NAME> ${starting_model}" $PWD/INFILE
csh $PWD/mother.sh ${ID} &
echo "--------------> Working on ID $ID RecordNUM: $NUM starting model: $starting_model"
sleep 5s

set ID = T302
set starting_model = SEMUCB_WM1
cat $PWD/back/eventinfo.data |head -n ${NUM} > $PWD/LSM_record_input
sed -i "/MODEL_NAME/c\<MODEL_NAME> ${starting_model}" $PWD/INFILE
csh $PWD/mother.sh ${ID} &
echo "--------------> Working on ID $ID RecordNUM: $NUM starting model: $starting_model"
sleep 5s



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

