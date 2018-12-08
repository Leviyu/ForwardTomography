#!/bin/tcsh


set PWD = `pwd`
#########################################################
## Test S Wave
#########################################################

set ID = K14
set NUM = 100
set iteration = 3
set starting_model = GYPSUM_S
##cat $PWD/back/eventinfo.S_Sdiff |head -n ${NUM}  > $PWD/LSM_record_input
cat $PWD/back/eventinfo.S40RTS |grep -w Sdiff |head -n ${NUM}  > $PWD/LSM_record_input
sed -i "/MODEL_NAME/c\<MODEL_NAME> ${starting_model}" $PWD/INFILE
sed -i "/Iteration_M/c\<Iteration_M> ${iteration}" $PWD/INFILE
csh $PWD/mother.sh ${ID} &
echo "--------------> Working on ID $ID RecordNUM: $NUM starting model: $starting_model"
sleep 1s
exit 0



#########################################################
## Test P Wave
#########################################################

set ID = J122
set NUM = 15000
set iteration = 5
set starting_model = GYPSUM_P
cat $PWD/back/eventinfo.P.G123 |awk '$3>80 {print $0}'|head -n ${NUM} > $PWD/LSM_record_input
sed -i "/MODEL_NAME/c\<MODEL_NAME> ${starting_model}" $PWD/INFILE
sed -i "/Iteration_M/c\<Iteration_M> ${iteration}" $PWD/INFILE
sed -i "/dv_type/c\<dv_type> vp" $PWD/INFILE
csh $PWD/mother.sh ${ID} &
echo "--------------> Working on ID $ID RecordNUM: $NUM starting model: $starting_model"
sed -i "/dv_type/c\<dv_type> vs" $PWD/INFILE
sleep 1s
exit 0

