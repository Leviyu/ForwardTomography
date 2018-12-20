#!/bin/tcsh


set PWD = `pwd`
#########################################################
## Test P Wave
#########################################################

set ID = H82
set NUM = 15000
set iteration = 3
set starting_model = GYPSUM_P
cat $PWD/back/eventinfo.P.G123 |awk '$3>80 {print $0}'|head -n ${NUM} > $PWD/LSM_record_input
sed -i "/MODEL_NAME/c\<MODEL_NAME> ${starting_model}" $PWD/INFILE
sed -i "/Iteration_M/c\<Iteration_M> ${iteration}" $PWD/INFILE
sed -i "/dv_type/c\<dv_type> vp" $PWD/INFILE
csh $PWD/mother.sh ${ID} &
echo "--------------> Working on ID $ID RecordNUM: $NUM starting model: $starting_model"
sleep 5s

#########################################################
## Test S Wave
#########################################################

set ID = H83
set NUM = 150000
set iteration = 3
set starting_model = GYPSUM_S
cat $PWD/back/eventinfo.S_Sdiff |head -n ${NUM}  > $PWD/LSM_record_input
##cat $PWD/back/eventinfo.S40RTS |grep -w Sdiff|awk '$3>120 {print $0}' |head -n ${NUM}  > $PWD/LSM_record_input
##cat $PWD/back/eventinfo.S40RTS |grep -w S|awk '$3>90 {print $0}' |head -n ${NUM}  > $PWD/LSM_record_input
sed -i "/MODEL_NAME/c\<MODEL_NAME> ${starting_model}" $PWD/INFILE
sed -i "/Iteration_M/c\<Iteration_M> ${iteration}" $PWD/INFILE
sed -i "/dv_type/c\<dv_type> vs" $PWD/INFILE
csh $PWD/mother.sh ${ID} &
echo "--------------> Working on ID $ID RecordNUM: $NUM starting model: $starting_model"
sleep 1s
exit 0



