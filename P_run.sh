#!/bin/tcsh


set PWD = `pwd`
#########################################################
## Test S Wave
#########################################################

set ID = J32
set NUM = 2000
set iteration = 5
set starting_model = S40RTS
cat $PWD/back/eventinfo.S_Sdiff |head -n ${NUM} |awk '{$19=25; print $0}' > $PWD/LSM_record_input
sed -i "/MODEL_NAME/c\<MODEL_NAME> ${starting_model}" $PWD/INFILE
sed -i "/Iteration_M/c\<Iteration_M> ${iteration}" $PWD/INFILE
csh $PWD/mother.sh ${ID} &
echo "--------------> Working on ID $ID RecordNUM: $NUM starting model: $starting_model"
sleep 1s
exit 0

#########################################################
## Test P Wave
#########################################################

set ID = J28
set NUM = 1500
set iteration = 1
set starting_model = S40RTS
cat $PWD/back/eventinfo.P.G123 |grep Pdiff |head -n ${NUM} > $PWD/LSM_record_input
sed -i "/MODEL_NAME/c\<MODEL_NAME> ${starting_model}" $PWD/INFILE
sed -i "/Iteration_M/c\<Iteration_M> ${iteration}" $PWD/INFILE
sed -i "/dv_type/c\<dv_type> vp" $PWD/INFILE
csh $PWD/mother.sh ${ID} &
echo "--------------> Working on ID $ID RecordNUM: $NUM starting model: $starting_model"
sed -i "/dv_type/c\<dv_type> vs" $PWD/INFILE
sleep 1s
exit 0



