#!/bin/tcsh


set PWD = `pwd`

#set deg_list = (10 20 30 )
set deg_list = (20)

foreach deg (`echo $deg_list`)

####################################################3
#sleep 5h
#echo "--> run for checker $deg "
# run1 use CHECK10time to run for all
#cat $PWD/back/event.CHECK20time |head -n 30000 > $PWD/LSM_record_input
#sed -i "/MODEL_NAME/c\<MODEL_NAME> EMPTY" $PWD/INFILE
#csh $PWD/mother.sh CHECK20run30k



############# For ALL ##################
#cat $PWD/back/event.CHECK20time  > $PWD/LSM_record_input
#sed -i "/MODEL_NAME/c\<MODEL_NAME> EMPTY" $PWD/INFILE
#csh $PWD/mother.sh CHECK${deg}runall

end  #############


set deg = S40RTS
############# For ALL ##################
cat $PWD/back/eventinfo.S40RTS |head -n 20000 > $PWD/LSM_record_input
sed -i "/MODEL_NAME/c\<MODEL_NAME> EMPTY" $PWD/INFILE
csh $PWD/mother.sh CHECKS40RTS20k

#sleep 10h
############# For ALL ##################
#cat $PWD/back/eventinfo.S40RTS  > $PWD/LSM_record_input
#sed -i "/MODEL_NAME/c\<MODEL_NAME> EMPTY" $PWD/INFILE
#csh $PWD/mother.sh CHECKS40RTSall

