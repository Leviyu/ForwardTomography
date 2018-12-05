#!/bin/tcsh



set PWD = $1

set ID = `cat $PWD/INFILE |grep TASK_NAME |awk '{print $2}'`
set WORKDIR = `cat $PWD/INFILE |grep WORKDIR |awk '{print $2}'`
if($ID == "") then
	echo "ERROR ID is empty!"
exit 0
endif

mkdir -p $WORKDIR

set DATADIR = $WORKDIR/$ID
set PLOTDIR = $PWD/PLOTDIR/$ID
/bin/rm -rf $DATADIR
/bin/rm -rf $PLOTDIR
mkdir -p $DATADIR
mkdir -p $PLOTDIR


## copy INFILE and INPUT into WORKDIR
cp $PWD/INFILE $DATADIR/
cp $PWD/LSM_record_input $DATADIR/


set logfile = $PWD/LOG/logfile.${ID}
cat /dev/null >! $logfile

#echo " --> running task for $ID"
#echo "WORKDIR is $WORKDIR"
#csh $PWD/code_dir/c00.forward_tomography.sh $PWD $ID
csh $PWD/code_dir/c00.forward_tomography.sh $PWD $ID >> & $logfile 
##csh $PWD/code_dir/c00.forward_tomography.sh $PWD $ID



# Convert Tomo into standard format
#csh $PWD/code_dir/c02.convert_model_into_standard_format.sh $PWD $ID >>& $logfile &

# to hongyu
##to_hongyu $PLOTDIR/$ID

sleep 3s


