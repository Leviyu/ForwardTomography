#!/bin/tcsh


set PWD = `pwd`
set INFILE = $PWD/INFILE

set ID = $1

if($ID == "" ) then
echo " ERROR ID is empty!"
exit 0
endif

sed -i '/TASK_NAME/c\<TASK_NAME> '$ID'' $INFILE


csh $PWD/code_dir/work.sh $PWD 
