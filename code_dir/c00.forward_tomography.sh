#!/bin/tcsh

# ==================================================
#	This Script is the main scrip for layer stripping method (LSM)
#	
#	
#	Hongyu DATE: June 25 2014
#	Key words: Layer Strippint Method LSM
# ==================================================

# ===========================================================
#					Read in Input
set PWD = $1
set ID  = $2
if($ID == "" ) then
echo "---> input ID "
exit 0
endif


#cd $PWD/DATADIR/${ID}

set SRCDIR = `grep "SRCDIR" INFILE |awk '{print $2}'`
set WORKDIR = `grep "WORKDIR" INFILE |awk '{print $2}'`
set PLOTDIR = `grep "PLOTDIR" INFILE |awk '{print $2}'`
set MODEL = `grep "MODEL_NAME" INFILE |awk 'NR==1 {print $2}'`
set MODEL_DIR = `grep "MODEL_DIR" INFILE |awk '{print $2}'`

set PLOTDIR = $PLOTDIR/$ID
set WORKDIR = $WORKDIR/$ID

set INFILE = $WORKDIR/INFILE
set INPUT = $WORKDIR/LSM_record_input


# ========================================================
# 3. compile c code and get it run
echo "---> Compile c code and get it run" 
cd $PWD/Maligaro/cpp_lib/forward_tomography
make 
#make > & /dev/null 
cd $WORKDIR
cp $PWD/Maligaro/cpp_lib/forward_tomography/forward_tomography .



echo "======================================================="
echo "======================================================="
echo "++++++++++++++++++++ Begin time +++++++++++++++++ `date`"
echo "$WORKDIR/LSM_record_input $MODEL_DIR/$MODEL $INFILE $MODEL ">! input

./forward_tomography
echo "++++++++++++++++++++ End time +++++++++++++++++++ `date`"
echo "======================================================="
echo "======================================================="

## Now begin to plot
csh $PWD/code_dir/c01.plot.sh $PWD $ID



exit 0
