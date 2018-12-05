#!/bin/csh

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

cd $WORKDIR
set script = $WORKDIR/plot.sh
cat /dev/null >! $script

set INFILE = $WORKDIR/INFILE
set INPUT = $WORKDIR/LSM_record_input


### ===================================================
# plot histogram info of current iteration
### ===================================================
#echo " Plot histogram info of current iteration"
#cd $WORKDIR
#foreach histogram_iteration (`ls timeinfo*`)
#python $SRCDIR/c01.plot_histogram.py $histogram_iteration $PLOTDIR
#end


## plot diffferent depth shell for tomo
foreach STEP (1 2 3 )
echo "csh $SRCDIR/c12.plot_different_depth_shell.sh $PLOTDIR $WORKDIR $ID $STEP" >> $script 
end
echo "csh $SRCDIR/c31.plot_dt_hist_per_update.sh $PLOTDIR $WORKDIR $ID" >> $script
echo "csh $SRCDIR/c41.plot_variance_const_per_iteration.sh $PLOTDIR $WORKDIR $ID" >> $script
# plot cross-section
set event = $WORKDIR/LSM_record_input
#echo "csh $SRCDIR/c21.plot_vertical_slice_for_all_record.sh $SRCDIR  $event $PLOTDIR $WORKDIR" >> $script

#set DEPTH = 2800
#foreach FLAG ( NE_PACIFIC ALASKA PERM E_SAmerica E_NAmerica W_NMerica M_SAmerica E_Africa  )
#csh $SRCDIR/c05.zoom_in_pacific_horizontal_slice.sh  $DEPTH $PLOTDIR $WORKDIR $SRCDIR $FLAG
#end





#cd $WORKDIR
#set MODEL = `grep "MODEL_NAME" INFILE |awk 'NR==1 {print $2}'`
#### convert  hit count / tomo dvs and tomo delta dvs  into vtk format
##foreach flag ( output_STD_tomo hit_count output_tomo output_delta_tomo output_gradient output_vertical_gradient)
#foreach flag (output_tomo )
#set CONVERT = $PWD/convert.py
#set tmpfile = $WORKDIR/convert.tmp.file.${flag}
#set file = ${flag}.${MODEL}
#if(! -e $file ) then
#continue
#endif
#set count_file = $file
#echo "--> Convert $flag into vtk"
#set out_file = $PLOTDIR/out_${flag}
#set xx = `cat $count_file |awk 'NR==1 {print $1}'`
#set yy = `cat $count_file |awk 'NR==2 {print $1}'`
#set zz = `cat $count_file |awk 'NR==3 {print $1}'`
#cat $count_file |awk 'NR>3 {print $0}' > $tmpfile
#python $CONVERT $tmpfile  $out_file $xx $yy $zz &
#end

### ===================================================
# plot horizontal slice at CMB
### ===================================================
# NOTES: the depth here should be exactaly as what is provided in INFILE
##foreach depth (2500 2600 2700 2800 2890)
#foreach depth ( 2800 )
#csh $SRCDIR/c05.plot_horizontal_slice_for_starting_residual_final_tomo.sh  $depth $PLOTDIR $WORKDIR $SRCDIR
#end


### ===================================================
### convert pathinfo cross-section into vtk format
#echo "--> Convert pathinfo cross-section into vtk"
#cd $WORKDIR
#set CONVERT = $PWD/convert.py
#set tmpfile = $WORKDIR/convert.tmp.file
#foreach pathinfo (`ls *pathinfo`)
	#set xx = `cat $pathinfo |awk 'NR==1 {print $1}'`
	#set yy = `cat $pathinfo |awk 'NR==2 {print $1}'`
	#set zz = `cat $pathinfo |awk 'NR==3 {print $1}'`
	#cat $pathinfo |awk 'NR>3 {print $0}' > $tmpfile
#python $CONVERT $tmpfile  $pathinfo $xx $yy $zz
#end
### ===================================================


### ===================================================
### convert pathinfo  into vtk format
#echo "--> Convert pathinfo cross-section into vtk"
#cd $WORKDIR
#set CONVERT = $PWD/convert.py
#set tmpfile = $WORKDIR/convert.tmp.file
#foreach pathinfo (`ls *pathline`)
	#set xx = `cat $pathinfo |awk 'NR==1 {print $1}'`
	#set yy = `cat $pathinfo |awk 'NR==2 {print $1}'`
	#set zz = `cat $pathinfo |awk 'NR==3 {print $1}'`
	#cat $pathinfo |awk 'NR>3 {print $0}' > $tmpfile
#python $CONVERT $tmpfile  $pathinfo $xx $yy $zz
#end
### ===================================================







# plot hit count map specific depth 
##csh $SRCDIR/c06.plot_hit_count_for_multi_depth.sh $SRCDIR $PLOTDIR  $PWD 

# =======================================================
# Plotting 3models together for different depth   starting + delta + final
##csh $SRCDIR/c05.plot_horizontal_slice_for_starting_residual_final_tomo_multi_depth.sh $SRCDIR $PLOTDIR $PWD $WORKDIR

# =======================================================


# =======================================================
# 4. plot horizontal slice map for review
	# for whole mantle loop
##echo "---> Plotting For Whole Mantle Run"
# =======================================================
##set flag = whole
##csh $SRCDIR/c02.plot_horizontal_slice_for_residual_tomo_for_multy_depth.sh $SRCDIR $PLOTDIR $flag $WORKDIR
##csh $SRCDIR/c03.plot_timeinfo_of_LSM_single_depth.sh $SRCDIR $PLOTDIR $WORKDIR   $flag							 
	
# =======================================================
# 5. plot horizontal slice map for review
	# for lower mantle loop
##echo "---> Plotting For Whole Mantle Run"
# =======================================================
##set flag = lower
##csh $SRCDIR/c02.plot_horizontal_slice_for_residual_tomo_for_multy_depth.sh $SRCDIR $PLOTDIR $flag $WORKDIR
##csh $SRCDIR/c03.plot_timeinfo_of_LSM_single_depth.sh $SRCDIR $PLOTDIR $WORKDIR   $flag							 



echo "--> Working on scripts"
##csh $script  > & /dev/null
csh $script  
#csh $script  

##to_hongyu $PLOTDIR
##
exit 0
