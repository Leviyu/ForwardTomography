#!/bin/csh

set PLOTDIR = $1
set WORKDIR = $2
set ID = $3

set OUT = $WORKDIR/COST.ps
set OUT_pdf = $PLOTDIR/COST.pdf

pstext -Y10.5i -JX6i/1i -R0/1/0/1 -K -N -P << EOF >! $OUT
0 0 15 0 0 LB Cost and Variance Reduction
EOF


foreach step (1 2 3 4 5 6 7 8 9)
set file2 = $WORKDIR/timeinfo.step.${step}.iteration.1.info
if(! -e $file2 ) then
continue
endif
	echo "--> worling on step $step"

set cost_file = $WORKDIR/.cost
set variance_file = $WORKDIR/.variance
cat /dev/null >! $variance_file

set file3 = $WORKDIR/timeinfo.step.${step}.iteration.1.info
set variance_tmp = `cat $file3 |awk 'NR==1 {print $10}'`
echo "0 $variance_tmp " >> $variance_file

	foreach iteration (1 2 3 4 5)
	set file = $WORKDIR/timeinfo.step.${step}.iteration.${iteration}.info
	set meta = $WORKDIR/timeinfo.step.${step}.iteration.${iteration}.meta

set current_num = `cat $meta |awk '{print $6}'`
set total_num = `cat $meta |awk '{print $8}'`
set orig_variance = `cat $meta |awk '{print $10}'`
set variance = `cat $meta |awk '{print $12}'`
set cost = `cat $meta |awk '{print $14}'`

##set dt_list = $WORKDIR/tmp.dt.${step}.${iteration}
##cp $file $dt_list
##set XMIN = -25    # min histog x-axis (here, degrees)
##set XMAX = 25  # max histog x-axis (here, degrees)
##set XINC = 1    # histog bin width
##set XNUM = 5   # increment for x-axis number

##pshistogram  $dt_list  -W$XINC -IO >! $WORKDIR/histo.tmp
##set sta_num = `cat $dt_list |wc -l`

##set YMAX = `minmax -C $WORKDIR/histo.tmp |awk '{print $4*1.2}'`
##set YNUM = ` echo $YMAX | awk '{print 20*int(1.0*$1/100.0) }' `
# have twice as many ticks as numbers on y-axis:
##set YTICK = ` echo $YNUM | awk '{print int(1.0*$1 / 2.0) }' `

##set XLABEL = ' '
##set YLABEL = ''
##echo "YMAX is $YMAX"
# change the GMT label font size default:
##gmtset LABEL_FONT_SIZE = 13p

# add some text info
##pstext -R0/1/0/1 -JX6i/2i -N  -Y-1.5i  -K -O -P << EOF >> $OUT
##0 0.5 13 0 0 LB Step: $step Iter: $iteration NUM: $current_num/$total_num Variance: $variance Cost: $cost
##EOF

##pshistogram $dt_list -R${XMIN}/${XMAX}/0/$YMAX -Ba${XNUM}f${XINC}:"${XLABEL}":/a${YNUM}f${YTICK}:"${YLABEL}":WSne -JX6.0i/0.6i -W$XINC -L0.5p -G50/50/250 -V  -Y-0.3i  -P -K -O  >> $OUT

echo "$iteration $variance " >> $variance_file
	end # iter

	
set cost_max = `cat $cost_file |head -n 1 |awk '{print $2}'`
set cost_min = `cat $cost_file |tail -n 1 |awk '{print $2}'`
echo " cost min max $cost_min $cost_max"
set inc = `echo "($cost_max - $cost_min) / 3"|bc -l`
set inc = `printf "%.0f" $inc`

psxy $cost_file -JX3i/0.5i -R0/5/$cost_min/$cost_max -Ba1f1:"Step ${step} Iteration":/a${inc}f{$inc}:"Cost":SW -O -K -Gred -Y-1.3i -Sc0.3 -N<< EOF >> $OUT
EOF

end # layer



ps2pdf $OUT $OUT_pdf
#to_hongyu $OUT_pdf
