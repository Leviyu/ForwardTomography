#!/bin/csh


# ============================================================
#	This Script plot the residual horizontal slice tomography model
# 	and the histogram of it
#
# Hongyu DATE: June 26 2014
# Key words: tomography horizontal slice plot residual
# ============================================================
##set DEPTH = 2871
##set PLOTDIR = /home/hongyu_ubuntu/PLOTDIR/TC03/100_forward_tomography/try_all_Sdiff
##//set WORKDIR = /DATA1/100_forward_tomography/try_all_Sdiff


set PLOTDIR = $1
set WORKDIR = $2
set ID = $3
set STEP = $4


echo "--> Plotting for Step $STEP"


set DEPTH_RANGE = (500 1000  2000  2800)
#set DEPTH_RANGE = (300 500 1100 1500 1900 2200 2500 2800)

set INFILE = $WORKDIR/INFILE
set MODEL = `grep -w MODEL_NAME $INFILE |awk 'NR==1 {print $2}'`


set OUT = $PLOTDIR/STEP_${STEP}_04_different_depth_shell_${MODEL}.ps
set OUT_pdf = $PLOTDIR/STEP_${STEP}_04_different_depth_shell_${MODEL}.pdf




set REG1 = -R0/10/0/10
set PROJ1 = -JX3i/3i
set REG = -R-180/180/-90/90
set PROJ = -JR180/2i
gmtset COLOR_FOREGROUND = 255/0/0
gmtset COLOR_BACKGROUND = 0/0/255
gmtset COLOR_NAN = 255/255/255
##gmtset COLOR_FOREGROUND = 0/29/223
##gmtset COLOR_BACKGROUND = 200/0/0
#gmtset COLOR_NAN = 0/0/0
gmtset TICK_PEN = 0.2p
gmtset ANNOT_FONT_SIZE_PRIMARY = 5p
gmtset ANNOT_OFFSET_PRIMARY = 0.1c

pstext -JX3i/3i -R0/10/0/10 -K  -Y6i -N << EOF > $OUT
0 8 10 0 0 LB Model $MODEL
EOF


set grid_space = 1
set CPT_MAX = 2.5
set CPT_MIN = `echo "$CPT_MAX * -1"|bc -l`
set CPT_DELTA = `echo "$CPT_MAX /5"|bc -l`
set CPT = $WORKDIR/tomo.cpt
makecpt  -Cpolar -I -T$CPT_MIN/$CPT_MAX/$CPT_DELTA -Z >! $CPT
set CPT_residual = $WORKDIR/tomo_residual.cpt
makecpt  -Cpolar -I -T-1/1/0.2 -Z >! $CPT_residual
##set E_option = -E210/20
set E_option = -E300

psxy -JX -R -O -K -N  << EOF >> $OUT
EOF

set NUM = 1
foreach DEPTH (  `echo $DEPTH_RANGE` )
echo "================= plot residual horizontal slice tomography model for depth $DEPTH============"

	set MODEL_FILE = $WORKDIR/output_tomo.${MODEL}.final.${STEP}
	#set orig_MODEL_FILE = $WORKDIR/output_tomo.${MODEL}.start
	set current_modle_file = $WORKDIR/data.${DEPTH}.${MODEL}.slice.${DEPTH}
	set orig_current_modle_file = $WORKDIR/orig_data.${DEPTH}.${MODEL}.slice.${DEPTH}
	set delta_model = $WORKDIR/output_tomo.${MODEL}.delta.${STEP}
	set delta_model_file = $WORKDIR/output_tomo.${MODEL}.slice.${DEPTH}
	cat $MODEL_FILE | awk -v dd=${DEPTH} '{if($1 == dd ) print $3,$2,$4}' >! $current_modle_file
	cat $delta_model | awk -v dd=${DEPTH} '{if($1 == dd && $4!= 0 ) print $3,$2,$4}' >! $delta_model_file

#set model_file = $WORKDIR/data.${DEPTH}.${MODEL}.slice.${DEPTH}
#set orig_model_file = $WORKDIR/orig_data.${DEPTH}.${MODEL}.slice.${DEPTH}
	set grd_file = $WORKDIR/data.${DEPTH}.${MODEL}.${DEPTH}.grd
	set orig_grd_file = $WORKDIR/orig_data.${DEPTH}.${MODEL}.${DEPTH}.grd
	set delta_grd_file = $WORKDIR/delta.${DEPTH}.${MODEL}.${DEPTH}.grd

pstext $PROJ1 $REG1 -O -K -N  <<EOF >>$OUT
0 4.5 7 0 0 LB  DEPTH: $DEPTH
EOF
xyz2grd $current_modle_file -G$grd_file -I$grid_space $REG 
#xyz2grd $orig_current_modle_file -G$orig_grd_file -I$grid_space $REG 
xyz2grd $delta_model_file -G$delta_grd_file -I$grid_space $REG 
grdimage $delta_grd_file -C$CPT_residual $REG $PROJ $E_option   -O -K -Sc  >>$OUT
pscoast $REG $PROJ -Ba0g45/a0g45/wsne -Dh -A4000 -W2   -O  -K >>$OUT
if( $NUM == 1) then
psscale -C$CPT -D2i/-0.2i/1.5i/0.1ih -B0.5/:"tomo @~\144@~Vs (%)":  -O  -K -N300 >>$OUT
endif
if( $NUM == 2) then
psscale -C$CPT_residual -D2i/-0.2i/1.5i/0.1ih -B0.3/:"residual tomo @~\144@~Vs (%)":  -O  -K -N300 >>$OUT
endif

psxy $REG $PROJ -O -K -X2i -N << EOF >> $OUT
EOF
grdimage $grd_file -C$CPT $REG $PROJ $E_option   -O -K -Sc   >>$OUT
pscoast $REG $PROJ -Ba0g45/a0g45/wsne -Dh -A4000 -W2   -O  -K >>$OUT

psxy -JX -R -O -K -N -X2.5i << EOF >> $OUT
EOF

# if num ==3 change line
if($NUM == 2 || $NUM == 4 || $NUM == 6 ||$NUM == 8  ) then
psxy -JX -R -O -K -N -X-9i -Y-1.8i << EOF >> $OUT
EOF

endif


@ NUM ++
end 



# add text
pstext $PROJ1 $REG1 -O  -N  -Y-2.3i  <<EOF>>$OUT
EOF

ps2pdf $OUT $OUT_pdf
rm $OUT
to_hongyu $OUT_pdf
