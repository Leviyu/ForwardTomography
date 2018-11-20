#!/bin/csh

# This script plot the cross-section of EQ STA geometry


set EQ = $1
set STA = $2
set PHASE = $3
set eventinfo = $4
set PLOTDIR = $5
set SRCDIR = $6
set WORKDIR = $7


set loop = 1
set flag = whole
set GCP = $WORKDIR/GCP.${EQ}.${STA}.${PHASE}.gcp
#set theta_r = vertical_theta_r_dvs.data.${EQ}.${STA}.${PHASE}
set theta_r = cross_section.${EQ}.${STA}.${PHASE}
set cross_points = cross_path.${EQ}.${STA}.${PHASE}
set hit_count = hit_count.${EQ}.${STA}.${PHASE}
#set cross_points = plot.output.cross_points.${EQ}.${STA}.${PHASE}
#set TIME_info = loop_${loop}_${flag}/out.record_info
set path_lat_lon = plot.output.cross_points_horizontal.${EQ}.${STA}.${PHASE}				# horizontal cross-points

set OUT = $PLOTDIR/${EQ}_${STA}_${PHASE}_vertical_slice.ps

echo " ============ plot cross-section plot for Station $STA ==================" 
set INPUT = vertical.tmp
awk 'NR>1 {print $0}' $theta_r > $INPUT

set tmp = `minmax -C $INPUT`
set theta_min = $tmp[1]
set theta_max = $tmp[2]
set r_min = $tmp[3]
set r_max = $tmp[4]
set MAX_ang = `echo "$theta_max - $theta_min"|bc -l`
set ang_offset = `echo "-1*(90 - $MAX_ang/2)"|bc -l`

set GRIDFILE = input.grd
set CPT = vertical.cpt
set INC_THETA = 2
set INC_DEP = 100
set MAP = ${theta_min}/$theta_max/3481/6371

#surface $INPUT -G$GRIDFILE -T0.25  -I$INC_THETA/$INC_DEP -R$theta_min/$theta_max/$r_min/$r_max  
##makecpt -Cpolar -I -T-2.5/2.5/0.5 -Z > $CPT
##gmtset COLOR_FOREGROUND = 255/255/255
gmtset COLOR_FOREGROUND = 0/29/223
gmtset COLOR_BACKGROUND = 200/0/0
##gmtset COLOR_BACKGROUND = 255/255/255
gmtset COLOR_NAN = 255/255/255


###################################################
# Plot for dvs cross section
xyz2grd $INPUT -G$GRIDFILE  -I$INC_THETA/$INC_DEP -R$theta_min/$theta_max/$r_min/$r_max  
makecpt -Cseis   -T-2/2/0.5 -Z > $CPT
set CPT = vertical.cpt

pstext -JX9i/1i -R0/1/0/1 -Y7i -X2i -K  << EOF >$OUT
EOF

set scale = `echo "($MAX_ang -20)/20 * 1.5 + 1.5"|bc -l`
set scale = `printf "%.0f" $scale`
set JP = -JP${scale}i/${ang_offset}z
grdimage $GRIDFILE -C$CPT $JP -Ba20f30/a1000f1000SWNE  -R$MAP -E300  -N -Y-3.5i -O -K  <<EOF>>$OUT
EOF
set vertical_GCP = ${STA}.${EQ}.${PHASE}.path
# add vertical great cicle path 
psxy $vertical_GCP  -JP -R   -O -N -K -P >>$OUT
# add cross-points 
psxy $cross_points -JP -Sx0.04 -N -R -O -K -P >> $OUT
# add EQ and STA mark
pstext -JP -R$MAP -O -K -P -N << EOF >>$OUT
$theta_min 6800 8 0 0 LB EQ
$theta_max 6800 8 0 0 LB Station
EOF
psxy -JP -R -Sc0.2 -O -K -N -G255/128/0 -P << EOF>>$OUT
$theta_min 6400
$theta_max 6400
EOF

#grdcontour $GRIDFILE -C$CPT -A- -JP -R$MAP -O -K -P >>$OUT
psscale -C$CPT -D-1i/1i/1.5i/0.15i -B0.5/:"dvs": -Y-0.1i -O -K -N -P <<EOF>>$OUT
EOF

gmtset COLOR_FOREGROUND = 255/255/255
gmtset COLOR_BACKGROUND = 0/0/0
##gmtset COLOR_BACKGROUND = 255/255/255
gmtset COLOR_NAN = 255/255/255
########################################################
# Plot for hit count cross-section

xyz2grd $hit_count  -G$GRIDFILE  -I$INC_THETA/$INC_DEP -R$theta_min/$theta_max/$r_min/$r_max  
makecpt -Chot  -T-0.1/10/3  -Z -I> $CPT
set CPT = vertical.cpt


grdimage $GRIDFILE -C$CPT $JP -Ba20f30/a1000f1000SWNE  -R$MAP -E300  -N -Y-2.5i -O -K  <<EOF>>$OUT
EOF
set vertical_GCP = ${STA}.${EQ}.${PHASE}.path
# add vertical great cicle path 
psxy $vertical_GCP  -JP -R   -O -N -K -P >>$OUT
# add cross-points 
psxy $cross_points -JP -Sx0.04 -N -R -O -K -P >> $OUT
# add EQ and STA mark
pstext -JP -R$MAP -O -K -P -N << EOF >>$OUT
$theta_min 6800 8 0 0 LB EQ
$theta_max 6800 8 0 0 LB Station
EOF
psxy -JP -R -Sc0.2 -O -K -N -G255/128/0 -P << EOF>>$OUT
$theta_min 6400
$theta_max 6400
EOF

#grdcontour $GRIDFILE -C$CPT -A- -JP -R$MAP -O -K -P >>$OUT
psscale -C$CPT -D-1i/1i/1.5i/0.15i -B2/:"hit count": -Y-0.1i -O -K -N -P <<EOF>>$OUT
EOF

# add the small global map showing the gcp of this cross-section
set LAT_LON = `awk 'NR==1 {print $0}' $theta_r`
set EQ_lat = $LAT_LON[7]
set EQ_lon = $LAT_LON[8]
set STA_lat = $LAT_LON[9]
set STA_lon = $LAT_LON[10]
set MID_lat = $LAT_LON[5]
set MID_lon = $LAT_LON[6]


set MAP = -R0/360/-90/90
set PROJ = -JG${MID_lon}/${MID_lat}/2i
set land = "255/225/160"
set sea = 103/204/0

pscoast $MAP $PROJ -Dc -A40000 -B45g45 -W1 -G$land -S$sea -Y5.3i  -X0i   -O -K -P >>$OUT

psxy $GCP $MAP $PROJ -:   -W0.001/"red" -O -P -K  >>$OUT
# add EQ and STA location
echo $EQ_lat $EQ_lon
psxy $MAP $PROJ -: -Sa0.4	-G255/0/0 -O -P -K << EOF>>$OUT
$EQ_lat $EQ_lon
EOF
psxy $MAP $PROJ -: -Si0.4  -G0/0/0	-O -K -P  << EOF>>$OUT
$STA_lat $STA_lon
EOF


# add horizontal cross-point 
##set horizontal_cross_tmp = horizontal_cross.tmp
##awk '{print $1,$2}' $path_lat_lon > $horizontal_cross_tmp
##psxy $horizontal_cross_tmp $MAP $PROJ -: -S+0.04  -O -K -P >> $OUT


# ======= add horizontal slice cross-points ======
##set MAP = -R${EQ_lon}/${STA_lon}/${EQ_lat}/${STA_lat}
##//set PROJ = -JG${MID_lon}/${MID_lat}/4i
##set land = "255/225/160"
##set sea = 103/204/0
##pscoast $MAP $PROJ -Dc -A4000 -B40g40 -W1 -Y-1i -O -K -P >>$OUT




# add text info for TC
#set tmp = `grep -w $STA $TIME_info`
#set PREM = $tmp[12]
#set PREM_TOMO = $tmp[14]
#set tomo_correction = $tmp[10]
#set residual = $tmp[16]
set obs_prem = `cat $eventinfo |grep -w $EQ |grep -w $STA |grep -w $PHASE |awk '{print $19}'`
set MODEL = `grep -w MODEL_NAME $WORKDIR/INFILE |awk 'NR==1 {print $2}'`


set LSM_INPUT = $WORKDIR/LSM_record_input
set TMP = `cat $LSM_INPUT |grep -w $EQ |grep -w $STA |grep -w $PHASE |awk 'NR==1 {print $0}'`
set EQ_lat = $TMP[8]
set EQ_lon = $TMP[9]
set EQ_dep = $TMP[10]
set STA_lat = $TMP[6]
set STA_lon = $TMP[7]
set DIST = $TMP[3]
#
#set PREM_taup_time_sta_EQ =  `taup_time -mod prem -ph $PHASE -h $EQ_dep -sta $STA_lat $STA_lon -evt $EQ_lat $EQ_lon |awk 'NR==6 {print $4}'`
#set time_deg = `taup_time -mod prem -ph $PHASE -h $EQ_dep -deg $DIST| awk 'NR==6 {print $4}'`
##echo "taup_time -mod prem -ph $PHASE -h $EQ_dep -sta $STA_lat $STA_lon -evt $EQ_lat $EQ_lon"
##echo "$EQ $STA $PHASE $EQ_dep $DIST"
#set delta_taup = `echo "$PREM_taup_time_sta_EQ - $time_deg"|bc -l`

#set delta_PREM2 = `echo "$PREM_taup_time_sta_EQ - $PREM"|bc -l`

#0.5 1.2 10 0 0 LB loop: $loop  flag: $flag
#0.5 0.2 10 0 0 LB Tomography Correction(Tomo): $tomo_correction
#0.5 0.0 10 0 0 LB Residual(Relative to Tomo): $residual
#0.5 0.6 10 0 0 LB PREM(Tomo): $PREM PREM(taup sta_EQ): $PREM_taup_time_sta_EQ dt: $delta_PREM2
pstext -JX7i/1i -R0/1/0/1 -O -K -P -Y-3.5i -X3i -N << EOF >>$OUT
0.5 1.0 10 0 0 LB Cross-Section Tomography of Station: $STA
0.5 0.8 10 0 0 LB Tomography Model: $MODEL 
0.5 0.4 10 0 0 LB Observation(Relative to PREM taup): $obs_prem
EOF

pstext -JX -R -O << EOF >>$OUT
EOF


exit 0

