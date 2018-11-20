#!/bin/tcsh

# This script plot the cross-section of EQ STA geometry for all record

set SRCDIR = $1
set event = $2
set PLOTDIR = $3
set WORKDIR = $4

set OUT = $PLOTDIR/vertical_slice.ps
set OUT_pdf = $PLOTDIR/vertical_slice_with_observation.pdf
/bin/rm -r $OUT >& /dev/null

foreach NR (`cat $event | awk '{print NR}'`)
set TMP = `awk -v nr=$NR 'NR==nr {print $0}' $event`
set EQ = $TMP[12]
set STA = $TMP[1]
set PHASE = $TMP[20]
echo "--> Plotting cross-section for EQ $EQ STA $STA PHASE: $PHASE"

csh $SRCDIR/c20.plot_vertical_slice_EQ_STA_cross_section_for_one_record.sh  $EQ $STA $PHASE $event $PLOTDIR $SRCDIR $WORKDIR

set ps_file = $PLOTDIR/${EQ}_${STA}_${PHASE}_vertical_slice.ps
cat $ps_file >> $OUT
rm $ps_file > & /dev/null
end

ps2pdf $OUT $OUT_pdf
rm $OUT
#to_hongyu $OUT_pdf
