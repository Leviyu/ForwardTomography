#!/bin/tcsh


set PWD = `pwd`
set MODEL_FILE = /DATA1/ForwardTomography/WORKDIR/T2/output_tomo2.S40RTS



set MODEL = `echo $MODEL_FILE |awk -F"." '{print $2}'`

set outfile = $PWD/MODEL_DIR/${MODEL}.data



set convert_py = $PWD/convert.py


set x = `cat $MODEL_FILE |awk 'NR==1 {print $1}'`
set y = `cat $MODEL_FILE |awk 'NR==2 {print $1}'`
set z = `cat $MODEL_FILE |awk 'NR==3 {print $1}'`
set tmp_model_file = $PWD/MODEL_DIR/tmp.data
cat $MODEL_FILE|awk 'NR>3 {print $0}' > $tmp_model_file

echo $tmp_model_file $outfile
echo $x $y $z
py $convert_py $tmp_model_file $outfile  $x $y $z


echo "--> Convert complete"


