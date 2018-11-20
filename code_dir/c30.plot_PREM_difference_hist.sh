#!/bin/tcsh




set PWD = `pwd`
set out = $PWD/PLOTDIR/hist.pdf


set file = ~/ForwardTomography/LOG/logfile.M20


set tmp = $PWD/PLOTDIR/tmp.data

#//cat $file |grep EQ_STA_PHASE
#cat $file |grep EQ_STA_PHASE|awk '{print $11}'
cat $file |grep EQ_STA_PHASE|awk '{print $11}' >! $tmp
#cat $file |grep EQ_STA_PHASE|awk '$4=="ScS" {print $11","$13}' >! $tmp
set script = $PWD/PLOTDIR/.script

cat << EOF >! $script
import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt

infile = '$tmp'
df = pd.read_csv(infile)
bb = np.arange(-2,2,0.1)
df.hist(bins=bb)
plt.plot()
plt.savefig('$out')
EOF
py $script
to_hongyu $out
