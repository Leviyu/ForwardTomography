#!/usr/bin/python2.7

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt




incs = ["0.005" ,"0.01", "0.02", "0.05", "0.1", "0.2", "0.3"]
#incs = ["0.005" ,"0.01", "0.02"]


df = pd.DataFrame( columns=incs)

for inc in incs:
    inputfile = "./WORKDIR/tmp."+str(inc)
    print inputfile

    # read into file
    dd = pd.read_csv(inputfile,names=[inc])
    #print dd
    df[inc] = dd
    #df[inc].read_csv(inputfile,sep='')

    #print df[inc]

#print df
dd = df["0.005"].copy()

index = 0
for inc in incs:
    df[inc] = df[inc] - dd
    #df[inc] = df[inc].astype(float)
    #//df.a.astype(float)
    #mask = (df[inc]<1 & df[inc]>-1)
    #df[inc] = df[ mask ]][inc]
    #plt.figure(index)
    #df[inc].hist(bins=100)
    #plt.title(inc)
    index +=1
df = df[ ((df["0.02"] > -1.0) & (df["0.02"]< 1.0) )]

#//print df.describe()
#plt.show()
#exit(0)

#print df

df.mean().plot()
plt.gca().invert_yaxis()
plt.xlabel("taup_path increment")
plt.ylabel("average of dt for 1000records () relative to inc=0.005 ")
plt.show()
