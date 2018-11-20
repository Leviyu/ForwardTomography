#!/usr/bin/python2.7

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt




colum = ["dt_anomaly","model prediction"]


df = pd.read_csv("./WORKDIR/out.c03",sep=' ',names=colum)


print (df["dt_anomaly"] - df["model prediction"]).mean()
print (df["dt_anomaly"] - df["model prediction"]).std()

df.plot()
plt.gca().invert_yaxis()
plt.xlabel("taup_path increment")
plt.ylabel("average of dt for 1000records () relative to inc=0.1 ")
plt.show()
