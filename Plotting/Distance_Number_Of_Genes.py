import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random
import glob


########################################################################
######### LOAD DATA


Simulated_Subfolders = [ f.path for f in os.scandir('./Simulated') if f.is_dir() ]


Simulated_Entropies=[]

NUMBER_TO_METRICS=[]

for FOLDER in Simulated_Subfolders:
    
    Distances_FILE=open(FOLDER+'/Tree_Distances.txt')
    Number_Of_Genes = len(glob.glob1(FOLDER+'/',"*.entr"))
    
    
    for line in Distances_FILE:
        line=line.strip().split()
        line=[float(x) for x in line]
        NUMBER_TO_METRICS.append([int(Number_Of_Genes),line])



print(NUMBER_TO_METRICS)




########################################################################

### PLOT DATA

fig = plt.figure()
ax = fig.add_subplot()

fig.suptitle('Number of Genes and Distance to original tree', fontsize=14, fontweight='bold')


ax.set_xlabel('Number of Genes/Proteins', fontsize=11, fontweight='bold')
ax.set_ylabel('Distance to Original Tree', fontsize=11, fontweight='bold')

ax.scatter([x[0] for x in NUMBER_TO_METRICS], [x[1][0] for x in NUMBER_TO_METRICS])

# for Z in range(0,len(REAL_DATA_Y)):
    # ax.annotate(REAL_DATA_TEXT[Z], (REAL_DATA_X[Z], REAL_DATA_Y[Z]))



plt.show()
