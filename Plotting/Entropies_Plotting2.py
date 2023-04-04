import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random


########################################################################
######### LOAD DATA


Simulated_Subfolders = [ f.path for f in os.scandir('./Simulated') if f.is_dir() ]


Simulated_Entropies=[]


for FOLDER in Simulated_Subfolders:
    
    Entropies=open(FOLDER+'/Entropies.txt')
    for line in Entropies:
        line=float(line.strip())
        Simulated_Entropies.append(line)



Real_Entropies={}

RL_ENTR=open('Entropies_Hominids/Combined_Entropies','r')

for line in RL_ENTR:
    line=line.strip().split(' : ')
    PROTEIN=line[0]
    PROT_ENTROPY=float(line[1])
    Real_Entropies[PROTEIN]=PROT_ENTROPY
    
print(Real_Entropies)




REAL_DATA_X=[]
REAL_DATA_Y=[]
REAL_DATA_LABELS=[]
REAL_DATA_TEXT=[]


for x,j in Real_Entropies.items():

    REAL_DATA_Y.append(j)
    # REAL_DATA_X.append(random.uniform(-0.1, 0.1))
    REAL_DATA_X.append(0)
    REAL_DATA_TEXT.append(x)
    






########################################################################

### PLOT DATA

fig = plt.figure()
ax = fig.add_subplot()

fig.suptitle('Real Proteins - Entropies', fontsize=14, fontweight='bold')


ax.set_xlabel('Proteins', fontsize=11, fontweight='bold')
ax.set_ylabel('Calculated Entropy', fontsize=11, fontweight='bold')

ax.scatter(REAL_DATA_X, REAL_DATA_Y)

for Z in range(0,len(REAL_DATA_Y)):
    ax.annotate(REAL_DATA_TEXT[Z], (REAL_DATA_X[Z], REAL_DATA_Y[Z]))



plt.show()
