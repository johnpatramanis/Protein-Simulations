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



COMBINED=[]
COMBINED_LABELS=[]

REAL_DATA=[]
REAL_DATA_LABELS=[]
REAL_DATA_TEXT=[]


for x,j in Real_Entropies.items():
    COMBINED.append(j)
    COMBINED_LABELS.append('Real Data')
    
    REAL_DATA.append(j)
    REAL_DATA_LABELS.append('Real Data')
    REAL_DATA_TEXT.append(F'${x}$')
    

for x in Simulated_Entropies:
    COMBINED.append(x)
    COMBINED_LABELS.append('Simulated Data')








########################################################################

### PLOT DATA

fig = plt.figure()
ax = fig.add_subplot()

fig.suptitle('Generated proteins vs Real proteins', fontsize=14, fontweight='bold')


ax.set_xlabel('Origin of Data', fontsize=11, fontweight='bold')
ax.set_ylabel('Calculated Entropy', fontsize=11, fontweight='bold')

sns.stripplot(x=COMBINED_LABELS, y=COMBINED)

# for X in range(0,len(REAL_DATA)):
    # ax.text(x=0.1+random.uniform(0,0.1),y=REAL_DATA[X],s=REAL_DATA_TEXT[X],transform=ax.transData,size=5)




plt.show()
