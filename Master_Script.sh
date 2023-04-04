#!/bin/bash
#### Pre install both conda environments using conda -f

### How many cores for snakemake to use
CORES=64



###### load first conda env, generate trees through Slim-Tree

eval "$(conda shell.bash hook)"
conda activate SlimTree

snakemake -F -j$CORES

conda deactivate




###### load second conda env, infere trees through Paleoprophyler


eval "$(conda shell.bash hook)"
conda activate Analyser

cp Newick_Files/For_PaleoProPhyler/Datasets.txt Dataset_Analysis/Datasets.txt
cp Newick_Files/For_PaleoProPhyler/*_Dataset.fasta Dataset_Analysis/Workspace/1_OG_Dataset/

rm -rf Dataset_Analysis/Workspace/2_DATASETS/*

cd Dataset_Analysis
snakemake -F -j$CORES
cd ..
conda deactivate






###### load third conda env Comapare input tree with output tree


eval "$(conda shell.bash hook)"
conda activate Entropy
cd Tree_Comparisons
snakemake -F -j$CORES
cd ..
conda deactivate