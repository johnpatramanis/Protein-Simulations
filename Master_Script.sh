#!/bin/bash
#### Pre install both conda environments using conda -f

###### load first conda env, generate trees through Slim-Tree
eval "$(conda shell.bash hook)"
conda activate SlimTree

snakemake -F -j10

cp Newick_Files/For_PaleoProPhyler/Datasets.txt Dataset_Analysis/Workspace
cp Newick_Files/For_PaleoProPhyler/*_Dataset.fasta Dataset_Analysis/Workspace/1_OG_Dataset/

conda deactivate



###### load second conda env, infere trees through Paleoprophyler


eval "$(conda shell.bash hook)"
conda activate Analyser
cd Dataset_Analysis
snakemake -F -j24
conda deactivate



###### load thrid conda env Comapare input tree with output tree
