args = commandArgs(trailingOnly=TRUE)
library('TreeDist')


### Run with Rguments: Rscript ../Newick_Files/EXAMP1/Newick_File.txt ../Dataset_Analysis/Workspace/2_DATASETS/EXAMP1_Dataset/CONCATINATED/CONCATINATED_aln_e.phy_phyml_tree.txt ./EXAMP1/Tree_Distances.txt

### Get paths for files

INITIAL_TREE_PATH=args[1]
INFERED_TREE_PATH=args[2]
OUTPUT_FILES=args[3]


###################### Import Trees ####################

Tree1=ape::read.tree(INITIAL_TREE_PATH) #### Input Tree

Tips=length(Tree1$tip.label)
counter=0
for (T in 1:Tips){
	Tree1$tip.label[T]=paste('sample',as.character(Tips-(counter+1)),sep='')
	counter=counter+1
	}




Tree2=ape::read.tree(INFERED_TREE_PATH) #### Infered Tree

# Tips=length(Tree2$tip.label)
# counter=0
# for (T in 1:Tips){
	# Tree2$tip.label[T]=paste('sample_',as.character(counter),sep='')
	# counter=counter+1
	# }
	
print('Trees loaded')	
##################### Calculate Distances #################


##### Normal Tree Distances	
distance_Generic <- TreeDistance(Tree1, Tree2)
distance_Nye <- NyeSimilarity(Tree1, Tree2, normalize = TRUE)
distance_JRF <- JaccardRobinsonFoulds(Tree1,Tree2)
distance_Bod_Giaro <- MatchingSplitDistance(Tree1,Tree2)
distance_KC <- KendallColijn(Tree1,Tree2)




##### For comparisons of tree with not - equal number of leaves
Clustering_Entropy1=ClusteringEntropy(Tree1)
Clustering_Entropy2=ClusteringEntropy(Tree2)
Mutual_Clustering_Info=MutualClusteringInfo(Tree1, Tree2)



#### Output Distances into file
METRICS=c(distance_Generic,distance_Nye,distance_JRF,distance_Bod_Giaro,distance_KC,Mutual_Clustering_Info)
write(paste(METRICS,collapse='\t'), file = OUTPUT_FILES)