library(ape)
library(kmerDistance)
args <- commandArgs(trailingOnly = TRUE)
klen <- as.integer(args[1])
input_dir <- args[2]
output_dir <- args[3]
kc_dist <- as.dist(kmerDistance.dif(klen,input_dir))
kc_tree <- nj(kc_dist)
#cophenetic(kc_tree)
kc_td <- dist.nodes(kc_tree)
kc_branch <- node.depth(kc_tree)
write.tree(kc_tree,paste(output_dir,"tree.csv",sep=""))
#print(as.matrix(cophenetic(kc_tree)))
write.csv(cophenetic(kc_tree),paste(output_dir,"dist.csv",sep=""))
