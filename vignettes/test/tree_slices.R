library(ape)
library(picante)
library(plotrix)
library(maptpx)
library(CountClust)
library(cluster)
library(FactoMineR)

treeSlice<-function(tree,slice){
  if(class(tree)!="phylo") stop("tree should be object of class 'phylo.'")
  tree<-reorder(tree) # reorder cladewise
  # compute node heights
  root<-length(tree$tip)+1
  node.height<-matrix(NA,nrow(tree$edge),2)
  for(i in 1:nrow(tree$edge)){
    if(tree$edge[i,1]==root){
      node.height[i,1]<-0.0
      node.height[i,2]<-tree$edge.length[i]
    } else {
      node.height[i,1]<-node.height[match(tree$edge[i,1],tree$edge[,2]),2]
      node.height[i,2]<-node.height[i,1]+tree$edge.length[i]
    }
  }
  edges<-which(node.height[,2]>slice&node.height[,1]<slice)
  nodes<-tree$edge[edges,2]; nodes<-nodes[nodes>length(tree$tip)]
  trees<-list(); class(trees)<-"multiPhylo"
  for(i in 1:length(nodes)){ 
    trees[[i]]<-extract.clade(tree,node=nodes[i])
    trees[[i]]$root.edge<-node.height[which(tree$edge[,2]==nodes[i]),2]-slice
  }
  return(trees)
}

new_counts <- as.matrix(read.csv("Counts_9_15_16.csv",header=TRUE,row.names=1))
grids <- new_counts

all.him.tree <- read.tree(file="AllHim_Mar_27_2015.tre")
grids.tree <- drop.tip(all.him.tree, drop.tip(all.him.tree, colnames(new_counts))$tip.label)
write.tree(grids.tree, file="grids_tree_9_21_16.tre")

large.tree <- read.tree(file="grids_tree_9_21_16.tre")

70.85267-5
65.85267

trees5<-treeSlice(large.tree,65.85267)

grids5 <- grids
for( i in 1:length(trees5)){
  print(trees5[[i]]$tip.label)
  new.column <- as.data.frame(rowSums(grids5[,trees5[[i]]$tip.label]))
  colnames(new.column) <- trees5[[i]]$tip.label[1]
  drops <- trees5[[i]]$tip.label
  grids5 <- grids5[,!(names(grids5) %in% drops)]
  grids5 <- cbind(grids5, new.column)
}


trees10 <-treeSlice(large.tree,60.85267)

grids10 <- as.data.frame(grids)
for( i in 1:length(trees10)){
  print(i)
  new.column <- as.data.frame(rowSums(grids10[,trees10[[i]]$tip.label]))
  colnames(new.column) <- trees10[[i]]$tip.label[1]
  drops <- trees10[[i]]$tip.label
  grids10 <- grids10[,!(names(grids10) %in% drops)]
  grids10 <- cbind(grids10, new.column)
}

trees15 <-treeSlice(large.tree,55.85267)

grids15 <-as.data.frame(grids)
for( i in 1:length(trees15)){
  new.column <- as.data.frame(rowSums(grids15[,trees15[[i]]$tip.label]))
  colnames(new.column) <- trees15[[i]]$tip.label[1]
  drops <- trees15[[i]]$tip.label
  grids15 <- grids15[,!(names(grids15) %in% drops)]
  grids15 <- cbind(grids15, new.column)
}

trees20 <-treeSlice(large.tree,50.85267)

grids20 <- as.data.frame(grids)
for( i in 1:length(trees20)){
  new.column <- as.data.frame(rowSums(grids20[,trees20[[i]]$tip.label]))
  colnames(new.column) <- trees20[[i]]$tip.label[1]
  drops <- trees20[[i]]$tip.label
  grids20 <- grids20[,!(names(grids20) %in% drops)]
  grids20 <- cbind(grids20, new.column)
}

trees25 <-treeSlice(large.tree,45.85267)

grids25 <- as.data.frame(grids)
for( i in 1:length(trees25)){
  new.column <- as.data.frame(rowSums(grids25[,trees25[[i]]$tip.label]))
  colnames(new.column) <- trees25[[i]]$tip.label[1]
  drops <- trees25[[i]]$tip.label
  grids25 <- grids25[,!(names(grids25) %in% drops)]
  grids25 <- cbind(grids25, new.column)
}

trees30 <-treeSlice(large.tree,40.85267)

grids30 <- as.data.frame(grids)
for( i in 1:length(trees30)){
  new.column <- as.data.frame(rowSums(grids30[,trees30[[i]]$tip.label]))
  colnames(new.column) <- trees30[[i]]$tip.label[1]
  drops <- trees30[[i]]$tip.label
  grids30 <- grids30[,!(names(grids30) %in% drops)]
  grids30 <- cbind(grids30, new.column)
}

trees35 <-treeSlice(large.tree,35.85267)

grids35 <- as.data.frame(grids)
for( i in 1:length(trees35)){
  new.column <- as.data.frame(rowSums(grids35[,trees35[[i]]$tip.label]))
  colnames(new.column) <- trees35[[i]]$tip.label[1]
  drops <- trees35[[i]]$tip.label
  grids35 <- grids35[,!(names(grids35) %in% drops)]
  grids35 <- cbind(grids35, new.column)
}

trees40 <-treeSlice(large.tree,30.85267)

grids40 <- as.data.frame(grids)
for( i in 1:length(trees40)){
  new.column <- as.data.frame(rowSums(grids40[,trees40[[i]]$tip.label]))
  colnames(new.column) <- trees40[[i]]$tip.label[1]
  drops <- trees40[[i]]$tip.label
  grids40 <- grids40[,!(names(grids40) %in% drops)]
  grids40 <- cbind(grids40, new.column)
}

trees45 <-treeSlice(large.tree,25.85267)

grids45 <- as.data.frame(grids)
for( i in 1:length(trees45)){
  new.column <- as.data.frame(rowSums(grids45[,trees45[[i]]$tip.label]))
  colnames(new.column) <- trees45[[i]]$tip.label[1]
  drops <- trees45[[i]]$tip.label
  grids45 <- grids45[,!(names(grids45) %in% drops)]
  grids45 <- cbind(grids45, new.column)
}

write.csv(grids5, file="grids5.csv", fileEncoding = 'macroman')
write.csv(grids10, file="grids10_9_21_16.csv", fileEncoding = 'macroman')
write.csv(grids15, file="grids15_9_21_16.csv", fileEncoding = 'macroman')
write.csv(grids20, file="grids20_9_21_16.csv", fileEncoding = 'macroman')
write.csv(grids25, file="grids25.csv", fileEncoding = 'macroman')
write.csv(grids30, file="grids30.csv", fileEncoding = 'macroman')
write.csv(grids35, file="grids35.csv", fileEncoding = 'macroman')
write.csv(grids40, file="grids40.csv", fileEncoding = 'macroman')
write.csv(grids45, file="grids45.csv", fileEncoding = 'macroman')

pruned.tree10<-drop.tip(large.tree,large.tree$tip.label[-match(colnames(grids10), large.tree$tip.label)])
pruned.tree15<-drop.tip(large.tree,large.tree$tip.label[-match(colnames(grids15), large.tree$tip.label)])
pruned.tree20<-drop.tip(large.tree,large.tree$tip.label[-match(colnames(grids20), large.tree$tip.label)])

plot.phylo(pruned.tree10, x.lim = c(0,60.85267), show.tip.label = F)
plot.phylo(pruned.tree20, x.lim = c(0,50.85267), show.tip.label = F)

plot.phylo(large.tree, show.tip.label = F)
lines(x=c(60.85267,60.85267),y=c(0,305),lwd=3)
lines(x=c(50.85267,50.85267),y=c(0,305),lwd=3)



