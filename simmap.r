library(phytools)
library(corHMM)
library(geiger)
library(viridis)
library(beepr)
library(tidyverse)

load('MCC_tree_sm.rdata') ## analysis with consensus tree
load('sm_trees_incertidumbre.rdata') ## objects of 1000 phy trees analysis

# Consensus tree
tree <- read.nexus('Gomphrena_MCCT_Santi.nex')

## photosynthesis trait

photo_data <- read.csv('Fotosíntesis.csv', header = T)
rownames(photo_data) <- photo_data$sp
head(photo_data)

chk <- name.check(tree,photo_data)
chk
## different transitions matrix
ARD_p <- corHMM(tree, photo_data, rate.cat = 1, model = 'ARD', node.states = 'marginal',
                 root.p = 'yang', n.cores = 9, get.tip.states = T)
ARD_p

SYM_p <- corHMM(tree, photo_data, rate.cat = 1, model = 'SYM', node.states = 'marginal',
                 root.p = 'yang', n.cores = 9, get.tip.states = T)
SYM_p

ER_p <- corHMM(tree, photo_data, rate.cat = 1, model = 'ER', node.states = 'marginal',
                root.p = 'yang', n.cores = 9, get.tip.states = T)
ER_p

## AICc weights to calculate the model-averaged marginal probabilities at the nodes
aicc_p<-setNames(
  c(ARD_p$AICc, SYM_p$AICc, ER_p$AICc),
  c("ARD","SYM", 'ER'))
aicc_p

AICC.W_p <- aic.w(aicc_p)
AICC.W_p

## Akaike weights to generate stochastic maps under all of the alternative models 
## under consideration, in proportion to the weight of evidence in support of that model
nsim<-1000 ### there will be 1000 simmap
## arbitrarily add (or subtract) the deficit (or surplus) to models chosen at random
Nsim<-round(nsim*AICC.W_p)

d<-if(sum(Nsim)>nsim) -1 else 1
nsim<-Nsim+d*sample(c(rep(1,abs(nsim-sum(Nsim))),
                      rep(0,length(Nsim)-abs(nsim-sum(Nsim)))))
nsim

sm_p_ARD <- makeSimmap(tree = tree, data = photo_data,
                        model =  ARD_p$solution,
                        rate.cat = 1, nSim = nsim['ARD'],
                        nCores = 9)

sm_p_SYM <- makeSimmap(tree = tree, data = photo_data,
                        model =  SYM_p$solution,
                        rate.cat = 1, nSim = nsim['SYM'],
                        nCores = 9)

sm_p_ER <- makeSimmap(tree = tree, data = photo_data,
                       model =  ER_p$solution,
                       rate.cat = 1, nSim = nsim['ER'],
                       nCores = 9)

smtrees <- c(sm_p_ARD, sm_p_SYM, sm_p_ER)

sum_smtrees <- describe.simmap(smtrees);beep('fanfare')
sum_smtrees$ace

colr_p<- setNames(viridis(2, end = 1, direction = -1) , c(0,1))
#pdf('mcct_simmap_photo.pdf')
#svg('mcct_simmap_photo.svg')
plot(sum_smtrees, fsize=0.4, ftype="i", colors=colr_p,
     ylim = c(-2, Ntip(tree)), cex = c(0.4,0.4))
add.simmap.legend(leg=c('c3', 'c4'), colors=colr_p,prompt=FALSE,x=0,
                  y=20,fsize=0.8)
dev.off()
## for supplementary material - probabilities at nodes, for each model
#pdf('nodelabelstree.pdf')
plot(tree,no.margin=TRUE,edge.width=2,cex=0.7)
nodelabels(text=1:tree$Nnode,node=1:tree$Nnode+Ntip(tree), cex = 0.6)
dev.off()

ard200trees <- describe.simmap(sm_p_ARD)
sym400trees <- describe.simmap(sm_p_SYM)
er400trees <- describe.simmap(sm_p_ER)

rn = as.character(c(1:tree$Nnode,rownames(ard200trees$ace)[42:83]))
rownames(ard200trees$ace) = rn
#write.csv(ard200trees$ace,'ard200.csv')
rownames(sym400trees$ace) = rn
#write.csv(sym400trees$ace,'sym400.csv')
rownames(er400trees$ace) = rn
#write.csv(er400trees$ace,'er400.csv')

## habit trait

h_data <- read.csv('Hábito_anual(0)_perenne(1).csv', header = T)
rownames(h_data) <- h_data$sp

chk <- name.check(tree,h_data)
chk

## transition matrix

ARD_h <- corHMM(tree, h_data, rate.cat = 1, model = 'ARD', node.states = 'marginal',
                root.p = 'yang', n.cores = 9, get.tip.states = T)

SYM_h <- corHMM(tree, h_data, rate.cat = 1, model = 'SYM', node.states = 'marginal',
                root.p = 'yang', n.cores = 9, get.tip.states = T)

ER_h <- corHMM(tree, h_data, rate.cat = 1, model = 'ER', node.states = 'marginal',
               root.p = 'yang', n.cores = 9, get.tip.states = T)



## AICc weights to calculate the model-averaged marginal probabilities at the nodes
aicc_h<-setNames(
  c(ARD_h$AICc, SYM_h$AICc, ER_h$AICc),
  c("ARD","SYM", 'ER'))
aicc_h

AICC.W_h <- aic.w(aicc_h)
AICC.W_h

## Akaike weights to generate stochastic maps under all of the alternative models 
## under consideration, in proportion to the weight of evidence in support of that model
nsim<-1000 ### there will be 1000 simmap
## arbitrarily add (or subtract) the deficit (or surplus) to models chosen at random
nsim<-round(nsim*AICC.W_h)



sm_h_ARD <- makeSimmap(tree = tree, data = h_data,
                       model =  ARD_h$solution,
                       rate.cat = 1, nSim = nsim['ARD'],
                       nCores = 9)

sm_h_SYM <- makeSimmap(tree = tree, data = h_data,
                       model =  SYM_h$solution,
                       rate.cat = 1, nSim = nsim['SYM'],
                       nCores = 9)

sm_h_ER <- makeSimmap(tree = tree, data = h_data,
                      model =  ER_h$solution,
                      rate.cat = 1, nSim = nsim['ER'],
                      nCores = 9)

smtrees_h <- c(sm_h_ARD, sm_h_SYM, sm_h_ER)

sum_smtrees_h <- describe.simmap(smtrees_h);beep('fanfare')

colr_h<- setNames(magma(2, begin = 0.3, end = 0.9, direction = -1) , c(0,1))
#pdf('mcct_simmap_h.pdf')
#svg('mcct_simmap_h.svg')
plot(sum_smtrees_h, fsize=0.4, ftype="i", colors=colr_h,
     ylim = c(-2, Ntip(tree)), cex = c(0.4,0.4))
add.simmap.legend(leg=c('annual','perennial'), colors=colr_h,prompt=FALSE,x=0,
                  y=20,fsize=0.8)
dev.off()

#save.image('MCC_tree_sm.rdata')
## for suplemmentary material - habitat trait probabilities per model
ard188trees_h <- describe.simmap(sm_h_ARD)
sym406trees_h <- describe.simmap(sm_h_SYM)
er406trees_h <- describe.simmap(sm_h_ER)

rownames(ard188trees_h$ace) = rn
#write.csv(ard188trees_h$ace,'ard188h.csv')
rownames(sym406trees_h$ace) = rn
#write.csv(sym406trees_h$ace,'sym406h.csv')
rownames(er406trees_h$ace) = rn
#write.csv(er406trees_h$ace,'er406h.csv')
#finish consensus tree



### 1000 sampled trees
trees <- read.nexus("Gomphrena_1000_trees_Santi.nex")

simmap_trees <- function(tre,x, sim){ 
  ARD <- corHMM(tre, x, rate.cat = 1, model = 'ARD', node.states = 'marginal',
                root.p = 'yang', n.cores = 9, get.tip.states = T)
  SYM <- corHMM(tre, x, rate.cat = 1, model = 'SYM', node.states = 'marginal',
                root.p = 'yang', n.cores = 9, get.tip.states = T)
  ER <- corHMM(tre, x, rate.cat = 1, model = 'ER', node.states = 'marginal',
               root.p = 'yang', n.cores = 9, get.tip.states = T)
  aicc<-setNames(
    c(ARD$AICc, SYM$AICc, ER$AICc),
    c("ARD","SYM", 'ER'))
  AICC.W <- aic.w(aicc)
  nsim <- sim
  Nsim<-round(nsim*AICC.W)
  if (sum(Nsim) != nsim){ 
    d<-if(sum(Nsim)<nsim) -1 else 1
    nsim<-Nsim+d*sample(c(rep(1,abs(nsim-sum(Nsim))),
                          rep(0,length(Nsim)-abs(nsim-sum(Nsim)))))
  }
  else {
    nsim <- Nsim
  }
  for (i in 1:length(nsim)) {
    if (nsim[i] == 0) { nsim[i] <- 1} ## necesary when nsim[i] = 0 to run.
  }
  sm_ARD <- makeSimmap(tree = tre, data = x,
                       model =  ARD$solution,
                       rate.cat = 1, nSim = nsim['ARD'],
                       nCores = 9)
  class(sm_ARD) <- 'multiPhylo'
  sm_SYM <- makeSimmap(tree = tre, data = x,
                       model =  SYM$solution,
                       rate.cat = 1, nSim = nsim['SYM'],
                       nCores = 9)
  class(sm_SYM) <- 'multiPhylo'
  sm_ER <- makeSimmap(tree = tre, data = x,
                      model =  ER$solution,
                      rate.cat = 1, nSim = nsim['ER'],
                      nCores = 9)
  class(sm_ER) <- 'multiPhylo'
  
  sm_trees <- c(sm_ARD, sm_SYM, sm_ER)
  
  if (length(sm_trees) > sim) {sm_trees <- sample(sm_trees, sim, replace = F)} ## when is not sim , aleatory select sim trees
  class(sm_trees) <- 'multiPhylo'
  return(sm_trees)
}

photo_data <- read.csv('Fotosíntesis.csv', header = T)
rownames(photo_data) <- photo_data$sp
# 1000 trees photo
## empty vector to start
sm_trees_p <- simmap_trees(trees[[1]], photo_data, 10) ## starting vector tree 1
sep_smtrees_p <- rep(list(0),1000)
sep_smtrees_p[[1]] <- sm_trees_p ## first tree added to the first item
for (i in 2:1000) {
  sm <- simmap_trees(trees[[i]], photo_data, 10)
  sm_trees_p <- c(sm_trees_p, sm) ## integrating across all sampled trees (1000)
  sep_smtrees_p[[i]] <- sm ## 10 simmaps per tree separated in each list index
}

sep_smtrees_p
# 1000 trees habit
h_data <- read.csv('Hábito_anual(0)_perenne(1).csv', header = T)

## vector to start
sm_trees_h <- simmap_trees(trees[[1]], h_data, 10) ## starting vector tree 1
sep_smtrees_h <- rep(list(0),1000)
sep_smtrees_h[[1]] <- sm_trees_h ## first tree added to the first item
for (i in 2:1000) {
  print(i)
  sm <- simmap_trees(trees[[i]], h_data, 10)
  sm_trees_h <- c(sm_trees_h, sm) ## integrating across all sampled trees (1000)
  sep_smtrees_h[[i]] <- sm ## 10 simmaps per tree separated in each list index
}

## nodes of interest to use above
BD <- c('Gomphrena_perennis', 'Gomphrena_celosioides') 
ED <- c('Gomphrena_macrocephala', 'Gomphrena_pulchella')
wHA <- c('Gomphrena_mendocina', 'Gomphrena_meyeniana')
wLA <- c('Gomphrena_boliviana', 'Gomphrena_martiana')
mr <- c('Gomphrena_mollis', 'Gomphrena_rupestris')
nodes <- list(BD,ED,wHA,wLA,mr)
names(nodes) <- c('BD', 'ED', 'wHA', 'wLA','mr')

## summary per tree of ancestral states - triait habit
annual <- rep(0,1000)
perennial <- rep(0,1000)
nodes.p_h <- array(data = c(annual, perennial), dim = c(1000,2,5))
for (i in 1:1000) {
  sum_h <- describe.simmap(sep_smtrees_h[[i]])
  for(c in 1:length(nodes)) {
    mrca<-findMRCA(trees[[i]], nodes[[c]])
    nodes.p_h[i,,c] <- sum_h$ace[mrca-42,]
  }
}
dev.off()

#pdf('boxplot-probabilities_nodes.pdf')
## to print the boxplots
par(mfrow = c(2,3))
for (i in 1:length(nodes)){
  df_probabilitiesbd_h <- data.frame(
    state = as.factor(c(rep('annual',1000),rep('perennial',1000))),
    probabilities = c(as.vector(nodes.p_h[,1,i]), as.vector(nodes.p_h[,2,i]))
  )
    boxplot(probabilities ~ state , data = df_probabilitiesbd_h, main = names(nodes[i]))
}
dev.off()


## summary per tree of ancestral states - trait photo
photo_0 <- rep(0,1000)
photo_1 <- rep(0,1000)
nodes.p_photo <- array(data = c(photo_0, photo_1), dim = c(1000,2,5)) ## 1000 trees x 2 states x 5 clades to get probabilities at clade nodes per tree per state

for (i in 1:length(sep_smtrees_p)) {
  sum <- describe.simmap(sep_smtrees_p[[i]]) ## summary of tree(i)
  for(c in 1:length(nodes)) {
    mrca <-findMRCA(trees[[i]], nodes[[c]]) ## find the most recent anc of each node of int
    nodes.p_photo[i,,c] <- sum$ace[mrca-42,]
  }
}

#pdf('boxplot-probabilities_nodes_photo.pdf')
par(mfrow = c(2,3))
for (i in 1:length(nodes)){
  df_probabilitiesbd_p <- data.frame(
    state = as.factor(c(rep('c3',1000),rep('c4',1000))),
    probabilities = c(as.vector(nodes.p_photo[,1,i]), as.vector(nodes.p_photo[,2,i]))
  )
  boxplot(probabilities ~ state , data = df_probabilitiesbd_p, main = names(nodes[i]))
}
dev.off()


nodes.p_h ## it is an array where the lines are de trees, the columns are the probability of states (0,1) and the third dim are de caldes
nodes.p_photo ## same for photo
#save.image('sm_trees_incertidumbre.rdata')



