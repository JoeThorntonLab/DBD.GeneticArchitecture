#####################################################
##Explore sequence space with and without epistasis## 
#####################################################
####################################
##1. Packages, data, and functions##
####################################
####Working directory
setwd("/Users/Chimera/Documents/ThorntonLab/Epistasis/2020")

####Packages
library(data.table)
library(Matrix)
library(ordinalNet)
library(Hmisc)
library(ggplot2)
library(boot)
library(gtools)
library(gplots)
library(glmnetcr)
library(dotCall64)
library(gridExtra)
library(grid)
library(seqinr)
library(stringr)
library(vioplot)
library(rgexf)
library(igraph)
library(nloptr)
library(pROC)
library(MatrixModels)
library(lamW)
library(plotrix)

####Custom functions
source("Functions.R")

####Data
#Experimental data
load("DT.JOINT.rda")
load("PREDICT.TE.CLASS.rda")
load("PREDICT.TE.LINK.rda")
load("PREDICT.PE.CLASS.rda")
load("PREDICT.PE.LINK.rda")
load("PREDICT.ME.CLASS.rda")
load("PREDICT.ME.LINK.rda")

##Amino acid states
AAs <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")

##Split S into two groups (S and Z) for analyses that take account of the genetic code
AAs.extend <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","Z")

##Amino acid neighbors given genetic code
AA.NEIGHBORS.S <- read.table(file="AA.neighbors.S.txt",header=T,row.names=1,quote="")
AA.DISTANCES.S <- read.table(file="AA.distance.S.txt",header=T,row.names=1,quote="")
GEN.CODE.NET.S <- graph.adjacency(as.matrix(AA.NEIGHBORS.S), mode="undirected")
GC.ROUTES.S    <- num.paths(GEN.CODE.NET.S,AAs,AAs)

##Amino acid neighbors accounting for unconnected serine codons
AA.NEIGHBORS.Z <- read.table(file="AA.neighbors.Z.txt",header=T,row.names=1,quote="")
AA.DISTANCES.Z <- read.table(file="AA.distance.Z.txt",header=T,row.names=1,quote="")
GEN.CODE.NET.Z <- graph.adjacency(as.matrix(AA.NEIGHBORS.Z), mode="undirected")
GC.ROUTES.Z    <- num.paths(GEN.CODE.NET.Z,AAs.extend,AAs.extend)

##Amino acid neighbors with no genetic code (i.e. hamming distance)
AA.NEIGHBORS.N <- read.table(file="AA.neighbors.N.txt",header=T,row.names=1,quote="")
AA.DISTANCES.N <- read.table(file="AA.distance.N.txt",header=T,row.names=1,quote="")
GEN.CODE.NET.N <- graph.adjacency(as.matrix(AA.NEIGHBORS.N), mode="undirected")
GC.ROUTES.N    <- num.paths(GEN.CODE.NET.N,AAs,AAs)

##Add predicted classes, link values, and joint classification to data table
#DT.JOINT[,PRED.TE.CLASS := PREDICT.TE.CLASS$class] 
#DT.JOINT[,PRED.TE.LINK  := PREDICT.TE.LINK] 
#DT.JOINT[,TE.JOINT.CLASS := 'null']
#DT.JOINT[which(DT.JOINT[1:160000,PRED.TE.CLASS] == "strong" & DT.JOINT[160001:320000,PRED.TE.CLASS] == "strong"),         'TE.JOINT.CLASS'] <- "promiscuous"
#DT.JOINT[which(DT.JOINT[1:160000,PRED.TE.CLASS] == "strong" & DT.JOINT[160001:320000,PRED.TE.CLASS] != "strong"),         'TE.JOINT.CLASS'] <- "ERE-specific"
#DT.JOINT[which(DT.JOINT[1:160000,PRED.TE.CLASS] != "strong" & DT.JOINT[160001:320000,PRED.TE.CLASS] == "strong"),         'TE.JOINT.CLASS'] <- "SRE-specific"
#DT.JOINT[which(DT.JOINT[1:160000,PRED.TE.CLASS] == "strong" & DT.JOINT[160001:320000,PRED.TE.CLASS] == "strong") + 160000,'TE.JOINT.CLASS'] <- "promiscuous"
#DT.JOINT[which(DT.JOINT[1:160000,PRED.TE.CLASS] == "strong" & DT.JOINT[160001:320000,PRED.TE.CLASS] != "strong") + 160000,'TE.JOINT.CLASS'] <- "ERE-specific"
#DT.JOINT[which(DT.JOINT[1:160000,PRED.TE.CLASS] != "strong" & DT.JOINT[160001:320000,PRED.TE.CLASS] == "strong") + 160000,'TE.JOINT.CLASS'] <- "SRE-specific"
#
#DT.JOINT[,PRED.PE.CLASS := PREDICT.PE.CLASS$class] 
#DT.JOINT[,PRED.PE.LINK  := PREDICT.PE.LINK] 
#DT.JOINT[,PE.JOINT.CLASS := 'null']
#DT.JOINT[which(DT.JOINT[1:160000,PRED.PE.CLASS] == "strong" & DT.JOINT[160001:320000,PRED.PE.CLASS] == "strong"),         'PE.JOINT.CLASS'] <- "promiscuous"
#DT.JOINT[which(DT.JOINT[1:160000,PRED.PE.CLASS] == "strong" & DT.JOINT[160001:320000,PRED.PE.CLASS] != "strong"),         'PE.JOINT.CLASS'] <- "ERE-specific"
#DT.JOINT[which(DT.JOINT[1:160000,PRED.PE.CLASS] != "strong" & DT.JOINT[160001:320000,PRED.PE.CLASS] == "strong"),         'PE.JOINT.CLASS'] <- "SRE-specific"
#DT.JOINT[which(DT.JOINT[1:160000,PRED.PE.CLASS] == "strong" & DT.JOINT[160001:320000,PRED.PE.CLASS] == "strong") + 160000,'PE.JOINT.CLASS'] <- "promiscuous"
#DT.JOINT[which(DT.JOINT[1:160000,PRED.PE.CLASS] == "strong" & DT.JOINT[160001:320000,PRED.PE.CLASS] != "strong") + 160000,'PE.JOINT.CLASS'] <- "ERE-specific"
#DT.JOINT[which(DT.JOINT[1:160000,PRED.PE.CLASS] != "strong" & DT.JOINT[160001:320000,PRED.PE.CLASS] == "strong") + 160000,'PE.JOINT.CLASS'] <- "SRE-specific"
#
#DT.JOINT[,PRED.ME.CLASS := PREDICT.ME.CLASS$class] 
#DT.JOINT[,PRED.ME.LINK  := PREDICT.ME.LINK] 
#DT.JOINT[,ME.JOINT.CLASS := 'null']
#DT.JOINT[which(DT.JOINT[1:160000,PRED.ME.CLASS] == "strong" & DT.JOINT[160001:320000,PRED.ME.CLASS] == "strong"),         'ME.JOINT.CLASS'] <- "promiscuous"
#DT.JOINT[which(DT.JOINT[1:160000,PRED.ME.CLASS] == "strong" & DT.JOINT[160001:320000,PRED.ME.CLASS] != "strong"),         'ME.JOINT.CLASS'] <- "ERE-specific"
#DT.JOINT[which(DT.JOINT[1:160000,PRED.ME.CLASS] != "strong" & DT.JOINT[160001:320000,PRED.ME.CLASS] == "strong"),         'ME.JOINT.CLASS'] <- "SRE-specific"
#DT.JOINT[which(DT.JOINT[1:160000,PRED.ME.CLASS] == "strong" & DT.JOINT[160001:320000,PRED.ME.CLASS] == "strong") + 160000,'ME.JOINT.CLASS'] <- "promiscuous"
#DT.JOINT[which(DT.JOINT[1:160000,PRED.ME.CLASS] == "strong" & DT.JOINT[160001:320000,PRED.ME.CLASS] != "strong") + 160000,'ME.JOINT.CLASS'] <- "ERE-specific"
#DT.JOINT[which(DT.JOINT[1:160000,PRED.ME.CLASS] != "strong" & DT.JOINT[160001:320000,PRED.ME.CLASS] == "strong") + 160000,'ME.JOINT.CLASS'] <- "SRE-specific"

##Add amino acid state "Z" to account for unconnected serine codons in genetic code
#DT.JOINT.APPEND <- DT.JOINT[1:320000,][AA1 == "S"]
#DT.JOINT.APPEND[,'AA1'] <- "Z"
#DT.JOINT <- rbind(DT.JOINT,DT.JOINT.APPEND)
#DT.JOINT.APPEND <- DT.JOINT[1:320000,][AA2 == "S"]
#DT.JOINT.APPEND[,'AA2'] <- "Z"
#DT.JOINT <- rbind(DT.JOINT,DT.JOINT.APPEND)
#DT.JOINT.APPEND <- DT.JOINT[1:320000,][AA3 == "S"]
#DT.JOINT.APPEND[,'AA3'] <- "Z"
#DT.JOINT <- rbind(DT.JOINT,DT.JOINT.APPEND)
#DT.JOINT.APPEND <- DT.JOINT[1:320000,][AA4 == "S"]
#DT.JOINT.APPEND[,'AA4'] <- "Z"
#DT.JOINT <- rbind(DT.JOINT,DT.JOINT.APPEND)
#DT.JOINT.APPEND <- DT.JOINT[1:320000,][AA1 == "S" & AA2 == "S"]
#DT.JOINT.APPEND[,'AA1'] <- "Z"
#DT.JOINT.APPEND[,'AA2'] <- "Z"
#DT.JOINT <- rbind(DT.JOINT,DT.JOINT.APPEND)
#DT.JOINT.APPEND <- DT.JOINT[1:320000,][AA1 == "S" & AA3 == "S"]
#DT.JOINT.APPEND[,'AA1'] <- "Z"
#DT.JOINT.APPEND[,'AA3'] <- "Z"
#DT.JOINT <- rbind(DT.JOINT,DT.JOINT.APPEND)
#DT.JOINT.APPEND <- DT.JOINT[1:320000,][AA1 == "S" & AA4 == "S"]
#DT.JOINT.APPEND[,'AA1'] <- "Z"
#DT.JOINT.APPEND[,'AA4'] <- "Z"
#DT.JOINT <- rbind(DT.JOINT,DT.JOINT.APPEND)
#DT.JOINT.APPEND <- DT.JOINT[1:320000,][AA2 == "S" & AA3 == "S"]
#DT.JOINT.APPEND[,'AA2'] <- "Z"
#DT.JOINT.APPEND[,'AA3'] <- "Z"
#DT.JOINT <- rbind(DT.JOINT,DT.JOINT.APPEND)
#DT.JOINT.APPEND <- DT.JOINT[1:320000,][AA2 == "S" & AA4 == "S"]
#DT.JOINT.APPEND[,'AA2'] <- "Z"
#DT.JOINT.APPEND[,'AA4'] <- "Z"
#DT.JOINT <- rbind(DT.JOINT,DT.JOINT.APPEND)
#DT.JOINT.APPEND <- DT.JOINT[1:320000,][AA3 == "S" & AA4 == "S"]
#DT.JOINT.APPEND[,'AA3'] <- "Z"
#DT.JOINT.APPEND[,'AA4'] <- "Z"
#DT.JOINT <- rbind(DT.JOINT,DT.JOINT.APPEND)
#DT.JOINT.APPEND <- DT.JOINT[1:320000,][AA1 == "S" & AA2 == "S" & AA3 == "S"]
#DT.JOINT.APPEND[,'AA1'] <- "Z"
#DT.JOINT.APPEND[,'AA2'] <- "Z"
#DT.JOINT.APPEND[,'AA3'] <- "Z"
#DT.JOINT <- rbind(DT.JOINT,DT.JOINT.APPEND)
#DT.JOINT.APPEND <- DT.JOINT[1:320000,][AA1 == "S" & AA2 == "S" & AA4 == "S"]
#DT.JOINT.APPEND[,'AA1'] <- "Z"
#DT.JOINT.APPEND[,'AA2'] <- "Z"
#DT.JOINT.APPEND[,'AA4'] <- "Z"
#DT.JOINT <- rbind(DT.JOINT,DT.JOINT.APPEND)
#DT.JOINT.APPEND <- DT.JOINT[1:320000,][AA1 == "S" & AA3 == "S" & AA4 == "S"]
#DT.JOINT.APPEND[,'AA1'] <- "Z"
#DT.JOINT.APPEND[,'AA3'] <- "Z"
#DT.JOINT.APPEND[,'AA4'] <- "Z"
#DT.JOINT <- rbind(DT.JOINT,DT.JOINT.APPEND)
#DT.JOINT.APPEND <- DT.JOINT[1:320000,][AA2 == "S" & AA3 == "S" & AA4 == "S"]
#DT.JOINT.APPEND[,'AA2'] <- "Z"
#DT.JOINT.APPEND[,'AA3'] <- "Z"
#DT.JOINT.APPEND[,'AA4'] <- "Z"
#DT.JOINT <- rbind(DT.JOINT,DT.JOINT.APPEND)
#DT.JOINT.APPEND <- DT.JOINT[1:320000,][AA1 == "S" & AA2 == "S" & AA3 == "S" & AA4 == "S"]
#DT.JOINT.APPEND[,'AA1'] <- "Z"
#DT.JOINT.APPEND[,'AA2'] <- "Z"
#DT.JOINT.APPEND[,'AA3'] <- "Z"
#DT.JOINT.APPEND[,'AA4'] <- "Z"
#DT.JOINT <- rbind(DT.JOINT,DT.JOINT.APPEND)
#DT.JOINT[,'AAseq'] <- paste0(DT.JOINT$RE,DT.JOINT$AA1,DT.JOINT$AA2,DT.JOINT$AA3,DT.JOINT$AA4)
#
#DT.JOINT$PRED.TE.CLASS  <- as.factor(DT.JOINT$PRED.TE.CLASS)
#DT.JOINT$TE.JOINT.CLASS <- as.factor(DT.JOINT$TE.JOINT.CLASS)
#DT.JOINT$PRED.PE.CLASS  <- as.factor(DT.JOINT$PRED.PE.CLASS)
#DT.JOINT$PE.JOINT.CLASS <- as.factor(DT.JOINT$PE.JOINT.CLASS)
#DT.JOINT$PRED.ME.CLASS  <- as.factor(DT.JOINT$PRED.ME.CLASS)
#DT.JOINT$ME.JOINT.CLASS <- as.factor(DT.JOINT$ME.JOINT.CLASS)
#save(DT.JOINT,file="DT.JOINT.APPEND.rda")
load("DT.JOINT.APPEND.rda")

###Model Matrixes
load("TE.MATRIX.rda")
load("PE.MATRIX.rda")
load("ME.MATRIX.rda")

###Model Coefficients
load("ME.COEFS.ADJ.rda")
load("PE.COEFS.ADJ.rda")
load("TE.COEFS.ADJ.rda")

###Model Thresholds
load("ME.THRESH.NULL.rda"); load("ME.THRESH.WEAK.rda")
load("PE.THRESH.NULL.rda"); load("PE.THRESH.WEAK.rda")
load("TE.THRESH.NULL.rda"); load("TE.THRESH.WEAK.rda")

###Model intercepts (B0 and S0 effects)
load("ME.COEFS.EFFECT.ADJ.B0.rda"); load("ME.COEFS.EFFECT.ADJ.S0.rda")
load("PE.COEFS.EFFECT.ADJ.B0.rda"); load("PE.COEFS.EFFECT.ADJ.S0.rda")
load("TE.COEFS.EFFECT.ADJ.B0.rda"); load("TE.COEFS.EFFECT.ADJ.S0.rda")

###################################
##2. Single-step changes in class##
###################################
##Sets of genotypes in each class (N=null, W=weak, S=strong) for each RE
TE.ERE.N.LIST <- substr(DT.JOINT$AAseq[(DT.JOINT$RE == "E" & DT.JOINT$PRED.TE.CLASS == "null")],2,5)
TE.ERE.W.LIST <- substr(DT.JOINT$AAseq[(DT.JOINT$RE == "E" & DT.JOINT$PRED.TE.CLASS == "weak")],2,5)
TE.ERE.S.LIST <- substr(DT.JOINT$AAseq[(DT.JOINT$RE == "E" & DT.JOINT$PRED.TE.CLASS == "strong")],2,5)
TE.SRE.N.LIST <- substr(DT.JOINT$AAseq[(DT.JOINT$RE == "S" & DT.JOINT$PRED.TE.CLASS == "null")],2,5)
TE.SRE.W.LIST <- substr(DT.JOINT$AAseq[(DT.JOINT$RE == "S" & DT.JOINT$PRED.TE.CLASS == "weak")],2,5)
TE.SRE.S.LIST <- substr(DT.JOINT$AAseq[(DT.JOINT$RE == "S" & DT.JOINT$PRED.TE.CLASS == "strong")],2,5)

#Sets of genotypes in each class for each RE without separating serine codons
TE.ERE.N.LIST.ZS <- unique(str_replace_all(TE.ERE.N.LIST,"Z","S"))
TE.ERE.W.LIST.ZS <- unique(str_replace_all(TE.ERE.W.LIST,"Z","S"))
TE.ERE.S.LIST.ZS <- unique(str_replace_all(TE.ERE.S.LIST,"Z","S"))
TE.SRE.N.LIST.ZS <- unique(str_replace_all(TE.SRE.N.LIST,"Z","S"))
TE.SRE.W.LIST.ZS <- unique(str_replace_all(TE.SRE.W.LIST,"Z","S"))
TE.SRE.S.LIST.ZS <- unique(str_replace_all(TE.SRE.S.LIST,"Z","S"))

##ERE
#ERE NULL to WEAK
ERE.NW.MUTS <- data.frame(matrix(0,nrow=0,ncol=5)); colnames(ERE.NW.MUTS) <- c("SOURCE","TARGET","SOURCE.AA","SITE","TARGET.AA")
I <- 1:length(TE.ERE.W.LIST.ZS)
for(i in I) {
  NEIGHBORS <- get.neighbors(TE.ERE.W.LIST.ZS[i],code = "N")
  NULL.NEIGHBORS <- which(NEIGHBORS %in% TE.ERE.N.LIST.ZS)
  if(length(NULL.NEIGHBORS) > 0) {
    J <- 1:length(NULL.NEIGHBORS)
    ADD.SET <- data.frame(matrix(0,nrow=length(J),ncol=5))
    colnames(ADD.SET) <- c("SOURCE","TARGET","SOURCE.AA","SITE","TARGET.AA")
    for(j in J) {
      ADD.SET$SOURCE[j] <- NEIGHBORS[NULL.NEIGHBORS[j]]
      ADD.SET$TARGET[j]   <- TE.ERE.W.LIST.ZS[i]
      SOURCE.SPLIT <- unlist(strsplit(ADD.SET$SOURCE[j],split = ""))
      TARGET.SPLIT   <- unlist(strsplit(ADD.SET$TARGET[j],  split = ""))
      ADD.SET$SITE[j] <- which(SOURCE.SPLIT != TARGET.SPLIT)
      ADD.SET$SOURCE.AA[j] <- SOURCE.SPLIT[ADD.SET$SITE[j]]
      ADD.SET$TARGET.AA[j]   <- TARGET.SPLIT[ADD.SET$SITE[j]]
    }
    ERE.NW.MUTS <- rbind(ERE.NW.MUTS,ADD.SET)
  }
}
#Number of WEAK genotypes
length(TE.ERE.W.LIST.ZS)
#Total number of mutations to WEAK
nrow(ERE.NW.MUTS)
#Number of unique NULL SOURCE genotypes
length(unique(ERE.NW.MUTS$SOURCE))
#Number of unique WEAK TARGET genotypes
length(unique(ERE.NW.MUTS$TARGET))
#Mutation types
ERE.NW.MUT.TYPES <- apply(ERE.NW.MUTS,1,function(x) { paste0(x[4],x[3],x[5])}) 
ERE.NW.MUT.TYPES.TABLE <- table(ERE.NW.MUT.TYPES)
ERE.NW.MUT.TYPES.TABLE <- ERE.NW.MUT.TYPES.TABLE[order(ERE.NW.MUT.TYPES.TABLE)]
length(ERE.NW.MUT.TYPES.TABLE)

#ERE NULL to STRONG
ERE.NS.MUTS <- data.frame(matrix(0,nrow=0,ncol=5)); colnames(ERE.NS.MUTS) <- c("SOURCE","TARGET","SOURCE.AA","SITE","TARGET.AA")
I <- 1:length(TE.ERE.S.LIST.ZS)
for(i in I) {
  NEIGHBORS <- get.neighbors(TE.ERE.S.LIST.ZS[i],code = "N")
  NULL.NEIGHBORS <- which(NEIGHBORS %in% TE.ERE.N.LIST.ZS)
  if(length(NULL.NEIGHBORS) > 0) {
    J <- 1:length(NULL.NEIGHBORS)
    ADD.SET <- data.frame(matrix(0,nrow=length(J),ncol=5))
    colnames(ADD.SET) <- c("SOURCE","TARGET","SOURCE.AA","SITE","TARGET.AA")
    for(j in J) {
      ADD.SET$SOURCE[j] <- NEIGHBORS[NULL.NEIGHBORS[j]]
      ADD.SET$TARGET[j]   <- TE.ERE.S.LIST.ZS[i]
      SOURCE.SPLIT <- unlist(strsplit(ADD.SET$SOURCE[j],split = ""))
      TARGET.SPLIT   <- unlist(strsplit(ADD.SET$TARGET[j],  split = ""))
      ADD.SET$SITE[j] <- which(SOURCE.SPLIT != TARGET.SPLIT)
      ADD.SET$SOURCE.AA[j] <- SOURCE.SPLIT[ADD.SET$SITE[j]]
      ADD.SET$TARGET.AA[j]   <- TARGET.SPLIT[ADD.SET$SITE[j]]
    }
    ERE.NS.MUTS <- rbind(ERE.NS.MUTS,ADD.SET)
  }
}
#Number of STRONG genotypes
length(TE.ERE.S.LIST.ZS)
#Total number of mutations to STRONG
nrow(ERE.NS.MUTS)
#Number of unique NULL SOURCE genotypes
length(unique(ERE.NS.MUTS$SOURCE))
#Number of unique STRONG TARGET genotypes
length(unique(ERE.NS.MUTS$TARGET))
#Mutation types
ERE.NS.MUT.TYPES <- apply(ERE.NS.MUTS,1,function(x) { paste0(x[4],x[3],x[5])}) 
ERE.NS.MUT.TYPES.TABLE <- table(ERE.NS.MUT.TYPES)
ERE.NS.MUT.TYPES.TABLE <- ERE.NS.MUT.TYPES.TABLE[order(ERE.NS.MUT.TYPES.TABLE)]
length(ERE.NS.MUT.TYPES.TABLE)

#ERE WEAK to STRONG
ERE.WS.MUTS <- data.frame(matrix(0,nrow=0,ncol=5)); colnames(ERE.WS.MUTS) <- c("SOURCE","TARGET","SOURCE.AA","SITE","TARGET.AA")
I <- 1:length(TE.ERE.S.LIST.ZS)
for(i in I) {
  NEIGHBORS <- get.neighbors(TE.ERE.S.LIST.ZS[i],code = "N")
  WEAK.NEIGHBORS <- which(NEIGHBORS %in% TE.ERE.W.LIST.ZS)
  if(length(WEAK.NEIGHBORS) > 0) {
    J <- 1:length(WEAK.NEIGHBORS)
    ADD.SET <- data.frame(matrix(0,nrow=length(J),ncol=5))
    colnames(ADD.SET) <- c("SOURCE","TARGET","SOURCE.AA","SITE","TARGET.AA")
    for(j in J) {
      ADD.SET$SOURCE[j] <- NEIGHBORS[WEAK.NEIGHBORS[j]]
      ADD.SET$TARGET[j]   <- TE.ERE.S.LIST.ZS[i]
      SOURCE.SPLIT <- unlist(strsplit(ADD.SET$SOURCE[j],split = ""))
      TARGET.SPLIT   <- unlist(strsplit(ADD.SET$TARGET[j],  split = ""))
      ADD.SET$SITE[j] <- which(SOURCE.SPLIT != TARGET.SPLIT)
      ADD.SET$SOURCE.AA[j] <- SOURCE.SPLIT[ADD.SET$SITE[j]]
      ADD.SET$TARGET.AA[j]   <- TARGET.SPLIT[ADD.SET$SITE[j]]
    }
    ERE.WS.MUTS <- rbind(ERE.WS.MUTS,ADD.SET)
  }
}
#Number of STRONG genotypes
length(TE.ERE.S.LIST.ZS)
#Total number of mutations to STRONG
nrow(ERE.WS.MUTS)
#Number of unique WEAK SOURCE genotypes
length(unique(ERE.WS.MUTS$SOURCE))
#Number of unique STRONG TARGET genotypes
length(unique(ERE.WS.MUTS$TARGET))
#Mutation types
ERE.WS.MUT.TYPES <- apply(ERE.WS.MUTS,1,function(x) { paste0(x[4],x[3],x[5])}) 
ERE.WS.MUT.TYPES.TABLE <- table(ERE.WS.MUT.TYPES)
ERE.WS.MUT.TYPES.TABLE <- ERE.WS.MUT.TYPES.TABLE[order(ERE.WS.MUT.TYPES.TABLE)]
length(ERE.WS.MUT.TYPES.TABLE)

##SRE
#SRE NULL to WEAK
SRE.NW.MUTS <- data.frame(matrix(0,nrow=0,ncol=5)); colnames(SRE.NW.MUTS) <- c("SOURCE","TARGET","SOURCE.AA","SITE","TARGET.AA")
I <- 1:length(TE.SRE.W.LIST.ZS)
for(i in I) {
  NEIGHBORS <- get.neighbors(TE.SRE.W.LIST.ZS[i],code = "N")
  NULL.NEIGHBORS <- which(NEIGHBORS %in% TE.SRE.N.LIST.ZS)
  if(length(NULL.NEIGHBORS) > 0) {
    J <- 1:length(NULL.NEIGHBORS)
    ADD.SET <- data.frame(matrix(0,nrow=length(J),ncol=5))
    colnames(ADD.SET) <- c("SOURCE","TARGET","SOURCE.AA","SITE","TARGET.AA")
    for(j in J) {
      ADD.SET$SOURCE[j] <- NEIGHBORS[NULL.NEIGHBORS[j]]
      ADD.SET$TARGET[j]   <- TE.SRE.W.LIST.ZS[i]
      SOURCE.SPLIT <- unlist(strsplit(ADD.SET$SOURCE[j],split = ""))
      TARGET.SPLIT   <- unlist(strsplit(ADD.SET$TARGET[j],  split = ""))
      ADD.SET$SITE[j] <- which(SOURCE.SPLIT != TARGET.SPLIT)
      ADD.SET$SOURCE.AA[j] <- SOURCE.SPLIT[ADD.SET$SITE[j]]
      ADD.SET$TARGET.AA[j]   <- TARGET.SPLIT[ADD.SET$SITE[j]]
    }
    SRE.NW.MUTS <- rbind(SRE.NW.MUTS,ADD.SET)
  }
}
#Number of WEAK genotypes
length(TE.SRE.W.LIST.ZS)
#Total number of mutations to WEAK
nrow(SRE.NW.MUTS)
#Number of unique NULL SOURCE genotypes
length(unique(SRE.NW.MUTS$SOURCE))
#Number of unique WEAK TARGET genotypes
length(unique(SRE.NW.MUTS$TARGET))
#Mutation types
SRE.NW.MUT.TYPES <- apply(SRE.NW.MUTS,1,function(x) { paste0(x[4],x[3],x[5])}) 
SRE.NW.MUT.TYPES.TABLE <- table(SRE.NW.MUT.TYPES)
SRE.NW.MUT.TYPES.TABLE <- SRE.NW.MUT.TYPES.TABLE[order(SRE.NW.MUT.TYPES.TABLE)]
length(SRE.NW.MUT.TYPES.TABLE)

#SRE NULL to STRONG
SRE.NS.MUTS <- data.frame(matrix(0,nrow=0,ncol=5)); colnames(SRE.NS.MUTS) <- c("SOURCE","TARGET","SOURCE.AA","SITE","TARGET.AA")
I <- 1:length(TE.SRE.S.LIST.ZS)
for(i in I) {
  NEIGHBORS <- get.neighbors(TE.SRE.S.LIST.ZS[i],code = "N")
  NULL.NEIGHBORS <- which(NEIGHBORS %in% TE.SRE.N.LIST.ZS)
  if(length(NULL.NEIGHBORS) > 0) {
    J <- 1:length(NULL.NEIGHBORS)
    ADD.SET <- data.frame(matrix(0,nrow=length(J),ncol=5))
    colnames(ADD.SET) <- c("SOURCE","TARGET","SOURCE.AA","SITE","TARGET.AA")
    for(j in J) {
      ADD.SET$SOURCE[j] <- NEIGHBORS[NULL.NEIGHBORS[j]]
      ADD.SET$TARGET[j]   <- TE.SRE.S.LIST.ZS[i]
      SOURCE.SPLIT <- unlist(strsplit(ADD.SET$SOURCE[j],split = ""))
      TARGET.SPLIT   <- unlist(strsplit(ADD.SET$TARGET[j],  split = ""))
      ADD.SET$SITE[j] <- which(SOURCE.SPLIT != TARGET.SPLIT)
      ADD.SET$SOURCE.AA[j] <- SOURCE.SPLIT[ADD.SET$SITE[j]]
      ADD.SET$TARGET.AA[j]   <- TARGET.SPLIT[ADD.SET$SITE[j]]
    }
    SRE.NS.MUTS <- rbind(SRE.NS.MUTS,ADD.SET)
  }
}
#Number of STRONG genotypes
length(TE.SRE.S.LIST.ZS)
#Total number of mutations to STRONG
nrow(SRE.NS.MUTS)
#Number of unique NULL SOURCE genotypes
length(unique(SRE.NS.MUTS$SOURCE))
#Number of unique STRONG TARGET genotypes
length(unique(SRE.NS.MUTS$TARGET))
#Mutation types
SRE.NS.MUT.TYPES <- apply(SRE.NS.MUTS,1,function(x) { paste0(x[4],x[3],x[5])}) 
SRE.NS.MUT.TYPES.TABLE <- table(SRE.NS.MUT.TYPES)
SRE.NS.MUT.TYPES.TABLE <- SRE.NS.MUT.TYPES.TABLE[order(SRE.NS.MUT.TYPES.TABLE)]
length(SRE.NS.MUT.TYPES.TABLE)

#SRE WEAK to STRONG
SRE.WS.MUTS <- data.frame(matrix(0,nrow=0,ncol=5)); colnames(SRE.WS.MUTS) <- c("SOURCE","TARGET","SOURCE.AA","SITE","TARGET.AA")
I <- 1:length(TE.SRE.S.LIST.ZS)
for(i in I) {
  NEIGHBORS <- get.neighbors(TE.SRE.S.LIST.ZS[i],code = "N")
  WEAK.NEIGHBORS <- which(NEIGHBORS %in% TE.SRE.W.LIST.ZS)
  if(length(WEAK.NEIGHBORS) > 0) {
    J <- 1:length(WEAK.NEIGHBORS)
    ADD.SET <- data.frame(matrix(0,nrow=length(J),ncol=5))
    colnames(ADD.SET) <- c("SOURCE","TARGET","SOURCE.AA","SITE","TARGET.AA")
    for(j in J) {
      ADD.SET$SOURCE[j] <- NEIGHBORS[WEAK.NEIGHBORS[j]]
      ADD.SET$TARGET[j]   <- TE.SRE.S.LIST.ZS[i]
      SOURCE.SPLIT <- unlist(strsplit(ADD.SET$SOURCE[j],split = ""))
      TARGET.SPLIT   <- unlist(strsplit(ADD.SET$TARGET[j],  split = ""))
      ADD.SET$SITE[j] <- which(SOURCE.SPLIT != TARGET.SPLIT)
      ADD.SET$SOURCE.AA[j] <- SOURCE.SPLIT[ADD.SET$SITE[j]]
      ADD.SET$TARGET.AA[j]   <- TARGET.SPLIT[ADD.SET$SITE[j]]
    }
    SRE.WS.MUTS <- rbind(SRE.WS.MUTS,ADD.SET)
  }
}
#Number of STRONG genotypes
length(TE.SRE.S.LIST.ZS)
#Total number of mutations to STRONG
nrow(SRE.WS.MUTS)
#Number of unique WEAK SOURCE genotypes
length(unique(SRE.WS.MUTS$SOURCE))
#Number of unique STRONG TARGET genotypes
length(unique(SRE.WS.MUTS$TARGET))
#Mutation types
SRE.WS.MUT.TYPES <- apply(SRE.WS.MUTS,1,function(x) { paste0(x[4],x[3],x[5])}) 
SRE.WS.MUT.TYPES.TABLE <- table(SRE.WS.MUT.TYPES)
SRE.WS.MUT.TYPES.TABLE <- SRE.WS.MUT.TYPES.TABLE[order(SRE.WS.MUT.TYPES.TABLE)]
length(SRE.WS.MUT.TYPES.TABLE)

##Mutation types, regardless of change in activator class or RE
TOTAL.MUT.TYPES <- c(ERE.NW.MUT.TYPES,SRE.NW.MUT.TYPES,ERE.NS.MUT.TYPES,SRE.NS.MUT.TYPES,ERE.WS.MUT.TYPES,SRE.WS.MUT.TYPES)
TOTAL.MUT.TYPES.TABLE <- table(TOTAL.MUT.TYPES)

##Mutation types for each change in activator class, regardless of RE
TOTAL.NW.MUT.TYPES <- c(ERE.NW.MUT.TYPES,SRE.NW.MUT.TYPES)
TOTAL.NS.MUT.TYPES <- c(ERE.NS.MUT.TYPES,SRE.NS.MUT.TYPES)
TOTAL.WS.MUT.TYPES <- c(ERE.WS.MUT.TYPES,SRE.WS.MUT.TYPES)
TOTAL.NW.MUT.TYPES.TABLE <- table(TOTAL.NW.MUT.TYPES)
TOTAL.NS.MUT.TYPES.TABLE <- table(TOTAL.NS.MUT.TYPES)
TOTAL.WS.MUT.TYPES.TABLE <- table(TOTAL.WS.MUT.TYPES)

#Mutation types for all three transitions
sum(names(TOTAL.NW.MUT.TYPES.TABLE) %in% names(TOTAL.NS.MUT.TYPES.TABLE) & names(TOTAL.NW.MUT.TYPES.TABLE) %in% names(TOTAL.WS.MUT.TYPES.TABLE))
#Mutations types for two transitions
sum(  names(TOTAL.NW.MUT.TYPES.TABLE) %in% names(TOTAL.NS.MUT.TYPES.TABLE)  & !(names(TOTAL.NW.MUT.TYPES.TABLE) %in% names(TOTAL.WS.MUT.TYPES.TABLE)))
sum(!(names(TOTAL.NW.MUT.TYPES.TABLE) %in% names(TOTAL.NS.MUT.TYPES.TABLE)) &  (names(TOTAL.NW.MUT.TYPES.TABLE) %in% names(TOTAL.WS.MUT.TYPES.TABLE)))
sum(  names(TOTAL.NS.MUT.TYPES.TABLE) %in% names(TOTAL.WS.MUT.TYPES.TABLE)  & !(names(TOTAL.NS.MUT.TYPES.TABLE) %in% names(TOTAL.NW.MUT.TYPES.TABLE)))
#Mutations types for one transition
sum(!(names(TOTAL.NW.MUT.TYPES.TABLE) %in% names(TOTAL.NS.MUT.TYPES.TABLE)) & !(names(TOTAL.NW.MUT.TYPES.TABLE) %in% names(TOTAL.WS.MUT.TYPES.TABLE)))
sum(!(names(TOTAL.NS.MUT.TYPES.TABLE) %in% names(TOTAL.NW.MUT.TYPES.TABLE)) & !(names(TOTAL.NS.MUT.TYPES.TABLE) %in% names(TOTAL.WS.MUT.TYPES.TABLE)))
sum(!(names(TOTAL.WS.MUT.TYPES.TABLE) %in% names(TOTAL.NW.MUT.TYPES.TABLE)) & !(names(TOTAL.WS.MUT.TYPES.TABLE) %in% names(TOTAL.NS.MUT.TYPES.TABLE)))

##Mutation types for changes on one RE, regardless of activator class
TOTAL.ERE.MUT.TYPES <- c(ERE.NW.MUT.TYPES,ERE.NS.MUT.TYPES,ERE.WS.MUT.TYPES)
TOTAL.SRE.MUT.TYPES <- c(SRE.NW.MUT.TYPES,SRE.NS.MUT.TYPES,SRE.WS.MUT.TYPES)
TOTAL.ERE.MUT.TYPES.TABLE <- table(TOTAL.ERE.MUT.TYPES)
TOTAL.SRE.MUT.TYPES.TABLE <- table(TOTAL.SRE.MUT.TYPES)

#Mutation types on both REs
sum(names(TOTAL.ERE.MUT.TYPES.TABLE) %in% names(TOTAL.SRE.MUT.TYPES.TABLE))
#Mutation types on one RE
sum(!(names(TOTAL.ERE.MUT.TYPES.TABLE) %in% names(TOTAL.SRE.MUT.TYPES.TABLE)))
sum(!(names(TOTAL.SRE.MUT.TYPES.TABLE) %in% names(TOTAL.ERE.MUT.TYPES.TABLE)))

###Decompose contributions to identify necessity and sufficientcy
load("EFFECT.TABLE.TE.rda")

##ERE NULL to WEAK
EFFECT.TABLE.TE.ERE <- EFFECT.TABLE.TE[1:160000,]
EFFECT.TABLE.TE.ERE$SEQ <- substr(EFFECT.TABLE.TE.ERE$SEQ,2,5)

ERE.NW.MUTS.OUT <- data.frame(matrix(0,nrow=nrow(ERE.NW.MUTS),ncol=5))
colnames(ERE.NW.MUTS.OUT) <- c("ERE.SOURCE","ERE.TARGET","ERE.1","ERE.2","ERE.3")

ERE.NW.MUTS.OUT$ERE.SOURCE <- rowSums(EFFECT.TABLE.TE.ERE[apply(ERE.NW.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.ERE$SEQ == x[1]) }),3:32])
ERE.NW.MUTS.OUT$ERE.TARGET <- rowSums(EFFECT.TABLE.TE.ERE[apply(ERE.NW.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.ERE$SEQ == x[2]) }),3:32])
ERE.NW.MUTS.OUT$ERE.1      <- rowSums(EFFECT.TABLE.TE.ERE[apply(ERE.NW.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.ERE$SEQ == x[2]) }),5:12])  - rowSums(EFFECT.TABLE.TE.ERE[apply(ERE.NW.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.ERE$SEQ == x[1]) }),5:12])
ERE.NW.MUTS.OUT$ERE.2      <- rowSums(EFFECT.TABLE.TE.ERE[apply(ERE.NW.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.ERE$SEQ == x[2]) }),13:24]) - rowSums(EFFECT.TABLE.TE.ERE[apply(ERE.NW.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.ERE$SEQ == x[1]) }),13:24])
ERE.NW.MUTS.OUT$ERE.3      <- rowSums(EFFECT.TABLE.TE.ERE[apply(ERE.NW.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.ERE$SEQ == x[2]) }),25:32]) - rowSums(EFFECT.TABLE.TE.ERE[apply(ERE.NW.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.ERE$SEQ == x[1]) }),25:32])

#Epistasis necessary or sufficient
ERE.NW.EPI.NEC <- which((ERE.NW.MUTS.OUT$ERE.TARGET - ERE.NW.MUTS.OUT$ERE.2 - ERE.NW.MUTS.OUT$ERE.3) < (-1*TE.THRESH.NULL))
ERE.NW.EPI.SUF <- which((ERE.NW.MUTS.OUT$ERE.TARGET - ERE.NW.MUTS.OUT$ERE.1) > (-1*TE.THRESH.NULL))


##ERE NULL to STRONG
EFFECT.TABLE.TE.ERE <- EFFECT.TABLE.TE[1:160000,]
EFFECT.TABLE.TE.ERE$SEQ <- substr(EFFECT.TABLE.TE.ERE$SEQ,2,5)

ERE.NS.MUTS.OUT <- data.frame(matrix(0,nrow=nrow(ERE.NS.MUTS),ncol=5))
colnames(ERE.NS.MUTS.OUT) <- c("ERE.SOURCE","ERE.TARGET","ERE.1","ERE.2","ERE.3")

ERE.NS.MUTS.OUT$ERE.SOURCE <- rowSums(EFFECT.TABLE.TE.ERE[apply(ERE.NS.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.ERE$SEQ == x[1]) }),3:32])
ERE.NS.MUTS.OUT$ERE.TARGET <- rowSums(EFFECT.TABLE.TE.ERE[apply(ERE.NS.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.ERE$SEQ == x[2]) }),3:32])
ERE.NS.MUTS.OUT$ERE.1      <- rowSums(EFFECT.TABLE.TE.ERE[apply(ERE.NS.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.ERE$SEQ == x[2]) }),5:12])  - rowSums(EFFECT.TABLE.TE.ERE[apply(ERE.NS.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.ERE$SEQ == x[1]) }),5:12])
ERE.NS.MUTS.OUT$ERE.2      <- rowSums(EFFECT.TABLE.TE.ERE[apply(ERE.NS.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.ERE$SEQ == x[2]) }),13:24]) - rowSums(EFFECT.TABLE.TE.ERE[apply(ERE.NS.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.ERE$SEQ == x[1]) }),13:24])
ERE.NS.MUTS.OUT$ERE.3      <- rowSums(EFFECT.TABLE.TE.ERE[apply(ERE.NS.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.ERE$SEQ == x[2]) }),25:32]) - rowSums(EFFECT.TABLE.TE.ERE[apply(ERE.NS.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.ERE$SEQ == x[1]) }),25:32])

#Epistasis necessary or sufficient
ERE.NS.EPI.NEC <- which((ERE.NS.MUTS.OUT$ERE.TARGET - ERE.NS.MUTS.OUT$ERE.2 - ERE.NS.MUTS.OUT$ERE.3) < (-1*TE.THRESH.WEAK))
ERE.NS.EPI.SUF <- which((ERE.NS.MUTS.OUT$ERE.TARGET - ERE.NS.MUTS.OUT$ERE.1) > (-1*TE.THRESH.WEAK))


##ERE WEAK to STRONG
EFFECT.TABLE.TE.ERE <- EFFECT.TABLE.TE[1:160000,]
EFFECT.TABLE.TE.ERE$SEQ <- substr(EFFECT.TABLE.TE.ERE$SEQ,2,5)

ERE.WS.MUTS.OUT <- data.frame(matrix(0,nrow=nrow(ERE.WS.MUTS),ncol=5))
colnames(ERE.WS.MUTS.OUT) <- c("ERE.SOURCE","ERE.TARGET","ERE.1","ERE.2","ERE.3")

ERE.WS.MUTS.OUT$ERE.SOURCE <- rowSums(EFFECT.TABLE.TE.ERE[apply(ERE.WS.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.ERE$SEQ == x[1]) }),3:32])
ERE.WS.MUTS.OUT$ERE.TARGET <- rowSums(EFFECT.TABLE.TE.ERE[apply(ERE.WS.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.ERE$SEQ == x[2]) }),3:32])
ERE.WS.MUTS.OUT$ERE.1      <- rowSums(EFFECT.TABLE.TE.ERE[apply(ERE.WS.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.ERE$SEQ == x[2]) }),5:12])  - rowSums(EFFECT.TABLE.TE.ERE[apply(ERE.WS.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.ERE$SEQ == x[1]) }),5:12])
ERE.WS.MUTS.OUT$ERE.2      <- rowSums(EFFECT.TABLE.TE.ERE[apply(ERE.WS.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.ERE$SEQ == x[2]) }),13:24]) - rowSums(EFFECT.TABLE.TE.ERE[apply(ERE.WS.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.ERE$SEQ == x[1]) }),13:24])
ERE.WS.MUTS.OUT$ERE.3      <- rowSums(EFFECT.TABLE.TE.ERE[apply(ERE.WS.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.ERE$SEQ == x[2]) }),25:32]) - rowSums(EFFECT.TABLE.TE.ERE[apply(ERE.WS.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.ERE$SEQ == x[1]) }),25:32])

#Epistasis necessary or sufficient
ERE.WS.EPI.NEC <- which((ERE.WS.MUTS.OUT$ERE.TARGET - ERE.WS.MUTS.OUT$ERE.2 - ERE.WS.MUTS.OUT$ERE.3) < (-1*TE.THRESH.WEAK))
ERE.WS.EPI.SUF <- which((ERE.WS.MUTS.OUT$ERE.TARGET - ERE.WS.MUTS.OUT$ERE.1) > (-1*TE.THRESH.WEAK))


##SRE NULL to WEAK
EFFECT.TABLE.TE.SRE <- EFFECT.TABLE.TE[160001:320000,]
EFFECT.TABLE.TE.SRE$SEQ <- substr(EFFECT.TABLE.TE.SRE$SEQ,2,5)

SRE.NW.MUTS.OUT <- data.frame(matrix(0,nrow=nrow(SRE.NW.MUTS),ncol=5))
colnames(SRE.NW.MUTS.OUT) <- c("SRE.SOURCE","SRE.TARGET","SRE.1","SRE.2","SRE.3")

SRE.NW.MUTS.OUT$SRE.SOURCE <- rowSums(EFFECT.TABLE.TE.SRE[apply(SRE.NW.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.SRE$SEQ == x[1]) }),3:32])
SRE.NW.MUTS.OUT$SRE.TARGET <- rowSums(EFFECT.TABLE.TE.SRE[apply(SRE.NW.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.SRE$SEQ == x[2]) }),3:32])
SRE.NW.MUTS.OUT$SRE.1      <- rowSums(EFFECT.TABLE.TE.SRE[apply(SRE.NW.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.SRE$SEQ == x[2]) }),5:12])  - rowSums(EFFECT.TABLE.TE.SRE[apply(SRE.NW.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.SRE$SEQ == x[1]) }),5:12])
SRE.NW.MUTS.OUT$SRE.2      <- rowSums(EFFECT.TABLE.TE.SRE[apply(SRE.NW.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.SRE$SEQ == x[2]) }),13:24]) - rowSums(EFFECT.TABLE.TE.SRE[apply(SRE.NW.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.SRE$SEQ == x[1]) }),13:24])
SRE.NW.MUTS.OUT$SRE.3      <- rowSums(EFFECT.TABLE.TE.SRE[apply(SRE.NW.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.SRE$SEQ == x[2]) }),25:32]) - rowSums(EFFECT.TABLE.TE.SRE[apply(SRE.NW.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.SRE$SEQ == x[1]) }),25:32])

#Epistasis necessary or sufficient
SRE.NW.EPI.NEC <- which((SRE.NW.MUTS.OUT$SRE.TARGET - SRE.NW.MUTS.OUT$SRE.2 - SRE.NW.MUTS.OUT$SRE.3) < (-1*TE.THRESH.NULL))
SRE.NW.EPI.SUF <- which((SRE.NW.MUTS.OUT$SRE.TARGET - SRE.NW.MUTS.OUT$SRE.1) > (-1*TE.THRESH.NULL))


##SRE NULL to STRONG
SRE.NS.MUTS.OUT <- data.frame(matrix(0,nrow=nrow(SRE.NS.MUTS),ncol=5))
colnames(SRE.NS.MUTS.OUT) <- c("SRE.SOURCE","SRE.TARGET","SRE.1","SRE.2","SRE.3")

SRE.NS.MUTS.OUT$SRE.SOURCE <- rowSums(EFFECT.TABLE.TE.SRE[apply(SRE.NS.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.SRE$SEQ == x[1]) }),3:32])
SRE.NS.MUTS.OUT$SRE.TARGET <- rowSums(EFFECT.TABLE.TE.SRE[apply(SRE.NS.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.SRE$SEQ == x[2]) }),3:32])
SRE.NS.MUTS.OUT$SRE.1      <- rowSums(EFFECT.TABLE.TE.SRE[apply(SRE.NS.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.SRE$SEQ == x[2]) }),5:12])  - rowSums(EFFECT.TABLE.TE.SRE[apply(SRE.NS.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.SRE$SEQ == x[1]) }),5:12])
SRE.NS.MUTS.OUT$SRE.2      <- rowSums(EFFECT.TABLE.TE.SRE[apply(SRE.NS.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.SRE$SEQ == x[2]) }),13:24]) - rowSums(EFFECT.TABLE.TE.SRE[apply(SRE.NS.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.SRE$SEQ == x[1]) }),13:24])
SRE.NS.MUTS.OUT$SRE.3      <- rowSums(EFFECT.TABLE.TE.SRE[apply(SRE.NS.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.SRE$SEQ == x[2]) }),25:32]) - rowSums(EFFECT.TABLE.TE.SRE[apply(SRE.NS.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.SRE$SEQ == x[1]) }),25:32])

#Epistasis necessary or sufficient
SRE.NS.EPI.NEC <- which((SRE.NS.MUTS.OUT$SRE.TARGET - SRE.NS.MUTS.OUT$SRE.2 - SRE.NS.MUTS.OUT$SRE.3) < (-1*TE.THRESH.WEAK))
SRE.NS.EPI.SUF <- which((SRE.NS.MUTS.OUT$SRE.TARGET - SRE.NS.MUTS.OUT$SRE.1) > (-1*TE.THRESH.WEAK))


##SRE WEAK to STRONG
SRE.WS.MUTS.OUT <- data.frame(matrix(0,nrow=nrow(SRE.WS.MUTS),ncol=5))
colnames(SRE.WS.MUTS.OUT) <- c("SRE.SOURCE","SRE.TARGET","SRE.1","SRE.2","SRE.3")

SRE.WS.MUTS.OUT$SRE.SOURCE <- rowSums(EFFECT.TABLE.TE.SRE[apply(SRE.WS.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.SRE$SEQ == x[1]) }),3:32])
SRE.WS.MUTS.OUT$SRE.TARGET <- rowSums(EFFECT.TABLE.TE.SRE[apply(SRE.WS.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.SRE$SEQ == x[2]) }),3:32])
SRE.WS.MUTS.OUT$SRE.1      <- rowSums(EFFECT.TABLE.TE.SRE[apply(SRE.WS.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.SRE$SEQ == x[2]) }),5:12])  - rowSums(EFFECT.TABLE.TE.SRE[apply(SRE.WS.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.SRE$SEQ == x[1]) }),5:12])
SRE.WS.MUTS.OUT$SRE.2      <- rowSums(EFFECT.TABLE.TE.SRE[apply(SRE.WS.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.SRE$SEQ == x[2]) }),13:24]) - rowSums(EFFECT.TABLE.TE.SRE[apply(SRE.WS.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.SRE$SEQ == x[1]) }),13:24])
SRE.WS.MUTS.OUT$SRE.3      <- rowSums(EFFECT.TABLE.TE.SRE[apply(SRE.WS.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.SRE$SEQ == x[2]) }),25:32]) - rowSums(EFFECT.TABLE.TE.SRE[apply(SRE.WS.MUTS,1,FUN=function(x) { which(EFFECT.TABLE.TE.SRE$SEQ == x[1]) }),25:32])

#Epistasis necessary or sufficient
SRE.WS.EPI.NEC <- which((SRE.WS.MUTS.OUT$SRE.TARGET - SRE.WS.MUTS.OUT$SRE.2 - SRE.WS.MUTS.OUT$SRE.3) < (-1*TE.THRESH.WEAK))
SRE.WS.EPI.SUF <- which((SRE.WS.MUTS.OUT$SRE.TARGET - SRE.WS.MUTS.OUT$SRE.1) > (-1*TE.THRESH.WEAK))



##Epistasis contributes positively
sum(rowSums(ERE.NS.MUTS.OUT[,4:5]) > 0)/nrow(ERE.NS.MUTS.OUT)
sum(rowSums(ERE.NW.MUTS.OUT[,4:5]) > 0)/nrow(ERE.NW.MUTS.OUT)
sum(rowSums(ERE.WS.MUTS.OUT[,4:5]) > 0)/nrow(ERE.WS.MUTS.OUT)
sum(rowSums(SRE.NS.MUTS.OUT[,4:5]) > 0)/nrow(SRE.NS.MUTS.OUT)
sum(rowSums(SRE.NW.MUTS.OUT[,4:5]) > 0)/nrow(SRE.NW.MUTS.OUT)
sum(rowSums(SRE.WS.MUTS.OUT[,4:5]) > 0)/nrow(SRE.WS.MUTS.OUT)

(sum(rowSums(ERE.NS.MUTS.OUT[,4:5]) > 0) + sum(rowSums(ERE.NW.MUTS.OUT[,4:5]) > 0) + sum(rowSums(ERE.WS.MUTS.OUT[,4:5]) > 0) + sum(rowSums(SRE.NS.MUTS.OUT[,4:5]) > 0) + sum(rowSums(SRE.NW.MUTS.OUT[,4:5]) > 0) + sum(rowSums(SRE.WS.MUTS.OUT[,4:5]) > 0))/
  (nrow(ERE.NS.MUTS.OUT)  + nrow(ERE.NW.MUTS.OUT)  + nrow(ERE.WS.MUTS.OUT)  + nrow(SRE.NS.MUTS.OUT)  + nrow(SRE.NW.MUTS.OUT)  + nrow(SRE.WS.MUTS.OUT))


##Epistasis necessary
length(ERE.NS.EPI.NEC)/nrow(ERE.NS.MUTS.OUT)
length(ERE.NW.EPI.NEC)/nrow(ERE.NW.MUTS.OUT)
length(ERE.WS.EPI.NEC)/nrow(ERE.WS.MUTS.OUT)
length(SRE.NS.EPI.NEC)/nrow(SRE.NS.MUTS.OUT)
length(SRE.NW.EPI.NEC)/nrow(SRE.NW.MUTS.OUT)
length(SRE.WS.EPI.NEC)/nrow(SRE.WS.MUTS.OUT)

(length(ERE.NS.EPI.NEC) + length(ERE.NW.EPI.NEC) + length(ERE.WS.EPI.NEC) + length(SRE.NS.EPI.NEC) + length(SRE.NW.EPI.NEC) + length(SRE.WS.EPI.NEC))/
(nrow(ERE.NS.MUTS.OUT)  + nrow(ERE.NW.MUTS.OUT)  + nrow(ERE.WS.MUTS.OUT)  + nrow(SRE.NS.MUTS.OUT)  + nrow(SRE.NW.MUTS.OUT)  + nrow(SRE.WS.MUTS.OUT))


##Epistasis sufficient
length(ERE.NS.EPI.SUF)/nrow(ERE.NS.MUTS.OUT)
length(ERE.NW.EPI.SUF)/nrow(ERE.NW.MUTS.OUT)
length(ERE.WS.EPI.SUF)/nrow(ERE.WS.MUTS.OUT)
length(SRE.NS.EPI.SUF)/nrow(SRE.NS.MUTS.OUT)
length(SRE.NW.EPI.SUF)/nrow(SRE.NW.MUTS.OUT)
length(SRE.WS.EPI.SUF)/nrow(SRE.WS.MUTS.OUT)

(length(ERE.NS.EPI.SUF) + length(ERE.NW.EPI.SUF) + length(ERE.WS.EPI.SUF) + length(SRE.NS.EPI.SUF) + length(SRE.NW.EPI.SUF) + length(SRE.WS.EPI.SUF))/
  (nrow(ERE.NS.MUTS.OUT)  + nrow(ERE.NW.MUTS.OUT)  + nrow(ERE.WS.MUTS.OUT)  + nrow(SRE.NS.MUTS.OUT)  + nrow(SRE.NW.MUTS.OUT)  + nrow(SRE.WS.MUTS.OUT))

########################
##3. Network distances##
########################
##First letter is the order of epistasis in the model used to assign function; N=0, M=1, P=2, T=3
##Second letter is type of distance used; G=Genetic code, H=Hamming Distance

####All activators
###TE Model
##Get list of activators
TG.ACT.LIST <- unique(substr(DT.JOINT[TE.JOINT.CLASS != 'null',]$AAseq,2,5))
TH.ACT.LIST <- TG.ACT.LIST[-grep("Z",TG.ACT.LIST)]

##Make adjacency matrixes
#TG.ACT.INDEX <- index.adjaceny(TG.ACT.LIST,code="Z")
#TG.ACT.ADJM  <- make.adjaceny(TG.ACT.LIST,TG.ACT.INDEX)
#colnames(TG.ACT.ADJM) <- TG.ACT.LIST; rownames(TG.ACT.ADJM) <- TG.ACT.LIST
#TH.ACT.INDEX <- index.adjaceny(TH.ACT.LIST,code="N")
#TH.ACT.ADJM  <- make.adjaceny(TH.ACT.LIST,TH.ACT.INDEX)
#colnames(TH.ACT.ADJM) <- TH.ACT.LIST; rownames(TH.ACT.ADJM) <- TH.ACT.LIST

##Create graphs from adjaceny matrixes
#TG.ACT.NET <- graph.adjacency(TG.ACT.ADJM, mode="undirected")
#TH.ACT.NET <- graph.adjacency(TH.ACT.ADJM, mode="undirected")

##Calculate distances
#TG.ACT.DIST <- distances(TG.ACT.NET,v=TG.ACT.LIST,to=TG.ACT.LIST,mode="all"); is.na(TG.ACT.DIST) <- sapply(TG.ACT.DIST,is.infinite)
#TH.ACT.DIST <- distances(TH.ACT.NET,v=TH.ACT.LIST,to=TH.ACT.LIST,mode="all"); is.na(TH.ACT.DIST) <- sapply(TH.ACT.DIST,is.infinite)

##Save networks
#save(TG.ACT.ADJM, file="TG.ACT.ADJM.rda"); save(TG.ACT.NET,  file="TG.ACT.NET.rda"); save(TG.ACT.DIST, file="TG.ACT.DIST.rda")
#save(TH.ACT.ADJM, file="TH.ACT.ADJM.rda"); save(TH.ACT.NET,  file="TH.ACT.NET.rda"); save(TH.ACT.DIST, file="TH.ACT.DIST.rda")

###PE Model
##Get list of activators
PG.ACT.LIST <- unique(substr(DT.JOINT[PE.JOINT.CLASS != 'null',]$AAseq,2,5))
PH.ACT.LIST <- PG.ACT.LIST[-grep("Z",PG.ACT.LIST)]

##Make adjacency matrixes
#PG.ACT.INDEX <- index.adjaceny(PG.ACT.LIST,code="Z")
#PG.ACT.ADJM  <- make.adjaceny(PG.ACT.LIST,PG.ACT.INDEX)
#colnames(PG.ACT.ADJM) <- PG.ACT.LIST; rownames(PG.ACT.ADJM) <- PG.ACT.LIST
#PH.ACT.INDEX <- index.adjaceny(PH.ACT.LIST,code="N")
#PH.ACT.ADJM  <- make.adjaceny(PH.ACT.LIST,PH.ACT.INDEX)
#colnames(PH.ACT.ADJM) <- PH.ACT.LIST; rownames(PH.ACT.ADJM) <- PH.ACT.LIST

##Create graphs from adjaceny matrixes
#PG.ACT.NET <- graph.adjacency(PG.ACT.ADJM, mode="undirected")
#PH.ACT.NET <- graph.adjacency(PH.ACT.ADJM, mode="undirected")

##Calculate distances
#PG.ACT.DIST <- distances(PG.ACT.NET,v=PG.ACT.LIST,to=PG.ACT.LIST,mode="all"); is.na(PG.ACT.DIST) <- sapply(PG.ACT.DIST,is.infinite)
#PH.ACT.DIST <- distances(PH.ACT.NET,v=PH.ACT.LIST,to=PH.ACT.LIST,mode="all"); is.na(PH.ACT.DIST) <- sapply(PH.ACT.DIST,is.infinite)

##Save networks
#save(PG.ACT.ADJM, file="PG.ACT.ADJM.rda"); save(PG.ACT.NET,  file="PG.ACT.NET.rda"); save(PG.ACT.DIST, file="PG.ACT.DIST.rda")
#save(PH.ACT.ADJM, file="PH.ACT.ADJM.rda"); save(PH.ACT.NET,  file="PH.ACT.NET.rda"); save(PH.ACT.DIST, file="PH.ACT.DIST.rda")

###ME Model
##Get list of activators
MG.ACT.LIST <- unique(substr(DT.JOINT[ME.JOINT.CLASS != 'null',]$AAseq,2,5))
MH.ACT.LIST <- MG.ACT.LIST[-grep("Z",MG.ACT.LIST)]

##Make adjacency matrixes
#MG.ACT.INDEX <- index.adjaceny(MG.ACT.LIST,code="Z")
#MG.ACT.ADJM  <- make.adjaceny(MG.ACT.LIST,MG.ACT.INDEX)
#colnames(MG.ACT.ADJM) <- MG.ACT.LIST; rownames(MG.ACT.ADJM) <- MG.ACT.LIST
#MH.ACT.INDEX <- index.adjaceny(MH.ACT.LIST,code="N")
#MH.ACT.ADJM  <- make.adjaceny(MH.ACT.LIST,MH.ACT.INDEX)
#colnames(MH.ACT.ADJM) <- MH.ACT.LIST; rownames(MH.ACT.ADJM) <- MH.ACT.LIST

##Create graphs from adjaceny matrixes
#MG.ACT.NET <- graph.adjacency(MG.ACT.ADJM, mode="undirected")
#MH.ACT.NET <- graph.adjacency(MH.ACT.ADJM, mode="undirected")

##Calculate distances
#MG.ACT.DIST <- distances(MG.ACT.NET,v=MG.ACT.LIST,to=MG.ACT.LIST,mode="all"); is.na(MG.ACT.DIST) <- sapply(MG.ACT.DIST,is.infinite)
#MH.ACT.DIST <- distances(MH.ACT.NET,v=MH.ACT.LIST,to=MH.ACT.LIST,mode="all"); is.na(MH.ACT.DIST) <- sapply(MH.ACT.DIST,is.infinite)

##Save networks
#save(MG.ACT.ADJM, file="MG.ACT.ADJM.rda"); save(MG.ACT.NET,  file="MG.ACT.NET.rda"); save(MG.ACT.DIST, file="MG.ACT.DIST.rda")
#save(MH.ACT.ADJM, file="MH.ACT.ADJM.rda"); save(MH.ACT.NET,  file="MH.ACT.NET.rda"); save(MH.ACT.DIST, file="MH.ACT.DIST.rda")

###Null Model
##Get list of 'activators'
NG.ACT.LIST <- unique(substr(DT.JOINT$AAseq,2,5))
NH.ACT.LIST <- NG.ACT.LIST[-grep("Z",NG.ACT.LIST)]

##Make adjacency matrixes
#NG.ACT.INDEX <- index.adjaceny(NG.ACT.LIST,code="Z")
#NG.ACT.ADJM  <- make.adjaceny(NG.ACT.LIST,NG.ACT.INDEX)
#colnames(NG.ACT.ADJM) <- NG.ACT.LIST; rownames(NG.ACT.ADJM) <- NG.ACT.LIST
#NH.ACT.INDEX <- index.adjaceny(NH.ACT.LIST,code="N")
#NH.ACT.ADJM  <- make.adjaceny(NH.ACT.LIST,NH.ACT.INDEX)
#colnames(NH.ACT.ADJM) <- NH.ACT.LIST; rownames(NH.ACT.ADJM) <- NH.ACT.LIST

##Create graphs from adjaceny matrixes
#NG.ACT.NET <- graph.adjacency(NG.ACT.ADJM, mode="undirected")
#NH.ACT.NET <- graph.adjacency(NH.ACT.ADJM, mode="undirected")

##Subset list to consider only activators in other models
NG.ACT.LIST <- unique(substr(DT.JOINT[ME.JOINT.CLASS != 'null' | PE.JOINT.CLASS != 'null' | TE.JOINT.CLASS != 'null',]$AAseq,2,5))
NH.ACT.LIST <- NG.ACT.LIST[-grep("Z",NG.ACT.LIST)]

##Calculate distances
#NG.ACT.DIST <- distances(NG.ACT.NET,v=NG.ACT.LIST,to=NG.ACT.LIST,mode="all"); is.na(NG.ACT.DIST) <- sapply(NG.ACT.DIST,is.infinite)
#NH.ACT.DIST <- distances(NH.ACT.NET,v=NH.ACT.LIST,to=NH.ACT.LIST,mode="all"); is.na(NH.ACT.DIST) <- sapply(NH.ACT.DIST,is.infinite)

##Save networks
#save(NG.ACT.ADJM, file="NG.ACT.ADJM.rda"); save(NG.ACT.NET,  file="NG.ACT.NET.rda"); save(NG.ACT.DIST, file="NG.ACT.DIST.rda")
#save(NH.ACT.ADJM, file="NH.ACT.ADJM.rda"); save(NH.ACT.NET,  file="NH.ACT.NET.rda"); save(NH.ACT.DIST, file="NH.ACT.DIST.rda")

###Load networks
load("TG.ACT.ADJM.rda"); load("TG.ACT.NET.rda"); load("TG.ACT.DIST.rda")
load("TH.ACT.ADJM.rda"); load("TH.ACT.NET.rda"); load("TH.ACT.DIST.rda")
load("PG.ACT.ADJM.rda"); load("PG.ACT.NET.rda"); load("PG.ACT.DIST.rda")
load("PH.ACT.ADJM.rda"); load("PH.ACT.NET.rda"); load("PH.ACT.DIST.rda")
load("MG.ACT.ADJM.rda"); load("MG.ACT.NET.rda"); load("MG.ACT.DIST.rda")
load("MH.ACT.ADJM.rda"); load("MH.ACT.NET.rda"); load("MH.ACT.DIST.rda")
load("NG.ACT.ADJM.rda"); load("NG.ACT.NET.rda"); load("NG.ACT.DIST.rda")
load("NH.ACT.ADJM.rda"); load("NH.ACT.NET.rda"); load("NH.ACT.DIST.rda")

####ERE-specific to all SRE-specific pairs
###TE Model
##Get lists of ERE-specific and SRE-specific activators
TG.ERE.LIST <- unique(substr(DT.JOINT[TE.JOINT.CLASS == 'ERE-specific',]$AAseq,2,5))
TG.SRE.LIST <- unique(substr(DT.JOINT[TE.JOINT.CLASS == 'SRE-specific',]$AAseq,2,5))
TH.ERE.LIST <- TG.ERE.LIST[-grep("Z",TG.ERE.LIST)]
TH.SRE.LIST <- TG.SRE.LIST[-grep("Z",TG.SRE.LIST)]

##Calculate distances
#TG.ERE.SRE.DIST <- distances(TG.ACT.NET,v=TG.ERE.LIST,to=TG.SRE.LIST,mode="all"); is.na(TG.ERE.SRE.DIST) <- sapply(TG.ERE.SRE.DIST,is.infinite)
#TH.ERE.SRE.DIST <- distances(TH.ACT.NET,v=TH.ERE.LIST,to=TH.SRE.LIST,mode="all"); is.na(TH.ERE.SRE.DIST) <- sapply(TH.ERE.SRE.DIST,is.infinite)

##Save networks
#save(TG.ERE.SRE.DIST, file="TG.ERE.SRE.DIST.rda")
#save(TH.ERE.SRE.DIST, file="TH.ERE.SRE.DIST.rda")

###PE Model
##Get lists of ERE-specific and SRE-specific activators
PG.ERE.LIST <- unique(substr(DT.JOINT[PE.JOINT.CLASS == 'ERE-specific',]$AAseq,2,5))
PG.SRE.LIST <- unique(substr(DT.JOINT[PE.JOINT.CLASS == 'SRE-specific',]$AAseq,2,5))
PH.ERE.LIST <- PG.ERE.LIST[-grep("Z",PG.ERE.LIST)]
PH.SRE.LIST <- PG.SRE.LIST[-grep("Z",PG.SRE.LIST)]

##Calculate distances
#PG.ERE.SRE.DIST <- distances(PG.ACT.NET,v=PG.ERE.LIST,to=PG.SRE.LIST,mode="all"); is.na(PG.ERE.SRE.DIST) <- sapply(PG.ERE.SRE.DIST,is.infinite)
#PH.ERE.SRE.DIST <- distances(PH.ACT.NET,v=PH.ERE.LIST,to=PH.SRE.LIST,mode="all"); is.na(PH.ERE.SRE.DIST) <- sapply(PH.ERE.SRE.DIST,is.infinite)

##Save networks
#save(PG.ERE.SRE.DIST, file="PG.ERE.SRE.DIST.rda")
#save(PH.ERE.SRE.DIST, file="PH.ERE.SRE.DIST.rda")

###ME Model
##Get lists of ERE-specific and SRE-specific activators
MG.ERE.LIST <- unique(substr(DT.JOINT[ME.JOINT.CLASS == 'ERE-specific',]$AAseq,2,5))
MG.SRE.LIST <- unique(substr(DT.JOINT[ME.JOINT.CLASS == 'SRE-specific',]$AAseq,2,5))
MH.ERE.LIST <- MG.ERE.LIST[-grep("Z",MG.ERE.LIST)]
MH.SRE.LIST <- MG.SRE.LIST[-grep("Z",MG.SRE.LIST)]

##Calculate distances
#MG.ERE.SRE.DIST <- distances(MG.ACT.NET,v=MG.ERE.LIST,to=MG.SRE.LIST,mode="all"); is.na(MG.ERE.SRE.DIST) <- sapply(MG.ERE.SRE.DIST,is.infinite)
#MH.ERE.SRE.DIST <- distances(MH.ACT.NET,v=MH.ERE.LIST,to=MH.SRE.LIST,mode="all"); is.na(MH.ERE.SRE.DIST) <- sapply(MH.ERE.SRE.DIST,is.infinite)

##Save networks
#save(MG.ERE.SRE.DIST, file="MG.ERE.SRE.DIST.rda")
#save(MH.ERE.SRE.DIST, file="MH.ERE.SRE.DIST.rda")

###Null Model
##Get lists of ERE-specific and SRE-specific activators
NG.ERE.LIST <- unique(substr(DT.JOINT[ME.JOINT.CLASS == 'ERE-specific' ,]$AAseq,2,5))
NG.SRE.LIST <- unique(substr(DT.JOINT[ME.JOINT.CLASS == 'SRE-specific' ,]$AAseq,2,5))
NH.ERE.LIST <- NG.ERE.LIST[-grep("Z",NG.ERE.LIST)]
NH.SRE.LIST <- NG.SRE.LIST[-grep("Z",NG.SRE.LIST)]

##Calculate distances
#NG.ERE.SRE.DIST <- distances(NG.ACT.NET,v=NG.ERE.LIST,to=NG.SRE.LIST,mode="all"); is.na(NG.ERE.SRE.DIST) <- sapply(NG.ERE.SRE.DIST,is.infinite)
#NH.ERE.SRE.DIST <- distances(NH.ACT.NET,v=NH.ERE.LIST,to=NH.SRE.LIST,mode="all"); is.na(NH.ERE.SRE.DIST) <- sapply(NH.ERE.SRE.DIST,is.infinite)

##Save networks
#save(NG.ERE.SRE.DIST, file="NG.ERE.SRE.DIST.rda")
#save(NH.ERE.SRE.DIST, file="NH.ERE.SRE.DIST.rda")

###Load networks
load("TG.ERE.SRE.DIST.rda")
load("TH.ERE.SRE.DIST.rda")
load("PG.ERE.SRE.DIST.rda")
load("PH.ERE.SRE.DIST.rda")
load("MG.ERE.SRE.DIST.rda")
load("MH.ERE.SRE.DIST.rda")
load("NG.ERE.SRE.DIST.rda")
load("NH.ERE.SRE.DIST.rda")

##Lists of promiscuous genotypes in each network
TG.PRO.LIST <- TG.ACT.LIST[(!TG.ACT.LIST %in% c(TG.ERE.LIST,TG.SRE.LIST))]
TH.PRO.LIST <- TH.ACT.LIST[(!TH.ACT.LIST %in% c(TH.ERE.LIST,TH.SRE.LIST))]
PG.PRO.LIST <- PG.ACT.LIST[(!PG.ACT.LIST %in% c(PG.ERE.LIST,PG.SRE.LIST))]
PH.PRO.LIST <- PH.ACT.LIST[(!PH.ACT.LIST %in% c(PH.ERE.LIST,PH.SRE.LIST))]
MG.PRO.LIST <- MG.ACT.LIST[(!MG.ACT.LIST %in% c(MG.ERE.LIST,MG.SRE.LIST))]
MH.PRO.LIST <- MH.ACT.LIST[(!MH.ACT.LIST %in% c(MH.ERE.LIST,MH.SRE.LIST))]
NG.PRO.LIST <- NG.ACT.LIST[(!NG.ACT.LIST %in% c(NG.ERE.LIST,NG.SRE.LIST))]
NH.PRO.LIST <- NH.ACT.LIST[(!NH.ACT.LIST %in% c(NH.ERE.LIST,NH.SRE.LIST))]

##Find unconnected nodes
TG.NA <- which(colSums(TG.ACT.DIST,na.rm=TRUE) == 0)
PG.NA <- which(colSums(PG.ACT.DIST,na.rm=TRUE) == 0)
MG.NA <- which(colSums(MG.ACT.DIST,na.rm=TRUE) == 0)
NG.NA <- which(colSums(NG.ACT.DIST,na.rm=TRUE) == 0)
TH.NA <- which(colSums(TH.ACT.DIST,na.rm=TRUE) == 0)
PH.NA <- which(colSums(PH.ACT.DIST,na.rm=TRUE) == 0)
MH.NA <- which(colSums(MH.ACT.DIST,na.rm=TRUE) == 0)
NH.NA <- which(colSums(NH.ACT.DIST,na.rm=TRUE) == 0)

##Identify a common set of genotypes that are in the M,P, and T models for activators, ERE-specific, SRE-specific, and promiscuous genotypes
TG.ACT.INTERSECT <- which(TG.ACT.LIST %in% TH.ACT.LIST & TG.ACT.LIST %in% PH.ACT.LIST & TG.ACT.LIST %in% MH.ACT.LIST & !(TG.ACT.LIST %in% names(c(TG.NA,PG.NA,MG.NA,NG.NA,TH.NA,PH.NA,MH.NA,NH.NA))))
PG.ACT.INTERSECT <- which(PG.ACT.LIST %in% TH.ACT.LIST & PG.ACT.LIST %in% PH.ACT.LIST & PG.ACT.LIST %in% MH.ACT.LIST & !(PG.ACT.LIST %in% names(c(TG.NA,PG.NA,MG.NA,NG.NA,TH.NA,PH.NA,MH.NA,NH.NA))))
MG.ACT.INTERSECT <- which(MG.ACT.LIST %in% TH.ACT.LIST & MG.ACT.LIST %in% PH.ACT.LIST & MG.ACT.LIST %in% MH.ACT.LIST & !(MG.ACT.LIST %in% names(c(TG.NA,PG.NA,MG.NA,NG.NA,TH.NA,PH.NA,MH.NA,NH.NA))))
NG.ACT.INTERSECT <- which(NG.ACT.LIST %in% TH.ACT.LIST & NG.ACT.LIST %in% PH.ACT.LIST & NG.ACT.LIST %in% MH.ACT.LIST & !(NG.ACT.LIST %in% names(c(TG.NA,PG.NA,MG.NA,NG.NA,TH.NA,PH.NA,MH.NA,NH.NA))))
TH.ACT.INTERSECT <- which(TH.ACT.LIST %in% TH.ACT.LIST & TH.ACT.LIST %in% PH.ACT.LIST & TH.ACT.LIST %in% MH.ACT.LIST & !(TH.ACT.LIST %in% names(c(TG.NA,PG.NA,MG.NA,NG.NA,TH.NA,PH.NA,MH.NA,NH.NA))))
PH.ACT.INTERSECT <- which(PH.ACT.LIST %in% TH.ACT.LIST & PH.ACT.LIST %in% PH.ACT.LIST & PH.ACT.LIST %in% MH.ACT.LIST & !(PH.ACT.LIST %in% names(c(TG.NA,PG.NA,MG.NA,NG.NA,TH.NA,PH.NA,MH.NA,NH.NA))))
MH.ACT.INTERSECT <- which(MH.ACT.LIST %in% TH.ACT.LIST & MH.ACT.LIST %in% PH.ACT.LIST & MH.ACT.LIST %in% MH.ACT.LIST & !(MH.ACT.LIST %in% names(c(TG.NA,PG.NA,MG.NA,NG.NA,TH.NA,PH.NA,MH.NA,NH.NA))))
NH.ACT.INTERSECT <- which(NH.ACT.LIST %in% TH.ACT.LIST & NH.ACT.LIST %in% PH.ACT.LIST & NH.ACT.LIST %in% MH.ACT.LIST & !(NH.ACT.LIST %in% names(c(TG.NA,PG.NA,MG.NA,NG.NA,TH.NA,PH.NA,MH.NA,NH.NA))))

TG.ERE.INTERSECT <- which(TG.ERE.LIST %in% TH.ERE.LIST & TG.ERE.LIST %in% PH.ERE.LIST & TG.ERE.LIST %in% MH.ERE.LIST & !(TG.ERE.LIST %in% names(c(TG.NA,PG.NA,MG.NA,NG.NA,TH.NA,PH.NA,MH.NA,NH.NA))))
PG.ERE.INTERSECT <- which(PG.ERE.LIST %in% TH.ERE.LIST & PG.ERE.LIST %in% PH.ERE.LIST & PG.ERE.LIST %in% MH.ERE.LIST & !(PG.ERE.LIST %in% names(c(TG.NA,PG.NA,MG.NA,NG.NA,TH.NA,PH.NA,MH.NA,NH.NA))))
MG.ERE.INTERSECT <- which(MG.ERE.LIST %in% TH.ERE.LIST & MG.ERE.LIST %in% PH.ERE.LIST & MG.ERE.LIST %in% MH.ERE.LIST & !(MG.ERE.LIST %in% names(c(TG.NA,PG.NA,MG.NA,NG.NA,TH.NA,PH.NA,MH.NA,NH.NA))))
NG.ERE.INTERSECT <- which(NG.ERE.LIST %in% TH.ERE.LIST & NG.ERE.LIST %in% PH.ERE.LIST & NG.ERE.LIST %in% MH.ERE.LIST & !(NG.ERE.LIST %in% names(c(TG.NA,PG.NA,MG.NA,NG.NA,TH.NA,PH.NA,MH.NA,NH.NA))))
TH.ERE.INTERSECT <- which(TH.ERE.LIST %in% TH.ERE.LIST & TH.ERE.LIST %in% PH.ERE.LIST & TH.ERE.LIST %in% MH.ERE.LIST & !(TH.ERE.LIST %in% names(c(TG.NA,PG.NA,MG.NA,NG.NA,TH.NA,PH.NA,MH.NA,NH.NA))))
PH.ERE.INTERSECT <- which(PH.ERE.LIST %in% TH.ERE.LIST & PH.ERE.LIST %in% PH.ERE.LIST & PH.ERE.LIST %in% MH.ERE.LIST & !(PH.ERE.LIST %in% names(c(TG.NA,PG.NA,MG.NA,NG.NA,TH.NA,PH.NA,MH.NA,NH.NA))))
MH.ERE.INTERSECT <- which(MH.ERE.LIST %in% TH.ERE.LIST & MH.ERE.LIST %in% PH.ERE.LIST & MH.ERE.LIST %in% MH.ERE.LIST & !(MH.ERE.LIST %in% names(c(TG.NA,PG.NA,MG.NA,NG.NA,TH.NA,PH.NA,MH.NA,NH.NA))))
NH.ERE.INTERSECT <- which(NH.ERE.LIST %in% TH.ERE.LIST & NH.ERE.LIST %in% PH.ERE.LIST & NH.ERE.LIST %in% MH.ERE.LIST & !(NH.ERE.LIST %in% names(c(TG.NA,PG.NA,MG.NA,NG.NA,TH.NA,PH.NA,MH.NA,NH.NA))))

TG.SRE.INTERSECT <- which(TG.SRE.LIST %in% TH.SRE.LIST & TG.SRE.LIST %in% PH.SRE.LIST & TG.SRE.LIST %in% MH.SRE.LIST & !(TG.SRE.LIST %in% names(c(TG.NA,PG.NA,MG.NA,NG.NA,TH.NA,PH.NA,MH.NA,NH.NA))))
PG.SRE.INTERSECT <- which(PG.SRE.LIST %in% TH.SRE.LIST & PG.SRE.LIST %in% PH.SRE.LIST & PG.SRE.LIST %in% MH.SRE.LIST & !(PG.SRE.LIST %in% names(c(TG.NA,PG.NA,MG.NA,NG.NA,TH.NA,PH.NA,MH.NA,NH.NA))))
MG.SRE.INTERSECT <- which(MG.SRE.LIST %in% TH.SRE.LIST & MG.SRE.LIST %in% PH.SRE.LIST & MG.SRE.LIST %in% MH.SRE.LIST & !(MG.SRE.LIST %in% names(c(TG.NA,PG.NA,MG.NA,NG.NA,TH.NA,PH.NA,MH.NA,NH.NA))))
NG.SRE.INTERSECT <- which(NG.SRE.LIST %in% TH.SRE.LIST & NG.SRE.LIST %in% PH.SRE.LIST & NG.SRE.LIST %in% MH.SRE.LIST & !(NG.SRE.LIST %in% names(c(TG.NA,PG.NA,MG.NA,NG.NA,TH.NA,PH.NA,MH.NA,NH.NA))))
TH.SRE.INTERSECT <- which(TH.SRE.LIST %in% TH.SRE.LIST & TH.SRE.LIST %in% PH.SRE.LIST & TH.SRE.LIST %in% MH.SRE.LIST & !(TH.SRE.LIST %in% names(c(TG.NA,PG.NA,MG.NA,NG.NA,TH.NA,PH.NA,MH.NA,NH.NA))))
PH.SRE.INTERSECT <- which(PH.SRE.LIST %in% TH.SRE.LIST & PH.SRE.LIST %in% PH.SRE.LIST & PH.SRE.LIST %in% MH.SRE.LIST & !(PH.SRE.LIST %in% names(c(TG.NA,PG.NA,MG.NA,NG.NA,TH.NA,PH.NA,MH.NA,NH.NA))))
MH.SRE.INTERSECT <- which(MH.SRE.LIST %in% TH.SRE.LIST & MH.SRE.LIST %in% PH.SRE.LIST & MH.SRE.LIST %in% MH.SRE.LIST & !(MH.SRE.LIST %in% names(c(TG.NA,PG.NA,MG.NA,NG.NA,TH.NA,PH.NA,MH.NA,NH.NA))))
NH.SRE.INTERSECT <- which(NH.SRE.LIST %in% TH.SRE.LIST & NH.SRE.LIST %in% PH.SRE.LIST & NH.SRE.LIST %in% MH.SRE.LIST & !(NH.SRE.LIST %in% names(c(TG.NA,PG.NA,MG.NA,NG.NA,TH.NA,PH.NA,MH.NA,NH.NA))))

TG.PRO.INTERSECT <- which(TG.PRO.LIST %in% TH.PRO.LIST & TG.PRO.LIST %in% PH.PRO.LIST & TG.PRO.LIST %in% MH.PRO.LIST & !(TG.PRO.LIST %in% names(c(TG.NA,PG.NA,MG.NA,NG.NA,TH.NA,PH.NA,MH.NA,NH.NA))))
PG.PRO.INTERSECT <- which(PG.PRO.LIST %in% TH.PRO.LIST & PG.PRO.LIST %in% PH.PRO.LIST & PG.PRO.LIST %in% MH.PRO.LIST & !(PG.PRO.LIST %in% names(c(TG.NA,PG.NA,MG.NA,NG.NA,TH.NA,PH.NA,MH.NA,NH.NA))))
MG.PRO.INTERSECT <- which(MG.PRO.LIST %in% TH.PRO.LIST & MG.PRO.LIST %in% PH.PRO.LIST & MG.PRO.LIST %in% MH.PRO.LIST & !(MG.PRO.LIST %in% names(c(TG.NA,PG.NA,MG.NA,NG.NA,TH.NA,PH.NA,MH.NA,NH.NA))))
NG.PRO.INTERSECT <- which(NG.PRO.LIST %in% TH.PRO.LIST & NG.PRO.LIST %in% PH.PRO.LIST & NG.PRO.LIST %in% MH.PRO.LIST & !(NG.PRO.LIST %in% names(c(TG.NA,PG.NA,MG.NA,NG.NA,TH.NA,PH.NA,MH.NA,NH.NA))))
TH.PRO.INTERSECT <- which(TH.PRO.LIST %in% TH.PRO.LIST & TH.PRO.LIST %in% PH.PRO.LIST & TH.PRO.LIST %in% MH.PRO.LIST & !(TH.PRO.LIST %in% names(c(TG.NA,PG.NA,MG.NA,NG.NA,TH.NA,PH.NA,MH.NA,NH.NA))))
PH.PRO.INTERSECT <- which(PH.PRO.LIST %in% TH.PRO.LIST & PH.PRO.LIST %in% PH.PRO.LIST & PH.PRO.LIST %in% MH.PRO.LIST & !(PH.PRO.LIST %in% names(c(TG.NA,PG.NA,MG.NA,NG.NA,TH.NA,PH.NA,MH.NA,NH.NA))))
MH.PRO.INTERSECT <- which(MH.PRO.LIST %in% TH.PRO.LIST & MH.PRO.LIST %in% PH.PRO.LIST & MH.PRO.LIST %in% MH.PRO.LIST & !(MH.PRO.LIST %in% names(c(TG.NA,PG.NA,MG.NA,NG.NA,TH.NA,PH.NA,MH.NA,NH.NA))))
NH.PRO.INTERSECT <- which(NH.PRO.LIST %in% TH.PRO.LIST & NH.PRO.LIST %in% PH.PRO.LIST & NH.PRO.LIST %in% MH.PRO.LIST & !(NH.PRO.LIST %in% names(c(TG.NA,PG.NA,MG.NA,NG.NA,TH.NA,PH.NA,MH.NA,NH.NA))))


####ERE-specific to closest SRE-specific pairs
##Calculate distances
#TG.ERE.SRE.MIN.DIST <- apply(TG.ERE.SRE.DIST, 1, min, na.rm=TRUE)
#PG.ERE.SRE.MIN.DIST <- apply(PG.ERE.SRE.DIST, 1, min, na.rm=TRUE)
#MG.ERE.SRE.MIN.DIST <- apply(MG.ERE.SRE.DIST, 1, min, na.rm=TRUE)
#NG.ERE.SRE.MIN.DIST <- apply(NG.ERE.SRE.DIST, 1, min, na.rm=TRUE)
#TH.ERE.SRE.MIN.DIST <- apply(TH.ERE.SRE.DIST, 1, min, na.rm=TRUE)
#PH.ERE.SRE.MIN.DIST <- apply(PH.ERE.SRE.DIST, 1, min, na.rm=TRUE)
#MH.ERE.SRE.MIN.DIST <- apply(MH.ERE.SRE.DIST, 1, min, na.rm=TRUE)
#NH.ERE.SRE.MIN.DIST <- apply(NH.ERE.SRE.DIST, 1, min, na.rm=TRUE)

#save(TG.ERE.SRE.MIN.DIST, file="TG.ERE.SRE.MIN.DIST.rda")
#save(PG.ERE.SRE.MIN.DIST, file="PG.ERE.SRE.MIN.DIST.rda")
#save(MG.ERE.SRE.MIN.DIST, file="MG.ERE.SRE.MIN.DIST.rda")
#save(NG.ERE.SRE.MIN.DIST, file="NG.ERE.SRE.MIN.DIST.rda")
#save(TH.ERE.SRE.MIN.DIST, file="TH.ERE.SRE.MIN.DIST.rda")
#save(PH.ERE.SRE.MIN.DIST, file="PH.ERE.SRE.MIN.DIST.rda")
#save(MH.ERE.SRE.MIN.DIST, file="MH.ERE.SRE.MIN.DIST.rda")
#save(NH.ERE.SRE.MIN.DIST, file="NH.ERE.SRE.MIN.DIST.rda")

###Load distances
load("TG.ERE.SRE.MIN.DIST.rda")
load("PG.ERE.SRE.MIN.DIST.rda")
load("MG.ERE.SRE.MIN.DIST.rda")
load("NG.ERE.SRE.MIN.DIST.rda")
load("TH.ERE.SRE.MIN.DIST.rda")
load("PH.ERE.SRE.MIN.DIST.rda")
load("MH.ERE.SRE.MIN.DIST.rda")
load("NH.ERE.SRE.MIN.DIST.rda")

####Plot results
##pdf("DISTANCE.DISTRIBUTION.pdf")
#par(mfrow=c(4,2))
###All Activators
#barplot(table(factor(TG.ACT.DIST[upper.tri(TG.ACT.DIST)], levels=1:13), useNA="always"), ylim=c(0,300000), ylab="", xlab="TG",main="All activators")
#barplot(table(factor(TH.ACT.DIST[upper.tri(TH.ACT.DIST)], levels=1:7 ), useNA="always"), ylim=c(0,300000), ylab="", xlab="TH",main="All activators")
#barplot(table(factor(PG.ACT.DIST[upper.tri(PG.ACT.DIST)], levels=1:13), useNA="always"), ylim=c(0,300000), ylab="", xlab="PG")
#barplot(table(factor(PH.ACT.DIST[upper.tri(PH.ACT.DIST)], levels=1:7 ), useNA="always"), ylim=c(0,300000), ylab="", xlab="PH")
#barplot(table(factor(MG.ACT.DIST[upper.tri(MG.ACT.DIST)], levels=1:13), useNA="always"), ylim=c(0,300000), ylab="", xlab="MG")
#barplot(table(factor(MH.ACT.DIST[upper.tri(MH.ACT.DIST)], levels=1:7 ), useNA="always"), ylim=c(0,300000), ylab="", xlab="MH")
#barplot(table(factor(NG.ACT.DIST[upper.tri(NG.ACT.DIST)], levels=1:13), useNA="always"), ylim=c(0,500000), ylab="", xlab="NG")
#barplot(table(factor(NH.ACT.DIST[upper.tri(NH.ACT.DIST)], levels=1:7 ), useNA="always"), ylim=c(0,500000), ylab="", xlab="NH")
#
###ERE-specific to SRE-specific activators
#barplot(table(factor(TG.ERE.SRE.DIST[upper.tri(TG.ERE.SRE.DIST)], levels=1:13), useNA="always"), ylim=c(0,100000), ylab="", xlab="TG",main="All ERE to SRE")
#barplot(table(factor(TH.ERE.SRE.DIST[upper.tri(TH.ERE.SRE.DIST)], levels=1:7 ), useNA="always"), ylim=c(0,100000), ylab="", xlab="TH",main="All ERE to SRE")
#barplot(table(factor(PG.ERE.SRE.DIST[upper.tri(PG.ERE.SRE.DIST)], levels=1:13), useNA="always"), ylim=c(0,100000), ylab="", xlab="PG")
#barplot(table(factor(PH.ERE.SRE.DIST[upper.tri(PH.ERE.SRE.DIST)], levels=1:7 ), useNA="always"), ylim=c(0,100000), ylab="", xlab="PH")
#barplot(table(factor(MG.ERE.SRE.DIST[upper.tri(MG.ERE.SRE.DIST)], levels=1:13), useNA="always"), ylim=c(0,100000), ylab="", xlab="MG")
#barplot(table(factor(MH.ERE.SRE.DIST[upper.tri(MH.ERE.SRE.DIST)], levels=1:7 ), useNA="always"), ylim=c(0,100000), ylab="", xlab="MH")
#barplot(table(factor(NG.ERE.SRE.DIST[upper.tri(NG.ERE.SRE.DIST)], levels=1:13), useNA="always"), ylim=c(0,50000 ), ylab="", xlab="NG")
#barplot(table(factor(NH.ERE.SRE.DIST[upper.tri(NH.ERE.SRE.DIST)], levels=1:7 ), useNA="always"), ylim=c(0,50000 ), ylab="", xlab="NH")
#
###ERE-specific to closest SRE-specific
#barplot(table(factor(TG.ERE.SRE.MIN.DIST, levels=1:4), useNA="always"), ylim=c(0,150), ylab="", xlab="TG",main="ERE to closest SRE")
#barplot(table(factor(TH.ERE.SRE.MIN.DIST, levels=1:4), useNA="always"), ylim=c(0,150), ylab="", xlab="TH",main="ERE to closest SRE")
#barplot(table(factor(PG.ERE.SRE.MIN.DIST, levels=1:4), useNA="always"), ylim=c(0,150), ylab="", xlab="PG")
#barplot(table(factor(PH.ERE.SRE.MIN.DIST, levels=1:4), useNA="always"), ylim=c(0,150), ylab="", xlab="PH")
#barplot(table(factor(MG.ERE.SRE.MIN.DIST, levels=1:4), useNA="always"), ylim=c(0,150), ylab="", xlab="MG")
#barplot(table(factor(MH.ERE.SRE.MIN.DIST, levels=1:4), useNA="always"), ylim=c(0,150), ylab="", xlab="MH")
#barplot(table(factor(NG.ERE.SRE.MIN.DIST, levels=1:4), useNA="always"), ylim=c(0,150), ylab="", xlab="NG")
#barplot(table(factor(NH.ERE.SRE.MIN.DIST, levels=1:4), useNA="always"), ylim=c(0,150), ylab="", xlab="NH")
##dev.off()

###Average distances
##Collect averages
ALL.DIST.AVG     <- c(mean(TG.ERE.SRE.DIST,na.rm=TRUE),mean(PG.ERE.SRE.DIST,na.rm=TRUE),mean(MG.ERE.SRE.DIST,na.rm=TRUE),mean(NG.ERE.SRE.DIST,na.rm=TRUE),
                      mean(TH.ERE.SRE.DIST,na.rm=TRUE),mean(PH.ERE.SRE.DIST,na.rm=TRUE),mean(MH.ERE.SRE.DIST,na.rm=TRUE),mean(NH.ERE.SRE.DIST,na.rm=TRUE))
ALL.MIN.DIST.AVG <- c(mean(TG.ERE.SRE.MIN.DIST[is.finite(TG.ERE.SRE.MIN.DIST)]),mean(PG.ERE.SRE.MIN.DIST[is.finite(PG.ERE.SRE.MIN.DIST)]),mean(MG.ERE.SRE.MIN.DIST[is.finite(MG.ERE.SRE.MIN.DIST)]),mean(NG.ERE.SRE.MIN.DIST[is.finite(NG.ERE.SRE.MIN.DIST)]),
                      mean(TH.ERE.SRE.MIN.DIST[is.finite(TH.ERE.SRE.MIN.DIST)]),mean(PH.ERE.SRE.MIN.DIST[is.finite(PH.ERE.SRE.MIN.DIST)]),mean(MH.ERE.SRE.MIN.DIST[is.finite(MH.ERE.SRE.MIN.DIST)]),mean(NH.ERE.SRE.MIN.DIST[is.finite(NH.ERE.SRE.MIN.DIST)]))
INT.DIST.AVG     <- c(mean(TG.ERE.SRE.DIST[TG.ERE.INTERSECT,TG.SRE.INTERSECT],na.rm=TRUE),mean(PG.ERE.SRE.DIST[PG.ERE.INTERSECT,PG.SRE.INTERSECT],na.rm=TRUE),mean(MG.ERE.SRE.DIST[MG.ERE.INTERSECT,MG.SRE.INTERSECT],na.rm=TRUE),mean(NG.ERE.SRE.DIST[NG.ERE.INTERSECT,NG.SRE.INTERSECT],na.rm=TRUE),
                      mean(TH.ERE.SRE.DIST[TH.ERE.INTERSECT,TH.SRE.INTERSECT],na.rm=TRUE),mean(PH.ERE.SRE.DIST[PH.ERE.INTERSECT,PH.SRE.INTERSECT],na.rm=TRUE),mean(MH.ERE.SRE.DIST[MH.ERE.INTERSECT,MH.SRE.INTERSECT],na.rm=TRUE),mean(NH.ERE.SRE.DIST[NH.ERE.INTERSECT,NH.SRE.INTERSECT],na.rm=TRUE))
INT.MIN.DIST.AVG <- c(mean(TG.ERE.SRE.MIN.DIST[TG.ERE.INTERSECT],na.rm=TRUE),mean(PG.ERE.SRE.MIN.DIST[PG.ERE.INTERSECT],na.rm=TRUE),mean(MG.ERE.SRE.MIN.DIST[MG.ERE.INTERSECT],na.rm=TRUE),mean(NG.ERE.SRE.MIN.DIST[NG.ERE.INTERSECT],na.rm=TRUE),
                      mean(TH.ERE.SRE.MIN.DIST[TH.ERE.INTERSECT],na.rm=TRUE),mean(PH.ERE.SRE.MIN.DIST[PH.ERE.INTERSECT],na.rm=TRUE),mean(MH.ERE.SRE.MIN.DIST[MH.ERE.INTERSECT],na.rm=TRUE),mean(NH.ERE.SRE.MIN.DIST[NH.ERE.INTERSECT],na.rm=TRUE))

##pdf("DISTANCE.pdf")
#par(mar=c(3,3,3,3))
#
##Distance from ERE to SRE for all genotypes in network
#par(mfrow=c(1,1))
#vioplot(c(TG.ERE.SRE.DIST),c(PG.ERE.SRE.DIST),c(MG.ERE.SRE.DIST),c(NG.ERE.SRE.DIST),names=c("3G","2G","1G","NG"),cex.axis=0.75,cex.main=0.5,col="darkred", rectCol="white",lineCol="white",pchMed=19,colMed="black",main="ERE to SRE Distance: All Genotypes")
#vioplot(c(TH.ERE.SRE.DIST),c(PH.ERE.SRE.DIST),c(MH.ERE.SRE.DIST),c(NH.ERE.SRE.DIST),names=c("3H","2H","1H","NH"),cex.axis=0.75,cex.main=0.5,col="skyblue2",rectCol="white",lineCol="white",pchMed=19,colMed="black",main="ERE to SRE Distance: All Genotypes")
#
#barplot(ALL.DIST.AVG[1:4],names=c("3G","2G","1G","NG"),cex.axis=0.75,cex.main=0.5,col="darkred", main="ERE to SRE Distance: All Genotypes")
#barplot(ALL.DIST.AVG[5:8],names=c("3H","2H","1H","NH"),cex.axis=0.75,cex.main=0.5,col="skyblue2",main="ERE to SRE Distance: All Genotypes")
#
##Distance from ERE to SRE for common genotypes among networks
#par(mfrow=c(1,1))
#vioplot(c(TG.ERE.SRE.DIST[TG.ERE.INTERSECT,TG.SRE.INTERSECT]),c(PG.ERE.SRE.DIST[PG.ERE.INTERSECT,PG.SRE.INTERSECT]),c(MG.ERE.SRE.DIST[MG.ERE.INTERSECT,MG.SRE.INTERSECT]),c(NG.ERE.SRE.DIST[NG.ERE.INTERSECT,NG.SRE.INTERSECT]),names=c("3G","2G","1G","NG"),cex.axis=0.75,cex.main=0.5,col="darkred", rectCol="white",lineCol="white",pchMed=19,colMed="black",main="ERE to SRE Distance: Common Genotypes")
#vioplot(c(TH.ERE.SRE.DIST[TH.ERE.INTERSECT,TH.SRE.INTERSECT]),c(PH.ERE.SRE.DIST[PH.ERE.INTERSECT,PH.SRE.INTERSECT]),c(MH.ERE.SRE.DIST[MH.ERE.INTERSECT,MH.SRE.INTERSECT]),c(NH.ERE.SRE.DIST[NH.ERE.INTERSECT,NH.SRE.INTERSECT]),names=c("3H","2H","1H","NH"),cex.axis=0.75,cex.main=0.5,col="skyblue2",rectCol="white",lineCol="white",pchMed=19,colMed="black",main="ERE to SRE Distance: Common Genotypes")
#
#barplot(INT.DIST.AVG[1:4],names=c("3G","2G","1G","NG"),cex.axis=0.75,cex.main=0.5,col="darkred", main="ERE to SRE Distance: Common Genotypes")
#barplot(INT.DIST.AVG[5:8],names=c("3H","2H","1H","NH"),cex.axis=0.75,cex.main=0.5,col="skyblue2",main="ERE to SRE Distance: Common Genotypes")

###Permutation tests
#par(mfrow=c(3,2))
##Average distance for all ERE to SRE 
#TG.PG.ALL.DIST.PERMUTATION.TEST <- permutation.test.matrix(TG.ERE.SRE.DIST,PG.ERE.SRE.DIST,10000)
#  hist(TG.PG.ALL.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.5,.5),xlab="Average Difference path distance",main="All TG vs PG"); abline(v=mean(TG.ERE.SRE.DIST,na.rm=TRUE) - mean(PG.ERE.SRE.DIST,na.rm=TRUE),col="red")
#  sum(mean(TG.ERE.SRE.DIST,na.rm=TRUE) - mean(PG.ERE.SRE.DIST,na.rm=TRUE) < TG.PG.ALL.DIST.PERMUTATION.TEST)/length(TG.PG.ALL.DIST.PERMUTATION.TEST)
#PG.MG.ALL.DIST.PERMUTATION.TEST <- permutation.test.matrix(PG.ERE.SRE.DIST,MG.ERE.SRE.DIST,10000)
#  hist(PG.MG.ALL.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.5,.5),xlab="Average Difference path distance",main="All PG vs MG"); abline(v=mean(PG.ERE.SRE.DIST,na.rm=TRUE) - mean(MG.ERE.SRE.DIST,na.rm=TRUE),col="red")
#  sum(mean(PG.ERE.SRE.DIST,na.rm=TRUE) - mean(MG.ERE.SRE.DIST,na.rm=TRUE) < PG.MG.ALL.DIST.PERMUTATION.TEST)/length(PG.MG.ALL.DIST.PERMUTATION.TEST)
#MG.NG.ALL.DIST.PERMUTATION.TEST <- permutation.test.matrix(MG.ERE.SRE.DIST,NG.ERE.SRE.DIST,10000)
#  hist(MG.NG.ALL.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.5,.5),xlab="Average Difference path distance",main="All MG vs NG"); abline(v=mean(MG.ERE.SRE.DIST,na.rm=TRUE) - mean(NG.ERE.SRE.DIST,na.rm=TRUE),col="red")
#  sum(mean(MG.ERE.SRE.DIST,na.rm=TRUE) - mean(NG.ERE.SRE.DIST,na.rm=TRUE) < MG.NG.ALL.DIST.PERMUTATION.TEST)/length(MG.NG.ALL.DIST.PERMUTATION.TEST)
#TG.MG.ALL.DIST.PERMUTATION.TEST <- permutation.test.matrix(TG.ERE.SRE.DIST,MG.ERE.SRE.DIST,10000)
#  hist(TG.MG.ALL.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.5,.5),xlab="Average Difference path distance",main="All TG vs MG"); abline(v=mean(TG.ERE.SRE.DIST,na.rm=TRUE) - mean(MG.ERE.SRE.DIST,na.rm=TRUE),col="red")
#  sum(mean(TG.ERE.SRE.DIST,na.rm=TRUE) - mean(MG.ERE.SRE.DIST,na.rm=TRUE) < TG.MG.ALL.DIST.PERMUTATION.TEST)/length(TG.MG.ALL.DIST.PERMUTATION.TEST)
#TG.NG.ALL.DIST.PERMUTATION.TEST <- permutation.test.matrix(TG.ERE.SRE.DIST,NG.ERE.SRE.DIST,10000)
#  hist(TG.NG.ALL.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.8,.8),xlab="Average Difference path distance",main="All TG vs NG"); abline(v=mean(TG.ERE.SRE.DIST,na.rm=TRUE) - mean(NG.ERE.SRE.DIST,na.rm=TRUE),col="red")
#  sum(mean(TG.ERE.SRE.DIST,na.rm=TRUE) - mean(NG.ERE.SRE.DIST,na.rm=TRUE) < TG.NG.ALL.DIST.PERMUTATION.TEST)/length(TG.NG.ALL.DIST.PERMUTATION.TEST)
#PG.NG.ALL.DIST.PERMUTATION.TEST <- permutation.test.matrix(PG.ERE.SRE.DIST,NG.ERE.SRE.DIST,10000)
#  hist(PG.NG.ALL.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.8,.8),xlab="Average Difference path distance",main="All PG vs NG"); abline(v=mean(PG.ERE.SRE.DIST,na.rm=TRUE) - mean(NG.ERE.SRE.DIST,na.rm=TRUE),col="red")
#  sum(mean(PG.ERE.SRE.DIST,na.rm=TRUE) - mean(NG.ERE.SRE.DIST,na.rm=TRUE) < PG.NG.ALL.DIST.PERMUTATION.TEST)/length(PG.NG.ALL.DIST.PERMUTATION.TEST)
#TH.PH.ALL.DIST.PERMUTATION.TEST <- permutation.test.matrix(TH.ERE.SRE.DIST,PH.ERE.SRE.DIST,10000)
#  hist(TH.PH.ALL.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.4,.4),xlab="Average Difference path distance",main="All TH vs PH"); abline(v=mean(TH.ERE.SRE.DIST,na.rm=TRUE) - mean(PH.ERE.SRE.DIST,na.rm=TRUE),col="red")
#  sum(mean(TH.ERE.SRE.DIST,na.rm=TRUE) - mean(PH.ERE.SRE.DIST,na.rm=TRUE) < TH.PH.ALL.DIST.PERMUTATION.TEST)/length(TH.PH.ALL.DIST.PERMUTATION.TEST)
#PH.MH.ALL.DIST.PERMUTATION.TEST <- permutation.test.matrix(PH.ERE.SRE.DIST,MH.ERE.SRE.DIST,10000)
#  hist(PH.MH.ALL.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.4,.4),xlab="Average Difference path distance",main="All PH vs MH"); abline(v=mean(PH.ERE.SRE.DIST,na.rm=TRUE) - mean(MH.ERE.SRE.DIST,na.rm=TRUE),col="red")
#  sum(mean(PH.ERE.SRE.DIST,na.rm=TRUE) - mean(MH.ERE.SRE.DIST,na.rm=TRUE) < PH.MH.ALL.DIST.PERMUTATION.TEST)/length(PH.MH.ALL.DIST.PERMUTATION.TEST)
#MH.NH.ALL.DIST.PERMUTATION.TEST <- permutation.test.matrix(MH.ERE.SRE.DIST,NH.ERE.SRE.DIST,10000)
#  hist(MH.NH.ALL.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.4,.4),xlab="Average Difference path distance",main="All MH vs NH"); abline(v=mean(MH.ERE.SRE.DIST,na.rm=TRUE) - mean(NH.ERE.SRE.DIST,na.rm=TRUE),col="red")
#  sum(mean(MH.ERE.SRE.DIST,na.rm=TRUE) - mean(NH.ERE.SRE.DIST,na.rm=TRUE) < MH.NH.ALL.DIST.PERMUTATION.TEST)/length(MH.NH.ALL.DIST.PERMUTATION.TEST)
#TH.MH.ALL.DIST.PERMUTATION.TEST <- permutation.test.matrix(TH.ERE.SRE.DIST,MH.ERE.SRE.DIST,10000)
#  hist(TH.MH.ALL.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.4,.4),xlab="Average Difference path distance",main="All TH vs MH"); abline(v=mean(TH.ERE.SRE.DIST,na.rm=TRUE) - mean(MH.ERE.SRE.DIST,na.rm=TRUE),col="red")
#  sum(mean(TH.ERE.SRE.DIST,na.rm=TRUE) - mean(MH.ERE.SRE.DIST,na.rm=TRUE) < TH.MH.ALL.DIST.PERMUTATION.TEST)/length(TH.MH.ALL.DIST.PERMUTATION.TEST)
#TH.NH.ALL.DIST.PERMUTATION.TEST <- permutation.test.matrix(TH.ERE.SRE.DIST,NH.ERE.SRE.DIST,10000)
#  hist(TH.NH.ALL.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.4,.4),xlab="Average Difference path distance",main="All TH vs NH"); abline(v=mean(TH.ERE.SRE.DIST,na.rm=TRUE) - mean(NH.ERE.SRE.DIST,na.rm=TRUE),col="red")
#  sum(mean(TH.ERE.SRE.DIST,na.rm=TRUE) - mean(NH.ERE.SRE.DIST,na.rm=TRUE) < TH.NH.ALL.DIST.PERMUTATION.TEST)/length(TH.NH.ALL.DIST.PERMUTATION.TEST)
#PH.NH.ALL.DIST.PERMUTATION.TEST <- permutation.test.matrix(PH.ERE.SRE.DIST,NH.ERE.SRE.DIST,10000)
#  hist(PH.NH.ALL.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.4,.4),xlab="Average Difference path distance",main="All PH vs NH"); abline(v=mean(PH.ERE.SRE.DIST,na.rm=TRUE) - mean(NH.ERE.SRE.DIST,na.rm=TRUE),col="red")
#  sum(mean(PH.ERE.SRE.DIST,na.rm=TRUE) - mean(NH.ERE.SRE.DIST,na.rm=TRUE) < PH.NH.ALL.DIST.PERMUTATION.TEST)/length(PH.NH.ALL.DIST.PERMUTATION.TEST)
#
##Average distance for common ERE to SRE 
#TG.PG.INT.DIST.PERMUTATION.TEST <- permutation.test.matrix(TG.ERE.SRE.DIST[TG.ERE.INTERSECT,TG.SRE.INTERSECT],PG.ERE.SRE.DIST[PG.ERE.INTERSECT,PG.SRE.INTERSECT],10000)
#  hist(TG.PG.INT.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.4,.4),xlab="Average Difference path distance",main="Common TG vs PG"); abline(v=mean(TG.ERE.SRE.DIST[TG.ERE.INTERSECT,TG.SRE.INTERSECT],na.rm=TRUE) - mean(PG.ERE.SRE.DIST[PG.ERE.INTERSECT,PG.SRE.INTERSECT],na.rm=TRUE),col="red")
#  sum(mean(TG.ERE.SRE.DIST[TG.ERE.INTERSECT,TG.SRE.INTERSECT],na.rm=TRUE) - mean(PG.ERE.SRE.DIST[PG.ERE.INTERSECT,PG.SRE.INTERSECT],na.rm=TRUE) < TG.PG.INT.DIST.PERMUTATION.TEST)/length(TG.PG.INT.DIST.PERMUTATION.TEST)
#PG.MG.INT.DIST.PERMUTATION.TEST <- permutation.test.matrix(PG.ERE.SRE.DIST[PG.ERE.INTERSECT,PG.SRE.INTERSECT],MG.ERE.SRE.DIST[MG.ERE.INTERSECT,MG.SRE.INTERSECT],10000)
#  hist(PG.MG.INT.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.4,.4),xlab="Average Difference path distance",main="Common PG vs MG"); abline(v=mean(PG.ERE.SRE.DIST[PG.ERE.INTERSECT,PG.SRE.INTERSECT],na.rm=TRUE) - mean(MG.ERE.SRE.DIST[MG.ERE.INTERSECT,MG.SRE.INTERSECT],na.rm=TRUE),col="red")
#  sum(mean(PG.ERE.SRE.DIST[PG.ERE.INTERSECT,PG.SRE.INTERSECT],na.rm=TRUE) - mean(MG.ERE.SRE.DIST[MG.ERE.INTERSECT,MG.SRE.INTERSECT],na.rm=TRUE) < PG.MG.INT.DIST.PERMUTATION.TEST)/length(PG.MG.INT.DIST.PERMUTATION.TEST)
#MG.NG.INT.DIST.PERMUTATION.TEST <- permutation.test.matrix(MG.ERE.SRE.DIST[MG.ERE.INTERSECT,MG.SRE.INTERSECT],NG.ERE.SRE.DIST[NG.ERE.INTERSECT,NG.SRE.INTERSECT],10000)
#  hist(MG.NG.INT.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.4,.4),xlab="Average Difference path distance",main="Common MG vs NG"); abline(v=mean(MG.ERE.SRE.DIST[MG.ERE.INTERSECT,MG.SRE.INTERSECT],na.rm=TRUE) - mean(NG.ERE.SRE.DIST[NG.ERE.INTERSECT,NG.SRE.INTERSECT],na.rm=TRUE),col="red")
#  sum(mean(MG.ERE.SRE.DIST[MG.ERE.INTERSECT,MG.SRE.INTERSECT],na.rm=TRUE) - mean(NG.ERE.SRE.DIST[NG.ERE.INTERSECT,NG.SRE.INTERSECT],na.rm=TRUE) < MG.NG.INT.DIST.PERMUTATION.TEST)/length(MG.NG.INT.DIST.PERMUTATION.TEST)
#TG.MG.INT.DIST.PERMUTATION.TEST <- permutation.test.matrix(TG.ERE.SRE.DIST[TG.ERE.INTERSECT,TG.SRE.INTERSECT],MG.ERE.SRE.DIST[MG.ERE.INTERSECT,MG.SRE.INTERSECT],10000)
#  hist(TG.MG.INT.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.4,.4),xlab="Average Difference path distance",main="Common TG vs MG"); abline(v=mean(TG.ERE.SRE.DIST[TG.ERE.INTERSECT,TG.SRE.INTERSECT],na.rm=TRUE) - mean(MG.ERE.SRE.DIST[MG.ERE.INTERSECT,MG.SRE.INTERSECT],na.rm=TRUE),col="red")
#  sum(mean(TG.ERE.SRE.DIST[TG.ERE.INTERSECT,TG.SRE.INTERSECT],na.rm=TRUE) - mean(MG.ERE.SRE.DIST[MG.ERE.INTERSECT,MG.SRE.INTERSECT],na.rm=TRUE) < TG.MG.INT.DIST.PERMUTATION.TEST)/length(TG.MG.INT.DIST.PERMUTATION.TEST)
#TG.NG.INT.DIST.PERMUTATION.TEST <- permutation.test.matrix(TG.ERE.SRE.DIST[TG.ERE.INTERSECT,TG.SRE.INTERSECT],NG.ERE.SRE.DIST[NG.ERE.INTERSECT,NG.SRE.INTERSECT],10000)
#  hist(TG.NG.INT.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.4,.4),xlab="Average Difference path distance",main="Common TG vs NG"); abline(v=mean(TG.ERE.SRE.DIST[TG.ERE.INTERSECT,TG.SRE.INTERSECT],na.rm=TRUE) - mean(NG.ERE.SRE.DIST[NG.ERE.INTERSECT,NG.SRE.INTERSECT],na.rm=TRUE),col="red")
#  sum(mean(TG.ERE.SRE.DIST[TG.ERE.INTERSECT,TG.SRE.INTERSECT],na.rm=TRUE) - mean(NG.ERE.SRE.DIST[NG.ERE.INTERSECT,NG.SRE.INTERSECT],na.rm=TRUE) < TG.NG.INT.DIST.PERMUTATION.TEST)/length(TG.NG.INT.DIST.PERMUTATION.TEST)
#PG.NG.INT.DIST.PERMUTATION.TEST <- permutation.test.matrix(PG.ERE.SRE.DIST[PG.ERE.INTERSECT,PG.SRE.INTERSECT],NG.ERE.SRE.DIST[NG.ERE.INTERSECT,NG.SRE.INTERSECT],10000)
#  hist(PG.NG.INT.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.4,.4),xlab="Average Difference path distance",main="Common PG vs NG"); abline(v=mean(PG.ERE.SRE.DIST[PG.ERE.INTERSECT,PG.SRE.INTERSECT],na.rm=TRUE) - mean(NG.ERE.SRE.DIST[NG.ERE.INTERSECT,NG.SRE.INTERSECT],na.rm=TRUE),col="red")
#  sum(mean(PG.ERE.SRE.DIST[PG.ERE.INTERSECT,PG.SRE.INTERSECT],na.rm=TRUE) - mean(NG.ERE.SRE.DIST[NG.ERE.INTERSECT,NG.SRE.INTERSECT],na.rm=TRUE) < PG.NG.INT.DIST.PERMUTATION.TEST)/length(PG.NG.INT.DIST.PERMUTATION.TEST)
#TH.PH.INT.DIST.PERMUTATION.TEST <- permutation.test.matrix(TH.ERE.SRE.DIST[TH.ERE.INTERSECT,TH.SRE.INTERSECT],PH.ERE.SRE.DIST[PH.ERE.INTERSECT,PH.SRE.INTERSECT],10000)
#  hist(TH.PH.INT.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.1,.1),xlab="Average Difference path distance",main="Common TH vs PH"); abline(v=mean(TH.ERE.SRE.DIST[TH.ERE.INTERSECT,TH.SRE.INTERSECT],na.rm=TRUE) - mean(PH.ERE.SRE.DIST[PH.ERE.INTERSECT,PH.SRE.INTERSECT],na.rm=TRUE),col="red")
#  sum(mean(TH.ERE.SRE.DIST[TH.ERE.INTERSECT,TH.SRE.INTERSECT],na.rm=TRUE) - mean(PH.ERE.SRE.DIST[PH.ERE.INTERSECT,PH.SRE.INTERSECT],na.rm=TRUE) < TH.PH.INT.DIST.PERMUTATION.TEST)/length(TH.PH.INT.DIST.PERMUTATION.TEST)
#PH.MH.INT.DIST.PERMUTATION.TEST <- permutation.test.matrix(PH.ERE.SRE.DIST[PH.ERE.INTERSECT,PH.SRE.INTERSECT],MH.ERE.SRE.DIST[MH.ERE.INTERSECT,MH.SRE.INTERSECT],10000)
#  hist(PH.MH.INT.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.1,.1),xlab="Average Difference path distance",main="Common PH vs MH"); abline(v=mean(PH.ERE.SRE.DIST[PH.ERE.INTERSECT,PH.SRE.INTERSECT],na.rm=TRUE) - mean(MH.ERE.SRE.DIST[MH.ERE.INTERSECT,MH.SRE.INTERSECT],na.rm=TRUE),col="red")
#  sum(mean(PH.ERE.SRE.DIST[PH.ERE.INTERSECT,PH.SRE.INTERSECT],na.rm=TRUE) - mean(MH.ERE.SRE.DIST[MH.ERE.INTERSECT,MH.SRE.INTERSECT],na.rm=TRUE) < PH.MH.INT.DIST.PERMUTATION.TEST)/length(PH.MH.INT.DIST.PERMUTATION.TEST)
#MH.NH.INT.DIST.PERMUTATION.TEST <- permutation.test.matrix(MH.ERE.SRE.DIST[MH.ERE.INTERSECT,MH.SRE.INTERSECT],NH.ERE.SRE.DIST[NH.ERE.INTERSECT,NH.SRE.INTERSECT],10000)
#  hist(MH.NH.INT.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.1,.1),xlab="Average Difference path distance",main="Common MH vs NH"); abline(v=mean(MH.ERE.SRE.DIST[MH.ERE.INTERSECT,MH.SRE.INTERSECT],na.rm=TRUE) - mean(NH.ERE.SRE.DIST[NH.ERE.INTERSECT,NH.SRE.INTERSECT],na.rm=TRUE),col="red")
#  sum(mean(MH.ERE.SRE.DIST[MH.ERE.INTERSECT,MH.SRE.INTERSECT],na.rm=TRUE) - mean(NH.ERE.SRE.DIST[NH.ERE.INTERSECT,NH.SRE.INTERSECT],na.rm=TRUE) < MH.NH.INT.DIST.PERMUTATION.TEST)/length(MH.NH.INT.DIST.PERMUTATION.TEST)
#TH.MH.INT.DIST.PERMUTATION.TEST <- permutation.test.matrix(TH.ERE.SRE.DIST[TH.ERE.INTERSECT,TH.SRE.INTERSECT],MH.ERE.SRE.DIST[MH.ERE.INTERSECT,MH.SRE.INTERSECT],10000)
#  hist(TH.MH.INT.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.1,.1),xlab="Average Difference path distance",main="Common TH vs MH"); abline(v=mean(TH.ERE.SRE.DIST[TH.ERE.INTERSECT,TH.SRE.INTERSECT],na.rm=TRUE) - mean(MH.ERE.SRE.DIST[MH.ERE.INTERSECT,MH.SRE.INTERSECT],na.rm=TRUE),col="red")
#  sum(mean(TH.ERE.SRE.DIST[TH.ERE.INTERSECT,TH.SRE.INTERSECT],na.rm=TRUE) - mean(MH.ERE.SRE.DIST[MH.ERE.INTERSECT,MH.SRE.INTERSECT],na.rm=TRUE) < TH.MH.INT.DIST.PERMUTATION.TEST)/length(TH.MH.INT.DIST.PERMUTATION.TEST)
#TH.NH.INT.DIST.PERMUTATION.TEST <- permutation.test.matrix(TH.ERE.SRE.DIST[TH.ERE.INTERSECT,TH.SRE.INTERSECT],NH.ERE.SRE.DIST[NH.ERE.INTERSECT,NH.SRE.INTERSECT],10000)
#  hist(TH.NH.INT.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.1,.1),xlab="Average Difference path distance",main="Common TH vs NH"); abline(v=mean(TH.ERE.SRE.DIST[TH.ERE.INTERSECT,TH.SRE.INTERSECT],na.rm=TRUE) - mean(NH.ERE.SRE.DIST[NH.ERE.INTERSECT,NH.SRE.INTERSECT],na.rm=TRUE),col="red")
#  sum(mean(TH.ERE.SRE.DIST[TH.ERE.INTERSECT,TH.SRE.INTERSECT],na.rm=TRUE) - mean(NH.ERE.SRE.DIST[NH.ERE.INTERSECT,NH.SRE.INTERSECT],na.rm=TRUE) < TH.NH.INT.DIST.PERMUTATION.TEST)/length(TH.NH.INT.DIST.PERMUTATION.TEST)
#PH.NH.INT.DIST.PERMUTATION.TEST <- permutation.test.matrix(PH.ERE.SRE.DIST[PH.ERE.INTERSECT,PH.SRE.INTERSECT],NH.ERE.SRE.DIST[NH.ERE.INTERSECT,NH.SRE.INTERSECT],10000)
#  hist(PH.NH.INT.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.1,.1),xlab="Average Difference path distance",main="Common PH vs NH"); abline(v=mean(PH.ERE.SRE.DIST[PH.ERE.INTERSECT,PH.SRE.INTERSECT],na.rm=TRUE) - mean(NH.ERE.SRE.DIST[NH.ERE.INTERSECT,NH.SRE.INTERSECT],na.rm=TRUE),col="red")
#  sum(mean(PH.ERE.SRE.DIST[PH.ERE.INTERSECT,PH.SRE.INTERSECT],na.rm=TRUE) - mean(NH.ERE.SRE.DIST[NH.ERE.INTERSECT,NH.SRE.INTERSECT],na.rm=TRUE) < PH.NH.INT.DIST.PERMUTATION.TEST)/length(PH.NH.INT.DIST.PERMUTATION.TEST)
##dev.off()

##pdf("MIN.DISTANCE.pdf")
#par(mar=c(3,3,3,3))
#  
##Distance to closest SRE for all genotypes in network
#par(mfrow=c(2,4))
#barplot(table(factor(TG.ERE.SRE.MIN.DIST,levels=1:4))/length(TG.ERE.SRE.MIN.DIST[is.finite(TG.ERE.SRE.MIN.DIST)]),horiz = TRUE,xlim=c(0,0.6),col="darkred", main="TG ERE to closest SRE Distance: All Genotypes")
#barplot(table(factor(PG.ERE.SRE.MIN.DIST,levels=1:4))/length(PG.ERE.SRE.MIN.DIST[is.finite(PG.ERE.SRE.MIN.DIST)]),horiz = TRUE,xlim=c(0,0.6),col="darkred", main="PG ERE to closest SRE Distance: All Genotypes")
#barplot(table(factor(MG.ERE.SRE.MIN.DIST,levels=1:4))/length(MG.ERE.SRE.MIN.DIST[is.finite(MG.ERE.SRE.MIN.DIST)]),horiz = TRUE,xlim=c(0,0.6),col="darkred", main="MG ERE to closest SRE Distance: All Genotypes")
#barplot(table(factor(NG.ERE.SRE.MIN.DIST,levels=1:4))/length(NG.ERE.SRE.MIN.DIST[is.finite(NG.ERE.SRE.MIN.DIST)]),horiz = TRUE,xlim=c(0,0.6),col="darkred", main="NG ERE to closest SRE Distance: All Genotypes")
#barplot(table(factor(TH.ERE.SRE.MIN.DIST,levels=1:4))/length(TH.ERE.SRE.MIN.DIST[is.finite(TH.ERE.SRE.MIN.DIST)]),horiz = TRUE,xlim=c(0,0.6),col="skyblue2",main="TH ERE to closest SRE Distance: All Genotypes")
#barplot(table(factor(PH.ERE.SRE.MIN.DIST,levels=1:4))/length(PH.ERE.SRE.MIN.DIST[is.finite(PH.ERE.SRE.MIN.DIST)]),horiz = TRUE,xlim=c(0,0.6),col="skyblue2",main="PH ERE to closest SRE Distance: All Genotypes")
#barplot(table(factor(MH.ERE.SRE.MIN.DIST,levels=1:4))/length(MH.ERE.SRE.MIN.DIST[is.finite(MH.ERE.SRE.MIN.DIST)]),horiz = TRUE,xlim=c(0,0.6),col="skyblue2",main="MH ERE to closest SRE Distance: All Genotypes")
#barplot(table(factor(NH.ERE.SRE.MIN.DIST,levels=1:4))/length(NH.ERE.SRE.MIN.DIST[is.finite(NH.ERE.SRE.MIN.DIST)]),horiz = TRUE,xlim=c(0,0.6),col="skyblue2",main="NH ERE to closest SRE Distance: All Genotypes")
#
#par(mfrow=c(1,1))
#barplot(ALL.MIN.DIST.AVG[1:4],names=c("3G","2G","1G","NG"),cex.axis=0.75,cex.main=0.5,col="darkred", main="ERE to closest SRE Distance: All Genotypes")
#barplot(ALL.MIN.DIST.AVG[5:8],names=c("3H","2H","1H","NH"),cex.axis=0.75,cex.main=0.5,col="skyblue2",main="ERE to closest SRE Distance: All Genotypes")
#  
##Distance to closest SRE for common genotypes in network
#par(mfrow=c(2,4))
#barplot(table(factor(TG.ERE.SRE.MIN.DIST[TG.ERE.INTERSECT],levels=1:4))/length(TG.ERE.SRE.MIN.DIST[TG.ERE.INTERSECT][is.finite(TG.ERE.SRE.MIN.DIST[TG.ERE.INTERSECT])]),horiz = TRUE,xlim=c(0,0.6),col="darkred", main="TG ERE to closest SRE Distance: Common Genotypes")
#barplot(table(factor(PG.ERE.SRE.MIN.DIST[PG.ERE.INTERSECT],levels=1:4))/length(PG.ERE.SRE.MIN.DIST[PG.ERE.INTERSECT][is.finite(PG.ERE.SRE.MIN.DIST[PG.ERE.INTERSECT])]),horiz = TRUE,xlim=c(0,0.6),col="darkred", main="PG ERE to closest SRE Distance: Common Genotypes")
#barplot(table(factor(MG.ERE.SRE.MIN.DIST[MG.ERE.INTERSECT],levels=1:4))/length(MG.ERE.SRE.MIN.DIST[MG.ERE.INTERSECT][is.finite(MG.ERE.SRE.MIN.DIST[MG.ERE.INTERSECT])]),horiz = TRUE,xlim=c(0,0.6),col="darkred", main="MG ERE to closest SRE Distance: Common Genotypes")
#barplot(table(factor(NG.ERE.SRE.MIN.DIST[NG.ERE.INTERSECT],levels=1:4))/length(NG.ERE.SRE.MIN.DIST[NG.ERE.INTERSECT][is.finite(NG.ERE.SRE.MIN.DIST[NG.ERE.INTERSECT])]),horiz = TRUE,xlim=c(0,0.6),col="darkred", main="NG ERE to closest SRE Distance: Common Genotypes")
#barplot(table(factor(TH.ERE.SRE.MIN.DIST[TH.ERE.INTERSECT],levels=1:4))/length(TH.ERE.SRE.MIN.DIST[TH.ERE.INTERSECT][is.finite(TH.ERE.SRE.MIN.DIST[TH.ERE.INTERSECT])]),horiz = TRUE,xlim=c(0,0.6),col="skyblue2",main="TH ERE to closest SRE Distance: Common Genotypes")
#barplot(table(factor(PH.ERE.SRE.MIN.DIST[PH.ERE.INTERSECT],levels=1:4))/length(PH.ERE.SRE.MIN.DIST[PH.ERE.INTERSECT][is.finite(PH.ERE.SRE.MIN.DIST[PH.ERE.INTERSECT])]),horiz = TRUE,xlim=c(0,0.6),col="skyblue2",main="PH ERE to closest SRE Distance: Common Genotypes")
#barplot(table(factor(MH.ERE.SRE.MIN.DIST[MH.ERE.INTERSECT],levels=1:4))/length(MH.ERE.SRE.MIN.DIST[MH.ERE.INTERSECT][is.finite(MH.ERE.SRE.MIN.DIST[MH.ERE.INTERSECT])]),horiz = TRUE,xlim=c(0,0.6),col="skyblue2",main="MH ERE to closest SRE Distance: Common Genotypes")
#barplot(table(factor(NH.ERE.SRE.MIN.DIST[NH.ERE.INTERSECT],levels=1:4))/length(NH.ERE.SRE.MIN.DIST[NH.ERE.INTERSECT][is.finite(NH.ERE.SRE.MIN.DIST[NH.ERE.INTERSECT])]),horiz = TRUE,xlim=c(0,0.6),col="skyblue2",main="NH ERE to closest SRE Distance: Common Genotypes")
#
#par(mfrow=c(1,1))
#barplot(INT.MIN.DIST.AVG[1:4],names=c("3G","2G","1G","NG"),cex.axis=0.75,cex.main=0.5,col="darkred", main="ERE to closest SRE Distance: Common Genotypes")
#barplot(INT.MIN.DIST.AVG[5:8],names=c("3H","2H","1H","NH"),cex.axis=0.75,cex.main=0.5,col="skyblue2",main="ERE to closest SRE Distance: Common Genotypes")
#
###Permutation tests
#par(mfrow=c(3,2))
##Minimum distance for all ERE to SRE
#TG.PG.ALL.MIN.DIST.PERMUTATION.TEST <- permutation.test.vector(TG.ERE.SRE.MIN.DIST,PG.ERE.SRE.MIN.DIST,10000)
#  hist(TG.PG.ALL.MIN.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.2,.2),xlab="Average Difference path min distance",main="All TG vs PG"); abline(v=mean(TG.ERE.SRE.MIN.DIST[is.finite(TG.ERE.SRE.MIN.DIST)]) - mean(PG.ERE.SRE.MIN.DIST[is.finite(PG.ERE.SRE.MIN.DIST)]),col="red")
#  sum(mean(TG.ERE.SRE.MIN.DIST[is.finite(TG.ERE.SRE.MIN.DIST)]) - mean(PG.ERE.SRE.MIN.DIST[is.finite(PG.ERE.SRE.MIN.DIST)]) < TG.PG.ALL.MIN.DIST.PERMUTATION.TEST)/length(TG.PG.ALL.MIN.DIST.PERMUTATION.TEST)
#PG.MG.ALL.MIN.DIST.PERMUTATION.TEST <- permutation.test.vector(PG.ERE.SRE.MIN.DIST,MG.ERE.SRE.MIN.DIST,10000)
#  hist(PG.MG.ALL.MIN.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.2,.2),xlab="Average Difference path min distance",main="All PG vs MG"); abline(v=mean(PG.ERE.SRE.MIN.DIST[is.finite(PG.ERE.SRE.MIN.DIST)]) - mean(MG.ERE.SRE.MIN.DIST[is.finite(MG.ERE.SRE.MIN.DIST)]),col="red")
#  sum(mean(PG.ERE.SRE.MIN.DIST[is.finite(PG.ERE.SRE.MIN.DIST)]) - mean(MG.ERE.SRE.MIN.DIST[is.finite(MG.ERE.SRE.MIN.DIST)]) < PG.MG.ALL.MIN.DIST.PERMUTATION.TEST)/length(PG.MG.ALL.MIN.DIST.PERMUTATION.TEST)
#MG.NG.ALL.MIN.DIST.PERMUTATION.TEST <- permutation.test.vector(MG.ERE.SRE.MIN.DIST,NG.ERE.SRE.MIN.DIST,10000)
#  hist(MG.NG.ALL.MIN.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.2,.2),xlab="Average Difference path min distance",main="All MG vs NG"); abline(v=mean(MG.ERE.SRE.MIN.DIST[is.finite(MG.ERE.SRE.MIN.DIST)]) - mean(NG.ERE.SRE.MIN.DIST[is.finite(NG.ERE.SRE.MIN.DIST)]),col="red")
#  sum(mean(MG.ERE.SRE.MIN.DIST[is.finite(MG.ERE.SRE.MIN.DIST)]) - mean(NG.ERE.SRE.MIN.DIST[is.finite(NG.ERE.SRE.MIN.DIST)]) < MG.NG.ALL.MIN.DIST.PERMUTATION.TEST)/length(MG.NG.ALL.MIN.DIST.PERMUTATION.TEST)
#TG.MG.ALL.MIN.DIST.PERMUTATION.TEST <- permutation.test.vector(TG.ERE.SRE.MIN.DIST,MG.ERE.SRE.MIN.DIST,10000)
#  hist(TG.MG.ALL.MIN.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.2,.2),xlab="Average Difference path min distance",main="All TG vs MG"); abline(v=mean(TG.ERE.SRE.MIN.DIST[is.finite(TG.ERE.SRE.MIN.DIST)]) - mean(MG.ERE.SRE.MIN.DIST[is.finite(MG.ERE.SRE.MIN.DIST)]),col="red")
#  sum(mean(TG.ERE.SRE.MIN.DIST[is.finite(TG.ERE.SRE.MIN.DIST)]) - mean(MG.ERE.SRE.MIN.DIST[is.finite(MG.ERE.SRE.MIN.DIST)]) < TG.MG.ALL.MIN.DIST.PERMUTATION.TEST)/length(TG.MG.ALL.MIN.DIST.PERMUTATION.TEST)
#TG.NG.ALL.MIN.DIST.PERMUTATION.TEST <- permutation.test.vector(TG.ERE.SRE.MIN.DIST,NG.ERE.SRE.MIN.DIST,10000)
#  hist(TG.NG.ALL.MIN.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.2,.2),xlab="Average Difference path min distance",main="All TG vs NG"); abline(v=mean(TG.ERE.SRE.MIN.DIST[is.finite(TG.ERE.SRE.MIN.DIST)]) - mean(NG.ERE.SRE.MIN.DIST[is.finite(NG.ERE.SRE.MIN.DIST)]),col="red")
#  sum(mean(TG.ERE.SRE.MIN.DIST[is.finite(TG.ERE.SRE.MIN.DIST)]) - mean(NG.ERE.SRE.MIN.DIST[is.finite(NG.ERE.SRE.MIN.DIST)]) < TG.NG.ALL.MIN.DIST.PERMUTATION.TEST)/length(TG.NG.ALL.MIN.DIST.PERMUTATION.TEST)
#PG.NG.ALL.MIN.DIST.PERMUTATION.TEST <- permutation.test.vector(PG.ERE.SRE.MIN.DIST,NG.ERE.SRE.MIN.DIST,10000)
#  hist(PG.NG.ALL.MIN.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.2,.2),xlab="Average Difference path min distance",main="All PG vs NG"); abline(v=mean(PG.ERE.SRE.MIN.DIST[is.finite(PG.ERE.SRE.MIN.DIST)]) - mean(NG.ERE.SRE.MIN.DIST[is.finite(NG.ERE.SRE.MIN.DIST)]),col="red")
#  sum(mean(PG.ERE.SRE.MIN.DIST[is.finite(PG.ERE.SRE.MIN.DIST)]) - mean(NG.ERE.SRE.MIN.DIST[is.finite(NG.ERE.SRE.MIN.DIST)]) < PG.NG.ALL.MIN.DIST.PERMUTATION.TEST)/length(PG.NG.ALL.MIN.DIST.PERMUTATION.TEST)
#TH.PH.ALL.MIN.DIST.PERMUTATION.TEST <- permutation.test.vector(TH.ERE.SRE.MIN.DIST,PH.ERE.SRE.MIN.DIST,10000)
#  hist(TH.PH.ALL.MIN.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.1,.1),xlab="Average Difference path min distance",main="All TH vs PH"); abline(v=mean(TH.ERE.SRE.MIN.DIST[is.finite(TH.ERE.SRE.MIN.DIST)]) - mean(PH.ERE.SRE.MIN.DIST[is.finite(PH.ERE.SRE.MIN.DIST)]),col="red")
#  sum(mean(TH.ERE.SRE.MIN.DIST[is.finite(TH.ERE.SRE.MIN.DIST)]) - mean(PH.ERE.SRE.MIN.DIST[is.finite(PH.ERE.SRE.MIN.DIST)]) < TH.PH.ALL.MIN.DIST.PERMUTATION.TEST)/length(TH.PH.ALL.MIN.DIST.PERMUTATION.TEST)
#PH.MH.ALL.MIN.DIST.PERMUTATION.TEST <- permutation.test.vector(PH.ERE.SRE.MIN.DIST,MH.ERE.SRE.MIN.DIST,10000)
#  hist(PH.MH.ALL.MIN.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.1,.1),xlab="Average Difference path min distance",main="All PH vs MH"); abline(v=mean(PH.ERE.SRE.MIN.DIST[is.finite(PH.ERE.SRE.MIN.DIST)]) - mean(MH.ERE.SRE.MIN.DIST[is.finite(MH.ERE.SRE.MIN.DIST)]),col="red")
#  sum(mean(PH.ERE.SRE.MIN.DIST[is.finite(PH.ERE.SRE.MIN.DIST)]) - mean(MH.ERE.SRE.MIN.DIST[is.finite(MH.ERE.SRE.MIN.DIST)]) < PH.MH.ALL.MIN.DIST.PERMUTATION.TEST)/length(PH.MH.ALL.MIN.DIST.PERMUTATION.TEST)
#MH.NH.ALL.MIN.DIST.PERMUTATION.TEST <- permutation.test.vector(MH.ERE.SRE.MIN.DIST,NH.ERE.SRE.MIN.DIST,10000)
#  hist(MH.NH.ALL.MIN.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.1,.1),xlab="Average Difference path min distance",main="All MH vs NH"); abline(v=mean(MH.ERE.SRE.MIN.DIST[is.finite(MH.ERE.SRE.MIN.DIST)]) - mean(NH.ERE.SRE.MIN.DIST[is.finite(NH.ERE.SRE.MIN.DIST)]),col="red")
#  sum(mean(MH.ERE.SRE.MIN.DIST[is.finite(MH.ERE.SRE.MIN.DIST)]) - mean(NH.ERE.SRE.MIN.DIST[is.finite(NH.ERE.SRE.MIN.DIST)]) < MH.NH.ALL.MIN.DIST.PERMUTATION.TEST)/length(MH.NH.ALL.MIN.DIST.PERMUTATION.TEST)
#TH.MH.ALL.MIN.DIST.PERMUTATION.TEST <- permutation.test.vector(TH.ERE.SRE.MIN.DIST,MH.ERE.SRE.MIN.DIST,10000)
#  hist(TH.MH.ALL.MIN.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.1,.1),xlab="Average Difference path min distance",main="All TH vs MH"); abline(v=mean(TH.ERE.SRE.MIN.DIST[is.finite(TH.ERE.SRE.MIN.DIST)]) - mean(MH.ERE.SRE.MIN.DIST[is.finite(MH.ERE.SRE.MIN.DIST)]),col="red")
#  sum(mean(TH.ERE.SRE.MIN.DIST[is.finite(TH.ERE.SRE.MIN.DIST)]) - mean(MH.ERE.SRE.MIN.DIST[is.finite(MH.ERE.SRE.MIN.DIST)]) < TH.MH.ALL.MIN.DIST.PERMUTATION.TEST)/length(TH.MH.ALL.MIN.DIST.PERMUTATION.TEST)
#TH.NH.ALL.MIN.DIST.PERMUTATION.TEST <- permutation.test.vector(TH.ERE.SRE.MIN.DIST,NH.ERE.SRE.MIN.DIST,10000)
#  hist(TH.NH.ALL.MIN.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.1,.1),xlab="Average Difference path min distance",main="All TH vs NH"); abline(v=mean(TH.ERE.SRE.MIN.DIST[is.finite(TH.ERE.SRE.MIN.DIST)]) - mean(NH.ERE.SRE.MIN.DIST[is.finite(NH.ERE.SRE.MIN.DIST)]),col="red")
#  sum(mean(TH.ERE.SRE.MIN.DIST[is.finite(TH.ERE.SRE.MIN.DIST)]) - mean(NH.ERE.SRE.MIN.DIST[is.finite(NH.ERE.SRE.MIN.DIST)]) < TH.NH.ALL.MIN.DIST.PERMUTATION.TEST)/length(TH.NH.ALL.MIN.DIST.PERMUTATION.TEST)
#PH.NH.ALL.MIN.DIST.PERMUTATION.TEST <- permutation.test.vector(PH.ERE.SRE.MIN.DIST,NH.ERE.SRE.MIN.DIST,10000)
#  hist(PH.NH.ALL.MIN.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.1,.1),xlab="Average Difference path min distance",main="All PH vs NH"); abline(v=mean(PH.ERE.SRE.MIN.DIST[is.finite(PH.ERE.SRE.MIN.DIST)]) - mean(NH.ERE.SRE.MIN.DIST[is.finite(NH.ERE.SRE.MIN.DIST)]),col="red")
#  sum(mean(PH.ERE.SRE.MIN.DIST[is.finite(PH.ERE.SRE.MIN.DIST)]) - mean(NH.ERE.SRE.MIN.DIST[is.finite(NH.ERE.SRE.MIN.DIST)]) < PH.NH.ALL.MIN.DIST.PERMUTATION.TEST)/length(PH.NH.ALL.MIN.DIST.PERMUTATION.TEST)
#
##Minimum distance for common ERE to SRE
#TG.PG.INT.MIN.DIST.PERMUTATION.TEST <- permutation.test.vector(TG.ERE.SRE.MIN.DIST[TG.ERE.INTERSECT],PG.ERE.SRE.MIN.DIST[PG.ERE.INTERSECT],10000)
#  hist(TG.PG.INT.MIN.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.3,.3),xlab="Average Difference path min distance",main="Common TG vs PG"); abline(v=mean(TG.ERE.SRE.MIN.DIST[TG.ERE.INTERSECT][is.finite(TG.ERE.SRE.MIN.DIST[TG.ERE.INTERSECT])]) - mean(PG.ERE.SRE.MIN.DIST[PG.ERE.INTERSECT][is.finite(PG.ERE.SRE.MIN.DIST[PG.ERE.INTERSECT])]),col="red")
#  sum(mean(TG.ERE.SRE.MIN.DIST[TG.ERE.INTERSECT][is.finite(TG.ERE.SRE.MIN.DIST[TG.ERE.INTERSECT])]) - mean(PG.ERE.SRE.MIN.DIST[PG.ERE.INTERSECT][is.finite(PG.ERE.SRE.MIN.DIST[PG.ERE.INTERSECT])]) < TG.PG.INT.MIN.DIST.PERMUTATION.TEST)/length(TG.PG.INT.MIN.DIST.PERMUTATION.TEST)
#PG.MG.INT.MIN.DIST.PERMUTATION.TEST <- permutation.test.vector(PG.ERE.SRE.MIN.DIST[PG.ERE.INTERSECT],MG.ERE.SRE.MIN.DIST[MG.ERE.INTERSECT],10000)
#  hist(PG.MG.INT.MIN.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.3,.3),xlab="Average Difference path min distance",main="Common PG vs MG"); abline(v=mean(PG.ERE.SRE.MIN.DIST[PG.ERE.INTERSECT][is.finite(PG.ERE.SRE.MIN.DIST[PG.ERE.INTERSECT])]) - mean(MG.ERE.SRE.MIN.DIST[MG.ERE.INTERSECT][is.finite(MG.ERE.SRE.MIN.DIST[MG.ERE.INTERSECT])]),col="red")
#  sum(mean(PG.ERE.SRE.MIN.DIST[PG.ERE.INTERSECT][is.finite(PG.ERE.SRE.MIN.DIST[PG.ERE.INTERSECT])]) - mean(MG.ERE.SRE.MIN.DIST[MG.ERE.INTERSECT][is.finite(MG.ERE.SRE.MIN.DIST[MG.ERE.INTERSECT])]) < PG.MG.INT.MIN.DIST.PERMUTATION.TEST)/length(PG.MG.INT.MIN.DIST.PERMUTATION.TEST)
#MG.NG.INT.MIN.DIST.PERMUTATION.TEST <- permutation.test.vector(MG.ERE.SRE.MIN.DIST[MG.ERE.INTERSECT],NG.ERE.SRE.MIN.DIST[NG.ERE.INTERSECT],10000)
#  hist(MG.NG.INT.MIN.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.3,.3),xlab="Average Difference path min distance",main="Common MG vs NG"); abline(v=mean(MG.ERE.SRE.MIN.DIST[MG.ERE.INTERSECT][is.finite(MG.ERE.SRE.MIN.DIST[MG.ERE.INTERSECT])]) - mean(NG.ERE.SRE.MIN.DIST[NG.ERE.INTERSECT][is.finite(NG.ERE.SRE.MIN.DIST[NG.ERE.INTERSECT])]),col="red")
#  sum(mean(MG.ERE.SRE.MIN.DIST[MG.ERE.INTERSECT][is.finite(MG.ERE.SRE.MIN.DIST[MG.ERE.INTERSECT])]) - mean(NG.ERE.SRE.MIN.DIST[NG.ERE.INTERSECT][is.finite(NG.ERE.SRE.MIN.DIST[NG.ERE.INTERSECT])]) < MG.NG.INT.MIN.DIST.PERMUTATION.TEST)/length(MG.NG.INT.MIN.DIST.PERMUTATION.TEST)
#TG.MG.INT.MIN.DIST.PERMUTATION.TEST <- permutation.test.vector(TG.ERE.SRE.MIN.DIST[TG.ERE.INTERSECT],MG.ERE.SRE.MIN.DIST[MG.ERE.INTERSECT],10000)
#  hist(TG.MG.INT.MIN.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.3,.3),xlab="Average Difference path min distance",main="Common TG vs MG"); abline(v=mean(TG.ERE.SRE.MIN.DIST[TG.ERE.INTERSECT][is.finite(TG.ERE.SRE.MIN.DIST[TG.ERE.INTERSECT])]) - mean(MG.ERE.SRE.MIN.DIST[MG.ERE.INTERSECT][is.finite(MG.ERE.SRE.MIN.DIST[MG.ERE.INTERSECT])]),col="red")
#  sum(mean(TG.ERE.SRE.MIN.DIST[TG.ERE.INTERSECT][is.finite(TG.ERE.SRE.MIN.DIST[TG.ERE.INTERSECT])]) - mean(MG.ERE.SRE.MIN.DIST[MG.ERE.INTERSECT][is.finite(MG.ERE.SRE.MIN.DIST[MG.ERE.INTERSECT])]) < TG.MG.INT.MIN.DIST.PERMUTATION.TEST)/length(TG.MG.INT.MIN.DIST.PERMUTATION.TEST)
#TG.NG.INT.MIN.DIST.PERMUTATION.TEST <- permutation.test.vector(TG.ERE.SRE.MIN.DIST[TG.ERE.INTERSECT],NG.ERE.SRE.MIN.DIST[NG.ERE.INTERSECT],10000)
#  hist(TG.NG.INT.MIN.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.3,.3),xlab="Average Difference path min distance",main="Common TG vs NG"); abline(v=mean(TG.ERE.SRE.MIN.DIST[TG.ERE.INTERSECT][is.finite(TG.ERE.SRE.MIN.DIST[TG.ERE.INTERSECT])]) - mean(NG.ERE.SRE.MIN.DIST[NG.ERE.INTERSECT][is.finite(NG.ERE.SRE.MIN.DIST[NG.ERE.INTERSECT])]),col="red")
#  sum(mean(TG.ERE.SRE.MIN.DIST[TG.ERE.INTERSECT][is.finite(TG.ERE.SRE.MIN.DIST[TG.ERE.INTERSECT])]) - mean(NG.ERE.SRE.MIN.DIST[NG.ERE.INTERSECT][is.finite(NG.ERE.SRE.MIN.DIST[NG.ERE.INTERSECT])]) < TG.NG.INT.MIN.DIST.PERMUTATION.TEST)/length(TG.NG.INT.MIN.DIST.PERMUTATION.TEST)
#PG.NG.INT.MIN.DIST.PERMUTATION.TEST <- permutation.test.vector(PG.ERE.SRE.MIN.DIST[PG.ERE.INTERSECT],NG.ERE.SRE.MIN.DIST[NG.ERE.INTERSECT],10000)
#  hist(PG.NG.INT.MIN.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.3,.3),xlab="Average Difference path min distance",main="Common PG vs NG"); abline(v=mean(PG.ERE.SRE.MIN.DIST[PG.ERE.INTERSECT][is.finite(PG.ERE.SRE.MIN.DIST[PG.ERE.INTERSECT])]) - mean(NG.ERE.SRE.MIN.DIST[NG.ERE.INTERSECT][is.finite(NG.ERE.SRE.MIN.DIST[NG.ERE.INTERSECT])]),col="red")
#  sum(mean(PG.ERE.SRE.MIN.DIST[PG.ERE.INTERSECT][is.finite(PG.ERE.SRE.MIN.DIST[PG.ERE.INTERSECT])]) - mean(NG.ERE.SRE.MIN.DIST[NG.ERE.INTERSECT][is.finite(NG.ERE.SRE.MIN.DIST[NG.ERE.INTERSECT])]) < PG.NG.INT.MIN.DIST.PERMUTATION.TEST)/length(PG.NG.INT.MIN.DIST.PERMUTATION.TEST)
#TH.PH.INT.MIN.DIST.PERMUTATION.TEST <- permutation.test.vector(TH.ERE.SRE.MIN.DIST[TH.ERE.INTERSECT],PH.ERE.SRE.MIN.DIST[PH.ERE.INTERSECT],10000)
#  hist(TH.PH.INT.MIN.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.15,.15),xlab="Average Difference path min distance",main="Common TH vs PH"); abline(v=mean(TH.ERE.SRE.MIN.DIST[TH.ERE.INTERSECT][is.finite(TH.ERE.SRE.MIN.DIST[TH.ERE.INTERSECT])]) - mean(PH.ERE.SRE.MIN.DIST[PH.ERE.INTERSECT][is.finite(PH.ERE.SRE.MIN.DIST[PH.ERE.INTERSECT])]),col="red")
#  sum(mean(TH.ERE.SRE.MIN.DIST[TH.ERE.INTERSECT][is.finite(TH.ERE.SRE.MIN.DIST[TH.ERE.INTERSECT])]) - mean(PH.ERE.SRE.MIN.DIST[PH.ERE.INTERSECT][is.finite(PH.ERE.SRE.MIN.DIST[PH.ERE.INTERSECT])]) < TH.PH.INT.MIN.DIST.PERMUTATION.TEST)/length(TH.PH.INT.MIN.DIST.PERMUTATION.TEST)
#PH.MH.INT.MIN.DIST.PERMUTATION.TEST <- permutation.test.vector(PH.ERE.SRE.MIN.DIST[PH.ERE.INTERSECT],MH.ERE.SRE.MIN.DIST[MH.ERE.INTERSECT],10000)
#  hist(PH.MH.INT.MIN.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.15,.15),xlab="Average Difference path min distance",main="Common PH vs MH"); abline(v=mean(PH.ERE.SRE.MIN.DIST[PH.ERE.INTERSECT][is.finite(PH.ERE.SRE.MIN.DIST[PH.ERE.INTERSECT])]) - mean(MH.ERE.SRE.MIN.DIST[MH.ERE.INTERSECT][is.finite(MH.ERE.SRE.MIN.DIST[MH.ERE.INTERSECT])]),col="red")
#  sum(mean(PH.ERE.SRE.MIN.DIST[PH.ERE.INTERSECT][is.finite(PH.ERE.SRE.MIN.DIST[PH.ERE.INTERSECT])]) - mean(MH.ERE.SRE.MIN.DIST[MH.ERE.INTERSECT][is.finite(MH.ERE.SRE.MIN.DIST[MH.ERE.INTERSECT])]) < PH.MH.INT.MIN.DIST.PERMUTATION.TEST)/length(PH.MH.INT.MIN.DIST.PERMUTATION.TEST)
#MH.NH.INT.MIN.DIST.PERMUTATION.TEST <- permutation.test.vector(MH.ERE.SRE.MIN.DIST[MH.ERE.INTERSECT],NH.ERE.SRE.MIN.DIST[NH.ERE.INTERSECT],10000)
#  hist(MH.NH.INT.MIN.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.15,.15),xlab="Average Difference path min distance",main="Common MH vs NH"); abline(v=mean(MH.ERE.SRE.MIN.DIST[MH.ERE.INTERSECT][is.finite(MH.ERE.SRE.MIN.DIST[MH.ERE.INTERSECT])]) - mean(NH.ERE.SRE.MIN.DIST[NH.ERE.INTERSECT][is.finite(NH.ERE.SRE.MIN.DIST[NH.ERE.INTERSECT])]),col="red")
#  sum(mean(MH.ERE.SRE.MIN.DIST[MH.ERE.INTERSECT][is.finite(MH.ERE.SRE.MIN.DIST[MH.ERE.INTERSECT])]) - mean(NH.ERE.SRE.MIN.DIST[NH.ERE.INTERSECT][is.finite(NH.ERE.SRE.MIN.DIST[NH.ERE.INTERSECT])]) < MH.NH.INT.MIN.DIST.PERMUTATION.TEST)/length(MH.NH.INT.MIN.DIST.PERMUTATION.TEST)
#TH.MH.INT.MIN.DIST.PERMUTATION.TEST <- permutation.test.vector(TH.ERE.SRE.MIN.DIST[TH.ERE.INTERSECT],MH.ERE.SRE.MIN.DIST[MH.ERE.INTERSECT],10000)
#  hist(TH.MH.INT.MIN.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.15,.15),xlab="Average Difference path min distance",main="Common TH vs MH"); abline(v=mean(TH.ERE.SRE.MIN.DIST[TH.ERE.INTERSECT][is.finite(TH.ERE.SRE.MIN.DIST[TH.ERE.INTERSECT])]) - mean(MH.ERE.SRE.MIN.DIST[MH.ERE.INTERSECT][is.finite(MH.ERE.SRE.MIN.DIST[MH.ERE.INTERSECT])]),col="red")
#  sum(mean(TH.ERE.SRE.MIN.DIST[TH.ERE.INTERSECT][is.finite(TH.ERE.SRE.MIN.DIST[TH.ERE.INTERSECT])]) - mean(MH.ERE.SRE.MIN.DIST[MH.ERE.INTERSECT][is.finite(MH.ERE.SRE.MIN.DIST[MH.ERE.INTERSECT])]) < TH.MH.INT.MIN.DIST.PERMUTATION.TEST)/length(TH.MH.INT.MIN.DIST.PERMUTATION.TEST)
#TH.NH.INT.MIN.DIST.PERMUTATION.TEST <- permutation.test.vector(TH.ERE.SRE.MIN.DIST[TH.ERE.INTERSECT],NH.ERE.SRE.MIN.DIST[NH.ERE.INTERSECT],10000)
#  hist(TH.NH.INT.MIN.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.15,.15),xlab="Average Difference path min distance",main="Common TH vs NH"); abline(v=mean(TH.ERE.SRE.MIN.DIST[TH.ERE.INTERSECT][is.finite(TH.ERE.SRE.MIN.DIST[TH.ERE.INTERSECT])]) - mean(NH.ERE.SRE.MIN.DIST[NH.ERE.INTERSECT][is.finite(NH.ERE.SRE.MIN.DIST[NH.ERE.INTERSECT])]),col="red")
#  sum(mean(TH.ERE.SRE.MIN.DIST[TH.ERE.INTERSECT][is.finite(TH.ERE.SRE.MIN.DIST[TH.ERE.INTERSECT])]) - mean(NH.ERE.SRE.MIN.DIST[NH.ERE.INTERSECT][is.finite(NH.ERE.SRE.MIN.DIST[NH.ERE.INTERSECT])]) < TH.NH.INT.MIN.DIST.PERMUTATION.TEST)/length(TH.NH.INT.MIN.DIST.PERMUTATION.TEST)
#PH.NH.INT.MIN.DIST.PERMUTATION.TEST <- permutation.test.vector(PH.ERE.SRE.MIN.DIST[PH.ERE.INTERSECT],NH.ERE.SRE.MIN.DIST[NH.ERE.INTERSECT],10000)
#  hist(PH.NH.INT.MIN.DIST.PERMUTATION.TEST,breaks=50,xlim=c(-.15,.15),xlab="Average Difference path min distance",main="Common PH vs NH"); abline(v=mean(PH.ERE.SRE.MIN.DIST[PH.ERE.INTERSECT][is.finite(PH.ERE.SRE.MIN.DIST[PH.ERE.INTERSECT])]) - mean(NH.ERE.SRE.MIN.DIST[NH.ERE.INTERSECT][is.finite(NH.ERE.SRE.MIN.DIST[NH.ERE.INTERSECT])]),col="red")
#  sum(mean(PH.ERE.SRE.MIN.DIST[PH.ERE.INTERSECT][is.finite(PH.ERE.SRE.MIN.DIST[PH.ERE.INTERSECT])]) - mean(NH.ERE.SRE.MIN.DIST[NH.ERE.INTERSECT][is.finite(NH.ERE.SRE.MIN.DIST[NH.ERE.INTERSECT])]) < PH.NH.INT.MIN.DIST.PERMUTATION.TEST)/length(PH.NH.INT.MIN.DIST.PERMUTATION.TEST)
##dev.off()
  
##################################
##4. Percent of additional steps##
##################################
##Activator Distances
TG.ACT.DIST.INTERSECT <- TG.ACT.DIST[TG.ACT.INTERSECT,TG.ACT.INTERSECT]
PG.ACT.DIST.INTERSECT <- PG.ACT.DIST[PG.ACT.INTERSECT,PG.ACT.INTERSECT]
MG.ACT.DIST.INTERSECT <- MG.ACT.DIST[MG.ACT.INTERSECT,MG.ACT.INTERSECT]
NG.ACT.DIST.INTERSECT <- NG.ACT.DIST[NG.ACT.INTERSECT,NG.ACT.INTERSECT]
TH.ACT.DIST.INTERSECT <- TH.ACT.DIST[TH.ACT.INTERSECT,TH.ACT.INTERSECT]
PH.ACT.DIST.INTERSECT <- PH.ACT.DIST[PH.ACT.INTERSECT,PH.ACT.INTERSECT]
MH.ACT.DIST.INTERSECT <- MH.ACT.DIST[MH.ACT.INTERSECT,MH.ACT.INTERSECT]
NH.ACT.DIST.INTERSECT <- NH.ACT.DIST[NH.ACT.INTERSECT,NH.ACT.INTERSECT]

TG.ACT.DIST.TOTAL <- sum(TG.ACT.DIST.INTERSECT, na.rm = TRUE)/2
PG.ACT.DIST.TOTAL <- sum(PG.ACT.DIST.INTERSECT, na.rm = TRUE)/2
MG.ACT.DIST.TOTAL <- sum(MG.ACT.DIST.INTERSECT, na.rm = TRUE)/2
NG.ACT.DIST.TOTAL <- sum(NG.ACT.DIST.INTERSECT, na.rm = TRUE)/2
TH.ACT.DIST.TOTAL <- sum(TH.ACT.DIST.INTERSECT, na.rm = TRUE)/2
PH.ACT.DIST.TOTAL <- sum(PH.ACT.DIST.INTERSECT, na.rm = TRUE)/2
MH.ACT.DIST.TOTAL <- sum(MH.ACT.DIST.INTERSECT, na.rm = TRUE)/2
NH.ACT.DIST.TOTAL <- sum(NH.ACT.DIST.INTERSECT, na.rm = TRUE)/2
ACT.DIST.TOTAL <- c(TG.ACT.DIST.TOTAL-PG.ACT.DIST.TOTAL,PG.ACT.DIST.TOTAL-MG.ACT.DIST.TOTAL,MG.ACT.DIST.TOTAL-NG.ACT.DIST.TOTAL,NG.ACT.DIST.TOTAL-NH.ACT.DIST.TOTAL,TG.ACT.DIST.TOTAL-TH.ACT.DIST.TOTAL)

#ERE-specific to SRE-specific distances
TG.ERE.SRE.DIST.INTERSECT <- TG.ERE.SRE.DIST[TG.ERE.INTERSECT,TG.SRE.INTERSECT]
PG.ERE.SRE.DIST.INTERSECT <- PG.ERE.SRE.DIST[PG.ERE.INTERSECT,PG.SRE.INTERSECT]
MG.ERE.SRE.DIST.INTERSECT <- MG.ERE.SRE.DIST[MG.ERE.INTERSECT,MG.SRE.INTERSECT]
NG.ERE.SRE.DIST.INTERSECT <- NG.ERE.SRE.DIST[NG.ERE.INTERSECT,NG.SRE.INTERSECT]
TH.ERE.SRE.DIST.INTERSECT <- TH.ERE.SRE.DIST[TH.ERE.INTERSECT,TH.SRE.INTERSECT]
PH.ERE.SRE.DIST.INTERSECT <- PH.ERE.SRE.DIST[PH.ERE.INTERSECT,PH.SRE.INTERSECT]
MH.ERE.SRE.DIST.INTERSECT <- MH.ERE.SRE.DIST[MH.ERE.INTERSECT,MH.SRE.INTERSECT]
NH.ERE.SRE.DIST.INTERSECT <- NH.ERE.SRE.DIST[NH.ERE.INTERSECT,NH.SRE.INTERSECT]

TG.ERE.SRE.DIST.TOTAL <- sum(TG.ERE.SRE.DIST.INTERSECT, na.rm = TRUE)
PG.ERE.SRE.DIST.TOTAL <- sum(PG.ERE.SRE.DIST.INTERSECT, na.rm = TRUE)
MG.ERE.SRE.DIST.TOTAL <- sum(MG.ERE.SRE.DIST.INTERSECT, na.rm = TRUE)
NG.ERE.SRE.DIST.TOTAL <- sum(NG.ERE.SRE.DIST.INTERSECT, na.rm = TRUE)
TH.ERE.SRE.DIST.TOTAL <- sum(TH.ERE.SRE.DIST.INTERSECT, na.rm = TRUE)
PH.ERE.SRE.DIST.TOTAL <- sum(PH.ERE.SRE.DIST.INTERSECT, na.rm = TRUE)
MH.ERE.SRE.DIST.TOTAL <- sum(MH.ERE.SRE.DIST.INTERSECT, na.rm = TRUE)
NH.ERE.SRE.DIST.TOTAL <- sum(NH.ERE.SRE.DIST.INTERSECT, na.rm = TRUE)
ERE.SRE.DIST.TOTAL <- c(TG.ERE.SRE.DIST.TOTAL-PG.ERE.SRE.DIST.TOTAL,PG.ERE.SRE.DIST.TOTAL-MG.ERE.SRE.DIST.TOTAL,MG.ERE.SRE.DIST.TOTAL-NG.ERE.SRE.DIST.TOTAL,NG.ERE.SRE.DIST.TOTAL-NH.ERE.SRE.DIST.TOTAL,TG.ERE.SRE.DIST.TOTAL-TH.ERE.SRE.DIST.TOTAL)

#Closest ERE-specific to SRE-specific distances for common starting genotypes
TG.ERE.SRE.MIN.DIST.INTERSECT <- apply(TG.ERE.SRE.DIST[TG.ERE.INTERSECT,TG.SRE.INTERSECT],1,min, na.rm=TRUE)
PG.ERE.SRE.MIN.DIST.INTERSECT <- apply(PG.ERE.SRE.DIST[PG.ERE.INTERSECT,PG.SRE.INTERSECT],1,min, na.rm=TRUE)
MG.ERE.SRE.MIN.DIST.INTERSECT <- apply(MG.ERE.SRE.DIST[MG.ERE.INTERSECT,MG.SRE.INTERSECT],1,min, na.rm=TRUE)
NG.ERE.SRE.MIN.DIST.INTERSECT <- apply(NG.ERE.SRE.DIST[NG.ERE.INTERSECT,NG.SRE.INTERSECT],1,min, na.rm=TRUE)
TH.ERE.SRE.MIN.DIST.INTERSECT <- apply(TH.ERE.SRE.DIST[TH.ERE.INTERSECT,TH.SRE.INTERSECT],1,min, na.rm=TRUE)
PH.ERE.SRE.MIN.DIST.INTERSECT <- apply(PH.ERE.SRE.DIST[PH.ERE.INTERSECT,PH.SRE.INTERSECT],1,min, na.rm=TRUE)
MH.ERE.SRE.MIN.DIST.INTERSECT <- apply(MH.ERE.SRE.DIST[MH.ERE.INTERSECT,MH.SRE.INTERSECT],1,min, na.rm=TRUE)
NH.ERE.SRE.MIN.DIST.INTERSECT <- apply(NH.ERE.SRE.DIST[NH.ERE.INTERSECT,NH.SRE.INTERSECT],1,min, na.rm=TRUE)

TG.ERE.SRE.MIN.DIST.TOTAL <- sum(TG.ERE.SRE.MIN.DIST.INTERSECT, na.rm = TRUE)
PG.ERE.SRE.MIN.DIST.TOTAL <- sum(PG.ERE.SRE.MIN.DIST.INTERSECT, na.rm = TRUE)
MG.ERE.SRE.MIN.DIST.TOTAL <- sum(MG.ERE.SRE.MIN.DIST.INTERSECT, na.rm = TRUE)
NG.ERE.SRE.MIN.DIST.TOTAL <- sum(NG.ERE.SRE.MIN.DIST.INTERSECT, na.rm = TRUE)
TH.ERE.SRE.MIN.DIST.TOTAL <- sum(TH.ERE.SRE.MIN.DIST.INTERSECT, na.rm = TRUE)
PH.ERE.SRE.MIN.DIST.TOTAL <- sum(PH.ERE.SRE.MIN.DIST.INTERSECT, na.rm = TRUE)
MH.ERE.SRE.MIN.DIST.TOTAL <- sum(MH.ERE.SRE.MIN.DIST.INTERSECT, na.rm = TRUE)
NH.ERE.SRE.MIN.DIST.TOTAL <- sum(NH.ERE.SRE.MIN.DIST.INTERSECT, na.rm = TRUE)
ERE.SRE.MIN.DIST.TOTAL <- c(TG.ERE.SRE.MIN.DIST.TOTAL-PG.ERE.SRE.MIN.DIST.TOTAL,PG.ERE.SRE.MIN.DIST.TOTAL-MG.ERE.SRE.MIN.DIST.TOTAL,MG.ERE.SRE.MIN.DIST.TOTAL-NG.ERE.SRE.MIN.DIST.TOTAL,NG.ERE.SRE.MIN.DIST.TOTAL-NH.ERE.SRE.MIN.DIST.TOTAL,TG.ERE.SRE.MIN.DIST.TOTAL-TH.ERE.SRE.MIN.DIST.TOTAL)

#Average closest ERE-specific to SRE-specific distances for all starting genotypes
ERE.SRE.MIN.DIST.ALL <- data.frame("DISTANCE" = c(
TG.ERE.SRE.MIN.DIST[is.finite(TG.ERE.SRE.MIN.DIST)],
PG.ERE.SRE.MIN.DIST[is.finite(PG.ERE.SRE.MIN.DIST)],
MG.ERE.SRE.MIN.DIST[is.finite(MG.ERE.SRE.MIN.DIST)],
NG.ERE.SRE.MIN.DIST[is.finite(NG.ERE.SRE.MIN.DIST)],
TH.ERE.SRE.MIN.DIST[is.finite(TH.ERE.SRE.MIN.DIST)],
PH.ERE.SRE.MIN.DIST[is.finite(PH.ERE.SRE.MIN.DIST)],
MH.ERE.SRE.MIN.DIST[is.finite(MH.ERE.SRE.MIN.DIST)],
NH.ERE.SRE.MIN.DIST[is.finite(NH.ERE.SRE.MIN.DIST)]),
"NETWORK" = c(rep("TG",length(TG.ERE.SRE.MIN.DIST[is.finite(TG.ERE.SRE.MIN.DIST)])),
rep("PG",length(PG.ERE.SRE.MIN.DIST[is.finite(PG.ERE.SRE.MIN.DIST)])),
rep("MG",length(MG.ERE.SRE.MIN.DIST[is.finite(MG.ERE.SRE.MIN.DIST)])),
rep("NG",length(NG.ERE.SRE.MIN.DIST[is.finite(NG.ERE.SRE.MIN.DIST)])),
rep("TH",length(TH.ERE.SRE.MIN.DIST[is.finite(TH.ERE.SRE.MIN.DIST)])),
rep("PH",length(PH.ERE.SRE.MIN.DIST[is.finite(PH.ERE.SRE.MIN.DIST)])),
rep("MH",length(MH.ERE.SRE.MIN.DIST[is.finite(MH.ERE.SRE.MIN.DIST)])),
rep("NH",length(NH.ERE.SRE.MIN.DIST[is.finite(NH.ERE.SRE.MIN.DIST)]))))
ERE.SRE.MIN.DIST.ALL$NETWORK <- factor(ERE.SRE.MIN.DIST.ALL$NETWORK, levels=c("TG","PG","MG","NG","TH","PH","MH","NH"))

ERE.SRE.MIN.DIST.ALL.INTERSECT <- data.frame("DISTANCE" = c(
TG.ERE.SRE.MIN.DIST.INTERSECT[is.finite(TG.ERE.SRE.MIN.DIST.INTERSECT)],
PG.ERE.SRE.MIN.DIST.INTERSECT[is.finite(PG.ERE.SRE.MIN.DIST.INTERSECT)],
MG.ERE.SRE.MIN.DIST.INTERSECT[is.finite(MG.ERE.SRE.MIN.DIST.INTERSECT)],
NG.ERE.SRE.MIN.DIST.INTERSECT[is.finite(NG.ERE.SRE.MIN.DIST.INTERSECT)],
TH.ERE.SRE.MIN.DIST.INTERSECT[is.finite(TH.ERE.SRE.MIN.DIST.INTERSECT)],
PH.ERE.SRE.MIN.DIST.INTERSECT[is.finite(PH.ERE.SRE.MIN.DIST.INTERSECT)],
MH.ERE.SRE.MIN.DIST.INTERSECT[is.finite(MH.ERE.SRE.MIN.DIST.INTERSECT)],
NH.ERE.SRE.MIN.DIST.INTERSECT[is.finite(NH.ERE.SRE.MIN.DIST.INTERSECT)]),
"NETWORK" = c(rep("TG",length(TG.ERE.SRE.MIN.DIST.INTERSECT[is.finite(TG.ERE.SRE.MIN.DIST.INTERSECT)])),
rep("PG",length(PG.ERE.SRE.MIN.DIST.INTERSECT[is.finite(PG.ERE.SRE.MIN.DIST.INTERSECT)])),
rep("MG",length(MG.ERE.SRE.MIN.DIST.INTERSECT[is.finite(MG.ERE.SRE.MIN.DIST.INTERSECT)])),
rep("NG",length(NG.ERE.SRE.MIN.DIST.INTERSECT[is.finite(NG.ERE.SRE.MIN.DIST.INTERSECT)])),
rep("TH",length(TH.ERE.SRE.MIN.DIST.INTERSECT[is.finite(TH.ERE.SRE.MIN.DIST.INTERSECT)])),
rep("PH",length(PH.ERE.SRE.MIN.DIST.INTERSECT[is.finite(PH.ERE.SRE.MIN.DIST.INTERSECT)])),
rep("MH",length(MH.ERE.SRE.MIN.DIST.INTERSECT[is.finite(MH.ERE.SRE.MIN.DIST.INTERSECT)])),
rep("NH",length(NH.ERE.SRE.MIN.DIST.INTERSECT[is.finite(NH.ERE.SRE.MIN.DIST.INTERSECT)]))))
ERE.SRE.MIN.DIST.ALL.INTERSECT$NETWORK <- factor(ERE.SRE.MIN.DIST.ALL.INTERSECT$NETWORK, levels=c("TG","PG","MG","NG","TH","PH","MH","NH"))

##pdf("ADDITIONAL.STEPS.pdf")
###Additional steps for intersected data to look at effects of epistasis
#par(mfrow=c(3,2))
#barplot(table(factor(NG.ACT.DIST.INTERSECT - NH.ACT.DIST.INTERSECT, levels=-3:4))/length(NH.ACT.DIST.INTERSECT), ylim=c(0,1), ylab="All Activators", xlab="Additional Steps-Genetic Code",main="All Activators")
#barplot(table(factor(MG.ACT.DIST.INTERSECT - NG.ACT.DIST.INTERSECT, levels=-3:4))/length(NH.ACT.DIST.INTERSECT), ylim=c(0,1), ylab="All Activators", xlab="Additional Steps-Main Effects")
#barplot(table(factor(PG.ACT.DIST.INTERSECT - MG.ACT.DIST.INTERSECT, levels=-3:4))/length(NH.ACT.DIST.INTERSECT), ylim=c(0,1), ylab="All Activators", xlab="Additional Steps-Pairwise Epistasis")
#barplot(table(factor(TG.ACT.DIST.INTERSECT - PG.ACT.DIST.INTERSECT, levels=-3:4))/length(NH.ACT.DIST.INTERSECT), ylim=c(0,1), ylab="All Activators", xlab="Additional Steps-Third Order Epistasis")
#barplot(table(factor(TG.ACT.DIST.INTERSECT - TH.ACT.DIST.INTERSECT, levels=-3:4))/length(NH.ACT.DIST.INTERSECT), ylim=c(0,1), ylab="All Activators", xlab="Additional Steps-Genetic Code with Epistasis")
#
#par(mfrow=c(3,2))
#barplot(table(factor(NG.ERE.SRE.DIST.INTERSECT - NH.ERE.SRE.DIST.INTERSECT, levels=-3:4))/length(NH.ERE.SRE.DIST.INTERSECT), ylim=c(0,1), ylab="ERE to SRE", xlab="Additional Steps-Genetic Code",main="ERE to SRE")
#barplot(table(factor(MG.ERE.SRE.DIST.INTERSECT - NG.ERE.SRE.DIST.INTERSECT, levels=-3:4))/length(NH.ERE.SRE.DIST.INTERSECT), ylim=c(0,1), ylab="ERE to SRE", xlab="Additional Steps-Main Effects")
#barplot(table(factor(PG.ERE.SRE.DIST.INTERSECT - MG.ERE.SRE.DIST.INTERSECT, levels=-3:4))/length(NH.ERE.SRE.DIST.INTERSECT), ylim=c(0,1), ylab="ERE to SRE", xlab="Additional Steps-Pairwise Epistasis")
#barplot(table(factor(TG.ERE.SRE.DIST.INTERSECT - PG.ERE.SRE.DIST.INTERSECT, levels=-3:4))/length(NH.ERE.SRE.DIST.INTERSECT), ylim=c(0,1), ylab="ERE to SRE", xlab="Additional Steps-Third Order Epistasis")
#barplot(table(factor(TG.ERE.SRE.DIST.INTERSECT - TH.ERE.SRE.DIST.INTERSECT, levels=-3:4))/length(NH.ERE.SRE.DIST.INTERSECT), ylim=c(0,1), ylab="ERE to SRE", xlab="Additional Steps-Genetic Code with Epistasis")
#
#par(mfrow=c(3,2))
#barplot(table(factor(NG.ERE.SRE.MIN.DIST.INTERSECT - NH.ERE.SRE.MIN.DIST.INTERSECT, levels=-3:4))/length(NH.ERE.SRE.MIN.DIST.INTERSECT), ylim=c(0,1), ylab="Min ERE to SRE", xlab="Additional Steps-Genetic Code",main="ERE to closest SRE")
#barplot(table(factor(MG.ERE.SRE.MIN.DIST.INTERSECT - NG.ERE.SRE.MIN.DIST.INTERSECT, levels=-3:4))/length(NH.ERE.SRE.MIN.DIST.INTERSECT), ylim=c(0,1), ylab="Min ERE to SRE", xlab="Additional Steps-Main Effects")
#barplot(table(factor(PG.ERE.SRE.MIN.DIST.INTERSECT - MG.ERE.SRE.MIN.DIST.INTERSECT, levels=-3:4))/length(NH.ERE.SRE.MIN.DIST.INTERSECT), ylim=c(0,1), ylab="Min ERE to SRE", xlab="Additional Steps-Pairwise Epistasis")
#barplot(table(factor(TG.ERE.SRE.MIN.DIST.INTERSECT - PG.ERE.SRE.MIN.DIST.INTERSECT, levels=-3:4))/length(NH.ERE.SRE.MIN.DIST.INTERSECT), ylim=c(0,1), ylab="Min ERE to SRE", xlab="Additional Steps-Third Order Epistasis")
#barplot(table(factor(TG.ERE.SRE.MIN.DIST.INTERSECT - TH.ERE.SRE.MIN.DIST.INTERSECT, levels=-3:4))/length(NH.ERE.SRE.MIN.DIST.INTERSECT), ylim=c(0,1), ylab="Min ERE to SRE", xlab="Additional Steps-Genetic Code with Epistasis")
#
####Percent additional steps 
#par(mfrow=c(1,2))
###Including genetic code
#barplot(cbind(ACT.DIST.TOTAL[1:4]/sum(ACT.DIST.TOTAL[1:4]),ERE.SRE.DIST.TOTAL[1:4]/sum(ERE.SRE.DIST.TOTAL[1:4]),ERE.SRE.MIN.DIST.TOTAL[1:4]/sum(ERE.SRE.MIN.DIST.TOTAL[1:4])),names.arg=c("All Activators","ERE to SRE","ERE to closest SRE"),ylab="% Additional Steps",main="Including Genetic Code",col=c("black","firebrick","skyblue1","darkgreen"))
#
##Excluding genetic code
#barplot(cbind(ACT.DIST.TOTAL[1:3]/sum(ACT.DIST.TOTAL[1:3]),ERE.SRE.DIST.TOTAL[1:3]/sum(ERE.SRE.DIST.TOTAL[1:3]),ERE.SRE.MIN.DIST.TOTAL[1:3]/sum(ERE.SRE.MIN.DIST.TOTAL[1:3])),names.arg=c("All Activators","ERE to SRE","ERE to closest SRE"),ylab="% Additional Steps",main="Excluding Genetic Code",col=c("black","firebrick","skyblue1","darkgreen"))

###Number of paths increased or decreased due to epistasis
#par(mfrow=c(2,2))
##All activators
#PG.MG.ACT.EPI.INCREASE <- sum((PG.ACT.DIST.INTERSECT - MG.ACT.DIST.INTERSECT) > 0)
#PG.MG.ACT.EPI.DECREASE <- sum((PG.ACT.DIST.INTERSECT - MG.ACT.DIST.INTERSECT) < 0)
#pie(c(PG.MG.ACT.EPI.INCREASE,PG.MG.ACT.EPI.DECREASE),labels=c("Increase","Decrease"),col=c("firebrick","skyblue1"),main="Pairwise Epistasis effect:Activators")
#TG.PG.ACT.EPI.INCREASE <- sum((TG.ACT.DIST.INTERSECT - PG.ACT.DIST.INTERSECT) > 0)
#TG.PG.ACT.EPI.DECREASE <- sum((TG.ACT.DIST.INTERSECT - PG.ACT.DIST.INTERSECT) < 0)
#pie(c(TG.PG.ACT.EPI.INCREASE,TG.PG.ACT.EPI.DECREASE),labels=c("Increase","Decrease"),col=c("firebrick","skyblue1"),main="Threeway Epistasis effect:Activators")
#
##ERE-specific to SRE-specific
#PG.MG.ERE.SRE.EPI.INCREASE <- sum((PG.ERE.SRE.DIST.INTERSECT - MG.ERE.SRE.DIST.INTERSECT) > 0)
#PG.MG.ERE.SRE.EPI.DECREASE <- sum((PG.ERE.SRE.DIST.INTERSECT - MG.ERE.SRE.DIST.INTERSECT) < 0)
#pie(c(PG.MG.ERE.SRE.EPI.INCREASE,PG.MG.ERE.SRE.EPI.DECREASE),labels=c("Increase","Decrease"),col=c("firebrick","skyblue1"),main="Pairwise Epistasis effect:ERE to SRE")
#TG.PG.ERE.SRE.EPI.INCREASE <- sum((TG.ERE.SRE.DIST.INTERSECT - PG.ERE.SRE.DIST.INTERSECT) > 0)
#TG.PG.ERE.SRE.EPI.DECREASE <- sum((TG.ERE.SRE.DIST.INTERSECT - PG.ERE.SRE.DIST.INTERSECT) < 0)
#pie(c(TG.PG.ERE.SRE.EPI.INCREASE,TG.PG.ERE.SRE.EPI.DECREASE),labels=c("Increase","Decrease"),col=c("firebrick","skyblue1"),main="Threeway Epistasis effect:ERE to SRE")
##dev.off()

####################
##5. Network plots##
####################
###Network where a genotype is functional in at least one model
DT.JOINT <- subset(DT.JOINT,DT.JOINT$RE=="E" & (DT.JOINT$TE.JOINT.CLASS != "null" | DT.JOINT$PE.JOINT.CLASS != "null" | DT.JOINT$ME.JOINT.CLASS != "null"))
NODES.ALL.ACT <- data.frame(id=1:length(DT.JOINT$AAseq),label=substr(DT.JOINT$AAseq,2,5),stringsAsFactors = F)
NODES.ALL.ACT$TE <- DT.JOINT$TE.JOINT.CLASS
NODES.ALL.ACT$PE <- DT.JOINT$PE.JOINT.CLASS
NODES.ALL.ACT$ME <- DT.JOINT$ME.JOINT.CLASS
NODES.ALL.ACT$X1 <- DT.JOINT$AA1
NODES.ALL.ACT$X2 <- DT.JOINT$AA2
NODES.ALL.ACT$X3 <- DT.JOINT$AA3
NODES.ALL.ACT$X4 <- DT.JOINT$AA4
SOURCE.ALL.ACT <- vector()
TARGET.ALL.ACT <- vector()
I <- 1:nrow(NODES.ALL.ACT)
for(i in I) {
  J <- which(NODES.ALL.ACT$label %in% get.neighbors(NODES.ALL.ACT$label[i],code="Z"))
  for(j in J){
    SOURCE.ALL.ACT <- c(SOURCE.ALL.ACT, i)
    TARGET.ALL.ACT <- c(TARGET.ALL.ACT, j)
  }
}
#Remove duplicate edges
EDGE.SET <- (SOURCE.ALL.ACT < TARGET.ALL.ACT)
EDGES.ALL.ACT <- data.frame(source=SOURCE.ALL.ACT[EDGE.SET], target=TARGET.ALL.ACT[EDGE.SET], stringsAsFactors = F)

#Label edges to nulls in a network
EDGES.ALL.ACT$ME.NULL <- apply(EDGES.ALL.ACT[,1:2],1,FUN = function(x) {
  if(NODES.ALL.ACT$ME[x[1]] == "null" | NODES.ALL.ACT$ME[x[2]] == "null") {
    return("null")
  } else {
    return("Activators")
  }
})
EDGES.ALL.ACT$PE.NULL <- apply(EDGES.ALL.ACT[,1:2],1,FUN = function(x) {
  if(NODES.ALL.ACT$PE[x[1]] == "null" | NODES.ALL.ACT$PE[x[2]] == "null") {
    return("null")
  } else {
    return("Activators")
  }
})
EDGES.ALL.ACT$TE.NULL <- apply(EDGES.ALL.ACT[,1:2],1,FUN = function(x) {
  if(NODES.ALL.ACT$TE[x[1]] == "null" | NODES.ALL.ACT$TE[x[2]] == "null") {
    return("null")
  } else {
    return("Activators")
  }
})
#write.gexf(nodes=NODES.ALL.ACT[,1:2],nodesAtt=NODES.ALL.ACT[,3:9],edges=EDGES.ALL.ACT[,1:2],edgesAtt=EDGES.ALL.ACT[,3:5],output="ALL.ACT.gexf")

##Only main and third order models
#DT.JOINT <- subset(DT.JOINT,DT.JOINT$RE=="E" & (DT.JOINT$TE.JOINT.CLASS != "null" | DT.JOINT$ME.JOINT.CLASS != "null"))
#NODES.ALL.ACT <- data.frame(id=1:length(DT.JOINT$AAseq),label=substr(DT.JOINT$AAseq,2,5),stringsAsFactors = F)
#NODES.ALL.ACT$TE <- DT.JOINT$TE.JOINT.CLASS
#NODES.ALL.ACT$ME <- DT.JOINT$ME.JOINT.CLASS
#NODES.ALL.ACT$X1 <- DT.JOINT$AA1
#NODES.ALL.ACT$X2 <- DT.JOINT$AA2
#NODES.ALL.ACT$X3 <- DT.JOINT$AA3
#NODES.ALL.ACT$X4 <- DT.JOINT$AA4
#SOURCE.ALL.ACT <- vector()
#TARGET.ALL.ACT <- vector()
#I <- 1:nrow(NODES.ALL.ACT)
#for(i in I) {
#  J <- which(NODES.ALL.ACT$label %in% get.neighbors(NODES.ALL.ACT$label[i],code="Z"))
#  for(j in J){
#    SOURCE.ALL.ACT <- c(SOURCE.ALL.ACT, i)
#    TARGET.ALL.ACT <- c(TARGET.ALL.ACT, j)
#  }
#}
##Remove duplicate edges
#EDGE.SET <- (SOURCE.ALL.ACT < TARGET.ALL.ACT)
#EDGES.ALL.ACT <- data.frame(source=SOURCE.ALL.ACT[EDGE.SET], target=TARGET.ALL.ACT[EDGE.SET], stringsAsFactors = F)
#
##Label edges to nulls in a network
#EDGES.ALL.ACT$ME.NULL <- apply(EDGES.ALL.ACT[,1:2],1,FUN = function(x) {
#  if(NODES.ALL.ACT$ME[x[1]] == "null" | NODES.ALL.ACT$ME[x[2]] == "null") {
#    return("null")
#  } else {
#    return("Activators")
#  }
#})
#EDGES.ALL.ACT$TE.NULL <- apply(EDGES.ALL.ACT[,1:2],1,FUN = function(x) {
#  if(NODES.ALL.ACT$TE[x[1]] == "null" | NODES.ALL.ACT$TE[x[2]] == "null") {
#    return("null")
#  } else {
#    return("Activators")
#  }
#})
#write.gexf(nodes=NODES.ALL.ACT[,1:2],nodesAtt=NODES.ALL.ACT[,3:8],edges=EDGES.ALL.ACT[,1:2],edgesAtt=EDGES.ALL.ACT[,3:4],output="ME.TE.ACT.gexf")


########################################
##6. Number of activators and overlaps##
########################################
NUM.ACTIVATORS <- matrix(0,nrow=4,ncol=3); rownames(NUM.ACTIVATORS) <- c("ERE-specific","SRE-specific","promiscuous","null"); colnames(NUM.ACTIVATORS) <- c("TE","PE","ME")
NUM.ACTIVATORS[c(1,4,3,2),1] <- table(NODES.ALL.ACT$TE)
NUM.ACTIVATORS[c(1,4,3,2),2] <- table(NODES.ALL.ACT$PE)
NUM.ACTIVATORS[c(1,4,3,2),3] <- table(NODES.ALL.ACT$ME)

NUM.ACT.CLASS <- matrix(0,nrow=16,ncol=3)
colnames(NUM.ACT.CLASS) <- c("ME.PE","PE.TE","ME.TE")
rownames(NUM.ACT.CLASS) <- c("NN","NP","NS","NE","PN","PP","PS","PE","SN","SP","SS","SE","EN","EP","ES","EE")

NUM.ACT.CLASS[1,1]  <- sum(DT.JOINT$PE.JOINT.CLASS == "null"         & DT.JOINT$ME.JOINT.CLASS == "null")
NUM.ACT.CLASS[2,1]  <- sum(DT.JOINT$PE.JOINT.CLASS == "null"         & DT.JOINT$ME.JOINT.CLASS == "promiscuous")
NUM.ACT.CLASS[3,1]  <- sum(DT.JOINT$PE.JOINT.CLASS == "null"         & DT.JOINT$ME.JOINT.CLASS == "SRE-specific")
NUM.ACT.CLASS[4,1]  <- sum(DT.JOINT$PE.JOINT.CLASS == "null"         & DT.JOINT$ME.JOINT.CLASS == "ERE-specific")
NUM.ACT.CLASS[5,1]  <- sum(DT.JOINT$PE.JOINT.CLASS == "promiscuous"  & DT.JOINT$ME.JOINT.CLASS == "null")
NUM.ACT.CLASS[6,1]  <- sum(DT.JOINT$PE.JOINT.CLASS == "promiscuous"  & DT.JOINT$ME.JOINT.CLASS == "promiscuous")
NUM.ACT.CLASS[7,1]  <- sum(DT.JOINT$PE.JOINT.CLASS == "promiscuous"  & DT.JOINT$ME.JOINT.CLASS == "SRE-specific")
NUM.ACT.CLASS[8,1]  <- sum(DT.JOINT$PE.JOINT.CLASS == "promiscuous"  & DT.JOINT$ME.JOINT.CLASS == "ERE-specific")
NUM.ACT.CLASS[9,1]  <- sum(DT.JOINT$PE.JOINT.CLASS == "SRE-specific" & DT.JOINT$ME.JOINT.CLASS == "null")
NUM.ACT.CLASS[10,1] <- sum(DT.JOINT$PE.JOINT.CLASS == "SRE-specific" & DT.JOINT$ME.JOINT.CLASS == "promiscuous")
NUM.ACT.CLASS[11,1] <- sum(DT.JOINT$PE.JOINT.CLASS == "SRE-specific" & DT.JOINT$ME.JOINT.CLASS == "SRE-specific")
NUM.ACT.CLASS[12,1] <- sum(DT.JOINT$PE.JOINT.CLASS == "SRE-specific" & DT.JOINT$ME.JOINT.CLASS == "ERE-specific")
NUM.ACT.CLASS[13,1] <- sum(DT.JOINT$PE.JOINT.CLASS == "ERE-specific" & DT.JOINT$ME.JOINT.CLASS == "null")
NUM.ACT.CLASS[14,1] <- sum(DT.JOINT$PE.JOINT.CLASS == "ERE-specific" & DT.JOINT$ME.JOINT.CLASS == "promiscuous")
NUM.ACT.CLASS[15,1] <- sum(DT.JOINT$PE.JOINT.CLASS == "ERE-specific" & DT.JOINT$ME.JOINT.CLASS == "SRE-specific")
NUM.ACT.CLASS[16,1] <- sum(DT.JOINT$PE.JOINT.CLASS == "ERE-specific" & DT.JOINT$ME.JOINT.CLASS == "ERE-specific")

NUM.ACT.CLASS[1,2]  <- sum(DT.JOINT$TE.JOINT.CLASS == "null"         & DT.JOINT$PE.JOINT.CLASS == "null")
NUM.ACT.CLASS[2,2]  <- sum(DT.JOINT$TE.JOINT.CLASS == "null"         & DT.JOINT$PE.JOINT.CLASS == "promiscuous")
NUM.ACT.CLASS[3,2]  <- sum(DT.JOINT$TE.JOINT.CLASS == "null"         & DT.JOINT$PE.JOINT.CLASS == "SRE-specific")
NUM.ACT.CLASS[4,2]  <- sum(DT.JOINT$TE.JOINT.CLASS == "null"         & DT.JOINT$PE.JOINT.CLASS == "ERE-specific")
NUM.ACT.CLASS[5,2]  <- sum(DT.JOINT$TE.JOINT.CLASS == "promiscuous"  & DT.JOINT$PE.JOINT.CLASS == "null")
NUM.ACT.CLASS[6,2]  <- sum(DT.JOINT$TE.JOINT.CLASS == "promiscuous"  & DT.JOINT$PE.JOINT.CLASS == "promiscuous")
NUM.ACT.CLASS[7,2]  <- sum(DT.JOINT$TE.JOINT.CLASS == "promiscuous"  & DT.JOINT$PE.JOINT.CLASS == "SRE-specific")
NUM.ACT.CLASS[8,2]  <- sum(DT.JOINT$TE.JOINT.CLASS == "promiscuous"  & DT.JOINT$PE.JOINT.CLASS == "ERE-specific")
NUM.ACT.CLASS[9,2]  <- sum(DT.JOINT$TE.JOINT.CLASS == "SRE-specific" & DT.JOINT$PE.JOINT.CLASS == "null")
NUM.ACT.CLASS[10,2] <- sum(DT.JOINT$TE.JOINT.CLASS == "SRE-specific" & DT.JOINT$PE.JOINT.CLASS == "promiscuous")
NUM.ACT.CLASS[11,2] <- sum(DT.JOINT$TE.JOINT.CLASS == "SRE-specific" & DT.JOINT$PE.JOINT.CLASS == "SRE-specific")
NUM.ACT.CLASS[12,2] <- sum(DT.JOINT$TE.JOINT.CLASS == "SRE-specific" & DT.JOINT$PE.JOINT.CLASS == "ERE-specific")
NUM.ACT.CLASS[13,2] <- sum(DT.JOINT$TE.JOINT.CLASS == "ERE-specific" & DT.JOINT$PE.JOINT.CLASS == "null")
NUM.ACT.CLASS[14,2] <- sum(DT.JOINT$TE.JOINT.CLASS == "ERE-specific" & DT.JOINT$PE.JOINT.CLASS == "promiscuous")
NUM.ACT.CLASS[15,2] <- sum(DT.JOINT$TE.JOINT.CLASS == "ERE-specific" & DT.JOINT$PE.JOINT.CLASS == "SRE-specific")
NUM.ACT.CLASS[16,2] <- sum(DT.JOINT$TE.JOINT.CLASS == "ERE-specific" & DT.JOINT$PE.JOINT.CLASS == "ERE-specific")

NUM.ACT.CLASS[1,3]  <- sum(DT.JOINT$TE.JOINT.CLASS == "null"         & DT.JOINT$ME.JOINT.CLASS == "null")
NUM.ACT.CLASS[2,3]  <- sum(DT.JOINT$TE.JOINT.CLASS == "null"         & DT.JOINT$ME.JOINT.CLASS == "promiscuous")
NUM.ACT.CLASS[3,3]  <- sum(DT.JOINT$TE.JOINT.CLASS == "null"         & DT.JOINT$ME.JOINT.CLASS == "SRE-specific")
NUM.ACT.CLASS[4,3]  <- sum(DT.JOINT$TE.JOINT.CLASS == "null"         & DT.JOINT$ME.JOINT.CLASS == "ERE-specific")
NUM.ACT.CLASS[5,3]  <- sum(DT.JOINT$TE.JOINT.CLASS == "promiscuous"  & DT.JOINT$ME.JOINT.CLASS == "null")
NUM.ACT.CLASS[6,3]  <- sum(DT.JOINT$TE.JOINT.CLASS == "promiscuous"  & DT.JOINT$ME.JOINT.CLASS == "promiscuous")
NUM.ACT.CLASS[7,3]  <- sum(DT.JOINT$TE.JOINT.CLASS == "promiscuous"  & DT.JOINT$ME.JOINT.CLASS == "SRE-specific")
NUM.ACT.CLASS[8,3]  <- sum(DT.JOINT$TE.JOINT.CLASS == "promiscuous"  & DT.JOINT$ME.JOINT.CLASS == "ERE-specific")
NUM.ACT.CLASS[9,3]  <- sum(DT.JOINT$TE.JOINT.CLASS == "SRE-specific" & DT.JOINT$ME.JOINT.CLASS == "null")
NUM.ACT.CLASS[10,3] <- sum(DT.JOINT$TE.JOINT.CLASS == "SRE-specific" & DT.JOINT$ME.JOINT.CLASS == "promiscuous")
NUM.ACT.CLASS[11,3] <- sum(DT.JOINT$TE.JOINT.CLASS == "SRE-specific" & DT.JOINT$ME.JOINT.CLASS == "SRE-specific")
NUM.ACT.CLASS[12,3] <- sum(DT.JOINT$TE.JOINT.CLASS == "SRE-specific" & DT.JOINT$ME.JOINT.CLASS == "ERE-specific")
NUM.ACT.CLASS[13,3] <- sum(DT.JOINT$TE.JOINT.CLASS == "ERE-specific" & DT.JOINT$ME.JOINT.CLASS == "null")
NUM.ACT.CLASS[14,3] <- sum(DT.JOINT$TE.JOINT.CLASS == "ERE-specific" & DT.JOINT$ME.JOINT.CLASS == "promiscuous")
NUM.ACT.CLASS[15,3] <- sum(DT.JOINT$TE.JOINT.CLASS == "ERE-specific" & DT.JOINT$ME.JOINT.CLASS == "SRE-specific")
NUM.ACT.CLASS[16,3] <- sum(DT.JOINT$TE.JOINT.CLASS == "ERE-specific" & DT.JOINT$ME.JOINT.CLASS == "ERE-specific")

#pdf("NUMBER.ACTIVATOR.pdf")
par(mfrow=c(1,1))
barplot(NUM.ACTIVATORS[4:1,3:1],col=c("black","cyan","forestgreen","purple"))
barplot(NUM.ACT.CLASS,col=c("black","cyan","forestgreen","purple"))
#dev.off()

TG.ONLY <- sum(!(TG.ACT.LIST %in% PG.ACT.LIST) & !(TG.ACT.LIST %in% MG.ACT.LIST))
PG.ONLY <- sum(!(PG.ACT.LIST %in% TG.ACT.LIST) & !(PG.ACT.LIST %in% MG.ACT.LIST))
MG.ONLY <- sum(!(MG.ACT.LIST %in% PG.ACT.LIST) & !(MG.ACT.LIST %in% TG.ACT.LIST))
TG.MG.INTERSECT <- sum(!(TG.ACT.LIST %in% PG.ACT.LIST) &  (TG.ACT.LIST %in% MG.ACT.LIST))
TG.PG.INTERSECT <- sum( (TG.ACT.LIST %in% PG.ACT.LIST) & !(TG.ACT.LIST %in% MG.ACT.LIST))
PG.MG.INTERSECT <- sum( (PG.ACT.LIST %in% MG.ACT.LIST) & !(PG.ACT.LIST %in% TG.ACT.LIST))
ALL.INTERSECT <-  sum((TG.ACT.LIST %in% PG.ACT.LIST) & (TG.ACT.LIST %in% MG.ACT.LIST))

c(TG.ONLY,PG.ONLY,MG.ONLY,TG.MG.INTERSECT,TG.PG.INTERSECT,PG.MG.INTERSECT,ALL.INTERSECT)

###################
##7. Path Example##
###################
##Find specific examples of paths for illustrative purposes
#GC +1, Function +1, Epistasis +1
EXAMPLE.1 <- which(NH.ERE.SRE.DIST.INTERSECT == 3 & NG.ERE.SRE.DIST.INTERSECT == 4 & MG.ERE.SRE.DIST.INTERSECT == 5 & PG.ERE.SRE.DIST.INTERSECT == 6,arr.ind = TRUE)

#GC +1, Function +1, Epistasis -1
EXAMPLE.2 <- which(NH.ERE.SRE.DIST.INTERSECT == 3 & NG.ERE.SRE.DIST.INTERSECT == 4 & MG.ERE.SRE.DIST.INTERSECT == 5 & PG.ERE.SRE.DIST.INTERSECT == 4,arr.ind = TRUE)

#Look for overlap in starting point of examples
EXAMPLE.1.INTERSECT <- which(EXAMPLE.1[,1] %in% EXAMPLE.2[,1])
EXAMPLE.2.INTERSECT <- which(EXAMPLE.2[,1] %in% EXAMPLE.1[,1])

subset(EXAMPLE.1, rownames(EXAMPLE.1) == "HGSC")
subset(EXAMPLE.2, rownames(EXAMPLE.2) == "HGSC")

NH.ERE.SRE.DIST.INTERSECT[ 13,c(189,163)]
NG.ERE.SRE.DIST.INTERSECT[ 13,c(189,163)]
MG.ERE.SRE.DIST.INTERSECT[ 13,c(189,163)]
PG.ERE.SRE.DIST.INTERSECT[ 13,c(189,163)]

all_shortest_paths(MG.ACT.NET,"HGSC","PASM")[[1]]
all_shortest_paths(PG.ACT.NET,"HGSC","PASM")[[1]]
all_shortest_paths(MG.ACT.NET,"HGSC","NASM")[[1]]
all_shortest_paths(PG.ACT.NET,"HGSC","NASM")[[1]]

EXAMPLE.SEQS <- c("NGSC","NGSS","NGSR","NGST","NGSM",
                  "NASC","NASS","NASR","NAST","NASM",
                  "HGSC","HGSS","HGSR","HGST","HGSM",
                  "HASC","HASS","HASR","HAST","HASM",
                  "PGSC","PGSS","PGSR","PGST","PGSM",
                  "PASC","PASS","PASR","PAST","PASM")

##Sets of ERE and SRE specific genotypes different at all 4 sites
#AA differences
TG.AA.DIFFERENCES <- adist(TG.ERE.LIST[TG.ERE.INTERSECT],TG.SRE.LIST[TG.SRE.INTERSECT])
rownames(TG.AA.DIFFERENCES) <- TG.ERE.LIST[TG.ERE.INTERSECT]
colnames(TG.AA.DIFFERENCES) <- TG.SRE.LIST[TG.SRE.INTERSECT]

TEST.ERE.E <- which(DT.JOINT$AAseq %in% paste0("E",TG.ERE.LIST[TG.ERE.INTERSECT]))
TEST.ERE.S <- which(DT.JOINT$AAseq %in% paste0("S",TG.ERE.LIST[TG.ERE.INTERSECT]))
SPEC.ERE <- DT.JOINT[TEST.ERE.E,PRED.TE.LINK] - DT.JOINT[TEST.ERE.S,PRED.TE.LINK]

TEST.SRE.E <- which(DT.JOINT$AAseq %in% paste0("E",TG.SRE.LIST[TG.SRE.INTERSECT]))
TEST.SRE.S <- which(DT.JOINT$AAseq %in% paste0("S",TG.SRE.LIST[TG.SRE.INTERSECT]))
SPEC.SRE <- DT.JOINT[TEST.SRE.E,PRED.TE.LINK] - DT.JOINT[TEST.SRE.S,PRED.TE.LINK]

TG.MAX.AA.DIFFERENCES <- which(TG.AA.DIFFERENCES == 4, arr.ind = TRUE)

I <- 1:nrow(TG.MAX.AA.DIFFERENCES)
for(i in I) {
  START <- DT.JOINT$AAseq[TEST.ERE.E][TG.MAX.AA.DIFFERENCES[i,1]]
  STOP  <- DT.JOINT$AAseq[TEST.SRE.S][TG.MAX.AA.DIFFERENCES[i,2]]
  
  SPLIT <- rbind(unlist(strsplit(START,""))[2:5],  unlist(strsplit(STOP ,""))[2:5])
  EXPAND <- as.matrix(expand.grid(SPLIT[,1],SPLIT[,2],SPLIT[,3],SPLIT[,4],stringsAsFactors = FALSE)) 
  PASTE <- paste0("E",EXPAND[,1],EXPAND[,2],EXPAND[,3],EXPAND[,4],sep="")
  
  DT.JOINT[DT.JOINT$AAseq %in% PASTE,]$TE.JOINT.CLASS
}  


######################################
##8. One step functional transitions##
######################################
###TE Network
##Find single step functional transitions
TE.ONE.STEP <- which(TH.ERE.SRE.DIST==1,arr.ind=TRUE)
TE.ONE.STEP.SOURCE <- rownames(TH.ERE.SRE.DIST)[TE.ONE.STEP[,1]]
TE.ONE.STEP.TARGET <- colnames(TH.ERE.SRE.DIST)[TE.ONE.STEP[,2]]

##Thresholds
TE.ERE.THRESH <- c(TE.THRESH.NULL,TE.THRESH.WEAK) - TE.COEFS.EFFECT.ADJ.B0 - TE.COEFS.EFFECT.ADJ.S0
TE.SRE.THRESH <- c(TE.THRESH.NULL,TE.THRESH.WEAK) - TE.COEFS.EFFECT.ADJ.B0 + TE.COEFS.EFFECT.ADJ.S0
TE.THRESH.DIFF <- (TE.ERE.THRESH-TE.SRE.THRESH)[1]

###Decompose contributions
#I <- 1:length(TE.ONE.STEP.SOURCE)
#TE.OUT <- matrix(0,nrow=length(I),ncol=10)
#colnames(TE.OUT) <- c("ERE_SOURCE","ERE_TARGET","ERE_MAIN","ERE_EPI","ERE_TE","SRE_SOURCE","SRE_TARGET","SRE_MAIN","SRE_EPI","SRE_TE")
#for(i in I) {
#  SOURCE <- TE.ONE.STEP.SOURCE[i]
#  E.SOURCE <- TE.COEFS.ADJ[which(TE.COEF.MATRIX.ADJ[which(rownames(TE.COEF.MATRIX.ADJ)==paste0("E",SOURCE)),] !=0)]
#  S.SOURCE <- TE.COEFS.ADJ[which(TE.COEF.MATRIX.ADJ[which(rownames(TE.COEF.MATRIX.ADJ)==paste0("S",SOURCE)),] !=0)]
#  E.SOURCE.DIRECTION <- TE.MATRIX[which(rownames(TE.MATRIX) == paste0("E",SOURCE)),which(TE.MATRIX[which(rownames(TE.MATRIX)==paste0("E",SOURCE)),] !=0)][-29]
#  S.SOURCE.DIRECTION <- TE.MATRIX[which(rownames(TE.MATRIX) == paste0("S",SOURCE)),which(TE.MATRIX[which(rownames(TE.MATRIX)==paste0("S",SOURCE)),] !=0)][-29]
#  
#  TARGET <- TE.ONE.STEP.TARGET[i]
#  E.TARGET <- TE.COEFS.ADJ[which(TE.COEF.MATRIX.ADJ[which(rownames(TE.COEF.MATRIX.ADJ)==paste0("E",TARGET)),] !=0)]
#  S.TARGET <- TE.COEFS.ADJ[which(TE.COEF.MATRIX.ADJ[which(rownames(TE.COEF.MATRIX.ADJ)==paste0("S",TARGET)),] !=0)]
#  E.TARGET.DIRECTION <- TE.MATRIX[which(rownames(TE.MATRIX) == paste0("E",TARGET)),which(TE.MATRIX[which(rownames(TE.MATRIX)==paste0("E",TARGET)),] !=0)][-29]
#  S.TARGET.DIRECTION <- TE.MATRIX[which(rownames(TE.MATRIX) == paste0("S",TARGET)),which(TE.MATRIX[which(rownames(TE.MATRIX)==paste0("S",TARGET)),] !=0)][-29]
#
#  x <- ((E.TARGET.DIRECTION*E.TARGET)[c(1:4,9:14,21:24)] + (E.TARGET.DIRECTION*E.TARGET)[c(5:8,15:20,25:28)]) - ((E.SOURCE.DIRECTION*E.SOURCE)[c(1:4,9:14,21:24)] + (E.SOURCE.DIRECTION*E.SOURCE)[c(5:8,15:20,25:28)])
#  y <- ((E.TARGET.DIRECTION*E.TARGET)[c(1:4,9:14,21:24)] - (E.TARGET.DIRECTION*E.TARGET)[c(5:8,15:20,25:28)]) - ((E.SOURCE.DIRECTION*E.SOURCE)[c(1:4,9:14,21:24)] - (E.SOURCE.DIRECTION*E.SOURCE)[c(5:8,15:20,25:28)])
#  
#  TE.OUT[i,1] <- sum(E.SOURCE.DIRECTION*E.SOURCE)
#  TE.OUT[i,2] <- sum(E.TARGET.DIRECTION*E.TARGET)
#  TE.OUT[i,3] <- sum(x[1:4])
#  TE.OUT[i,4] <- sum(x[5:10])
#  TE.OUT[i,5] <- sum(x[11:14])
#  TE.OUT[i,6] <- sum(S.SOURCE.DIRECTION*S.SOURCE)
#  TE.OUT[i,7] <- sum(S.TARGET.DIRECTION*S.TARGET)
#  TE.OUT[i,8] <- sum(y[1:4])
#  TE.OUT[i,9] <- sum(y[5:10])
#  TE.OUT[i,10]<- sum(y[11:14])
#}
#save(TE.OUT,file="TE.ONE.STEP.rda")
load("TE.ONE.STEP.rda")

#Require change in ERE and SRE binding to bridge the strong to null gap size
TE.ONE.STEP.SOURCE <- TE.ONE.STEP.SOURCE[which(TE.OUT[,1] - TE.OUT[,2] < diff(TE.ERE.THRESH) & TE.OUT[,5] - TE.OUT[,6] > -1*diff(TE.SRE.THRESH))]
TE.ONE.STEP.TARGET <- TE.ONE.STEP.TARGET[which(TE.OUT[,1] - TE.OUT[,2] < diff(TE.ERE.THRESH) & TE.OUT[,5] - TE.OUT[,6] > -1*diff(TE.SRE.THRESH))]
TE.OUT <- TE.OUT[which(TE.OUT[,1] - TE.OUT[,2] < diff(TE.ERE.THRESH) & TE.OUT[,5] - TE.OUT[,6] > -1*diff(TE.SRE.THRESH)),]

#Is epistasis sufficient for one step changes?
TE.EPI.ERE.SUFFICIENT.ALL <- length(which(-1*(TE.OUT[,1] + TE.OUT[,4] + TE.OUT[,5]) < -1*TE.ERE.THRESH[1]))/nrow(TE.OUT)
TE.EPI.ERE.SUFFICIENT.P   <- length(which(-1*(TE.OUT[,1] + TE.OUT[,4]) < -1*TE.ERE.THRESH[1]))/nrow(TE.OUT)
TE.EPI.ERE.SUFFICIENT.T   <- length(which(-1*(TE.OUT[,1] + TE.OUT[,5]) < -1*TE.ERE.THRESH[1]))/nrow(TE.OUT)

TE.EPI.SRE.SUFFICIENT.ALL <- length(which(-1*(TE.OUT[,6] + TE.OUT[,9] + TE.OUT[,10]) > -1*TE.SRE.THRESH[2]))/nrow(TE.OUT)
TE.EPI.SRE.SUFFICIENT.P   <- length(which(-1*(TE.OUT[,6] + TE.OUT[,9]) > -1*TE.SRE.THRESH[2]))/nrow(TE.OUT)
TE.EPI.SRE.SUFFICIENT.T   <- length(which(-1*(TE.OUT[,6] + TE.OUT[,10]) > -1*TE.SRE.THRESH[2]))/nrow(TE.OUT)

TE.EPI.BOTH.SUFFICIENT.ALL <- length(which(-1*(TE.OUT[,1] + TE.OUT[,4] + TE.OUT[,5]) < -1*TE.ERE.THRESH[1] & -1*(TE.OUT[,6] + TE.OUT[,9] + TE.OUT[,10]) > -1*TE.SRE.THRESH[2]))/nrow(TE.OUT)
TE.EPI.BOTH.SUFFICIENT.P   <- length(which(-1*(TE.OUT[,1] + TE.OUT[,4]) < -1*TE.ERE.THRESH[1] & -1*(TE.OUT[,6] + TE.OUT[,9]) > -1*TE.SRE.THRESH[2]))/nrow(TE.OUT)
TE.EPI.BOTH.SUFFICIENT.T   <- length(which(-1*(TE.OUT[,1] + TE.OUT[,5]) < -1*TE.ERE.THRESH[1] & -1*(TE.OUT[,6] + TE.OUT[,10]) > -1*TE.SRE.THRESH[2]))/nrow(TE.OUT)

#Is epistasis necessary for one step changes?
TE.EPI.ERE.NECESSARY.ALL <- length(which(-1*(TE.OUT[,1] + TE.OUT[,3]) > -1*TE.ERE.THRESH[1]))/nrow(TE.OUT)
TE.EPI.ERE.NECESSARY.P   <- length(which(-1*(TE.OUT[,1] + TE.OUT[,3] + TE.OUT[,5]) > -1*TE.ERE.THRESH[1]))/nrow(TE.OUT)
TE.EPI.ERE.NECESSARY.T   <- length(which(-1*(TE.OUT[,1] + TE.OUT[,3] + TE.OUT[,4]) > -1*TE.ERE.THRESH[1]))/nrow(TE.OUT)

TE.EPI.SRE.NECESSARY.ALL <- length(which(-1*(TE.OUT[,6] + TE.OUT[,8]) < -1*TE.SRE.THRESH[2]))/nrow(TE.OUT)
TE.EPI.SRE.NECESSARY.P   <- length(which(-1*(TE.OUT[,6] + TE.OUT[,8] + TE.OUT[,10]) < -1*TE.SRE.THRESH[2]))/nrow(TE.OUT)
TE.EPI.SRE.NECESSARY.T   <- length(which(-1*(TE.OUT[,6] + TE.OUT[,8] + TE.OUT[,9]) < -1*TE.SRE.THRESH[2]))/nrow(TE.OUT)

TE.EPI.BOTH.NECESSARY.ALL <- length(which(-1*(TE.OUT[,1] + TE.OUT[,3]) > -1*TE.ERE.THRESH[1] | -1*(TE.OUT[,6] + TE.OUT[,8]) < -1*TE.SRE.THRESH[2]))/nrow(TE.OUT)
TE.EPI.BOTH.NECESSARY.P   <- length(which(-1*(TE.OUT[,1] + TE.OUT[,3] + TE.OUT[,5]) > -1*TE.ERE.THRESH[1] | -1*(TE.OUT[,6] + TE.OUT[,8] + TE.OUT[,10]) < -1*TE.SRE.THRESH[2]))/nrow(TE.OUT)
TE.EPI.BOTH.NECESSARY.T   <- length(which(-1*(TE.OUT[,1] + TE.OUT[,3] + TE.OUT[,4]) > -1*TE.ERE.THRESH[1] | -1*(TE.OUT[,6] + TE.OUT[,8] + TE.OUT[,9]) < -1*TE.SRE.THRESH[2]))/nrow(TE.OUT)

#Identify types of changes
TE.ONE.STEP.TYPES <- data.frame(SITE=rep(1,nrow(TE.OUT),SOURCE=rep("A",nrow(TE.OUT)),TARGET=rep("A",nrow(TE.OUT))),stringsAsFactors = FALSE)

I <- 1:nrow(TE.OUT)
for(i in I) {
  TE.ONE.STEP.TYPES$SITE[i]   <- which(strsplit(TE.ONE.STEP.TARGET[i],"")[[1]] != strsplit(TE.ONE.STEP.SOURCE[i],"")[[1]])
  
  TE.ONE.STEP.TYPES$SOURCE[i] <- strsplit(TE.ONE.STEP.SOURCE[i],"")[[1]][TE.ONE.STEP.TYPES$SITE[i]]
  TE.ONE.STEP.TYPES$TARGET[i] <- strsplit(TE.ONE.STEP.TARGET[i],"")[[1]][TE.ONE.STEP.TYPES$SITE[i]]
}
unique(TE.ONE.STEP.TYPES)

###PE Network
##Find single step functional transitions
PE.ONE.STEP <- which(PH.ERE.SRE.DIST==1,arr.ind=TRUE)
PE.ONE.STEP.SOURCE <- rownames(PH.ERE.SRE.DIST)[PE.ONE.STEP[,1]]
PE.ONE.STEP.TARGET <- colnames(PH.ERE.SRE.DIST)[PE.ONE.STEP[,2]]

PE.MATRIX.NO.S0 <- PE.MATRIX[,-ncol(PE.MATRIX)]
PE.COEF.MATRIX.ADJ  <- t(t(PE.MATRIX.NO.S0)*PE.COEFS.ADJ)

#Thresholds
PE.ERE.THRESH <- c(PE.THRESH.NULL,PE.THRESH.WEAK) - PE.COEFS.EFFECT.ADJ.B0 - PE.COEFS.EFFECT.ADJ.S0
PE.SRE.THRESH <- c(PE.THRESH.NULL,PE.THRESH.WEAK) - PE.COEFS.EFFECT.ADJ.B0 + PE.COEFS.EFFECT.ADJ.S0
PE.THRESH.DIFF <- (PE.ERE.THRESH-PE.SRE.THRESH)[1]

##Decompose contributions
I <- 1:length(PE.ONE.STEP.SOURCE)
PE.OUT <- matrix(0,nrow=length(I),ncol=8)
colnames(PE.OUT) <- c("ERE_SOURCE","ERE_TARGET","ERE_MAIN","ERE_EPI","SRE_SOURCE","SRE_TARGET","SRE_MAIN","SRE_EPI" )
for(i in I) {
  SOURCE <- PE.ONE.STEP.SOURCE[i]
  E.SOURCE <- PE.COEFS.ADJ[which(PE.COEF.MATRIX.ADJ[which(rownames(PE.COEF.MATRIX.ADJ)==paste0("E",SOURCE)),] !=0)]
  S.SOURCE <- PE.COEFS.ADJ[which(PE.COEF.MATRIX.ADJ[which(rownames(PE.COEF.MATRIX.ADJ)==paste0("S",SOURCE)),] !=0)]
  E.SOURCE.DIRECTION <- PE.MATRIX.NO.S0[which(rownames(PE.MATRIX.NO.S0) == paste0("E",SOURCE)),which(PE.MATRIX.NO.S0[which(rownames(PE.MATRIX.NO.S0)==paste0("E",SOURCE)),] !=0)]
  S.SOURCE.DIRECTION <- PE.MATRIX.NO.S0[which(rownames(PE.MATRIX.NO.S0) == paste0("S",SOURCE)),which(PE.MATRIX.NO.S0[which(rownames(PE.MATRIX.NO.S0)==paste0("S",SOURCE)),] !=0)]
  
  TARGET <- PE.ONE.STEP.TARGET[i]
  E.TARGET <- PE.COEFS.ADJ[which(PE.COEF.MATRIX.ADJ[which(rownames(PE.COEF.MATRIX.ADJ)==paste0("E",TARGET)),] !=0)]
  S.TARGET <- PE.COEFS.ADJ[which(PE.COEF.MATRIX.ADJ[which(rownames(PE.COEF.MATRIX.ADJ)==paste0("S",TARGET)),] !=0)]
  E.TARGET.DIRECTION <- PE.MATRIX.NO.S0[which(rownames(PE.MATRIX.NO.S0) == paste0("E",TARGET)),which(PE.MATRIX.NO.S0[which(rownames(PE.MATRIX.NO.S0)==paste0("E",TARGET)),] !=0)]
  S.TARGET.DIRECTION <- PE.MATRIX.NO.S0[which(rownames(PE.MATRIX.NO.S0) == paste0("S",TARGET)),which(PE.MATRIX.NO.S0[which(rownames(PE.MATRIX.NO.S0)==paste0("S",TARGET)),] !=0)]
  
  x <- ((E.TARGET.DIRECTION*E.TARGET)[c(1:4,9:14)] + (E.TARGET.DIRECTION*E.TARGET)[c(5:8,15:20)]) - ((E.SOURCE.DIRECTION*E.SOURCE)[c(1:4,9:14)] + (E.SOURCE.DIRECTION*E.SOURCE)[c(5:8,15:20)])
  y <- ((E.TARGET.DIRECTION*E.TARGET)[c(1:4,9:14)] - (E.TARGET.DIRECTION*E.TARGET)[c(5:8,15:20)]) - ((E.SOURCE.DIRECTION*E.SOURCE)[c(1:4,9:14)] - (E.SOURCE.DIRECTION*E.SOURCE)[c(5:8,15:20)])

  PE.OUT[i,1] <- sum(E.SOURCE.DIRECTION*E.SOURCE)
  PE.OUT[i,2] <- sum(E.TARGET.DIRECTION*E.TARGET)
  PE.OUT[i,3] <- sum(x[1:4])
  PE.OUT[i,4] <- sum(x[5:10])
  PE.OUT[i,5] <- sum(S.SOURCE.DIRECTION*S.SOURCE)
  PE.OUT[i,6] <- sum(S.TARGET.DIRECTION*S.TARGET)
  PE.OUT[i,7] <- sum(y[1:4])
  PE.OUT[i,8] <- sum(y[5:10])
}
save(PE.OUT,file="PE.ONE.STEP.rda")
load("PE.ONE.STEP.rda")

#Require change in ERE and SRE binding to bridge the strong to null gap size
PE.ONE.STEP.SOURCE <- PE.ONE.STEP.SOURCE[which(PE.OUT[,1] - PE.OUT[,2] < diff(PE.ERE.THRESH) & PE.OUT[,5] - PE.OUT[,6] > -1*diff(PE.SRE.THRESH))]
PE.ONE.STEP.TARGET <- PE.ONE.STEP.TARGET[which(PE.OUT[,1] - PE.OUT[,2] < diff(PE.ERE.THRESH) & PE.OUT[,5] - PE.OUT[,6] > -1*diff(PE.SRE.THRESH))]
PE.OUT <- PE.OUT[which(PE.OUT[,1] - PE.OUT[,2] < diff(PE.ERE.THRESH) & PE.OUT[,5] - PE.OUT[,6] > -1*diff(PE.SRE.THRESH)),]

#Is epistasis sufficient for one step changes?
PE.EPI.ERE.SUFFICIENT.ALL  <- length(which(-1*(PE.OUT[,1] + PE.OUT[,4]) < -1*PE.ERE.THRESH[1]))/nrow(PE.OUT)
PE.EPI.SRE.SUFFICIENT.ALL  <- length(which(-1*(PE.OUT[,5] + PE.OUT[,8]) > -1*PE.SRE.THRESH[2]))/nrow(PE.OUT)
PE.EPI.BOTH.SUFFICIENT.ALL <- length(which(-1*(PE.OUT[,1] + PE.OUT[,4]) < -1*PE.ERE.THRESH[1] & -1*(PE.OUT[,5] + PE.OUT[,8]) > -1*PE.SRE.THRESH[2]))/nrow(PE.OUT)

#Is epistasis necssary for one step changes?
PE.EPI.ERE.NECESSARY.ALL  <- length(which(-1*(PE.OUT[,1] + PE.OUT[,3]) > -1*PE.ERE.THRESH[1]))/nrow(PE.OUT)
PE.EPI.SRE.NECESSARY.ALL  <- length(which(-1*(PE.OUT[,5] + PE.OUT[,7]) < -1*PE.SRE.THRESH[2]))/nrow(PE.OUT)
PE.EPI.BOTH.NECESSARY.ALL <- length(which(-1*(PE.OUT[,1] + PE.OUT[,3]) > -1*PE.ERE.THRESH[1] | -1*(PE.OUT[,5] + PE.OUT[,7]) < -1*PE.SRE.THRESH[2]))/nrow(PE.OUT)


#pdf("EPI.NECESSARY.SUFFICIENT.pdf")
par(mfrow=c(2,2))
barplot(c(PE.EPI.ERE.NECESSARY.ALL, PE.EPI.SRE.NECESSARY.ALL, PE.EPI.BOTH.NECESSARY.ALL), ylim=c(0,1),main="PE.EPI.NECESSARY.ALL", names.arg=c("ERE","SRE","Both"))
barplot(c(PE.EPI.ERE.SUFFICIENT.ALL,PE.EPI.SRE.SUFFICIENT.ALL,PE.EPI.BOTH.SUFFICIENT.ALL),ylim=c(0,1),main="PE.EPI.SUFFICIENT.ALL",names.arg=c("ERE","SRE","Both"))
barplot(c(TE.EPI.ERE.NECESSARY.ALL, TE.EPI.SRE.NECESSARY.ALL, TE.EPI.BOTH.NECESSARY.ALL), ylim=c(0,1),main="TE.EPI.NECESSARY.ALL", names.arg=c("ERE","SRE","Both"))
barplot(c(TE.EPI.ERE.SUFFICIENT.ALL,TE.EPI.SRE.SUFFICIENT.ALL,TE.EPI.BOTH.SUFFICIENT.ALL),ylim=c(0,1),main="TE.EPI.SUFFICIENT.ALL",names.arg=c("ERE","SRE","Both"))
barplot(c(TE.EPI.ERE.NECESSARY.P,   TE.EPI.SRE.NECESSARY.P,   TE.EPI.BOTH.NECESSARY.P),   ylim=c(0,1),main="TE.EPI.NECESSARY.P",   names.arg=c("ERE","SRE","Both"))
barplot(c(TE.EPI.ERE.SUFFICIENT.P,  TE.EPI.SRE.SUFFICIENT.P,  TE.EPI.BOTH.SUFFICIENT.P),  ylim=c(0,1),main="TE.EPI.SUFFICIENT.P",  names.arg=c("ERE","SRE","Both"))
barplot(c(TE.EPI.ERE.NECESSARY.T,   TE.EPI.SRE.NECESSARY.T,   TE.EPI.BOTH.NECESSARY.T),   ylim=c(0,1),main="TE.EPI.NECESSARY.T",   names.arg=c("ERE","SRE","Both"))
barplot(c(TE.EPI.ERE.SUFFICIENT.T,  TE.EPI.SRE.SUFFICIENT.T,  TE.EPI.BOTH.SUFFICIENT.T),  ylim=c(0,1),main="TE.EPI.SUFFICIENT.T",  names.arg=c("ERE","SRE","Both"))
#dev.off()

#Identify types of changes
PE.ONE.STEP.TYPES <- data.frame(SITE=rep(1,nrow(PE.OUT),SOURCE=rep("A",nrow(PE.OUT)),TARGET=rep("A",nrow(PE.OUT))),stringsAsFactors = FALSE)

I <- 1:nrow(PE.OUT)
for(i in I) {
  PE.ONE.STEP.TYPES$SITE[i]   <- which(strsplit(PE.ONE.STEP.TARGET[i],"")[[1]] != strsplit(PE.ONE.STEP.SOURCE[i],"")[[1]])
  
  PE.ONE.STEP.TYPES$SOURCE[i] <- strsplit(PE.ONE.STEP.SOURCE[i],"")[[1]][PE.ONE.STEP.TYPES$SITE[i]]
  PE.ONE.STEP.TYPES$TARGET[i] <- strsplit(PE.ONE.STEP.TARGET[i],"")[[1]][PE.ONE.STEP.TYPES$SITE[i]]
}
unique(PE.ONE.STEP.TYPES)

###ME Network
##Find single step functional transitions
ME.ONE.STEP <- which(MH.ERE.SRE.DIST==1,arr.ind=TRUE)
ME.ONE.STEP.SOURCE <- rownames(MH.ERE.SRE.DIST)[ME.ONE.STEP[,1]]
ME.ONE.STEP.TARGET <- colnames(MH.ERE.SRE.DIST)[ME.ONE.STEP[,2]]

ME.MATRIX.NO.S0 <- ME.MATRIX[,-ncol(ME.MATRIX)]
ME.COEF.MATRIX.ADJ    <- t(t(ME.MATRIX.NO.S0)*ME.COEFS.ADJ)

#Thresholds
ME.ERE.THRESH <- c(ME.THRESH.NULL,ME.THRESH.WEAK) - ME.COEFS.EFFECT.ADJ.B0 - ME.COEFS.EFFECT.ADJ.S0
ME.SRE.THRESH <- c(ME.THRESH.NULL,ME.THRESH.WEAK) - ME.COEFS.EFFECT.ADJ.B0 + ME.COEFS.EFFECT.ADJ.S0
ME.THRESH.DIFF <- (ME.ERE.THRESH-ME.SRE.THRESH)[1]

###Decompose contributions
#I <- 1:length(ME.ONE.STEP.SOURCE)
#ME.OUT <- matrix(0,nrow=length(I),ncol=8)
#colnames(ME.OUT) <- c("ERE_SOURCE","ERE_TARGET","ERE_MAIN","ERE_EPI","SRE_SOURCE","SRE_TARGET","SRE_MAIN","SRE_EPI" )
#for(i in I) {
#  SOURCE <- ME.ONE.STEP.SOURCE[i]
#  E.SOURCE <- ME.COEFS.ADJ[which(ME.COEF.MATRIX.ADJ[which(rownames(ME.COEF.MATRIX.ADJ)==paste0("E",SOURCE)),] !=0)]
#  S.SOURCE <- ME.COEFS.ADJ[which(ME.COEF.MATRIX.ADJ[which(rownames(ME.COEF.MATRIX.ADJ)==paste0("S",SOURCE)),] !=0)]
#  E.SOURCE.DIRECTION <- ME.MATRIX.NO.S0[which(rownames(ME.MATRIX.NO.S0) == paste0("E",SOURCE)),which(ME.MATRIX.NO.S0[which(rownames(ME.MATRIX.NO.S0)==paste0("E",SOURCE)),] !=0)]
#  S.SOURCE.DIRECTION <- ME.MATRIX.NO.S0[which(rownames(ME.MATRIX.NO.S0) == paste0("S",SOURCE)),which(ME.MATRIX.NO.S0[which(rownames(ME.MATRIX.NO.S0)==paste0("S",SOURCE)),] !=0)]
#  
#  TARGET <- ME.ONE.STEP.TARGET[i]
#  E.TARGET <- ME.COEFS.ADJ[which(ME.COEF.MATRIX.ADJ[which(rownames(ME.COEF.MATRIX.ADJ)==paste0("E",TARGET)),] !=0)]
#  S.TARGET <- ME.COEFS.ADJ[which(ME.COEF.MATRIX.ADJ[which(rownames(ME.COEF.MATRIX.ADJ)==paste0("S",TARGET)),] !=0)]
#  E.TARGET.DIRECTION <- ME.MATRIX.NO.S0[which(rownames(ME.MATRIX.NO.S0) == paste0("E",TARGET)),which(ME.MATRIX.NO.S0[which(rownames(ME.MATRIX.NO.S0)==paste0("E",TARGET)),] !=0)]
#  S.TARGET.DIRECTION <- ME.MATRIX.NO.S0[which(rownames(ME.MATRIX.NO.S0) == paste0("S",TARGET)),which(ME.MATRIX.NO.S0[which(rownames(ME.MATRIX.NO.S0)==paste0("S",TARGET)),] !=0)]
#  
#  x <- ((E.TARGET.DIRECTION*E.TARGET)[1:4] + (E.TARGET.DIRECTION*E.TARGET)[5:8]) - ((E.SOURCE.DIRECTION*E.SOURCE)[1:4] + (E.SOURCE.DIRECTION*E.SOURCE)[5:8])
#  y <- ((E.TARGET.DIRECTION*E.TARGET)[1:4] - (E.TARGET.DIRECTION*E.TARGET)[5:8]) - ((E.SOURCE.DIRECTION*E.SOURCE)[1:4] - (E.SOURCE.DIRECTION*E.SOURCE)[5:8])
#  
#  ME.OUT[i,1] <- sum(E.SOURCE.DIRECTION*E.SOURCE)
#  ME.OUT[i,2] <- sum(E.TARGET.DIRECTION*E.TARGET)
#  ME.OUT[i,3] <- sum(x[1:4])
#  ME.OUT[i,4] <- 0
#  ME.OUT[i,5] <- sum(S.SOURCE.DIRECTION*S.SOURCE)
#  ME.OUT[i,6] <- sum(S.TARGET.DIRECTION*S.TARGET)
#  ME.OUT[i,7] <- sum(y[1:4])
#  ME.OUT[i,8] <- 0
#}
#save(ME.OUT,file="ME.ONE.STEP.rda")
load("ME.ONE.STEP.rda")

#Require change in ERE and SRE binding to bridge the strong to null gap size
ME.ONE.STEP.SOURCE <- ME.ONE.STEP.SOURCE[which(ME.OUT[,1] - ME.OUT[,2] < diff(ME.ERE.THRESH) & ME.OUT[,5] - ME.OUT[,6] > -1*diff(ME.SRE.THRESH))]
ME.ONE.STEP.TARGET <- ME.ONE.STEP.TARGET[which(ME.OUT[,1] - ME.OUT[,2] < diff(ME.ERE.THRESH) & ME.OUT[,5] - ME.OUT[,6] > -1*diff(ME.SRE.THRESH))]
ME.OUT <- ME.OUT[which(ME.OUT[,1] - ME.OUT[,2] < diff(ME.ERE.THRESH) & ME.OUT[,5] - ME.OUT[,6] > -1*diff(ME.SRE.THRESH)),]

##Identify types of changes
ME.ONE.STEP.TYPES <- data.frame(SITE=rep(1,nrow(ME.OUT),SOURCE=rep("A",nrow(ME.OUT)),TARGET=rep("A",nrow(ME.OUT))),stringsAsFactors = FALSE)

I <- 1:nrow(ME.OUT)
for(i in I) {
  ME.ONE.STEP.TYPES$SITE[i]   <- which(strsplit(ME.ONE.STEP.TARGET[i],"")[[1]] != strsplit(ME.ONE.STEP.SOURCE[i],"")[[1]])
  
  ME.ONE.STEP.TYPES$SOURCE[i] <- strsplit(ME.ONE.STEP.SOURCE[i],"")[[1]][ME.ONE.STEP.TYPES$SITE[i]]
  ME.ONE.STEP.TYPES$TARGET[i] <- strsplit(ME.ONE.STEP.TARGET[i],"")[[1]][ME.ONE.STEP.TYPES$SITE[i]]
}
unique(ME.ONE.STEP.TYPES)

###############################
##9. Path Number and Identity##
###############################
###Number of paths in each network from each ERE-specific to each SRE-specific
#TG.NUM.PATHS <- t(num.paths(TG.ACT.NET,TG.ERE.LIST,TG.SRE.LIST))  
#PG.NUM.PATHS <- t(num.paths(PG.ACT.NET,PG.ERE.LIST,PG.SRE.LIST))  
#MG.NUM.PATHS <- t(num.paths(MG.ACT.NET,MG.ERE.LIST,MG.SRE.LIST))  
#NG.NUM.PATHS <- t(num.paths.GC(NG.ERE.LIST,NG.SRE.LIST,CODE = "Z"))  
#TH.NUM.PATHS <- t(num.paths(TH.ACT.NET,TH.ERE.LIST,TH.SRE.LIST))  
#PH.NUM.PATHS <- t(num.paths(PH.ACT.NET,PH.ERE.LIST,PH.SRE.LIST))  
#MH.NUM.PATHS <- t(num.paths(MH.ACT.NET,MH.ERE.LIST,MH.SRE.LIST))  
#NH.NUM.PATHS <- t(num.paths.HD(NH.ERE.LIST,NH.SRE.LIST))  

#save(TG.NUM.PATHS,file="TG.NUM.PATHS.rda")
#save(PG.NUM.PATHS,file="PG.NUM.PATHS.rda")
#save(MG.NUM.PATHS,file="MG.NUM.PATHS.rda")
#save(NG.NUM.PATHS,file="NG.NUM.PATHS.rda")
#save(TH.NUM.PATHS,file="TH.NUM.PATHS.rda")
#save(PH.NUM.PATHS,file="PH.NUM.PATHS.rda")
#save(MH.NUM.PATHS,file="MH.NUM.PATHS.rda")
#save(NH.NUM.PATHS,file="NH.NUM.PATHS.rda")

load("TG.NUM.PATHS.rda")
load("PG.NUM.PATHS.rda")
load("MG.NUM.PATHS.rda")
load("NG.NUM.PATHS.rda")
load("TH.NUM.PATHS.rda")
load("PH.NUM.PATHS.rda")
load("MH.NUM.PATHS.rda")
load("NH.NUM.PATHS.rda")

##Collect mean values
ALL.NUM.PATHS.AVG <- c(mean(TG.NUM.PATHS),mean(PG.NUM.PATHS),mean(MG.NUM.PATHS),mean(NG.NUM.PATHS),mean(TH.NUM.PATHS),mean(PH.NUM.PATHS),mean(MH.NUM.PATHS),mean(NH.NUM.PATHS))
INT.NUM.PATHS.AVG <- c(mean(TG.NUM.PATHS[TG.ERE.INTERSECT,TG.SRE.INTERSECT]),mean(PG.NUM.PATHS[PG.ERE.INTERSECT,PG.SRE.INTERSECT]),mean(MG.NUM.PATHS[MG.ERE.INTERSECT,MG.SRE.INTERSECT]),mean(NG.NUM.PATHS[NG.ERE.INTERSECT,NG.SRE.INTERSECT]),mean(TH.NUM.PATHS[TH.ERE.INTERSECT,TH.SRE.INTERSECT]),mean(PH.NUM.PATHS[PH.ERE.INTERSECT,PH.SRE.INTERSECT]),mean(MH.NUM.PATHS[MH.ERE.INTERSECT,MH.SRE.INTERSECT]),mean(NH.NUM.PATHS[NH.ERE.INTERSECT,NH.SRE.INTERSECT]))

#pdf("NUM.PATHS.pdf")
par(mfrow=c(1,1))

vioplot(c(log10(TG.NUM.PATHS[TG.NUM.PATHS != 0])),c(log10(PG.NUM.PATHS[PG.NUM.PATHS != 0])),c(log10(MG.NUM.PATHS[MG.NUM.PATHS != 0])),c(log10(NG.NUM.PATHS[NG.NUM.PATHS != 0])),names=c("3G","2G","1G","NG"),cex.axis=0.75,cex.main=0.5,col="darkred", rectCol="white",lineCol="white",pchMed=19,colMed="black",main="All Number of Paths")
vioplot(c(log10(TH.NUM.PATHS[TH.NUM.PATHS != 0])),c(log10(PH.NUM.PATHS[PH.NUM.PATHS != 0])),c(log10(MH.NUM.PATHS[MH.NUM.PATHS != 0])),c(log10(NH.NUM.PATHS[NH.NUM.PATHS != 0])),names=c("3H","2H","1H","NH"),cex.axis=0.75,cex.main=0.5,col="skyblue2",rectCol="white",lineCol="white",pchMed=19,colMed="black",main="All Number of Paths")

barplot(log10(ALL.NUM.PATHS.AVG[1:4]),names=c("3G","2G","1G","NG"),cex.axis=0.75,cex.main=0.5,col="darkred" ,main="All Number of Paths")
barplot(log10(ALL.NUM.PATHS.AVG[5:8]),names=c("3H","2H","1H","NH"),cex.axis=0.75,cex.main=0.5,col="skyblue2",main="All Number of Paths")

vioplot(c(log10(TG.NUM.PATHS[TG.ERE.INTERSECT,TG.SRE.INTERSECT])),c(log10(PG.NUM.PATHS[PG.ERE.INTERSECT,PG.SRE.INTERSECT])),c(log10(MG.NUM.PATHS[MG.ERE.INTERSECT,MG.SRE.INTERSECT])),c(log10(NG.NUM.PATHS[NG.ERE.INTERSECT,NG.SRE.INTERSECT])),names=c("3G","2G","1G","NG"),cex.axis=0.75,cex.main=0.5,col="darkred", rectCol="white",lineCol="white",pchMed=19,colMed="black",main="Common Number of Paths")
vioplot(c(log10(TH.NUM.PATHS[TH.ERE.INTERSECT,TH.SRE.INTERSECT])),c(log10(PH.NUM.PATHS[PH.ERE.INTERSECT,PH.SRE.INTERSECT])),c(log10(MH.NUM.PATHS[MH.ERE.INTERSECT,MH.SRE.INTERSECT])),c(log10(NH.NUM.PATHS[NH.ERE.INTERSECT,NH.SRE.INTERSECT])),names=c("3H","2H","1H","NH"),cex.axis=0.75,cex.main=0.5,col="skyblue2",rectCol="white",lineCol="white",pchMed=19,colMed="black",main="Common Number of Paths")

barplot(INT.NUM.PATHS.AVG[1:3],names=c("3G","2G","1G"),cex.axis=0.75,cex.main=0.5,ylim=c(0,70),col="darkred" ,main="Common Number of Paths")
barplot(INT.NUM.PATHS.AVG[5:7],names=c("3H","2H","1H"),cex.axis=0.75,cex.main=0.5,col="skyblue2",main="Common Number of Paths")

###Permutation tests
par(mfrow=c(3,2))
#All genotypes
TG.PG.ALL.NUM.PATH.PERMUTATION.TEST  <- permutation.test.matrix(TG.NUM.PATHS,PG.NUM.PATHS,10000)
  hist(TG.PG.ALL.NUM.PATH.PERMUTATION.TEST ,breaks=50,xlim=c(-10,10),xlab="Average Difference number of paths",main="All TG vs PG"); abline(v=mean(TG.NUM.PATHS,na.rm=TRUE) - mean(PG.NUM.PATHS,na.rm=TRUE),col="red")
  sum(mean(TG.NUM.PATHS,na.rm=TRUE) - mean(PG.NUM.PATHS,na.rm=TRUE) < TG.PG.ALL.NUM.PATH.PERMUTATION.TEST )/length(TG.PG.ALL.NUM.PATH.PERMUTATION.TEST )
PG.MG.ALL.NUM.PATH.PERMUTATION.TEST  <- permutation.test.matrix(PG.NUM.PATHS,MG.NUM.PATHS,10000)
  hist(PG.MG.ALL.NUM.PATH.PERMUTATION.TEST ,breaks=50,xlim=c(-10,10),xlab="Average Difference number of paths",main="All PG vs MG"); abline(v=mean(PG.NUM.PATHS,na.rm=TRUE) - mean(MG.NUM.PATHS,na.rm=TRUE),col="red")
  sum(mean(PG.NUM.PATHS,na.rm=TRUE) - mean(MG.NUM.PATHS,na.rm=TRUE) < PG.MG.ALL.NUM.PATH.PERMUTATION.TEST )/length(PG.MG.ALL.NUM.PATH.PERMUTATION.TEST )
MG.NG.ALL.NUM.PATH.PERMUTATION.TEST  <- permutation.test.matrix(MG.NUM.PATHS,NG.NUM.PATHS,10000)
  hist(MG.NG.ALL.NUM.PATH.PERMUTATION.TEST ,breaks=50,xlim=c(-2500,2500),xlab="Average Difference number of paths",main="All MG vs NG"); abline(v=mean(MG.NUM.PATHS,na.rm=TRUE) - mean(NG.NUM.PATHS,na.rm=TRUE),col="red")
  sum(mean(MG.NUM.PATHS,na.rm=TRUE) - mean(NG.NUM.PATHS,na.rm=TRUE) < MG.NG.ALL.NUM.PATH.PERMUTATION.TEST )/length(MG.NG.ALL.NUM.PATH.PERMUTATION.TEST )
TG.MG.ALL.NUM.PATH.PERMUTATION.TEST  <- permutation.test.matrix(TG.NUM.PATHS,MG.NUM.PATHS,10000)
  hist(TG.MG.ALL.NUM.PATH.PERMUTATION.TEST ,breaks=50,xlim=c(-10,10),xlab="Average Difference number of paths",main="All TG vs MG"); abline(v=mean(TG.NUM.PATHS,na.rm=TRUE) - mean(MG.NUM.PATHS,na.rm=TRUE),col="red")
  sum(mean(TG.NUM.PATHS,na.rm=TRUE) - mean(MG.NUM.PATHS,na.rm=TRUE) < TG.MG.ALL.NUM.PATH.PERMUTATION.TEST )/length(TG.MG.ALL.NUM.PATH.PERMUTATION.TEST )
TG.NG.ALL.NUM.PATH.PERMUTATION.TEST  <- permutation.test.matrix(TG.NUM.PATHS,NG.NUM.PATHS,10000)
  hist(TG.NG.ALL.NUM.PATH.PERMUTATION.TEST ,breaks=50,xlim=c(-2500,2500),xlab="Average Difference number of paths",main="All TG vs NG"); abline(v=mean(TG.NUM.PATHS,na.rm=TRUE) - mean(NG.NUM.PATHS,na.rm=TRUE),col="red")
  sum(mean(TG.NUM.PATHS,na.rm=TRUE) - mean(NG.NUM.PATHS,na.rm=TRUE) < TG.NG.ALL.NUM.PATH.PERMUTATION.TEST )/length(TG.NG.ALL.NUM.PATH.PERMUTATION.TEST )
PG.NG.ALL.NUM.PATH.PERMUTATION.TEST  <- permutation.test.matrix(PG.NUM.PATHS,NG.NUM.PATHS,10000)
  hist(PG.NG.ALL.NUM.PATH.PERMUTATION.TEST ,breaks=50,xlim=c(-2500,2500),xlab="Average Difference number of paths",main="All PG vs NG"); abline(v=mean(PG.NUM.PATHS,na.rm=TRUE) - mean(NG.NUM.PATHS,na.rm=TRUE),col="red")
  sum(mean(PG.NUM.PATHS,na.rm=TRUE) - mean(NG.NUM.PATHS,na.rm=TRUE) < PG.NG.ALL.NUM.PATH.PERMUTATION.TEST )/length(PG.NG.ALL.NUM.PATH.PERMUTATION.TEST )
TH.PH.ALL.NUM.PATH.PERMUTATION.TEST  <- permutation.test.matrix(TH.NUM.PATHS,PH.NUM.PATHS,10000)
  hist(TH.PH.ALL.NUM.PATH.PERMUTATION.TEST ,breaks=50,xlim=c(-.7,.7),xlab="Average Difference number of paths",main="All TH vs PH"); abline(v=mean(TH.NUM.PATHS,na.rm=TRUE) - mean(PH.NUM.PATHS,na.rm=TRUE),col="red")
  sum(mean(TH.NUM.PATHS,na.rm=TRUE) - mean(PH.NUM.PATHS,na.rm=TRUE) < TH.PH.ALL.NUM.PATH.PERMUTATION.TEST )/length(TH.PH.ALL.NUM.PATH.PERMUTATION.TEST )
PH.MH.ALL.NUM.PATH.PERMUTATION.TEST  <- permutation.test.matrix(PH.NUM.PATHS,MH.NUM.PATHS,10000)
  hist(PH.MH.ALL.NUM.PATH.PERMUTATION.TEST ,breaks=50,xlim=c(-.7,.7),xlab="Average Difference number of paths",main="All PH vs MH"); abline(v=mean(PH.NUM.PATHS,na.rm=TRUE) - mean(MH.NUM.PATHS,na.rm=TRUE),col="red")
  sum(mean(PH.NUM.PATHS,na.rm=TRUE) - mean(MH.NUM.PATHS,na.rm=TRUE) < PH.MH.ALL.NUM.PATH.PERMUTATION.TEST )/length(PH.MH.ALL.NUM.PATH.PERMUTATION.TEST )
MH.NH.ALL.NUM.PATH.PERMUTATION.TEST  <- permutation.test.matrix(MH.NUM.PATHS,NH.NUM.PATHS,10000)
  hist(MH.NH.ALL.NUM.PATH.PERMUTATION.TEST ,breaks=50,xlim=c(-10,10),xlab="Average Difference number of paths",main="All MH vs NH"); abline(v=mean(MH.NUM.PATHS,na.rm=TRUE) - mean(NH.NUM.PATHS,na.rm=TRUE),col="red")
  sum(mean(MH.NUM.PATHS,na.rm=TRUE) - mean(NH.NUM.PATHS,na.rm=TRUE) < MH.NH.ALL.NUM.PATH.PERMUTATION.TEST )/length(MH.NH.ALL.NUM.PATH.PERMUTATION.TEST )
TH.MH.ALL.NUM.PATH.PERMUTATION.TEST  <- permutation.test.matrix(TH.NUM.PATHS,MH.NUM.PATHS,10000)
  hist(TH.MH.ALL.NUM.PATH.PERMUTATION.TEST ,breaks=50,xlim=c(-.7,.7),xlab="Average Difference number of paths",main="All TH vs MH"); abline(v=mean(TH.NUM.PATHS,na.rm=TRUE) - mean(MH.NUM.PATHS,na.rm=TRUE),col="red")
  sum(mean(TH.NUM.PATHS,na.rm=TRUE) - mean(MH.NUM.PATHS,na.rm=TRUE) < TH.MH.ALL.NUM.PATH.PERMUTATION.TEST )/length(TH.MH.ALL.NUM.PATH.PERMUTATION.TEST )
TH.NH.ALL.NUM.PATH.PERMUTATION.TEST  <- permutation.test.matrix(TH.NUM.PATHS,NH.NUM.PATHS,10000)
  hist(TH.NH.ALL.NUM.PATH.PERMUTATION.TEST ,breaks=50,xlim=c(-10,10),xlab="Average Difference number of paths",main="All TH vs NH"); abline(v=mean(TH.NUM.PATHS,na.rm=TRUE) - mean(NH.NUM.PATHS,na.rm=TRUE),col="red")
  sum(mean(TH.NUM.PATHS,na.rm=TRUE) - mean(NH.NUM.PATHS,na.rm=TRUE) < TH.NH.ALL.NUM.PATH.PERMUTATION.TEST )/length(TH.NH.ALL.NUM.PATH.PERMUTATION.TEST )
PH.NH.ALL.NUM.PATH.PERMUTATION.TEST  <- permutation.test.matrix(PH.NUM.PATHS,NH.NUM.PATHS,10000)
  hist(PH.NH.ALL.NUM.PATH.PERMUTATION.TEST ,breaks=50,xlim=c(-10,10),xlab="Average Difference number of paths",main="All PH vs NH"); abline(v=mean(PH.NUM.PATHS,na.rm=TRUE) - mean(NH.NUM.PATHS,na.rm=TRUE),col="red")
  sum(mean(PH.NUM.PATHS,na.rm=TRUE) - mean(NH.NUM.PATHS,na.rm=TRUE) < PH.NH.ALL.NUM.PATH.PERMUTATION.TEST )/length(PH.NH.ALL.NUM.PATH.PERMUTATION.TEST )

#Common genotypes
TG.PG.INT.NUM.PATH.PERMUTATION.TEST  <- permutation.test.matrix(TG.NUM.PATHS[TG.ERE.INTERSECT,TG.SRE.INTERSECT],PG.NUM.PATHS[PG.ERE.INTERSECT,PG.SRE.INTERSECT],10000)
  hist(TG.PG.INT.NUM.PATH.PERMUTATION.TEST ,breaks=50,xlim=c(-15,15),xlab="Average Difference number of paths",main="Common TG vs PG"); abline(v=mean(TG.NUM.PATHS[TG.ERE.INTERSECT,TG.SRE.INTERSECT],na.rm=TRUE) - mean(PG.NUM.PATHS[PG.ERE.INTERSECT,PG.SRE.INTERSECT],na.rm=TRUE),col="red")
  sum(mean(TG.NUM.PATHS[TG.ERE.INTERSECT,TG.SRE.INTERSECT],na.rm=TRUE) - mean(PG.NUM.PATHS[PG.ERE.INTERSECT,PG.SRE.INTERSECT],na.rm=TRUE) < TG.PG.INT.NUM.PATH.PERMUTATION.TEST )/length(TG.PG.INT.NUM.PATH.PERMUTATION.TEST )
PG.MG.INT.NUM.PATH.PERMUTATION.TEST  <- permutation.test.matrix(PG.NUM.PATHS[PG.ERE.INTERSECT,PG.SRE.INTERSECT],MG.NUM.PATHS[MG.ERE.INTERSECT,MG.SRE.INTERSECT],10000)
  hist(PG.MG.INT.NUM.PATH.PERMUTATION.TEST ,breaks=50,xlim=c(-15,15),xlab="Average Difference number of paths",main="Common PG vs MG"); abline(v=mean(PG.NUM.PATHS[PG.ERE.INTERSECT,PG.SRE.INTERSECT],na.rm=TRUE) - mean(MG.NUM.PATHS[MG.ERE.INTERSECT,MG.SRE.INTERSECT],na.rm=TRUE),col="red")
  sum(mean(PG.NUM.PATHS[PG.ERE.INTERSECT,PG.SRE.INTERSECT],na.rm=TRUE) - mean(MG.NUM.PATHS[MG.ERE.INTERSECT,MG.SRE.INTERSECT],na.rm=TRUE) < PG.MG.INT.NUM.PATH.PERMUTATION.TEST )/length(PG.MG.INT.NUM.PATH.PERMUTATION.TEST )
MG.NG.INT.NUM.PATH.PERMUTATION.TEST  <- permutation.test.matrix(MG.NUM.PATHS[MG.ERE.INTERSECT,MG.SRE.INTERSECT],NG.NUM.PATHS[NG.ERE.INTERSECT,NG.SRE.INTERSECT],10000)
  hist(MG.NG.INT.NUM.PATH.PERMUTATION.TEST ,breaks=50,xlim=c(-2500,2500),xlab="Average Difference number of paths",main="Common MG vs NG"); abline(v=mean(MG.NUM.PATHS[MG.ERE.INTERSECT,MG.SRE.INTERSECT],na.rm=TRUE) - mean(NG.NUM.PATHS[NG.ERE.INTERSECT,NG.SRE.INTERSECT],na.rm=TRUE),col="red")
  sum(mean(MG.NUM.PATHS[MG.ERE.INTERSECT,MG.SRE.INTERSECT],na.rm=TRUE) - mean(NG.NUM.PATHS[NG.ERE.INTERSECT,NG.SRE.INTERSECT],na.rm=TRUE) < MG.NG.INT.NUM.PATH.PERMUTATION.TEST )/length(MG.NG.INT.NUM.PATH.PERMUTATION.TEST )
TG.MG.INT.NUM.PATH.PERMUTATION.TEST  <- permutation.test.matrix(TG.NUM.PATHS[TG.ERE.INTERSECT,TG.SRE.INTERSECT],MG.NUM.PATHS[MG.ERE.INTERSECT,MG.SRE.INTERSECT],10000)
  hist(TG.MG.INT.NUM.PATH.PERMUTATION.TEST ,breaks=50,xlim=c(-15,15),xlab="Average Difference number of paths",main="Common TG vs MG"); abline(v=mean(TG.NUM.PATHS[TG.ERE.INTERSECT,TG.SRE.INTERSECT],na.rm=TRUE) - mean(MG.NUM.PATHS[MG.ERE.INTERSECT,MG.SRE.INTERSECT],na.rm=TRUE),col="red")
  sum(mean(TG.NUM.PATHS[TG.ERE.INTERSECT,TG.SRE.INTERSECT],na.rm=TRUE) - mean(MG.NUM.PATHS[MG.ERE.INTERSECT,MG.SRE.INTERSECT],na.rm=TRUE) < TG.MG.INT.NUM.PATH.PERMUTATION.TEST )/length(TG.MG.INT.NUM.PATH.PERMUTATION.TEST )
TG.NG.INT.NUM.PATH.PERMUTATION.TEST  <- permutation.test.matrix(TG.NUM.PATHS[TG.ERE.INTERSECT,TG.SRE.INTERSECT],NG.NUM.PATHS[NG.ERE.INTERSECT,NG.SRE.INTERSECT],10000)
  hist(TG.NG.INT.NUM.PATH.PERMUTATION.TEST ,breaks=50,xlim=c(-2500,2500),xlab="Average Difference number of paths",main="Common TG vs NG"); abline(v=mean(TG.NUM.PATHS[TG.ERE.INTERSECT,TG.SRE.INTERSECT],na.rm=TRUE) - mean(NG.NUM.PATHS[NG.ERE.INTERSECT,NG.SRE.INTERSECT],na.rm=TRUE),col="red")
  sum(mean(TG.NUM.PATHS[TG.ERE.INTERSECT,TG.SRE.INTERSECT],na.rm=TRUE) - mean(NG.NUM.PATHS[NG.ERE.INTERSECT,NG.SRE.INTERSECT],na.rm=TRUE) < TG.NG.INT.NUM.PATH.PERMUTATION.TEST )/length(TG.NG.INT.NUM.PATH.PERMUTATION.TEST )
PG.NG.INT.NUM.PATH.PERMUTATION.TEST  <- permutation.test.matrix(PG.NUM.PATHS[PG.ERE.INTERSECT,PG.SRE.INTERSECT],NG.NUM.PATHS[NG.ERE.INTERSECT,NG.SRE.INTERSECT],10000)
  hist(PG.NG.INT.NUM.PATH.PERMUTATION.TEST ,breaks=50,xlim=c(-2500,2500),xlab="Average Difference number of paths",main="Common PG vs NG"); abline(v=mean(PG.NUM.PATHS[PG.ERE.INTERSECT,PG.SRE.INTERSECT],na.rm=TRUE) - mean(NG.NUM.PATHS[NG.ERE.INTERSECT,NG.SRE.INTERSECT],na.rm=TRUE),col="red")
  sum(mean(PG.NUM.PATHS[PG.ERE.INTERSECT,PG.SRE.INTERSECT],na.rm=TRUE) - mean(NG.NUM.PATHS[NG.ERE.INTERSECT,NG.SRE.INTERSECT],na.rm=TRUE) < PG.NG.INT.NUM.PATH.PERMUTATION.TEST )/length(PG.NG.INT.NUM.PATH.PERMUTATION.TEST )
TH.PH.INT.NUM.PATH.PERMUTATION.TEST  <- permutation.test.matrix(TH.NUM.PATHS[TH.ERE.INTERSECT,TH.SRE.INTERSECT],PH.NUM.PATHS[PH.ERE.INTERSECT,PH.SRE.INTERSECT],10000)
  hist(TH.PH.INT.NUM.PATH.PERMUTATION.TEST ,breaks=50,xlim=c(-.7,.7),xlab="Average Difference number of paths",main="Common TH vs PH"); abline(v=mean(TH.NUM.PATHS[TH.ERE.INTERSECT,TH.SRE.INTERSECT],na.rm=TRUE) - mean(PH.NUM.PATHS[PH.ERE.INTERSECT,PH.SRE.INTERSECT],na.rm=TRUE),col="red")
  sum(mean(TH.NUM.PATHS[TH.ERE.INTERSECT,TH.SRE.INTERSECT],na.rm=TRUE) - mean(PH.NUM.PATHS[PH.ERE.INTERSECT,PH.SRE.INTERSECT],na.rm=TRUE) < TH.PH.INT.NUM.PATH.PERMUTATION.TEST )/length(TH.PH.INT.NUM.PATH.PERMUTATION.TEST )
PH.MH.INT.NUM.PATH.PERMUTATION.TEST  <- permutation.test.matrix(PH.NUM.PATHS[PH.ERE.INTERSECT,PH.SRE.INTERSECT],MH.NUM.PATHS[MH.ERE.INTERSECT,MH.SRE.INTERSECT],10000)
  hist(PH.MH.INT.NUM.PATH.PERMUTATION.TEST ,breaks=50,xlim=c(-.7,.7),xlab="Average Difference number of paths",main="Common PH vs MH"); abline(v=mean(PH.NUM.PATHS[PH.ERE.INTERSECT,PH.SRE.INTERSECT],na.rm=TRUE) - mean(MH.NUM.PATHS[MH.ERE.INTERSECT,MH.SRE.INTERSECT],na.rm=TRUE),col="red")
  sum(mean(PH.NUM.PATHS[PH.ERE.INTERSECT,PH.SRE.INTERSECT],na.rm=TRUE) - mean(MH.NUM.PATHS[MH.ERE.INTERSECT,MH.SRE.INTERSECT],na.rm=TRUE) < PH.MH.INT.NUM.PATH.PERMUTATION.TEST )/length(PH.MH.INT.NUM.PATH.PERMUTATION.TEST )
MH.NH.INT.NUM.PATH.PERMUTATION.TEST  <- permutation.test.matrix(MH.NUM.PATHS[MH.ERE.INTERSECT,MH.SRE.INTERSECT],NH.NUM.PATHS[NH.ERE.INTERSECT,NH.SRE.INTERSECT],10000)
  hist(MH.NH.INT.NUM.PATH.PERMUTATION.TEST ,breaks=50,xlim=c(-10,10),xlab="Average Difference number of paths",main="Common MH vs NH"); abline(v=mean(MH.NUM.PATHS[MH.ERE.INTERSECT,MH.SRE.INTERSECT],na.rm=TRUE) - mean(NH.NUM.PATHS[NH.ERE.INTERSECT,NH.SRE.INTERSECT],na.rm=TRUE),col="red")
  sum(mean(MH.NUM.PATHS[MH.ERE.INTERSECT,MH.SRE.INTERSECT],na.rm=TRUE) - mean(NH.NUM.PATHS[NH.ERE.INTERSECT,NH.SRE.INTERSECT],na.rm=TRUE) < MH.NH.INT.NUM.PATH.PERMUTATION.TEST )/length(MH.NH.INT.NUM.PATH.PERMUTATION.TEST )
TH.MH.INT.NUM.PATH.PERMUTATION.TEST  <- permutation.test.matrix(TH.NUM.PATHS[TH.ERE.INTERSECT,TH.SRE.INTERSECT],MH.NUM.PATHS[MH.ERE.INTERSECT,MH.SRE.INTERSECT],10000)
  hist(TH.MH.INT.NUM.PATH.PERMUTATION.TEST ,breaks=50,xlim=c(-.7,.7),xlab="Average Difference number of paths",main="Common TH vs MH"); abline(v=mean(TH.NUM.PATHS[TH.ERE.INTERSECT,TH.SRE.INTERSECT],na.rm=TRUE) - mean(MH.NUM.PATHS[MH.ERE.INTERSECT,MH.SRE.INTERSECT],na.rm=TRUE),col="red")
  sum(mean(TH.NUM.PATHS[TH.ERE.INTERSECT,TH.SRE.INTERSECT],na.rm=TRUE) - mean(MH.NUM.PATHS[MH.ERE.INTERSECT,MH.SRE.INTERSECT],na.rm=TRUE) < TH.MH.INT.NUM.PATH.PERMUTATION.TEST )/length(TH.MH.INT.NUM.PATH.PERMUTATION.TEST )
TH.NH.INT.NUM.PATH.PERMUTATION.TEST  <- permutation.test.matrix(TH.NUM.PATHS[TH.ERE.INTERSECT,TH.SRE.INTERSECT],NH.NUM.PATHS[NH.ERE.INTERSECT,NH.SRE.INTERSECT],10000)
  hist(TH.NH.INT.NUM.PATH.PERMUTATION.TEST ,breaks=50,xlim=c(-10,10),xlab="Average Difference number of paths",main="Common TH vs NH"); abline(v=mean(TH.NUM.PATHS[TH.ERE.INTERSECT,TH.SRE.INTERSECT],na.rm=TRUE) - mean(NH.NUM.PATHS[NH.ERE.INTERSECT,NH.SRE.INTERSECT],na.rm=TRUE),col="red")
  sum(mean(TH.NUM.PATHS[TH.ERE.INTERSECT,TH.SRE.INTERSECT],na.rm=TRUE) - mean(NH.NUM.PATHS[NH.ERE.INTERSECT,NH.SRE.INTERSECT],na.rm=TRUE) < TH.NH.INT.NUM.PATH.PERMUTATION.TEST )/length(TH.NH.INT.NUM.PATH.PERMUTATION.TEST )
PH.NH.INT.NUM.PATH.PERMUTATION.TEST  <- permutation.test.matrix(PH.NUM.PATHS[PH.ERE.INTERSECT,PH.SRE.INTERSECT],NH.NUM.PATHS[NH.ERE.INTERSECT,NH.SRE.INTERSECT],10000)
  hist(PH.NH.INT.NUM.PATH.PERMUTATION.TEST ,breaks=50,xlim=c(-10,10),xlab="Average Difference number of paths",main="Common PH vs NH"); abline(v=mean(PH.NUM.PATHS[PH.ERE.INTERSECT,PH.SRE.INTERSECT],na.rm=TRUE) - mean(NH.NUM.PATHS[NH.ERE.INTERSECT,NH.SRE.INTERSECT],na.rm=TRUE),col="red")
  sum(mean(PH.NUM.PATHS[PH.ERE.INTERSECT,PH.SRE.INTERSECT],na.rm=TRUE) - mean(NH.NUM.PATHS[NH.ERE.INTERSECT,NH.SRE.INTERSECT],na.rm=TRUE) < PH.NH.INT.NUM.PATH.PERMUTATION.TEST )/length(PH.NH.INT.NUM.PATH.PERMUTATION.TEST )
#dev.off()

###Path distinctivness
##The more unique genotypes a set of paths goes through, the more distinct the paths are
##The more edges connect genotypes in a set of paths, the less distince the paths are 

#Read in the identity of amino acids along shortest routes between each amino acid pair
AA.INT.Z <- read.table("AA.INT.Z.txt",header=TRUE,colClasses = "character")

##Calculate effective number of paths between ERE and SRE
#TG.ERE.DISTINCT <- distinct.paths(TG.ACT.NET,TG.ERE.LIST,TG.SRE.LIST)
#save(TG.ERE.DISTINCT,file="TG.ERE.DISTINCT.rda")
#PG.ERE.DISTINCT <- distinct.paths(PG.ACT.NET,PG.ERE.LIST,PG.SRE.LIST)
#save(PG.ERE.DISTINCT,file="PG.ERE.DISTINCT.rda")
#MG.ERE.DISTINCT <- distinct.paths(MG.ACT.NET,MG.ERE.LIST,MG.SRE.LIST)
#save(MG.ERE.DISTINCT,file="MG.ERE.DISTINCT.rda")
#NG.ERE.DISTINCT <- distinct.paths.NG(NG.ERE.LIST,NG.SRE.LIST)
#save(NG.ERE.DISTINCT,file="NG.ERE.DISTINCT.rda")
#TH.ERE.DISTINCT <- distinct.paths(TH.ACT.NET,TH.ERE.LIST,TH.SRE.LIST)
#save(TH.ERE.DISTINCT,file="TH.ERE.DISTINCT.rda")
#PH.ERE.DISTINCT <- distinct.paths(PH.ACT.NET,PH.ERE.LIST,PH.SRE.LIST)
#save(PH.ERE.DISTINCT,file="PH.ERE.DISTINCT.rda")
#MH.ERE.DISTINCT <- distinct.paths(MH.ACT.NET,MH.ERE.LIST,MH.SRE.LIST)
#save(MH.ERE.DISTINCT,file="MH.ERE.DISTINCT.rda")
#NH.ERE.DISTINCT <- matrix(0,nrow=nrow(NH.NUM.PATHS),ncol=ncol(NH.NUM.PATHS))
#NH.ERE.DISTINCT[which(NH.NUM.PATHS == 1, arr.ind = TRUE)] <- distinct.paths(NH.ACT.NET,"AAAA","AAAC")
#NH.ERE.DISTINCT[which(NH.NUM.PATHS == 2, arr.ind = TRUE)] <- distinct.paths(NH.ACT.NET,"AAAA","AACC")
#NH.ERE.DISTINCT[which(NH.NUM.PATHS == 6, arr.ind = TRUE)] <- distinct.paths(NH.ACT.NET,"AAAA","ACCC")
#NH.ERE.DISTINCT[which(NH.NUM.PATHS == 24,arr.ind = TRUE)] <- distinct.paths(NH.ACT.NET,"AAAA","CCCC")
#save(NH.ERE.DISTINCT,file="NH.ERE.DISTINCT.rda")

load("TG.ERE.DISTINCT.rda")
load("PG.ERE.DISTINCT.rda")
load("MG.ERE.DISTINCT.rda")
load("NG.ERE.DISTINCT.rda")
load("TH.ERE.DISTINCT.rda")
load("PH.ERE.DISTINCT.rda")
load("MH.ERE.DISTINCT.rda")
load("NH.ERE.DISTINCT.rda")

##Collect mean values
ALL.ERE.DISTINCT.AVG <- c(mean(TG.ERE.DISTINCT,na.rm=TRUE),mean(PG.ERE.DISTINCT,na.rm=TRUE),mean(MG.ERE.DISTINCT,na.rm=TRUE),mean(NG.ERE.DISTINCT,na.rm=TRUE),mean(TH.ERE.DISTINCT,na.rm=TRUE),mean(PH.ERE.DISTINCT,na.rm=TRUE),mean(MH.ERE.DISTINCT,na.rm=TRUE),mean(NH.ERE.DISTINCT,na.rm=TRUE))
INT.ERE.DISTINCT.AVG <- c(mean(TG.ERE.DISTINCT[TG.ERE.INTERSECT,TG.SRE.INTERSECT]),mean(PG.ERE.DISTINCT[PG.ERE.INTERSECT,PG.SRE.INTERSECT]),mean(MG.ERE.DISTINCT[MG.ERE.INTERSECT,MG.SRE.INTERSECT]),mean(NG.ERE.DISTINCT[NG.ERE.INTERSECT,NG.SRE.INTERSECT]),mean(TH.ERE.DISTINCT[TH.ERE.INTERSECT,TH.SRE.INTERSECT]),mean(PH.ERE.DISTINCT[PH.ERE.INTERSECT,PH.SRE.INTERSECT]),mean(MH.ERE.DISTINCT[MH.ERE.INTERSECT,MH.SRE.INTERSECT]),mean(NH.ERE.DISTINCT[NH.ERE.INTERSECT,NH.SRE.INTERSECT]))

#pdf("DISTINCT.PATHS.pdf")
par(mfrow=c(1,1))

vioplot(c(TG.ERE.DISTINCT),c(PG.ERE.DISTINCT),c(MG.ERE.DISTINCT),c(NG.ERE.DISTINCT),names=c("3G","2G","1G","NG"),cex.axis=0.75,cex.main=0.5,col="darkred", rectCol="white",lineCol="white",pchMed=19,colMed="black",main="All Distinct Paths")
vioplot(c(TH.ERE.DISTINCT),c(PH.ERE.DISTINCT),c(MH.ERE.DISTINCT),c(NH.ERE.DISTINCT),names=c("3H","2H","1H","NH"),cex.axis=0.75,cex.main=0.5,col="skyblue2",rectCol="white",lineCol="white",pchMed=19,colMed="black",main="All Distinct Paths")

barplot(ALL.ERE.DISTINCT.AVG[1:4],names=c("3G","2G","1G","NG"),cex.axis=0.75,cex.main=0.5,col="darkred" ,main="All Distinct Paths")
barplot(ALL.ERE.DISTINCT.AVG[5:8],names=c("3H","2H","1H","NH"),cex.axis=0.75,cex.main=0.5,col="skyblue2",main="All Distinct Paths")

vioplot(c(TG.ERE.DISTINCT[TG.ERE.INTERSECT,TG.SRE.INTERSECT]),c(PG.ERE.DISTINCT[PG.ERE.INTERSECT,PG.SRE.INTERSECT]),c(MG.ERE.DISTINCT[MG.ERE.INTERSECT,MG.SRE.INTERSECT]),c(NG.ERE.DISTINCT[NG.ERE.INTERSECT,NG.SRE.INTERSECT]),names=c("3G","2G","1G","NG"),cex.axis=0.75,cex.main=0.5,col="darkred", rectCol="white",lineCol="white",pchMed=19,colMed="black",main="Common Distinct Paths")
vioplot(c(TH.ERE.DISTINCT[TH.ERE.INTERSECT,TH.SRE.INTERSECT]),c(PH.ERE.DISTINCT[PH.ERE.INTERSECT,PH.SRE.INTERSECT]),c(MH.ERE.DISTINCT[MH.ERE.INTERSECT,MH.SRE.INTERSECT]),c(NH.ERE.DISTINCT[NH.ERE.INTERSECT,NH.SRE.INTERSECT]),names=c("3H","2H","1H","NH"),cex.axis=0.75,cex.main=0.5,col="skyblue2",rectCol="white",lineCol="white",pchMed=19,colMed="black",main="Common Distinct Paths")

barplot(INT.ERE.DISTINCT.AVG[1:4],names=c("3G","2G","1G","NG"),cex.axis=0.75,cex.main=0.5,col="darkred" ,main="Common Distinct Paths")
barplot(INT.ERE.DISTINCT.AVG[5:8],names=c("3H","2H","1H","NH"),cex.axis=0.75,cex.main=0.5,col="skyblue2",main="Common Distinct Paths")

###Permutation tests
par(mfrow=c(3,2))
#All genotypes
TG.PG.ALL.DISTINCT.PATH.PERMUTATION.TEST <- permutation.test.matrix(TG.ERE.DISTINCT,PG.ERE.DISTINCT,10000)
  hist(TG.PG.ALL.DISTINCT.PATH.PERMUTATION.TEST,breaks=50,xlim=c(-.1,.1),xlab="Average Difference distinct paths",main="All TG vs PG"); abline(v=mean(TG.ERE.DISTINCT,na.rm=TRUE) - mean(PG.ERE.DISTINCT,na.rm=TRUE),col="red")
  sum(mean(TG.ERE.DISTINCT,na.rm=TRUE) - mean(PG.ERE.DISTINCT,na.rm=TRUE) < TG.PG.ALL.DISTINCT.PATH.PERMUTATION.TEST)/length(TG.PG.ALL.DISTINCT.PATH.PERMUTATION.TEST)
PG.MG.ALL.DISTINCT.PATH.PERMUTATION.TEST <- permutation.test.matrix(PG.ERE.DISTINCT,MG.ERE.DISTINCT,10000)
  hist(PG.MG.ALL.DISTINCT.PATH.PERMUTATION.TEST,breaks=50,xlim=c(-.1,.1),xlab="Average Difference distinct paths",main="All PG vs MG"); abline(v=mean(PG.ERE.DISTINCT,na.rm=TRUE) - mean(MG.ERE.DISTINCT,na.rm=TRUE),col="red")
  sum(mean(PG.ERE.DISTINCT,na.rm=TRUE) - mean(MG.ERE.DISTINCT,na.rm=TRUE) < PG.MG.ALL.DISTINCT.PATH.PERMUTATION.TEST)/length(PG.MG.ALL.DISTINCT.PATH.PERMUTATION.TEST)
MG.NG.ALL.DISTINCT.PATH.PERMUTATION.TEST <- permutation.test.matrix(MG.ERE.DISTINCT,NG.ERE.DISTINCT,10000)
  hist(MG.NG.ALL.DISTINCT.PATH.PERMUTATION.TEST,breaks=50,xlim=c(-3,3),xlab="Average Difference distinct paths",main="All MG vs NG"); abline(v=mean(MG.ERE.DISTINCT,na.rm=TRUE) - mean(NG.ERE.DISTINCT,na.rm=TRUE),col="red")
  sum(mean(MG.ERE.DISTINCT,na.rm=TRUE) - mean(NG.ERE.DISTINCT,na.rm=TRUE) < MG.NG.ALL.DISTINCT.PATH.PERMUTATION.TEST)/length(MG.NG.ALL.DISTINCT.PATH.PERMUTATION.TEST)
TG.MG.ALL.DISTINCT.PATH.PERMUTATION.TEST <- permutation.test.matrix(TG.ERE.DISTINCT,MG.ERE.DISTINCT,10000)
  hist(TG.MG.ALL.DISTINCT.PATH.PERMUTATION.TEST,breaks=50,xlim=c(-.1,.1),xlab="Average Difference distinct paths",main="All TG vs MG"); abline(v=mean(TG.ERE.DISTINCT,na.rm=TRUE) - mean(MG.ERE.DISTINCT,na.rm=TRUE),col="red")
  sum(mean(TG.ERE.DISTINCT,na.rm=TRUE) - mean(MG.ERE.DISTINCT,na.rm=TRUE) < TG.MG.ALL.DISTINCT.PATH.PERMUTATION.TEST)/length(TG.MG.ALL.DISTINCT.PATH.PERMUTATION.TEST)
TG.NG.ALL.DISTINCT.PATH.PERMUTATION.TEST <- permutation.test.matrix(TG.ERE.DISTINCT,NG.ERE.DISTINCT,10000)
  hist(TG.NG.ALL.DISTINCT.PATH.PERMUTATION.TEST,breaks=50,xlim=c(-3,3),xlab="Average Difference distinct paths",main="All TG vs NG"); abline(v=mean(TG.ERE.DISTINCT,na.rm=TRUE) - mean(NG.ERE.DISTINCT,na.rm=TRUE),col="red")
  sum(mean(TG.ERE.DISTINCT,na.rm=TRUE) - mean(NG.ERE.DISTINCT,na.rm=TRUE) < TG.NG.ALL.DISTINCT.PATH.PERMUTATION.TEST)/length(TG.NG.ALL.DISTINCT.PATH.PERMUTATION.TEST)
PG.NG.ALL.DISTINCT.PATH.PERMUTATION.TEST <- permutation.test.matrix(PG.ERE.DISTINCT,NG.ERE.DISTINCT,10000)
  hist(PG.NG.ALL.DISTINCT.PATH.PERMUTATION.TEST,breaks=50,xlim=c(-3,3),xlab="Average Difference distinct paths",main="All PG vs NG"); abline(v=mean(PG.ERE.DISTINCT,na.rm=TRUE) - mean(NG.ERE.DISTINCT,na.rm=TRUE),col="red")
  sum(mean(PG.ERE.DISTINCT,na.rm=TRUE) - mean(NG.ERE.DISTINCT,na.rm=TRUE) < PG.NG.ALL.DISTINCT.PATH.PERMUTATION.TEST)/length(PG.NG.ALL.DISTINCT.PATH.PERMUTATION.TEST)
TH.PH.ALL.DISTINCT.PATH.PERMUTATION.TEST <- permutation.test.matrix(TH.ERE.DISTINCT,PH.ERE.DISTINCT,10000)
  hist(TH.PH.ALL.DISTINCT.PATH.PERMUTATION.TEST,breaks=50,xlim=c(-.7,.7),xlab="Average Difference distinct paths",main="All TH vs PH"); abline(v=mean(TH.ERE.DISTINCT,na.rm=TRUE) - mean(PH.ERE.DISTINCT,na.rm=TRUE),col="red")
  sum(mean(TH.ERE.DISTINCT,na.rm=TRUE) - mean(PH.ERE.DISTINCT,na.rm=TRUE) < TH.PH.ALL.DISTINCT.PATH.PERMUTATION.TEST)/length(TH.PH.ALL.DISTINCT.PATH.PERMUTATION.TEST)
PH.MH.ALL.DISTINCT.PATH.PERMUTATION.TEST <- permutation.test.matrix(PH.ERE.DISTINCT,MH.ERE.DISTINCT,10000)
  hist(PH.MH.ALL.DISTINCT.PATH.PERMUTATION.TEST,breaks=50,xlim=c(-.7,.7),xlab="Average Difference distinct paths",main="All PH vs MH"); abline(v=mean(PH.ERE.DISTINCT,na.rm=TRUE) - mean(MH.ERE.DISTINCT,na.rm=TRUE),col="red")
  sum(mean(PH.ERE.DISTINCT,na.rm=TRUE) - mean(MH.ERE.DISTINCT,na.rm=TRUE) < PH.MH.ALL.DISTINCT.PATH.PERMUTATION.TEST)/length(PH.MH.ALL.DISTINCT.PATH.PERMUTATION.TEST)
MH.NH.ALL.DISTINCT.PATH.PERMUTATION.TEST <- permutation.test.matrix(MH.ERE.DISTINCT,NH.ERE.DISTINCT,10000)
  hist(MH.NH.ALL.DISTINCT.PATH.PERMUTATION.TEST,breaks=50,xlim=c(-.7,.7),xlab="Average Difference distinct paths",main="All MH vs NH"); abline(v=mean(MH.ERE.DISTINCT,na.rm=TRUE) - mean(NH.ERE.DISTINCT,na.rm=TRUE),col="red")
  sum(mean(MH.ERE.DISTINCT,na.rm=TRUE) - mean(NH.ERE.DISTINCT,na.rm=TRUE) < MH.NH.ALL.DISTINCT.PATH.PERMUTATION.TEST)/length(MH.NH.ALL.DISTINCT.PATH.PERMUTATION.TEST)
TH.MH.ALL.DISTINCT.PATH.PERMUTATION.TEST <- permutation.test.matrix(TH.ERE.DISTINCT,MH.ERE.DISTINCT,10000)
  hist(TH.MH.ALL.DISTINCT.PATH.PERMUTATION.TEST,breaks=50,xlim=c(-.7,.7),xlab="Average Difference distinct paths",main="All TH vs MH"); abline(v=mean(TH.ERE.DISTINCT,na.rm=TRUE) - mean(MH.ERE.DISTINCT,na.rm=TRUE),col="red")
  sum(mean(TH.ERE.DISTINCT,na.rm=TRUE) - mean(MH.ERE.DISTINCT,na.rm=TRUE) < TH.MH.ALL.DISTINCT.PATH.PERMUTATION.TEST)/length(TH.MH.ALL.DISTINCT.PATH.PERMUTATION.TEST)
TH.NH.ALL.DISTINCT.PATH.PERMUTATION.TEST <- permutation.test.matrix(TH.ERE.DISTINCT,NH.ERE.DISTINCT,10000)
  hist(TH.NH.ALL.DISTINCT.PATH.PERMUTATION.TEST,breaks=50,xlim=c(-.7,.7),xlab="Average Difference distinct paths",main="All TH vs NH"); abline(v=mean(TH.ERE.DISTINCT,na.rm=TRUE) - mean(NH.ERE.DISTINCT,na.rm=TRUE),col="red")
  sum(mean(TH.ERE.DISTINCT,na.rm=TRUE) - mean(NH.ERE.DISTINCT,na.rm=TRUE) < TH.NH.ALL.DISTINCT.PATH.PERMUTATION.TEST)/length(TH.NH.ALL.DISTINCT.PATH.PERMUTATION.TEST)
PH.NH.ALL.DISTINCT.PATH.PERMUTATION.TEST <- permutation.test.matrix(PH.ERE.DISTINCT,NH.ERE.DISTINCT,10000)
  hist(PH.NH.ALL.DISTINCT.PATH.PERMUTATION.TEST,breaks=50,xlim=c(-.7,.7),xlab="Average Difference distinct paths",main="All PH vs NH"); abline(v=mean(PH.ERE.DISTINCT,na.rm=TRUE) - mean(NH.ERE.DISTINCT,na.rm=TRUE),col="red")
  sum(mean(PH.ERE.DISTINCT,na.rm=TRUE) - mean(NH.ERE.DISTINCT,na.rm=TRUE) < PH.NH.ALL.DISTINCT.PATH.PERMUTATION.TEST)/length(PH.NH.ALL.DISTINCT.PATH.PERMUTATION.TEST)

#Common genotypes
TG.PG.INT.DISTINCT.PATH.PERMUTATION.TEST <- permutation.test.matrix(TG.ERE.DISTINCT[TG.ERE.INTERSECT,TG.SRE.INTERSECT],PG.ERE.DISTINCT[PG.ERE.INTERSECT,PG.SRE.INTERSECT],10000)
  hist(TG.PG.INT.DISTINCT.PATH.PERMUTATION.TEST,breaks=50,xlim=c(-.3,.3),xlab="Average Difference distinct paths",main="Common TG vs PG"); abline(v=mean(TG.ERE.DISTINCT[TG.ERE.INTERSECT,TG.SRE.INTERSECT],na.rm=TRUE) - mean(PG.ERE.DISTINCT[PG.ERE.INTERSECT,PG.SRE.INTERSECT],na.rm=TRUE),col="red")
  sum(mean(TG.ERE.DISTINCT[TG.ERE.INTERSECT,TG.SRE.INTERSECT],na.rm=TRUE) - mean(PG.ERE.DISTINCT[PG.ERE.INTERSECT,PG.SRE.INTERSECT],na.rm=TRUE) < TG.PG.INT.DISTINCT.PATH.PERMUTATION.TEST)/length(TG.PG.INT.DISTINCT.PATH.PERMUTATION.TEST)
PG.MG.INT.DISTINCT.PATH.PERMUTATION.TEST <- permutation.test.matrix(PG.ERE.DISTINCT[PG.ERE.INTERSECT,PG.SRE.INTERSECT],MG.ERE.DISTINCT[MG.ERE.INTERSECT,MG.SRE.INTERSECT],10000)
  hist(PG.MG.INT.DISTINCT.PATH.PERMUTATION.TEST,breaks=50,xlim=c(-.3,.3),xlab="Average Difference distinct paths",main="Common PG vs MG"); abline(v=mean(PG.ERE.DISTINCT[PG.ERE.INTERSECT,PG.SRE.INTERSECT],na.rm=TRUE) - mean(MG.ERE.DISTINCT[MG.ERE.INTERSECT,MG.SRE.INTERSECT],na.rm=TRUE),col="red")
  sum(mean(PG.ERE.DISTINCT[PG.ERE.INTERSECT,PG.SRE.INTERSECT],na.rm=TRUE) - mean(MG.ERE.DISTINCT[MG.ERE.INTERSECT,MG.SRE.INTERSECT],na.rm=TRUE) < PG.MG.INT.DISTINCT.PATH.PERMUTATION.TEST)/length(PG.MG.INT.DISTINCT.PATH.PERMUTATION.TEST)
MG.NG.INT.DISTINCT.PATH.PERMUTATION.TEST <- permutation.test.matrix(MG.ERE.DISTINCT[MG.ERE.INTERSECT,MG.SRE.INTERSECT],NG.ERE.DISTINCT[NG.ERE.INTERSECT,NG.SRE.INTERSECT],10000)
  hist(MG.NG.INT.DISTINCT.PATH.PERMUTATION.TEST,breaks=50,xlim=c(-2,2),xlab="Average Difference distinct paths",main="Common MG vs NG"); abline(v=mean(MG.ERE.DISTINCT[MG.ERE.INTERSECT,MG.SRE.INTERSECT],na.rm=TRUE) - mean(NG.ERE.DISTINCT[NG.ERE.INTERSECT,NG.SRE.INTERSECT],na.rm=TRUE),col="red")
  sum(mean(MG.ERE.DISTINCT[MG.ERE.INTERSECT,MG.SRE.INTERSECT],na.rm=TRUE) - mean(NG.ERE.DISTINCT[NG.ERE.INTERSECT,NG.SRE.INTERSECT],na.rm=TRUE) < MG.NG.INT.DISTINCT.PATH.PERMUTATION.TEST)/length(MG.NG.INT.DISTINCT.PATH.PERMUTATION.TEST)
TG.MG.INT.DISTINCT.PATH.PERMUTATION.TEST <- permutation.test.matrix(TG.ERE.DISTINCT[TG.ERE.INTERSECT,TG.SRE.INTERSECT],MG.ERE.DISTINCT[MG.ERE.INTERSECT,MG.SRE.INTERSECT],10000)
  hist(TG.MG.INT.DISTINCT.PATH.PERMUTATION.TEST,breaks=50,xlim=c(-.3,.3),xlab="Average Difference distinct paths",main="Common TG vs MG"); abline(v=mean(TG.ERE.DISTINCT[TG.ERE.INTERSECT,TG.SRE.INTERSECT],na.rm=TRUE) - mean(MG.ERE.DISTINCT[MG.ERE.INTERSECT,MG.SRE.INTERSECT],na.rm=TRUE),col="red")
  sum(mean(TG.ERE.DISTINCT[TG.ERE.INTERSECT,TG.SRE.INTERSECT],na.rm=TRUE) - mean(MG.ERE.DISTINCT[MG.ERE.INTERSECT,MG.SRE.INTERSECT],na.rm=TRUE) < TG.MG.INT.DISTINCT.PATH.PERMUTATION.TEST)/length(TG.MG.INT.DISTINCT.PATH.PERMUTATION.TEST)
TG.NG.INT.DISTINCT.PATH.PERMUTATION.TEST <- permutation.test.matrix(TG.ERE.DISTINCT[TG.ERE.INTERSECT,TG.SRE.INTERSECT],NG.ERE.DISTINCT[NG.ERE.INTERSECT,NG.SRE.INTERSECT],10000)
  hist(TG.NG.INT.DISTINCT.PATH.PERMUTATION.TEST,breaks=50,xlim=c(-2,2),xlab="Average Difference distinct paths",main="Common TG vs NG"); abline(v=mean(TG.ERE.DISTINCT[TG.ERE.INTERSECT,TG.SRE.INTERSECT],na.rm=TRUE) - mean(NG.ERE.DISTINCT[NG.ERE.INTERSECT,NG.SRE.INTERSECT],na.rm=TRUE),col="red")
  sum(mean(TG.ERE.DISTINCT[TG.ERE.INTERSECT,TG.SRE.INTERSECT],na.rm=TRUE) - mean(NG.ERE.DISTINCT[NG.ERE.INTERSECT,NG.SRE.INTERSECT],na.rm=TRUE) < TG.NG.INT.DISTINCT.PATH.PERMUTATION.TEST)/length(TG.NG.INT.DISTINCT.PATH.PERMUTATION.TEST)
PG.NG.INT.DISTINCT.PATH.PERMUTATION.TEST <- permutation.test.matrix(PG.ERE.DISTINCT[PG.ERE.INTERSECT,PG.SRE.INTERSECT],NG.ERE.DISTINCT[NG.ERE.INTERSECT,NG.SRE.INTERSECT],10000)
  hist(PG.NG.INT.DISTINCT.PATH.PERMUTATION.TEST,breaks=50,xlim=c(-2,2),xlab="Average Difference distinct paths",main="Common PG vs NG"); abline(v=mean(PG.ERE.DISTINCT[PG.ERE.INTERSECT,PG.SRE.INTERSECT],na.rm=TRUE) - mean(NG.ERE.DISTINCT[NG.ERE.INTERSECT,NG.SRE.INTERSECT],na.rm=TRUE),col="red")
  sum(mean(PG.ERE.DISTINCT[PG.ERE.INTERSECT,PG.SRE.INTERSECT],na.rm=TRUE) - mean(NG.ERE.DISTINCT[NG.ERE.INTERSECT,NG.SRE.INTERSECT],na.rm=TRUE) < PG.NG.INT.DISTINCT.PATH.PERMUTATION.TEST)/length(PG.NG.INT.DISTINCT.PATH.PERMUTATION.TEST)
TH.PH.INT.DISTINCT.PATH.PERMUTATION.TEST <- permutation.test.matrix(TH.ERE.DISTINCT[TH.ERE.INTERSECT,TH.SRE.INTERSECT],PH.ERE.DISTINCT[PH.ERE.INTERSECT,PH.SRE.INTERSECT],10000)
  hist(TH.PH.INT.DISTINCT.PATH.PERMUTATION.TEST,breaks=50,xlim=c(-.2,.2),xlab="Average Difference distinct paths",main="Common TH vs PH"); abline(v=mean(TH.ERE.DISTINCT[TH.ERE.INTERSECT,TH.SRE.INTERSECT],na.rm=TRUE) - mean(PH.ERE.DISTINCT[PH.ERE.INTERSECT,PH.SRE.INTERSECT],na.rm=TRUE),col="red")
  sum(mean(TH.ERE.DISTINCT[TH.ERE.INTERSECT,TH.SRE.INTERSECT],na.rm=TRUE) - mean(PH.ERE.DISTINCT[PH.ERE.INTERSECT,PH.SRE.INTERSECT],na.rm=TRUE) < TH.PH.INT.DISTINCT.PATH.PERMUTATION.TEST)/length(TH.PH.INT.DISTINCT.PATH.PERMUTATION.TEST)
PH.MH.INT.DISTINCT.PATH.PERMUTATION.TEST <- permutation.test.matrix(PH.ERE.DISTINCT[PH.ERE.INTERSECT,PH.SRE.INTERSECT],MH.ERE.DISTINCT[MH.ERE.INTERSECT,MH.SRE.INTERSECT],10000)
  hist(PH.MH.INT.DISTINCT.PATH.PERMUTATION.TEST,breaks=50,xlim=c(-.2,.2),xlab="Average Difference distinct paths",main="Common PH vs MH"); abline(v=mean(PH.ERE.DISTINCT[PH.ERE.INTERSECT,PH.SRE.INTERSECT],na.rm=TRUE) - mean(MH.ERE.DISTINCT[MH.ERE.INTERSECT,MH.SRE.INTERSECT],na.rm=TRUE),col="red")
  sum(mean(PH.ERE.DISTINCT[PH.ERE.INTERSECT,PH.SRE.INTERSECT],na.rm=TRUE) - mean(MH.ERE.DISTINCT[MH.ERE.INTERSECT,MH.SRE.INTERSECT],na.rm=TRUE) < PH.MH.INT.DISTINCT.PATH.PERMUTATION.TEST)/length(PH.MH.INT.DISTINCT.PATH.PERMUTATION.TEST)
MH.NH.INT.DISTINCT.PATH.PERMUTATION.TEST <- permutation.test.matrix(MH.ERE.DISTINCT[MH.ERE.INTERSECT,MH.SRE.INTERSECT],NH.ERE.DISTINCT[NH.ERE.INTERSECT,NH.SRE.INTERSECT],10000)
  hist(MH.NH.INT.DISTINCT.PATH.PERMUTATION.TEST,breaks=50,xlim=c(-.7,.7),xlab="Average Difference distinct paths",main="Common MH vs NH"); abline(v=mean(MH.ERE.DISTINCT[MH.ERE.INTERSECT,MH.SRE.INTERSECT],na.rm=TRUE) - mean(NH.ERE.DISTINCT[NH.ERE.INTERSECT,NH.SRE.INTERSECT],na.rm=TRUE),col="red")
  sum(mean(MH.ERE.DISTINCT[MH.ERE.INTERSECT,MH.SRE.INTERSECT],na.rm=TRUE) - mean(NH.ERE.DISTINCT[NH.ERE.INTERSECT,NH.SRE.INTERSECT],na.rm=TRUE) < MH.NH.INT.DISTINCT.PATH.PERMUTATION.TEST)/length(MH.NH.INT.DISTINCT.PATH.PERMUTATION.TEST)
TH.MH.INT.DISTINCT.PATH.PERMUTATION.TEST <- permutation.test.matrix(TH.ERE.DISTINCT[TH.ERE.INTERSECT,TH.SRE.INTERSECT],MH.ERE.DISTINCT[MH.ERE.INTERSECT,MH.SRE.INTERSECT],10000)
  hist(TH.MH.INT.DISTINCT.PATH.PERMUTATION.TEST,breaks=50,xlim=c(-.2,.2),xlab="Average Difference distinct paths",main="Common TH vs MH"); abline(v=mean(TH.ERE.DISTINCT[TH.ERE.INTERSECT,TH.SRE.INTERSECT],na.rm=TRUE) - mean(MH.ERE.DISTINCT[MH.ERE.INTERSECT,MH.SRE.INTERSECT],na.rm=TRUE),col="red")
  sum(mean(TH.ERE.DISTINCT[TH.ERE.INTERSECT,TH.SRE.INTERSECT],na.rm=TRUE) - mean(MH.ERE.DISTINCT[MH.ERE.INTERSECT,MH.SRE.INTERSECT],na.rm=TRUE) < TH.MH.INT.DISTINCT.PATH.PERMUTATION.TEST)/length(TH.MH.INT.DISTINCT.PATH.PERMUTATION.TEST)
TH.NH.INT.DISTINCT.PATH.PERMUTATION.TEST <- permutation.test.matrix(TH.ERE.DISTINCT[TH.ERE.INTERSECT,TH.SRE.INTERSECT],NH.ERE.DISTINCT[NH.ERE.INTERSECT,NH.SRE.INTERSECT],10000)
  hist(TH.NH.INT.DISTINCT.PATH.PERMUTATION.TEST,breaks=50,xlim=c(-.7,.7),xlab="Average Difference distinct paths",main="Common TH vs NH"); abline(v=mean(TH.ERE.DISTINCT[TH.ERE.INTERSECT,TH.SRE.INTERSECT],na.rm=TRUE) - mean(NH.ERE.DISTINCT[NH.ERE.INTERSECT,NH.SRE.INTERSECT],na.rm=TRUE),col="red")
  sum(mean(TH.ERE.DISTINCT[TH.ERE.INTERSECT,TH.SRE.INTERSECT],na.rm=TRUE) - mean(NH.ERE.DISTINCT[NH.ERE.INTERSECT,NH.SRE.INTERSECT],na.rm=TRUE) < TH.NH.INT.DISTINCT.PATH.PERMUTATION.TEST)/length(TH.NH.INT.DISTINCT.PATH.PERMUTATION.TEST)
PH.NH.INT.DISTINCT.PATH.PERMUTATION.TEST <- permutation.test.matrix(PH.ERE.DISTINCT[PH.ERE.INTERSECT,PH.SRE.INTERSECT],NH.ERE.DISTINCT[NH.ERE.INTERSECT,NH.SRE.INTERSECT],10000)
  hist(PH.NH.INT.DISTINCT.PATH.PERMUTATION.TEST,breaks=50,xlim=c(-.7,.7),xlab="Average Difference distinct paths",main="Common PH vs NH"); abline(v=mean(PH.ERE.DISTINCT[PH.ERE.INTERSECT,PH.SRE.INTERSECT],na.rm=TRUE) - mean(NH.ERE.DISTINCT[NH.ERE.INTERSECT,NH.SRE.INTERSECT],na.rm=TRUE),col="red")
  sum(mean(PH.ERE.DISTINCT[PH.ERE.INTERSECT,PH.SRE.INTERSECT],na.rm=TRUE) - mean(NH.ERE.DISTINCT[NH.ERE.INTERSECT,NH.SRE.INTERSECT],na.rm=TRUE) < PH.NH.INT.DISTINCT.PATH.PERMUTATION.TEST)/length(PH.NH.INT.DISTINCT.PATH.PERMUTATION.TEST)
#dev.off()
  
######################
##10. Neighbor Class###
######################

##Calculate class of ERE neighbors
MG.ALL_ERE.ALL_GENO.NEIGHBOR.CLASS <- lapply(1:length(MG.ERE.LIST),FUN=function(x) {
  NODE <- MG.ERE.LIST[x]
  neighbors(MG.ACT.NET,v=NODE)
})
MG.ALL_ERE.ALL_GENO.NEIGHBOR.CLASS <- names(unlist(MG.ALL_ERE.ALL_GENO.NEIGHBOR.CLASS))
MG.ALL_ERE.ALL_GENO.NEIGHBOR.CLASS <- sapply(MG.ALL_ERE.ALL_GENO.NEIGHBOR.CLASS,FUN=get.genotype.class,NETWORKTYPE="M")

MG.INT_ERE.ALL_GENO.NEIGHBOR.CLASS <- lapply(1:length(which(MG.ERE.LIST %in% c(PG.ERE.LIST,TG.ERE.LIST))),FUN=function(x) {
  NODE <- MG.ERE.LIST[which(MG.ERE.LIST %in% c(PG.ERE.LIST,TG.ERE.LIST))][x]
  neighbors(MG.ACT.NET,v=NODE)
})
MG.INT_ERE.ALL_GENO.NEIGHBOR.CLASS <- names(unlist(MG.INT_ERE.ALL_GENO.NEIGHBOR.CLASS))
MG.INT_ERE.ALL_GENO.NEIGHBOR.CLASS <- sapply(MG.INT_ERE.ALL_GENO.NEIGHBOR.CLASS,FUN=get.genotype.class,NETWORKTYPE="M")


PG.ALL_ERE.ALL_GENO.NEIGHBOR.CLASS <- lapply(1:length(PG.ERE.LIST),FUN=function(x) {
  NODE <- PG.ERE.LIST[x]
  neighbors(PG.ACT.NET,v=NODE)
})
PG.ALL_ERE.ALL_GENO.NEIGHBOR.CLASS <- names(unlist(PG.ALL_ERE.ALL_GENO.NEIGHBOR.CLASS))
PG.ALL_ERE.ALL_GENO.NEIGHBOR.CLASS <- sapply(PG.ALL_ERE.ALL_GENO.NEIGHBOR.CLASS,FUN=get.genotype.class,NETWORKTYPE="P")

PG.INT_ERE.ALL_GENO.NEIGHBOR.CLASS <- lapply(1:length(which(PG.ERE.LIST %in% c(MG.ERE.LIST,TG.ERE.LIST))),FUN=function(x) {
  NODE <- PG.ERE.LIST[which(PG.ERE.LIST %in% c(MG.ERE.LIST,TG.ERE.LIST))][x]
  neighbors(PG.ACT.NET,v=NODE)
})
PG.INT_ERE.ALL_GENO.NEIGHBOR.CLASS <- names(unlist(PG.INT_ERE.ALL_GENO.NEIGHBOR.CLASS))
PG.INT_ERE.ALL_GENO.NEIGHBOR.CLASS <- sapply(PG.INT_ERE.ALL_GENO.NEIGHBOR.CLASS,FUN=get.genotype.class,NETWORKTYPE="P")


TG.ALL_ERE.ALL_GENO.NEIGHBOR.CLASS <- lapply(1:length(TG.ERE.LIST),FUN=function(x) {
  NODE <- TG.ERE.LIST[x]
  neighbors(TG.ACT.NET,v=NODE)
})
TG.ALL_ERE.ALL_GENO.NEIGHBOR.CLASS <- names(unlist(TG.ALL_ERE.ALL_GENO.NEIGHBOR.CLASS))
TG.ALL_ERE.ALL_GENO.NEIGHBOR.CLASS <- sapply(TG.ALL_ERE.ALL_GENO.NEIGHBOR.CLASS,FUN=get.genotype.class,NETWORKTYPE="T")

TG.INT_ERE.ALL_GENO.NEIGHBOR.CLASS <- lapply(1:length(which(TG.ERE.LIST %in% c(MG.ERE.LIST,PG.ERE.LIST))),FUN=function(x) {
  NODE <- TG.ERE.LIST[which(TG.ERE.LIST %in% c(MG.ERE.LIST,PG.ERE.LIST))][x]
  neighbors(TG.ACT.NET,v=NODE)
})
TG.INT_ERE.ALL_GENO.NEIGHBOR.CLASS <- names(unlist(TG.INT_ERE.ALL_GENO.NEIGHBOR.CLASS))
TG.INT_ERE.ALL_GENO.NEIGHBOR.CLASS <- sapply(TG.INT_ERE.ALL_GENO.NEIGHBOR.CLASS,FUN=get.genotype.class,NETWORKTYPE="T")

##Calculate class of promiscuous neighbors
MG.ALL_PRO.ALL_GENO.NEIGHBOR.CLASS <- lapply(1:length(MG.PRO.LIST),FUN=function(x) {
  NODE <- MG.PRO.LIST[x]
  neighbors(MG.ACT.NET,v=NODE)
})
MG.ALL_PRO.ALL_GENO.NEIGHBOR.CLASS <- names(unlist(MG.ALL_PRO.ALL_GENO.NEIGHBOR.CLASS))
MG.ALL_PRO.ALL_GENO.NEIGHBOR.CLASS <- sapply(MG.ALL_PRO.ALL_GENO.NEIGHBOR.CLASS,FUN=get.genotype.class,NETWORKTYPE="M")

MG.INT_PRO.ALL_GENO.NEIGHBOR.CLASS <- lapply(1:length(which(MG.PRO.LIST %in% c(PG.PRO.LIST,TG.PRO.LIST))),FUN=function(x) {
  NODE <- MG.PRO.LIST[which(MG.PRO.LIST %in% c(PG.PRO.LIST,TG.PRO.LIST))][x]
  neighbors(MG.ACT.NET,v=NODE)
})
MG.INT_PRO.ALL_GENO.NEIGHBOR.CLASS <- names(unlist(MG.INT_PRO.ALL_GENO.NEIGHBOR.CLASS))
MG.INT_PRO.ALL_GENO.NEIGHBOR.CLASS <- sapply(MG.INT_PRO.ALL_GENO.NEIGHBOR.CLASS,FUN=get.genotype.class,NETWORKTYPE="M")


PG.ALL_PRO.ALL_GENO.NEIGHBOR.CLASS <- lapply(1:length(PG.PRO.LIST),FUN=function(x) {
  NODE <- PG.PRO.LIST[x]
  neighbors(PG.ACT.NET,v=NODE)
})
PG.ALL_PRO.ALL_GENO.NEIGHBOR.CLASS <- names(unlist(PG.ALL_PRO.ALL_GENO.NEIGHBOR.CLASS))
PG.ALL_PRO.ALL_GENO.NEIGHBOR.CLASS <- sapply(PG.ALL_PRO.ALL_GENO.NEIGHBOR.CLASS,FUN=get.genotype.class,NETWORKTYPE="P")

PG.INT_PRO.ALL_GENO.NEIGHBOR.CLASS <- lapply(1:length(which(PG.PRO.LIST %in% c(MG.PRO.LIST,TG.PRO.LIST))),FUN=function(x) {
  NODE <- PG.PRO.LIST[which(PG.PRO.LIST %in% c(MG.PRO.LIST,TG.PRO.LIST))][x]
  neighbors(PG.ACT.NET,v=NODE)
})
PG.INT_PRO.ALL_GENO.NEIGHBOR.CLASS <- names(unlist(PG.INT_PRO.ALL_GENO.NEIGHBOR.CLASS))
PG.INT_PRO.ALL_GENO.NEIGHBOR.CLASS <- sapply(PG.INT_PRO.ALL_GENO.NEIGHBOR.CLASS,FUN=get.genotype.class,NETWORKTYPE="P")


TG.ALL_PRO.ALL_GENO.NEIGHBOR.CLASS <- lapply(1:length(TG.PRO.LIST),FUN=function(x) {
  NODE <- TG.PRO.LIST[x]
  neighbors(TG.ACT.NET,v=NODE)
})
TG.ALL_PRO.ALL_GENO.NEIGHBOR.CLASS <- names(unlist(TG.ALL_PRO.ALL_GENO.NEIGHBOR.CLASS))
TG.ALL_PRO.ALL_GENO.NEIGHBOR.CLASS <- sapply(TG.ALL_PRO.ALL_GENO.NEIGHBOR.CLASS,FUN=get.genotype.class,NETWORKTYPE="T")

TG.INT_PRO.ALL_GENO.NEIGHBOR.CLASS <- lapply(1:length(which(TG.PRO.LIST %in% c(MG.PRO.LIST,PG.PRO.LIST))),FUN=function(x) {
  NODE <- TG.PRO.LIST[which(TG.PRO.LIST %in% c(MG.PRO.LIST,PG.PRO.LIST))][x]
  neighbors(TG.ACT.NET,v=NODE)
})
TG.INT_PRO.ALL_GENO.NEIGHBOR.CLASS <- names(unlist(TG.INT_PRO.ALL_GENO.NEIGHBOR.CLASS))
TG.INT_PRO.ALL_GENO.NEIGHBOR.CLASS <- sapply(TG.INT_PRO.ALL_GENO.NEIGHBOR.CLASS,FUN=get.genotype.class,NETWORKTYPE="T")

##Calculate class of SRE neighbors
MG.ALL_SRE.ALL_GENO.NEIGHBOR.CLASS <- lapply(1:length(MG.SRE.LIST),FUN=function(x) {
  NODE <- MG.SRE.LIST[x]
  neighbors(MG.ACT.NET,v=NODE)
})
MG.ALL_SRE.ALL_GENO.NEIGHBOR.CLASS <- names(unlist(MG.ALL_SRE.ALL_GENO.NEIGHBOR.CLASS))
MG.ALL_SRE.ALL_GENO.NEIGHBOR.CLASS <- sapply(MG.ALL_SRE.ALL_GENO.NEIGHBOR.CLASS,FUN=get.genotype.class,NETWORKTYPE="M")

MG.INT_SRE.ALL_GENO.NEIGHBOR.CLASS <- lapply(1:length(which(MG.SRE.LIST %in% c(PG.SRE.LIST,TG.SRE.LIST))),FUN=function(x) {
  NODE <- MG.SRE.LIST[which(MG.SRE.LIST %in% c(PG.SRE.LIST,TG.SRE.LIST))][x]
  neighbors(MG.ACT.NET,v=NODE)
})
MG.INT_SRE.ALL_GENO.NEIGHBOR.CLASS <- names(unlist(MG.INT_SRE.ALL_GENO.NEIGHBOR.CLASS))
MG.INT_SRE.ALL_GENO.NEIGHBOR.CLASS <- sapply(MG.INT_SRE.ALL_GENO.NEIGHBOR.CLASS,FUN=get.genotype.class,NETWORKTYPE="M")


PG.ALL_SRE.ALL_GENO.NEIGHBOR.CLASS <- lapply(1:length(PG.SRE.LIST),FUN=function(x) {
  NODE <- PG.SRE.LIST[x]
  neighbors(PG.ACT.NET,v=NODE)
})
PG.ALL_SRE.ALL_GENO.NEIGHBOR.CLASS <- names(unlist(PG.ALL_SRE.ALL_GENO.NEIGHBOR.CLASS))
PG.ALL_SRE.ALL_GENO.NEIGHBOR.CLASS <- sapply(PG.ALL_SRE.ALL_GENO.NEIGHBOR.CLASS,FUN=get.genotype.class,NETWORKTYPE="P")

PG.INT_SRE.ALL_GENO.NEIGHBOR.CLASS <- lapply(1:length(which(PG.SRE.LIST %in% c(MG.SRE.LIST,TG.SRE.LIST))),FUN=function(x) {
  NODE <- PG.SRE.LIST[which(PG.SRE.LIST %in% c(MG.SRE.LIST,TG.SRE.LIST))][x]
  neighbors(PG.ACT.NET,v=NODE)
})
PG.INT_SRE.ALL_GENO.NEIGHBOR.CLASS <- names(unlist(PG.INT_SRE.ALL_GENO.NEIGHBOR.CLASS))
PG.INT_SRE.ALL_GENO.NEIGHBOR.CLASS <- sapply(PG.INT_SRE.ALL_GENO.NEIGHBOR.CLASS,FUN=get.genotype.class,NETWORKTYPE="P")


TG.ALL_SRE.ALL_GENO.NEIGHBOR.CLASS <- lapply(1:length(TG.SRE.LIST),FUN=function(x) {
  NODE <- TG.SRE.LIST[x]
  neighbors(TG.ACT.NET,v=NODE)
})
TG.ALL_SRE.ALL_GENO.NEIGHBOR.CLASS <- names(unlist(TG.ALL_SRE.ALL_GENO.NEIGHBOR.CLASS))
TG.ALL_SRE.ALL_GENO.NEIGHBOR.CLASS <- sapply(TG.ALL_SRE.ALL_GENO.NEIGHBOR.CLASS,FUN=get.genotype.class,NETWORKTYPE="T")

TG.INT_SRE.ALL_GENO.NEIGHBOR.CLASS <- lapply(1:length(which(TG.SRE.LIST %in% c(MG.SRE.LIST,PG.SRE.LIST))),FUN=function(x) {
  NODE <- TG.SRE.LIST[which(TG.SRE.LIST %in% c(MG.SRE.LIST,PG.SRE.LIST))][x]
  neighbors(TG.ACT.NET,v=NODE)
})
TG.INT_SRE.ALL_GENO.NEIGHBOR.CLASS <- names(unlist(TG.INT_SRE.ALL_GENO.NEIGHBOR.CLASS))
TG.INT_SRE.ALL_GENO.NEIGHBOR.CLASS <- sapply(TG.INT_SRE.ALL_GENO.NEIGHBOR.CLASS,FUN=get.genotype.class,NETWORKTYPE="T")




#Collect results
ERE.NEIGHBOR.CLASS.COUNT <- matrix(0,nrow=9,ncol=9)
colnames(ERE.NEIGHBOR.CLASS.COUNT) <- c("E.COUNT","P.COUNT","S.COUNT","T.COUNT","N.GENO","E.DEGREE","P.DEGREE","S.DEGREE","T.DEGREE")
rownames(ERE.NEIGHBOR.CLASS.COUNT) <- c("MG.ALL","PG.ALL","TG.ALL","MG.INT","PG.INT","TG.INT","MG.UNI","PG.UNI","TG.UNI")

ERE.NEIGHBOR.CLASS.COUNT[1,1:3] <- table(MG.ALL_ERE.ALL_GENO.NEIGHBOR.CLASS)[c(1,3,4)]
ERE.NEIGHBOR.CLASS.COUNT[2,1:3] <- table(PG.ALL_ERE.ALL_GENO.NEIGHBOR.CLASS)[c(1,3,4)]
ERE.NEIGHBOR.CLASS.COUNT[3,1:3] <- table(TG.ALL_ERE.ALL_GENO.NEIGHBOR.CLASS)[c(1,3,4)]
ERE.NEIGHBOR.CLASS.COUNT[4,1:3] <- table(MG.INT_ERE.ALL_GENO.NEIGHBOR.CLASS)[c(1,3,4)]
ERE.NEIGHBOR.CLASS.COUNT[5,1:3] <- table(PG.INT_ERE.ALL_GENO.NEIGHBOR.CLASS)[c(1,3,4)]
ERE.NEIGHBOR.CLASS.COUNT[6,1:3] <- table(TG.INT_ERE.ALL_GENO.NEIGHBOR.CLASS)[c(1,3,4)]
ERE.NEIGHBOR.CLASS.COUNT[1,5] <- length(MG.ERE.LIST)
ERE.NEIGHBOR.CLASS.COUNT[2,5] <- length(PG.ERE.LIST)
ERE.NEIGHBOR.CLASS.COUNT[3,5] <- length(TG.ERE.LIST)
ERE.NEIGHBOR.CLASS.COUNT[4,5] <- sum(MG.ERE.LIST %in% c(PG.ERE.LIST,TG.ERE.LIST))
ERE.NEIGHBOR.CLASS.COUNT[5,5] <- sum(PG.ERE.LIST %in% c(MG.ERE.LIST,TG.ERE.LIST))
ERE.NEIGHBOR.CLASS.COUNT[6,5] <- sum(TG.ERE.LIST %in% c(MG.ERE.LIST,PG.ERE.LIST))
ERE.NEIGHBOR.CLASS.COUNT[7,c(1:3,5)] <- ERE.NEIGHBOR.CLASS.COUNT[1,c(1:3,5)] - ERE.NEIGHBOR.CLASS.COUNT[4,c(1:3,5)]
ERE.NEIGHBOR.CLASS.COUNT[8,c(1:3,5)] <- ERE.NEIGHBOR.CLASS.COUNT[2,c(1:3,5)] - ERE.NEIGHBOR.CLASS.COUNT[5,c(1:3,5)]
ERE.NEIGHBOR.CLASS.COUNT[9,c(1:3,5)] <- ERE.NEIGHBOR.CLASS.COUNT[3,c(1:3,5)] - ERE.NEIGHBOR.CLASS.COUNT[6,c(1:3,5)]
ERE.NEIGHBOR.CLASS.COUNT[,4] <- rowSums(ERE.NEIGHBOR.CLASS.COUNT[,1:3])
ERE.NEIGHBOR.CLASS.COUNT[,6:9] <- ERE.NEIGHBOR.CLASS.COUNT[,1:4]/ERE.NEIGHBOR.CLASS.COUNT[,5]


PRO.NEIGHBOR.CLASS.COUNT <- matrix(0,nrow=9,ncol=9)
colnames(PRO.NEIGHBOR.CLASS.COUNT) <- c("E.COUNT","P.COUNT","S.COUNT","T.COUNT","N.GENO","E.DEGREE","P.DEGREE","S.DEGREE","T.DEGREE")
rownames(PRO.NEIGHBOR.CLASS.COUNT) <- c("MG.ALL","PG.ALL","TG.ALL","MG.INT","PG.INT","TG.INT","MG.UNI","PG.UNI","TG.UNI")

PRO.NEIGHBOR.CLASS.COUNT[1,1:3] <- table(MG.ALL_PRO.ALL_GENO.NEIGHBOR.CLASS)[c(1,3,4)]
PRO.NEIGHBOR.CLASS.COUNT[2,1:3] <- table(PG.ALL_PRO.ALL_GENO.NEIGHBOR.CLASS)[c(1,3,4)]
PRO.NEIGHBOR.CLASS.COUNT[3,1:3] <- table(TG.ALL_PRO.ALL_GENO.NEIGHBOR.CLASS)[c(1,3,4)]
PRO.NEIGHBOR.CLASS.COUNT[4,1:3] <- table(MG.INT_PRO.ALL_GENO.NEIGHBOR.CLASS)[c(1,3,4)]
PRO.NEIGHBOR.CLASS.COUNT[5,1:3] <- table(PG.INT_PRO.ALL_GENO.NEIGHBOR.CLASS)[c(1,3,4)]
PRO.NEIGHBOR.CLASS.COUNT[6,1:3] <- table(TG.INT_PRO.ALL_GENO.NEIGHBOR.CLASS)[c(1,3,4)]
PRO.NEIGHBOR.CLASS.COUNT[1,5] <- length(MG.PRO.LIST)
PRO.NEIGHBOR.CLASS.COUNT[2,5] <- length(PG.PRO.LIST)
PRO.NEIGHBOR.CLASS.COUNT[3,5] <- length(TG.PRO.LIST)
PRO.NEIGHBOR.CLASS.COUNT[4,5] <- sum(MG.PRO.LIST %in% c(PG.PRO.LIST,TG.PRO.LIST))
PRO.NEIGHBOR.CLASS.COUNT[5,5] <- sum(PG.PRO.LIST %in% c(MG.PRO.LIST,TG.PRO.LIST))
PRO.NEIGHBOR.CLASS.COUNT[6,5] <- sum(TG.PRO.LIST %in% c(MG.PRO.LIST,PG.PRO.LIST))
PRO.NEIGHBOR.CLASS.COUNT[7,c(1:3,5)] <- PRO.NEIGHBOR.CLASS.COUNT[1,c(1:3,5)] - PRO.NEIGHBOR.CLASS.COUNT[4,c(1:3,5)]
PRO.NEIGHBOR.CLASS.COUNT[8,c(1:3,5)] <- PRO.NEIGHBOR.CLASS.COUNT[2,c(1:3,5)] - PRO.NEIGHBOR.CLASS.COUNT[5,c(1:3,5)]
PRO.NEIGHBOR.CLASS.COUNT[9,c(1:3,5)] <- PRO.NEIGHBOR.CLASS.COUNT[3,c(1:3,5)] - PRO.NEIGHBOR.CLASS.COUNT[6,c(1:3,5)]
PRO.NEIGHBOR.CLASS.COUNT[,4] <- rowSums(PRO.NEIGHBOR.CLASS.COUNT[,1:3])
PRO.NEIGHBOR.CLASS.COUNT[,6:9] <- PRO.NEIGHBOR.CLASS.COUNT[,1:4]/PRO.NEIGHBOR.CLASS.COUNT[,5]


SRE.NEIGHBOR.CLASS.COUNT <- matrix(0,nrow=9,ncol=9)
colnames(SRE.NEIGHBOR.CLASS.COUNT) <- c("E.COUNT","P.COUNT","S.COUNT","T.COUNT","N.GENO","E.DEGREE","P.DEGREE","S.DEGREE","T.DEGREE")
rownames(SRE.NEIGHBOR.CLASS.COUNT) <- c("MG.ALL","PG.ALL","TG.ALL","MG.INT","PG.INT","TG.INT","MG.UNI","PG.UNI","TG.UNI")

SRE.NEIGHBOR.CLASS.COUNT[1,1:3] <- table(MG.ALL_SRE.ALL_GENO.NEIGHBOR.CLASS)[c(1,3,4)]
SRE.NEIGHBOR.CLASS.COUNT[2,1:3] <- table(PG.ALL_SRE.ALL_GENO.NEIGHBOR.CLASS)[c(1,3,4)]
SRE.NEIGHBOR.CLASS.COUNT[3,1:3] <- table(TG.ALL_SRE.ALL_GENO.NEIGHBOR.CLASS)[c(1,3,4)]
SRE.NEIGHBOR.CLASS.COUNT[4,1:3] <- table(MG.INT_SRE.ALL_GENO.NEIGHBOR.CLASS)[c(1,3,4)]
SRE.NEIGHBOR.CLASS.COUNT[5,1:3] <- table(PG.INT_SRE.ALL_GENO.NEIGHBOR.CLASS)[c(1,3,4)]
SRE.NEIGHBOR.CLASS.COUNT[6,1:3] <- table(TG.INT_SRE.ALL_GENO.NEIGHBOR.CLASS)[c(1,3,4)]
SRE.NEIGHBOR.CLASS.COUNT[1,5] <- length(MG.SRE.LIST)
SRE.NEIGHBOR.CLASS.COUNT[2,5] <- length(PG.SRE.LIST)
SRE.NEIGHBOR.CLASS.COUNT[3,5] <- length(TG.SRE.LIST)
SRE.NEIGHBOR.CLASS.COUNT[4,5] <- sum(MG.SRE.LIST %in% c(PG.SRE.LIST,TG.SRE.LIST))
SRE.NEIGHBOR.CLASS.COUNT[5,5] <- sum(PG.SRE.LIST %in% c(MG.SRE.LIST,TG.SRE.LIST))
SRE.NEIGHBOR.CLASS.COUNT[6,5] <- sum(TG.SRE.LIST %in% c(MG.SRE.LIST,PG.SRE.LIST))
SRE.NEIGHBOR.CLASS.COUNT[7,c(1:3,5)] <- SRE.NEIGHBOR.CLASS.COUNT[1,c(1:3,5)] - SRE.NEIGHBOR.CLASS.COUNT[4,c(1:3,5)]
SRE.NEIGHBOR.CLASS.COUNT[8,c(1:3,5)] <- SRE.NEIGHBOR.CLASS.COUNT[2,c(1:3,5)] - SRE.NEIGHBOR.CLASS.COUNT[5,c(1:3,5)]
SRE.NEIGHBOR.CLASS.COUNT[9,c(1:3,5)] <- SRE.NEIGHBOR.CLASS.COUNT[3,c(1:3,5)] - SRE.NEIGHBOR.CLASS.COUNT[6,c(1:3,5)]
SRE.NEIGHBOR.CLASS.COUNT[,4] <- rowSums(SRE.NEIGHBOR.CLASS.COUNT[,1:3])
SRE.NEIGHBOR.CLASS.COUNT[,6:9] <- SRE.NEIGHBOR.CLASS.COUNT[,1:4]/SRE.NEIGHBOR.CLASS.COUNT[,5]

#pdf("Degree.pdf")
par(mfrow=c(3,2))
barplot(t(ERE.NEIGHBOR.CLASS.COUNT[,6:8]),col=c("purple","skyblue2","green"),main="ERE Neighbors",ylim=c(0,12),cex.names=.3)
barplot(t(ERE.NEIGHBOR.CLASS.COUNT[,6:8]/rowSums(ERE.NEIGHBOR.CLASS.COUNT[,6:8])),col=c("purple","skyblue2","green"),main="ERE Neighbors",ylim=c(0,1),cex.names=.3)
barplot(t(PRO.NEIGHBOR.CLASS.COUNT[,6:8]),col=c("purple","skyblue2","green"),main="PRO Neighbors",ylim=c(0,12),cex.names=.3)
barplot(t(PRO.NEIGHBOR.CLASS.COUNT[,6:8]/rowSums(PRO.NEIGHBOR.CLASS.COUNT[,6:8])),col=c("purple","skyblue2","green"),main="PRO Neighbors",ylim=c(0,1),cex.names=.3)
barplot(t(SRE.NEIGHBOR.CLASS.COUNT[,6:8]),col=c("purple","skyblue2","green"),main="SRE Neighbors",ylim=c(0,12),cex.names=.3)
barplot(t(SRE.NEIGHBOR.CLASS.COUNT[,6:8]/rowSums(SRE.NEIGHBOR.CLASS.COUNT[,6:8])),col=c("purple","skyblue2","green"),main="SRE Neighbors",ylim=c(0,1),cex.names=.3)

#dev.off()
