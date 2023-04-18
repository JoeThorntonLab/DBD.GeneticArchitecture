################################################################
##Compute penalized ordinal logistic regressions to decompose ## 
##SR seq-function landscape into its genetic determinants     ##
################################################################
##Single comment code was run and saved as an .Rdata object for faster loading or as a pdf

###########################
##1. Packages & Functions##
###########################
####Package #Version
library(data.table)   #1.11.8
library(ggseqlogo)    #0.1
library(Formula)      #1.2-3
library(survival)     #2.42-6
library(lattice)      #0.20-35
library(ggplot2)      #3.0.0
library(Hmisc)        #4.1-1
library(MatrixModels) #0.4-1
library(boot)         #1.3-20 
library(Matrix)       #1.2-15
library(ordinalNet)   #2.5 
library(pROC)         #1.16.2
library(glmnet)       #2.0-16
library(foreach)      #1.4.4
library(glmnetcr)     #1.0.4 
library(dotCall64)    #1.0-0


library(gtools)
library(gplots)
library(gridExtra)
library(grid)
library(seqinr)
library(stringr)
library(vioplot)
library(rgexf)
library(igraph)
library(nloptr)
library(egg)
library(moments)
library(arules)
library(ggExtra)
library(textTinyR)

####Custom Functions
setwd("/Users/Chimera/Documents/ThorntonLab/Epistasis/2020")
source("Functions.R")

############
##2. Data ##
############
####Amino Acid Data
##Amino acid states
AAs <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")

##Extended set of amino acids that splits S into two groups (S and Z) for analyses that take account of the genetic code
AAs.extend <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","Z")

##Amino acid neighbors and distances with (S or Z) and without (N) the genetic code
AA.NEIGHBORS.S <- read.table(file="AA.neighbors.S.txt",header=T,row.names=1,quote="")
AA.DISTANCES.S <- read.table(file="AA.distance.S.txt", header=T,row.names=1,quote="")
AA.NEIGHBORS.Z <- read.table(file="AA.neighbors.Z.txt",header=T,row.names=1,quote="")
AA.DISTANCES.Z <- read.table(file="AA.distance.Z.txt", header=T,row.names=1,quote="")
AA.NEIGHBORS.N <- read.table(file="AA.neighbors.N.txt",header=T,row.names=1,quote="")
AA.DISTANCES.N <- read.table(file="AA.distance.N.txt", header=T,row.names=1,quote="")

####Deep Mutational Scanning Data from Starr et al. 2017
##Original data
load("DT.11P.CODING.rda")
##Clean up data table columns
#DT <- dt.11P.coding[,.(AAseq,AA1,AA2,AA3,AA4,ERE.pooled.cfu,SRE.pooled.cfu,ERE.pooled.meanF,SRE.pooled.meanF,ERE.SE.meanF,SRE.SE.meanF,ERE.pooled.class,SRE.pooled.class,ERE.prediction,SRE.prediction,ERE.full.class,SRE.full.class,specificity)]; rm(dt.11P.coding)
##Add columns giving amino acid state subsets
#DT[,AA12:=as.factor(paste(AA1,AA2,sep=""))];DT[,AA13:=as.factor(paste(AA1,AA3,sep=""))];DT[,AA14:=as.factor(paste(AA1,AA4,sep=""))];DT[,AA23:=as.factor(paste(AA2,AA3,sep=""))];DT[,AA24:=as.factor(paste(AA2,AA4,sep=""))];DT[,AA34:=as.factor(paste(AA3,AA4,sep=""))];DT[,AA123:=as.factor(paste(AA1,AA2,AA3,sep=""))];DT[,AA124:=as.factor(paste(AA1,AA2,AA4,sep=""))];DT[,AA134:=as.factor(paste(AA1,AA3,AA4,sep=""))];DT[,AA234:=as.factor(paste(AA2,AA3,AA4,sep=""))]
##Set as NA for SRE.pooled.meanF and ERE.pooled.meanF if xRE.pooled.cfu < 15
#DT[ERE.pooled.cfu<15,ERE.pooled.meanF:=NA];DT[SRE.pooled.cfu<15,SRE.pooled.meanF:=NA]
##Remove columns containing prior class predictions
#KEEP <- which(!(names(DT) %in% c("ERE.prediction","SRE.prediction","ERE.full.class","SRE.full.class")))
#DT <- DT[,..KEEP]
##Recode specificity for observed data
#DT$specificity <- "null"
##Promiscuous sequences are strong on both REs
#DT$specificity[DT$ERE.pooled.class == "strong" | DT$SRE.pooled.class == "strong"] <- "promiscuous"
##ERE-specific sequences are strong on ERE and not strong on SRE
#DT$specificity[DT$ERE.pooled.class == "strong" & DT$SRE.pooled.class != "strong"] <- "ERE-specific"
##SRE-specific sequences are strong on SRE and not strong on ERE
#DT$specificity[DT$ERE.pooled.class != "strong" & DT$SRE.pooled.class == "strong"] <- "SRE-specific"
##Specificity is NA if binding to either RE is missing
#DT$specificity[is.na(DT$ERE.pooled.class) | is.na(DT$SRE.pooled.class)] <- NA
##Cleanup data sets and reorder columns
#DT[,ERE.pooled.class:=factor(ERE.pooled.class,levels=c("null","weak","strong"))]
#DT[,SRE.pooled.class:=factor(SRE.pooled.class,levels=c("null","weak","strong"))]
#DT <- droplevels(DT)
#setcolorder(DT, c("AAseq","AA1","AA2","AA3","AA4","AA12","AA13","AA14","AA23","AA24","AA34","AA123","AA124","AA134","AA234","ERE.pooled.cfu","SRE.pooled.cfu","ERE.pooled.meanF","SRE.pooled.meanF","ERE.SE.meanF","SRE.SE.meanF","ERE.pooled.class","SRE.pooled.class","specificity"))
##Create a joint dataset that contains the ERE and SRE data as separate entries (2 entries per genotype, one ERE one SRE)
##ERE Data
#DT.JOINT.ERE <- as.data.frame(DT[,c("AAseq","AA1","AA2","AA3","AA4","ERE.pooled.cfu","ERE.pooled.meanF","ERE.SE.meanF","ERE.pooled.class")])
#DT.JOINT.ERE$RE <- "E"
#DT.JOINT.ERE$AAseq <- as.character(DT.JOINT.ERE$AAseq)
#DT.JOINT.ERE$AAseq <- paste0("E",DT.JOINT.ERE$AAseq)
#colnames(DT.JOINT.ERE ) <- c("AAseq","AA1","AA2","AA3","AA4","pooled.cfu","pooled.meanF","SE.meanF","pooled.class","RE")
##SRE Data
#DT.JOINT.SRE <- as.data.frame(DT[,c("AAseq","AA1","AA2","AA3","AA4","SRE.pooled.cfu","SRE.pooled.meanF","SRE.SE.meanF","SRE.pooled.class")])
#DT.JOINT.SRE$RE <- "S" 
#DT.JOINT.SRE$AAseq <- as.character(DT.JOINT.SRE$AAseq)
#DT.JOINT.SRE$AAseq <- paste0("S",DT.JOINT.SRE$AAseq)
#colnames(DT.JOINT.SRE) <- c("AAseq","AA1","AA2","AA3","AA4","pooled.cfu","pooled.meanF","SE.meanF","pooled.class","RE")
##Joint ERE and SRE Data
#DT.JOINT <- data.table(rbind(DT.JOINT.ERE,DT.JOINT.SRE));setkey(DT,AAseq)
#DT.JOINT$RE <- as.factor(DT.JOINT$RE)
#rm(DT.JOINT.ERE );rm(DT.JOINT.SRE)
#save(DT.JOINT,file="DT.JOINT.rda")

####Load cleaned-up data
load("DT.JOINT.rda")

####Distribution of mean fluorescence 
#par(mfrow=c(2,3))
#hist(subset(DT.JOINT,DT.JOINT$RE=="E" & DT.JOINT$pooled.class=="null"  )$pooled.meanF,xlab="E null",  breaks=seq(1,10,0.1));abline(v=6.8,col="red")
#hist(subset(DT.JOINT,DT.JOINT$RE=="E" & DT.JOINT$pooled.class=="weak"  )$pooled.meanF,xlab="E weak",  breaks=seq(1,10,0.1));abline(v=6.8,col="red")
#hist(subset(DT.JOINT,DT.JOINT$RE=="E" & DT.JOINT$pooled.class=="strong")$pooled.meanF,xlab="E strong",breaks=seq(1,10,0.1));abline(v=6.8,col="red")
#hist(subset(DT.JOINT,DT.JOINT$RE=="S" & DT.JOINT$pooled.class=="null"  )$pooled.meanF,xlab="S null",  breaks=seq(1,10,0.1));abline(v=6.8,col="red")
#hist(subset(DT.JOINT,DT.JOINT$RE=="S" & DT.JOINT$pooled.class=="weak"  )$pooled.meanF,xlab="S weak",  breaks=seq(1,10,0.1));abline(v=6.8,col="red")
#hist(subset(DT.JOINT,DT.JOINT$RE=="S" & DT.JOINT$pooled.class=="strong")$pooled.meanF,xlab="S strong",breaks=seq(1,10,0.1));abline(v=6.8,col="red")

####Proportional odds assumption (Difference in effect across classes should be constant for all amino acids)
##Convert class to integer
#DT.PO.TEST <- copy(dt.11P.coding)
#DT.PO.TEST$ERE.pooled.class <- as.integer(DT.PO.TEST$ERE.pooled.class) - 1
#DT.PO.TEST$SRE.pooled.class <- as.integer(DT.PO.TEST$SRE.pooled.class) - 1

#pdf(file="RE.PO.ASSUMPTION.pdf",4,11,useDingbats = F)
##Fit main effect model on ERE across classes 
#ERE.PO.TEST <- summary(ERE.pooled.class ~ AA1 + AA2 + AA3 + AA4, fun=test.PO, data=DT.PO.TEST)
#plot(ERE.PO.TEST, which=1:2, pch=1:2, xlim=c(min(as.vector(ERE.PO.TEST)[is.finite(as.vector(ERE.PO.TEST))]), 0), xlab="logit", vnames="names", main="ERE", width.factor=1.5,cex=0.8)
##Fit main effect model on SRE across classes 
#SRE.PO.TEST  <- summary(SRE.pooled.class ~ AA1 + AA2 + AA3 + AA4, fun=test.PO, data=DT.PO.TEST)
#plot(SRE.PO.TEST, which=1:2, pch=1:2, xlim=c(min(as.vector(SRE.PO.TEST)[is.finite(as.vector(SRE.PO.TEST))]), 0), xlab="logit", vnames="names", main="SRE", width.factor=1.5,cex=0.8)
#dev.off()

#####################
##3. Model matrices##
#####################
####Set up amino acids for generating model matrices
##Individual amino acids
#A1 <- AAs; A2 <- AAs; A3 <- AAs; A4 <- AAs

##Amino acid present in each sequence
#AA.SEQ1 <- matrix(0,nrow=2*length(A1)*length(A2)*length(A3)*length(A4),ncol=5)
#AA.SEQ1[,1] <- rep(c("E","S"),each=length(A1)*length(A2)*length(A3)*length(A4))
#AA.SEQ1[,2] <- rep(rep(A1,each=length(A2)*length(A3)*length(A4)),2)
#AA.SEQ1[,3] <- rep(rep(rep(A2,each=length(A3)*length(A4)),length(A1)),2)
#AA.SEQ1[,4] <- rep(rep(rep(A3,each=length(A4)),length(A1)*length(A2)),2)
#AA.SEQ1[,5] <- rep(rep(A4,length(A1)*length(A2)*length(A3)),2)

##One-hot encoding of amino acid states at each site
#AA.SEQ2 <- matrix(0,nrow=2*length(A1)*length(A2)*length(A3)*length(A4),ncol=length(A1)+length(A2)+length(A3)+length(A4))
#I <- 1:(length(A1)+length(A2)+length(A3)+length(A4))
#for(i in I) {
#  if(i %in% 1:length(A1)) {
#    j <- i
#    AA.SEQ2[,i] <- rep(c(rep(0,(j-1)*length(A2)*length(A3)*length(A4)),rep(1,length(A2)*length(A3)*length(A4)),rep(0,(length(A1) - j)*length(A2)*length(A3)*length(A4))),2)
#  }
#  if((i-length(A1)) %in% 1:length(A2)) {
#    j <- i - length(A1)
#    AA.SEQ2[,i] <- rep(rep(c(rep(0,(j-1)*length(A3)*length(A4)),rep(1,length(A3)*length(A4)),rep(0,(length(A2) - j)*length(A3)*length(A4))),length(A1)),2)
#  }
#  if((i-(length(A1)+length(A2))) %in% 1:length(A3)) {
#    j <- i - length(A1) - length(A2)
#    AA.SEQ2[,i] <- rep(rep(c(rep(0,(j-1)*length(A4)),rep(1,length(A4)),rep(0,(length(A3) - j)*length(A4))),length(A1)*length(A2)),2)
#  }
#  if((i-(length(A1)+length(A2)+length(A3))) %in% 1:length(A4)) {
#    j <- i - length(A1) - length(A2) - length(A3)
#    AA.SEQ2[,i] <- rep(rep(c(rep(0,(j-1)),rep(1,1),rep(0,(length(A4) - j))),length(A1)*length(A2)*length(A3)),2)
#  }
#}
#AA.SEQ <- cbind(data.frame(AA.SEQ1),AA.SEQ2)
#rownames(AA.SEQ) <- paste(AA.SEQ[,1],AA.SEQ[,2],AA.SEQ[,3],AA.SEQ[,4],AA.SEQ[,5],sep="")
#colnames(AA.SEQ) <- c("RE","X1","X2","X3","X4",paste("X1",A1,sep=""),paste("X2",A2,sep=""),paste("X3",A3,sep=""),paste("X4",A4,sep=""))

##Set contrasts to center levels at average
#contrasts(AA.SEQ$RE) <- contr.sum
#attributes(AA.SEQ$X1)$contrasts <- contrasts(AA.SEQ$X1,contrasts=F)
#attributes(AA.SEQ$X2)$contrasts <- contrasts(AA.SEQ$X2,contrasts=F)
#attributes(AA.SEQ$X3)$contrasts <- contrasts(AA.SEQ$X3,contrasts=F)
#attributes(AA.SEQ$X4)$contrasts <- contrasts(AA.SEQ$X4,contrasts=F)
#save(AA.SEQ,file="AA.SEQ.rda")

load("AA.SEQ.rda")

####Generate model matrices
##Model matrix for main effects model without epistasis
#ME.MATRIX  <- model.Matrix(~ 0+(X1+X2+X3+X4)*RE,data=AA.SEQ,sparse=TRUE)
##Move global specificity term to last column
#ME.MATRIX <- ME.MATRIX[,c(1:80,82:161,81)]
#save(ME.MATRIX, file = "ME.MATRIX.rda")

##Model matrix for model with pairwise epistasis
#PE.MATRIX  <- model.Matrix(~ 0+(X1+X2+X3+X4)^2*RE,data=AA.SEQ,sparse=TRUE)
##Move global specificity term to last column
#PE.MATRIX <- PE.MATRIX[,c(1:80,2482:2561,82:2481,2562:4961,81)]
#save(PE.MATRIX, file = "PE.MATRIX.rda")
##Model matrix for model with third-order epistasis.
#TE.MATRIX.S1  <- model.Matrix(~ 0+(X1+X2+X3+X4)^2*RE + (X1+X2+X3+X4)^3,data=AA.SEQ,sparse=TRUE)
##Third order for ERE
#TE.MATRIX.S2 <- TE.MATRIX.S1[1:160000,4962:36961]
##Third order for SRE
#TE.MATRIX.S3 <- -1*TE.MATRIX.S1[160001:320000,4962:36961]
#TE.MATRIX.S4 <- rbind(TE.MATRIX.S2,TE.MATRIX.S3)
#colnames(TE.MATRIX.S4) <- paste0(colnames(TE.MATRIX.S4),":RE1")
##Combine binding and specificity
#TE.MATRIX <- cbind(TE.MATRIX.S1,TE.MATRIX.S4)
##Move global specificity term to last column
#TE.MATRIX <- TE.MATRIX[,c(1:80,2482:2561,82:2481,2562:4961,4962:36961,36962:68961,81)]
#save(TE.MATRIX, file = "TE.MATRIX.rda")

####Load model matrices
load("ME.MATRIX.rda")
load("PE.MATRIX.rda")
load("TE.MATRIX.rda")

#########################
##4. Ordinal regression##
#########################
####Fit of ordinal regression model. See functions for details.
####ME-Main Effects. PE-Pairwise Effects. TE-Third order Effects.

####L1 Fits (alpha = 1)
#ME.L1.MODEL <- glmnetpo64(x=ME.MATRIX[!is.na(DT.JOINT$pooled.class),],y=DT.JOINT$pooled.class[!is.na(DT.JOINT$pooled.class)],method="cumulative",nlambda=500,lambda.min.ratio=0.00001,alpha=1,maxit=5000,standardize=FALSE)
#save(ME.L1.MODEL, file = "ME.L1.MODEL.rda")
#PE.L1.MODEL <- glmnetpo64(x=PE.MATRIX[!is.na(DT.JOINT$pooled.class),],y=DT.JOINT$pooled.class[!is.na(DT.JOINT$pooled.class)],method="cumulative",nlambda=500,lambda.min.ratio=0.0001 ,alpha=1,maxit=5000,standardize=FALSE)
#save(PE.L1.MODEL, file = "PE.L1.MODEL.rda")
#TE.L1.MODEL <- glmnetpo64(x=TE.MATRIX[!is.na(DT.JOINT$pooled.class),],y=DT.JOINT$pooled.class[!is.na(DT.JOINT$pooled.class)],method="cumulative",nlambda=500,lambda.min.ratio=0.0001 ,alpha=1,maxit=5000,standardize=FALSE)
#save(TE.L1.MODEL, file = "TE.L1.MODEL.rda")

####L2 fits (alpha = 0)
#ME.L2.MODEL <- glmnetpo64(x=ME.MATRIX[!is.na(DT.JOINT$pooled.class),],y=DT.JOINT$pooled.class[!is.na(DT.JOINT$pooled.class)],method="cumulative",lambda=10^(seq(1,-8,length.out = 900)),alpha=0,maxit=5000,standardize=FALSE)
#save(ME.L2.MODEL, file = "ME.L2.MODEL.rda")
#PE.L2.MODEL <- glmnetpo64(x=PE.MATRIX[!is.na(DT.JOINT$pooled.class),],y=DT.JOINT$pooled.class[!is.na(DT.JOINT$pooled.class)],method="cumulative",lambda=10^(seq(1,-8,length.out = 900)),alpha=0,maxit=5000,standardize=FALSE)
#save(PE.L2.MODEL, file = "PE.L2.MODEL.rda")
#TE.L2.MODEL <- glmnetpo64(x=TE.MATRIX[!is.na(DT.JOINT$pooled.class),],y=DT.JOINT$pooled.class[!is.na(DT.JOINT$pooled.class)],method="cumulative",lambda=10^(seq(1,-8,length.out = 900)),alpha=0,maxit=5000,standardize=FALSE)
#save(TE.L2.MODEL, file = "TE.L2.MODEL.rda")

#######################
##5. Cross-Validation##
#######################
####Cross-validation to find best fit model. Uses 10-fold cross validation to estimate misclassification rate
####Misclassification to/from weak is a single error. Misclassification between null/strong is two errors.
##Set number of CV replicates
I <- 1:10

##L1 Fits (alpha = 1)
#ME.L1.CV <- list()
#for(i in I) {
#  ME.L1.CV[[i]] <- cv.glmnetpo64(x=ME.MATRIX[!is.na(DT.JOINT$pooled.class),],y=DT.JOINT$pooled.class[!is.na(DT.JOINT$pooled.class)],nlambda=500,lambda.min.ratio=0.00001,alpha=1,maxit=5000,standardize=FALSE,nfolds=10,type.measure="class",method="cumulative")
#}
#save(ME.L1.CV, file = "ME.L1.CV.rda")
#PE.L1.CV <- list()
#for(i in I) {
#  PE.L1.CV[[i]] <- cv.glmnetpo64(x=PE.MATRIX[!is.na(DT.JOINT$pooled.class),],y=DT.JOINT$pooled.class[!is.na(DT.JOINT$pooled.class)],nlambda=500,lambda.min.ratio=0.0001, alpha=1,maxit=5000,standardize=FALSE,nfolds=10,type.measure="class",method="cumulative")
#}
#save(PE.L1.CV, file = "PE.L1.CV.rda")
#TE.L1.CV <- list()
#for(i in I) {
#  TE.L1.CV[[i]] <- cv.glmnetpo64(x=TE.MATRIX[!is.na(DT.JOINT$pooled.class),],y=DT.JOINT$pooled.class[!is.na(DT.JOINT$pooled.class)],nlambda=500,lambda.min.ratio=0.0001, alpha=1,maxit=5000,standardize=FALSE,nfolds=10,type.measure="class",method="cumulative")
#}
#save(TE.L1.CV, file = "TE.L1.CV.rda")

##L2 Fits (alpha = 0)
#ME.L2.CV <- list()
#for(i in I) {
#  ME.L2.CV[[i]] <- cv.glmnetpo64(x=ME.MATRIX[!is.na(DT.JOINT$pooled.class),],y=DT.JOINT$pooled.class[!is.na(DT.JOINT$pooled.class)],lambda=10^(seq(1,-8,length.out = 900)),alpha=0,maxit=5000,standardize=FALSE,nfolds=10,type.measure="class",method="cumulative")
#}
#save(ME.L2.CV, file = "ME.L2.CV.rda")
#PE.L2.CV <- list()
#for(i in I) {
#  PE.L2.CV[[i]] <- cv.glmnetpo64(x=PE.MATRIX[!is.na(DT.JOINT$pooled.class),],y=DT.JOINT$pooled.class[!is.na(DT.JOINT$pooled.class)],lambda=10^(seq(1,-8,length.out = 900)),alpha=0,maxit=5000,standardize=FALSE,nfolds=10,type.measure="class",method="cumulative")
#}
#save(PE.L2.CV, file = "PE.L2.CV.rda")
#TE.L2.CV <- list()
#for(i in I) {
#  TE.L2.CV[[i]] <- cv.glmnetpo64(x=TE.MATRIX[!is.na(DT.JOINT$pooled.class),],y=DT.JOINT$pooled.class[!is.na(DT.JOINT$pooled.class)],lambda=10^(seq(1,-8,length.out = 900)),alpha=0,maxit=5000,standardize=FALSE,nfolds=10,type.measure="class",method="cumulative")
#}
#save(TE.L2.CV, file = "TE.L2.CV.rda")

#####################
##6. Penalty choice##
#####################
##Choose L1 or L2 penalty for regularized regression
MODEL <- "L2"
if(MODEL == "L1") {
  load("ME.L1.MODEL.rda")
  load("PE.L1.MODEL.rda")
  load("TE.L1.MODEL.rda")
  load("ME.L1.CV.rda")
  load("PE.L1.CV.rda")
  load("TE.L1.CV.rda")
  ME.MODEL <- ME.L1.MODEL
  PE.MODEL <- PE.L1.MODEL
  TE.MODEL <- TE.L1.MODEL
  ME.CV <- ME.L1.CV
  PE.CV <- PE.L1.CV
  TE.CV <- TE.L1.CV
  rm(ME.L1.CV)
  rm(PE.L1.CV)
  rm(TE.L1.CV)
}else if (MODEL == "L2") {
  load("ME.L2.MODEL.rda")
  load("PE.L2.MODEL.rda")
  load("TE.L2.MODEL.rda")
  load("ME.L2.CV.rda")
  load("PE.L2.CV.rda")
  load("TE.L2.CV.rda")
  ME.MODEL <- ME.L2.MODEL
  PE.MODEL <- PE.L2.MODEL
  TE.MODEL <- TE.L2.MODEL
  ME.CV <- ME.L2.CV
  PE.CV <- PE.L2.CV
  TE.CV <- TE.L2.CV
  rm(ME.L2.CV)
  rm(PE.L2.CV)
  rm(TE.L2.CV)
}

####################
##7. Lambda choice##
####################
####Choose lambda value from cross-validation. Dashed lines are for models with the minimum cross-validated misclassification rate 

#pdf("LAMBDA.CHOICE.PLOTS.pdf")
####ME
##Extract set of lambda values
ME.CV.LAMBDA <- unlist(lapply(ME.CV, function(x) x$lambda))
##Extract misclassification rate 
ME.CV.MU <- unlist(lapply(ME.CV, function(x) x$cvm))
##Average misclassification rate for each lambda
ME.CV.MU.MEAN <- aggregate(ME.CV.MU, list(ME.CV.LAMBDA),mean); colnames(ME.CV.MU.MEAN) <- c("lambda","mean")
##Order average misclassification rate by lambda
ME.CV.MU.MEAN <- ME.CV.MU.MEAN[order(ME.CV.MU.MEAN$lambda, decreasing=TRUE),]
##Standard deviation in misclassification rate for each lambda
ME.CV.MU.SD   <- aggregate(ME.CV.MU, list(ME.CV.LAMBDA),sd); colnames(ME.CV.MU.SD) <- c("lambda","sd")
##Order standard deviation by lambda
ME.CV.MU.SD <- ME.CV.MU.SD[order(ME.CV.MU.SD$lambda, decreasing=TRUE),]
##Remove extreme lambda values
ME.CV.MU.MEAN <- ME.CV.MU.MEAN[complete.cases(ME.CV.MU.SD),]
ME.CV.MU.SD <- ME.CV.MU.SD[complete.cases(ME.CV.MU.SD),]
##Identify lambda with minimum misclassification rate
ME.BEST.LAMBDA.INDEX <- which.min(ME.CV.MU.MEAN$mean)
ME.BEST.LAMBDA <- ME.CV.MU.MEAN$lambda[ME.BEST.LAMBDA.INDEX]

####PE
##Extract set of lambda values
PE.CV.LAMBDA <- unlist(lapply(PE.CV, function(x) x$lambda))
##Extract misclassification rate 
PE.CV.MU <- unlist(lapply(PE.CV, function(x) x$cvm))
##Average misclassification rate for each lambda
PE.CV.MU.MEAN <- aggregate(PE.CV.MU, list(PE.CV.LAMBDA),mean); colnames(PE.CV.MU.MEAN) <- c("lambda","mean")
##Order average misclassification rate by lambda
PE.CV.MU.MEAN <- PE.CV.MU.MEAN[order(PE.CV.MU.MEAN$lambda, decreasing=TRUE),]
##Standard deviation in misclassification rate for each lambda
PE.CV.MU.SD   <- aggregate(PE.CV.MU, list(PE.CV.LAMBDA),sd); colnames(PE.CV.MU.SD) <- c("lambda","sd")
##Order standard deviation by lambda
PE.CV.MU.SD <- PE.CV.MU.SD[order(PE.CV.MU.SD$lambda, decreasing=TRUE),]
##Remove extreme lambda values
PE.CV.MU.MEAN <- PE.CV.MU.MEAN[complete.cases(PE.CV.MU.SD),]
PE.CV.MU.SD <- PE.CV.MU.SD[complete.cases(PE.CV.MU.SD),]
##Identify lambda with minimum misclassification rate
PE.BEST.LAMBDA.INDEX <- which.min(PE.CV.MU.MEAN$mean)
PE.BEST.LAMBDA <- PE.CV.MU.MEAN$lambda[PE.BEST.LAMBDA.INDEX]

####TE
##Extract set of lambda values
TE.CV.LAMBDA <- unlist(lapply(TE.CV, function(x) x$lambda))
##Extract misclassification rate 
TE.CV.MU <- unlist(lapply(TE.CV, function(x) x$cvm))
##Average misclassification rate for each lambda
TE.CV.MU.MEAN <- aggregate(TE.CV.MU, list(TE.CV.LAMBDA),mean); colnames(TE.CV.MU.MEAN) <- c("lambda","mean")
##Order average misclassification rate by lambda
TE.CV.MU.MEAN <- TE.CV.MU.MEAN[order(TE.CV.MU.MEAN$lambda, decreasing=TRUE),]
##Standard deviation in misclassification rate for each lambda
TE.CV.MU.SD   <- aggregate(TE.CV.MU, list(TE.CV.LAMBDA),sd); colnames(TE.CV.MU.SD) <- c("lambda","sd")
##Order standard deviation by lambda
TE.CV.MU.SD <- TE.CV.MU.SD[order(TE.CV.MU.SD$lambda, decreasing=TRUE),]
##Remove extreme lambda values
TE.CV.MU.MEAN <- TE.CV.MU.MEAN[complete.cases(TE.CV.MU.SD),]
TE.CV.MU.SD <- TE.CV.MU.SD[complete.cases(TE.CV.MU.SD),]
##Identify lambda with minimum misclassification rate
TE.BEST.LAMBDA.INDEX <- which.min(TE.CV.MU.MEAN$mean)
TE.BEST.LAMBDA <- TE.CV.MU.MEAN$lambda[TE.BEST.LAMBDA.INDEX]

####Plot
##Cross validation results with error bars
#plot(-log10(TE.CV.MU.MEAN$lambda),TE.CV.MU.MEAN$mean,xlim=c(-1,8),col="purple",xlab="Lambda",ylab="Misclassification rate",type="l")
#error.bars(-log10(TE.CV.MU.SD$lambda),TE.CV.MU.MEAN$mean+1.96*TE.CV.MU.SD$sd/sqrt(length(I)),TE.CV.MU.MEAN$mean-1.96*TE.CV.MU.SD$sd/sqrt(length(I)),width=0.0005)
#error.bars(-log10(PE.CV.MU.SD$lambda),PE.CV.MU.MEAN$mean+1.96*PE.CV.MU.SD$sd/sqrt(length(I)),PE.CV.MU.MEAN$mean-1.96*PE.CV.MU.SD$sd/sqrt(length(I)),width=0.0005)
#error.bars(-log10(ME.CV.MU.SD$lambda),ME.CV.MU.MEAN$mean+1.96*ME.CV.MU.SD$sd/sqrt(length(I)),ME.CV.MU.MEAN$mean-1.96*ME.CV.MU.SD$sd/sqrt(length(I)),width=0.0005)
#points(-log10(TE.CV.MU.MEAN$lambda)[-1],TE.CV.MU.MEAN$mean[-1],pch=19,cex=0.6,col="green",type="l",lwd=2)
#points(-log10(PE.CV.MU.MEAN$lambda)[-1],PE.CV.MU.MEAN$mean[-1],pch=19,cex=0.6,col="blue",type="l",lwd=2)
#points(-log10(ME.CV.MU.MEAN$lambda)[-1],ME.CV.MU.MEAN$mean[-1],pch=19,cex=0.6,col="red",type="l",lwd=2)
#abline(v=-log10(TE.BEST.LAMBDA),lty=3,col="green")
#abline(v=-log10(PE.BEST.LAMBDA),lty=3,col="blue")
#abline(v=-log10(ME.BEST.LAMBDA),lty=3,col="red")

##Cross validation results with 95%CI regions
#plot(  -log10(TE.CV.MU.MEAN$lambda),TE.CV.MU.MEAN$mean,xlim=c(-1,8),col="purple",xlab="Lambda",ylab="Misclassification rate",type="l")
#points(-log10(TE.CV.MU.MEAN$lambda),TE.CV.MU.MEAN$mean,type="l",lwd=2,col="green")
#points(-log10(PE.CV.MU.MEAN$lambda),PE.CV.MU.MEAN$mean,type="l",lwd=2,col="blue")
#points(-log10(ME.CV.MU.MEAN$lambda),ME.CV.MU.MEAN$mean,type="l",lwd=2,col="red")
#points(-log10(TE.CV.MU.SD$lambda),TE.CV.MU.MEAN$mean+1.96*TE.CV.MU.SD$sd/sqrt(length(I)),type="l",lty=2,col="green")
#points(-log10(TE.CV.MU.SD$lambda),TE.CV.MU.MEAN$mean-1.96*TE.CV.MU.SD$sd/sqrt(length(I)),type="l",lty=2,col="green")
#points(-log10(PE.CV.MU.SD$lambda),PE.CV.MU.MEAN$mean+1.96*PE.CV.MU.SD$sd/sqrt(length(I)),type="l",lty=2,col="blue")
#points(-log10(PE.CV.MU.SD$lambda),PE.CV.MU.MEAN$mean-1.96*PE.CV.MU.SD$sd/sqrt(length(I)),type="l",lty=2,col="blue")
#points(-log10(ME.CV.MU.SD$lambda),ME.CV.MU.MEAN$mean+1.96*ME.CV.MU.SD$sd/sqrt(length(I)),type="l",lty=2,col="red")
#points(-log10(ME.CV.MU.SD$lambda),ME.CV.MU.MEAN$mean-1.96*ME.CV.MU.SD$sd/sqrt(length(I)),type="l",lty=2,col="red")
#abline(v=-log10(TE.BEST.LAMBDA),lty=3,col="green")
#abline(v=-log10(PE.BEST.LAMBDA),lty=3,col="blue")
#abline(v=-log10(ME.BEST.LAMBDA),lty=3,col="red")

##Deviance ratio (R^2 analog) for each lambda
#plot(-log10(TE.MODEL$lambda),TE.MODEL$dev.ratio,type="l",col="green",xlab="-log10(Lambda)",ylab="Deviance ratio")
#points(-log10(PE.MODEL$lambda),PE.MODEL$dev.ratio,type="l",col="blue")
#points(-log10(ME.MODEL$lambda),ME.MODEL$dev.ratio,type="l",col="red")
#abline(v=-log10(TE.BEST.LAMBDA),lty=3,col="green")
#abline(v=-log10(PE.BEST.LAMBDA),lty=3,col="blue")
#abline(v=-log10(ME.BEST.LAMBDA),lty=3,col="red")
#abline(h=1)
##dev.off()

#Remove large objects
rm(ME.CV)
rm(PE.CV)
rm(TE.CV)

#############################
##8. Model Predicted values##
#############################
####Predicated activator class
PREDICT.ME.CLASS <- predict.glmnetpo64(ME.MODEL,newx=ME.MATRIX , method="cumulative", s = ME.BEST.LAMBDA.INDEX)
PREDICT.PE.CLASS <- predict.glmnetpo64(PE.MODEL,newx=PE.MATRIX , method="cumulative", s = PE.BEST.LAMBDA.INDEX)
PREDICT.TE.CLASS <- predict.glmnetpo64(TE.MODEL,newx=TE.MATRIX , method="cumulative", s = TE.BEST.LAMBDA.INDEX)

#save(PREDICT.ME.CLASS,file="PREDICT.ME.CLASS.rda")
#save(PREDICT.PE.CLASS,file="PREDICT.PE.CLASS.rda")
#save(PREDICT.TE.CLASS,file="PREDICT.TE.CLASS.rda")

####Predicted Genetic Score
PREDICT.ME.LINK <- logit(PREDICT.ME.CLASS$probs[,1,])  - (ME.MODEL$beta[,ME.BEST.LAMBDA.INDEX][names(ME.MODEL$beta[,ME.BEST.LAMBDA.INDEX]) == "cp1"] + ME.MODEL$a0[ME.BEST.LAMBDA.INDEX])
PREDICT.PE.LINK <- logit(PREDICT.PE.CLASS$probs[,1,])  - (PE.MODEL$beta[,PE.BEST.LAMBDA.INDEX][names(PE.MODEL$beta[,PE.BEST.LAMBDA.INDEX]) == "cp1"] + PE.MODEL$a0[PE.BEST.LAMBDA.INDEX])
PREDICT.TE.LINK <- logit(PREDICT.TE.CLASS$probs[,1,])  - (TE.MODEL$beta[,TE.BEST.LAMBDA.INDEX][names(TE.MODEL$beta[,TE.BEST.LAMBDA.INDEX]) == "cp1"] + TE.MODEL$a0[TE.BEST.LAMBDA.INDEX])

#save(PREDICT.ME.LINK,file="PREDICT.ME.LINK.rda")
#save(PREDICT.PE.LINK,file="PREDICT.PE.LINK.rda")
#save(PREDICT.TE.LINK,file="PREDICT.TE.LINK.rda")

####Compare genetic score to class
#pdf("Linear_Approx.pdf")
##ME
ME.THRESH.NULL <- -1*(ME.MODEL$beta[,ME.BEST.LAMBDA.INDEX][names(ME.MODEL$beta[,ME.BEST.LAMBDA.INDEX]) == "cp1"] + ME.MODEL$a0[ME.BEST.LAMBDA.INDEX])
ME.THRESH.WEAK <- -1*(ME.MODEL$beta[,ME.BEST.LAMBDA.INDEX][names(ME.MODEL$beta[,ME.BEST.LAMBDA.INDEX]) == "cp2"] + ME.MODEL$a0[ME.BEST.LAMBDA.INDEX])
PREDICT.ME.PROBS.NULL   <- inv.logit(seq(-19,10,0.01) - ME.THRESH.NULL)
PREDICT.ME.PROBS.STRONG <- 1 - (inv.logit(seq(-19,10,0.01) - ME.THRESH.WEAK))
PREDICT.ME.PROBS.WEAK <- 1 - PREDICT.ME.PROBS.NULL - PREDICT.ME.PROBS.STRONG

#save(ME.THRESH.NULL,file="ME.THRESH.NULL.rda")
#save(ME.THRESH.WEAK,file="ME.THRESH.WEAK.rda")

#RANGE <- seq(-19,10,0.01)
#plot(  RANGE,PREDICT.ME.PROBS.NULL,  type="l",ylim=c(0,1),xlim=c(-19,-2),ylab="Probability",xlab="Genetic Score",main="ME")
#points(RANGE,PREDICT.ME.PROBS.WEAK,  type="l",col="blue")
#points(RANGE,PREDICT.ME.PROBS.STRONG,type="l",col="red")
#points(RANGE,PREDICT.ME.PROBS.STRONG + PREDICT.ME.PROBS.WEAK,type="l",col="purple")

#PE
PE.THRESH.NULL <- -1*(PE.MODEL$beta[,PE.BEST.LAMBDA.INDEX][names(PE.MODEL$beta[,PE.BEST.LAMBDA.INDEX]) == "cp1"] + PE.MODEL$a0[PE.BEST.LAMBDA.INDEX])
PE.THRESH.WEAK <- -1*(PE.MODEL$beta[,PE.BEST.LAMBDA.INDEX][names(PE.MODEL$beta[,PE.BEST.LAMBDA.INDEX]) == "cp2"] + PE.MODEL$a0[PE.BEST.LAMBDA.INDEX])
PREDICT.PE.PROBS.NULL   <- inv.logit(seq(-19,10,0.01) - PE.THRESH.NULL)
PREDICT.PE.PROBS.STRONG <- 1 - (inv.logit(seq(-19,10,0.01) - PE.THRESH.WEAK))
PREDICT.PE.PROBS.WEAK <- 1 - PREDICT.PE.PROBS.NULL - PREDICT.PE.PROBS.STRONG

#save(PE.THRESH.NULL,file="PE.THRESH.NULL.rda")
#save(PE.THRESH.WEAK,file="PE.THRESH.WEAK.rda")

#plot(  RANGE,PREDICT.PE.PROBS.NULL,  type="l",ylim=c(0,1),xlim=c(-19,-2),ylab="Probability",xlab="Link",main="PE")
#points(RANGE,PREDICT.PE.PROBS.WEAK,  type="l",col="blue")
#points(RANGE,PREDICT.PE.PROBS.STRONG,type="l",col="red")
#points(RANGE,PREDICT.PE.PROBS.STRONG + PREDICT.PE.PROBS.WEAK,type="l",col="purple")

#TE
TE.THRESH.NULL <- -1*(TE.MODEL$beta[,TE.BEST.LAMBDA.INDEX][names(TE.MODEL$beta[,TE.BEST.LAMBDA.INDEX]) == "cp1"] + TE.MODEL$a0[TE.BEST.LAMBDA.INDEX])
TE.THRESH.WEAK <- -1*(TE.MODEL$beta[,TE.BEST.LAMBDA.INDEX][names(TE.MODEL$beta[,TE.BEST.LAMBDA.INDEX]) == "cp2"] + TE.MODEL$a0[TE.BEST.LAMBDA.INDEX])
PREDICT.TE.PROBS.NULL   <- inv.logit(seq(-19,10,0.01) - TE.THRESH.NULL)
PREDICT.TE.PROBS.STRONG <- 1 - (inv.logit(seq(-19,10,0.01) - TE.THRESH.WEAK))
PREDICT.TE.PROBS.WEAK <- 1 - PREDICT.TE.PROBS.NULL - PREDICT.TE.PROBS.STRONG

#save(TE.THRESH.NULL,file="TE.THRESH.NULL.rda")
#save(TE.THRESH.WEAK,file="TE.THRESH.WEAK.rda")

#plot(  RANGE,PREDICT.TE.PROBS.NULL,  type="l",ylim=c(0,1),xlim=c(-19,-2),ylab="Probability",xlab="Link",main="TE")
#points(RANGE,PREDICT.TE.PROBS.WEAK,  type="l",col="blue")
#points(RANGE,PREDICT.TE.PROBS.STRONG,type="l",col="red")
#points(RANGE,PREDICT.TE.PROBS.STRONG + PREDICT.TE.PROBS.WEAK,type="l",col="purple")
#dev.off()

####################################
##9. Comparison with observed data##
####################################
##Vertical lines give cutoffs for the predicted classifications
##Horizontal lines show the 5th percentile for observed mean fluorescence of genotypes in a given class
##Orange = weak/strong cutoff
##Brown == null/weak cutoff

##ME
#png("Misclassification.ME.png",res=500,width=4,height=4,units="in")
#par(pty="s")
#plot(-1*PREDICT.ME.LINK,DT.JOINT$pooled.meanF,pch=19,cex=0.05,col="#000000",xlab="Genetic Score",ylab="Estimated Fluorescence",xlim=c(-20,25))
#abline(v = min(-1*PREDICT.ME.LINK[which(PREDICT.ME.CLASS$class=="weak")]), col="#A97C50",lwd=3)
#abline(v = min(-1*PREDICT.ME.LINK[which(PREDICT.ME.CLASS$class=="strong")]),  col="#F15A29",lwd=3)
#abline(h = quantile(DT.JOINT$pooled.meanF[which(DT.JOINT$pooled.class=="weak")],0.05), col="#A97C50",lwd=3)
#abline(h = quantile(DT.JOINT$pooled.meanF[which(DT.JOINT$pooled.class=="strong")],0.05), col="#F15A29",lwd=3)
#dev.off()

##PE
#png("Misclassification.PE.png",res=500,width=4,height=4,units="in")
#par(pty="s")
#plot(-1*PREDICT.PE.LINK,DT.JOINT$pooled.meanF,pch=19,cex=0.05,col="#000000",xlab="Genetic Score",ylab="Estimated Fluorescence",xlim=c(-20,25))
#abline(v = min(-1*PREDICT.PE.LINK[which(PREDICT.PE.CLASS$class=="weak")]),  col="#A97C50",lwd=3)
#abline(v = min(-1*PREDICT.PE.LINK[which(PREDICT.PE.CLASS$class=="strong")]), col="#F15A29",lwd=3)
#abline(h = quantile(DT.JOINT$pooled.meanF[which(DT.JOINT$pooled.class=="weak")],0.05),  col="#A97C50",lwd=3)
#abline(h = quantile(DT.JOINT$pooled.meanF[which(DT.JOINT$pooled.class=="strong")],0.05),  col="#F15A29",lwd=3)
#dev.off()

##TE
#png("Misclassification.TE.png",res=500,width=4,height=4,units="in")
#par(pty="s")
#plot(-1*PREDICT.TE.LINK,DT.JOINT$pooled.meanF,pch=19,cex=0.05,col="#000000",xlab="Genetic Score",ylab="Estimated Fluorescence",xlim=c(-10,20))
#abline(v = min(-1*PREDICT.TE.LINK[which(PREDICT.TE.CLASS$class=="weak")]), col="#A97C50",lwd=3)
#abline(v = min(-1*PREDICT.TE.LINK[which(PREDICT.TE.CLASS$class=="strong")]),  col="#F15A29",lwd=3)
#abline(h = quantile(DT.JOINT$pooled.meanF[which(DT.JOINT$pooled.class=="weak")],0.05), col="#A97C50",lwd=3)
#abline(h = quantile(DT.JOINT$pooled.meanF[which(DT.JOINT$pooled.class=="strong")],0.05),  col="#F15A29",lwd=3)
#dev.off()

###############################
##10. Misclassification rates##
###############################
####ME model
##Use null, weak, strong classes
##Columns: Observed
##Rows: Predicted
#M1.ME <- matrix(0,nrow=3,ncol=3)
#M1.ME[1,3] <- sum(PREDICT.ME.CLASS$class == "strong" & DT.JOINT$pooled.class == "strong", na.rm=TRUE)
#M1.ME[2,3] <- sum(PREDICT.ME.CLASS$class == "weak"   & DT.JOINT$pooled.class == "strong", na.rm=TRUE)
#M1.ME[3,3] <- sum(PREDICT.ME.CLASS$class == "null"   & DT.JOINT$pooled.class == "strong", na.rm=TRUE)
#M1.ME[1,2] <- sum(PREDICT.ME.CLASS$class == "strong" & DT.JOINT$pooled.class == "weak",   na.rm=TRUE)
#M1.ME[2,2] <- sum(PREDICT.ME.CLASS$class == "weak"   & DT.JOINT$pooled.class == "weak",   na.rm=TRUE)
#M1.ME[3,2] <- sum(PREDICT.ME.CLASS$class == "null"   & DT.JOINT$pooled.class == "weak",   na.rm=TRUE)
#M1.ME[1,1] <- sum(PREDICT.ME.CLASS$class == "strong" & DT.JOINT$pooled.class == "null",   na.rm=TRUE)
#M1.ME[2,1] <- sum(PREDICT.ME.CLASS$class == "weak"   & DT.JOINT$pooled.class == "null",   na.rm=TRUE)
#M1.ME[3,1] <- sum(PREDICT.ME.CLASS$class == "null"   & DT.JOINT$pooled.class == "null",   na.rm=TRUE)
#colnames(M1.ME) <- c("O.Null","O.Weak","O.Strong")
#rownames(M1.ME) <- c("P.Strong","P.Weak","P.Null")

##True/False Positives/Negatives
##Double counting null/strong errors to account for trinary classification
#M1.ME.TP <- 2*M1.ME[1,3] +   M1.ME[2,2] 
#M1.ME.TN <-   M1.ME[2,2] + 2*M1.ME[3,1]
#M1.ME.FP <- 2*M1.ME[1,1] +   M1.ME[1,2] +   M1.ME[2,1]
#M1.ME.FN <-   M1.ME[2,3] +   M1.ME[3,2] + 2*M1.ME[3,3]
#M1.ME.CM <- matrix(c(M1.ME.TP,M1.ME.FP,M1.ME.FN,M1.ME.TN),nrow=2,byrow=TRUE);colnames(M1.ME.CM) <- c("O.TRUE","O.FALSE");rownames(M1.ME.CM) <- c("P.TRUE","P.FALSE")

##Sensitivity/Specificity and Precision/Recall for each category
##Null
#M1.ME.NULL.TP   <- M1.ME[3,1]
#M1.ME.NULL.TN   <- M1.ME[1,2] +  M1.ME[1,3] + M1.ME[2,2] +  M1.ME[2,3] 
#M1.ME.NULL.FP   <- M1.ME[3,2] +  M1.ME[3,3]
#M1.ME.NULL.FN   <- M1.ME[1,1] +  M1.ME[2,1] 
#M1.ME.NULL.SENS <- M1.ME.NULL.TP/(M1.ME.NULL.TP + M1.ME.NULL.FN)
#M1.ME.NULL.SPEC <- M1.ME.NULL.TN/(M1.ME.NULL.TN + M1.ME.NULL.FP) 
#M1.ME.NULL.PR   <- M1.ME.NULL.TP/(M1.ME.NULL.TP + M1.ME.NULL.FP) 

##Weak
#M1.ME.WEAK.TP   <- M1.ME[2,2]
#M1.ME.WEAK.TN   <- M1.ME[1,1] +  M1.ME[1,3] + M1.ME[3,1] +  M1.ME[3,3] 
#M1.ME.WEAK.FP   <- M1.ME[2,1] +  M1.ME[2,3]
#M1.ME.WEAK.FN   <- M1.ME[1,2] +  M1.ME[3,2] 
#M1.ME.WEAK.SENS <- M1.ME.WEAK.TP/(M1.ME.WEAK.TP + M1.ME.WEAK.FN)
#M1.ME.WEAK.SPEC <- M1.ME.WEAK.TN/(M1.ME.WEAK.TN + M1.ME.WEAK.FP) 
#M1.ME.WEAK.PR   <- M1.ME.WEAK.TP/(M1.ME.WEAK.TP + M1.ME.WEAK.FP) 

##Strong
#M1.ME.STRONG.TP   <- M1.ME[1,3]
#M1.ME.STRONG.TN   <- M1.ME[2,1] +  M1.ME[2,2] + M1.ME[3,1] +  M1.ME[3,2] 
#M1.ME.STRONG.FP   <- M1.ME[1,1] +  M1.ME[1,2]
#M1.ME.STRONG.FN   <- M1.ME[2,3] +  M1.ME[3,3] 
#M1.ME.STRONG.SENS <- M1.ME.STRONG.TP/(M1.ME.STRONG.TP + M1.ME.STRONG.FN)
#M1.ME.STRONG.SPEC <- M1.ME.STRONG.TN/(M1.ME.STRONG.TN + M1.ME.STRONG.FP) 
#M1.ME.STRONG.PR   <- M1.ME.STRONG.TP/(M1.ME.STRONG.TP + M1.ME.STRONG.FP) 

##Activator (weak + strong)
#M1.ME.ACTIVATE.TP   <- M1.ME[1,2] +  M1.ME[1,3] + M1.ME[2,2] +  M1.ME[2,3] 
#M1.ME.ACTIVATE.TN   <- M1.ME[3,1] 
#M1.ME.ACTIVATE.FP   <- M1.ME[1,1] +  M1.ME[2,1]
#M1.ME.ACTIVATE.FN   <- M1.ME[3,2] +  M1.ME[3,3] 
#M1.ME.ACTIVATE.SENS <- M1.ME.ACTIVATE.TP/(M1.ME.ACTIVATE.TP + M1.ME.ACTIVATE.FN)
#M1.ME.ACTIVATE.SPEC <- M1.ME.ACTIVATE.TN/(M1.ME.ACTIVATE.TN + M1.ME.ACTIVATE.FP) 
#M1.ME.ACTIVATE.PR   <- M1.ME.ACTIVATE.TP/(M1.ME.ACTIVATE.TP + M1.ME.ACTIVATE.FP) 

##Not-Strong (null + weak)
#M1.ME.NSTRONG.TP   <- M1.ME[2,1] +  M1.ME[2,2] + M1.ME[3,1] +  M1.ME[3,2] 
#M1.ME.NSTRONG.TN   <- M1.ME[1,3] 
#M1.ME.NSTRONG.FP   <- M1.ME[2,3] +  M1.ME[3,3]
#M1.ME.NSTRONG.FN   <- M1.ME[1,1] +  M1.ME[1,2] 
#M1.ME.NSTRONG.SENS <- M1.ME.NSTRONG.TP/(M1.ME.NSTRONG.TP + M1.ME.NSTRONG.FN)
#M1.ME.NSTRONG.SPEC <- M1.ME.NSTRONG.TN/(M1.ME.NSTRONG.TN + M1.ME.NSTRONG.FP) 
#M1.ME.NSTRONG.PR   <- M1.ME.NSTRONG.TP/(M1.ME.NSTRONG.TP + M1.ME.NSTRONG.FP) 
#M1.ME.NSTRONG.RE   <- M1.ME.NSTRONG.TP/(M1.ME.NSTRONG.TP + M1.ME.NSTRONG.FN) 

##Rank data
#RANK.M1.ME <- data.frame(matrix(0,nrow=nrow(PREDICT.ME.CLASS$probs),ncol=10));colnames(RANK.M1.ME) = c("P.NULL","P.WEAK","P.STRONG","P.ACTIVATE","P.NSTRONG","O.NULL","O.WEAK","O.STRONG","O.ACTIVATE","O.NSTRONG")
#RANK.M1.ME$P.NULL     <- PREDICT.ME.CLASS$probs[,1,] 
#RANK.M1.ME$P.WEAK     <- PREDICT.ME.CLASS$probs[,2,] 
#RANK.M1.ME$P.STRONG   <- PREDICT.ME.CLASS$probs[,3,] 
#RANK.M1.ME$P.ACTIVATE <- PREDICT.ME.CLASS$probs[,2,] + PREDICT.ME.CLASS$probs[,3,] 
#RANK.M1.ME$P.NSTRONG  <- PREDICT.ME.CLASS$probs[,1,] + PREDICT.ME.CLASS$probs[,2,] 
#RANK.M1.ME[which(DT.JOINT$pooled.class %in% c("null")),6] <- 1
#RANK.M1.ME[which(DT.JOINT$pooled.class %in% c("weak")),7] <- 1
#RANK.M1.ME[which(DT.JOINT$pooled.class %in% c("strong")),8] <- 1
#RANK.M1.ME[which(DT.JOINT$pooled.class %in% c("weak","strong")),9] <- 1
#RANK.M1.ME[which(DT.JOINT$pooled.class %in% c("null","weak")),10] <- 1
#RANK.M1.ME[is.na(DT.JOINT$pooled.class),6:10] <- NA
#RANK.M1.ME <- RANK.M1.ME[order(RANK.M1.ME$P.STRONG),]
#RANK.M1.ME <- na.omit(RANK.M1.ME)

##ROC calculation
#ROC.M1.ME.NULL     <- roc(RANK.M1.ME$O.NULL,RANK.M1.ME$P.NULL,        ci=TRUE,plot=FALSE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,print.auc=TRUE,main="M1.ME",xlab="FP/(TN + FP)",ylab="TP/(TP + FN)")
#ROC.M1.ME.WEAK     <- roc(RANK.M1.ME$O.WEAK,RANK.M1.ME$P.WEAK,        ci=TRUE,plot=FALSE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,print.auc=TRUE,main="M1.ME",xlab="FP/(TN + FP)",ylab="TP/(TP + FN)")
#ROC.M1.ME.STRONG   <- roc(RANK.M1.ME$O.STRONG,RANK.M1.ME$P.STRONG,    ci=TRUE,plot=FALSE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,print.auc=TRUE,main="M1.ME",xlab="FP/(TN + FP)",ylab="TP/(TP + FN)")
#ROC.M1.ME.ACTIVATE <- roc(RANK.M1.ME$O.ACTIVATE,RANK.M1.ME$P.ACTIVATE,ci=TRUE,plot=FALSE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,print.auc=TRUE,main="M1.ME",xlab="FP/(TN + FP)",ylab="TP/(TP + FN)")
#ROC.M1.ME.NSTRONG  <- roc(RANK.M1.ME$O.NSTRONG,RANK.M1.ME$P.NSTRONG,  ci=TRUE,plot=FALSE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,print.auc=TRUE,main="M1.ME",xlab="FP/(TN + FP)",ylab="TP/(TP + FN)")

##Extract metrics
#METRICS.M1.ME.NULL     <- coords(ROC.M1.ME.NULL,     ret = "all", transpose = FALSE)
#METRICS.M1.ME.WEAK     <- coords(ROC.M1.ME.WEAK,     ret = "all", transpose = FALSE)
#METRICS.M1.ME.STRONG   <- coords(ROC.M1.ME.STRONG,   ret = "all", transpose = FALSE)
#METRICS.M1.ME.ACTIVATE <- coords(ROC.M1.ME.ACTIVATE, ret = "all", transpose = FALSE)
#METRICS.M1.ME.NSTRONG  <- coords(ROC.M1.ME.NSTRONG,  ret = "all", transpose = FALSE)

#pdf("Accuracy.ME.pdf")
#par(mfrow=c(2,2))
#par(pty="s")
##ROC
#plot(  1-rev(ROC.M1.ME.NULL$specificities),      rev(ROC.M1.ME.NULL$sensitivities),      type="l",col="red", xlab="FP/(TN + FP)",ylab="TP/(TP + FN)",main="M1.ME")
#points(1-rev(ROC.M1.ME.WEAK$specificities),      rev(ROC.M1.ME.WEAK$sensitivities),      type="l",col="blue")
#points(1-rev(ROC.M1.ME.STRONG$specificities),    rev(ROC.M1.ME.STRONG$sensitivities),    type="l",col="green")
#points(1-rev(ROC.M1.ME.ACTIVATE$specificities),  rev(ROC.M1.ME.ACTIVATE$sensitivities),  type="l",col="black")
#points(1-rev(ROC.M1.ME.NSTRONG$specificities),   rev(ROC.M1.ME.NSTRONG$sensitivities),   type="l",col="orange")
#points(1-M1.ME.NULL.SPEC,      M1.ME.NULL.SENS,      pch=19,cex=1,col="red")
#points(1-M1.ME.WEAK.SPEC,      M1.ME.WEAK.SENS,      pch=19,cex=1,col="blue")
#points(1-M1.ME.STRONG.SPEC,    M1.ME.STRONG.SENS,    pch=19,cex=1,col="green")
#points(1-M1.ME.ACTIVATE.SPEC,  M1.ME.ACTIVATE.SENS,  pch=19,cex=1,col="black")
#points(1-M1.ME.NSTRONG.SPEC,   M1.ME.NSTRONG.SENS,   pch=19,cex=1,col="orange")
#abline(a=0,b=1)
#legend(.7,.6,c("Null","Weak","Strong","Activator","Non-Strong"),fill=c("red","blue","green","black","orange"),cex=0.4)

##Precision/Recall
#plot(  METRICS.M1.ME.NULL$recall,     METRICS.M1.ME.NULL$precision,      type="l",col="red",xlab="TP/(TP+FN)",ylab="TP/(TP+FP)",main="M1.ME",ylim=c(0,1))
#points(METRICS.M1.ME.WEAK$recall,     METRICS.M1.ME.WEAK$precision,      type="l",col="blue")
#points(METRICS.M1.ME.STRONG$recall,   METRICS.M1.ME.STRONG$precision,    type="l",col="green")
#points(METRICS.M1.ME.ACTIVATE$recall, METRICS.M1.ME.ACTIVATE$precision,  type="l",col="black")
#points(METRICS.M1.ME.NSTRONG$recall,  METRICS.M1.ME.NSTRONG$precision,   type="l",col="orange")
#points(M1.ME.NULL.SENS,    M1.ME.NULL.PR,    pch=19,cex=0.6)
#points(M1.ME.WEAK.SENS,    M1.ME.WEAK.PR,    pch=19,cex=0.6)
#points(M1.ME.STRONG.SENS,  M1.ME.STRONG.PR,  pch=19,cex=0.6)
#points(M1.ME.ACTIVATE.SENS,M1.ME.ACTIVATE.PR,pch=19,cex=0.6)
#points(M1.ME.NSTRONG.SENS, M1.ME.NSTRONG.PR, pch=19,cex=0.6)
#legend(.3,.3,c("Null","Weak","Strong","Activator","Non-Strong"),fill=c("red","blue","green","black","orange"),cex=0.4)

##Lift
#AVG.RATE.M1.ME.NULL     <- sum(RANK.M1.ME$O.NULL    )/nrow(RANK.M1.ME)
#AVG.RATE.M1.ME.WEAK     <- sum(RANK.M1.ME$O.WEAK    )/nrow(RANK.M1.ME)
#AVG.RATE.M1.ME.STRONG   <- sum(RANK.M1.ME$O.STRONG  )/nrow(RANK.M1.ME)
#AVG.RATE.M1.ME.ACTIVATE <- sum(RANK.M1.ME$O.ACTIVATE)/nrow(RANK.M1.ME)
#AVG.RATE.M1.ME.NSTRONG  <- sum(RANK.M1.ME$O.NSTRONG )/nrow(RANK.M1.ME)

#I <- 1:100
#PER.RATE.M1.ME.NULL     <- numeric(length(I))
#PER.RATE.M1.ME.WEAK     <- numeric(length(I))
#PER.RATE.M1.ME.STRONG   <- numeric(length(I))
#PER.RATE.M1.ME.ACTIVATE <- numeric(length(I))
#PER.RATE.M1.ME.NSTRONG  <- numeric(length(I))
#for(i in I) {
#  LB <- floor(quantile(1:length(RANK.M1.ME$O.NULL),0.01*(i-1)))
#  UB <- floor(quantile(1:length(RANK.M1.ME$O.NULL),0.01*i))
#  PER.RATE.M1.ME.NULL[i] <- sum(RANK.M1.ME$O.NULL[LB:UB])/(UB-LB)
#  
#  LB <- floor(quantile(1:length(RANK.M1.ME$O.WEAK),0.01*(i-1)))
#  UB <- floor(quantile(1:length(RANK.M1.ME$O.WEAK),0.01*i))
#  PER.RATE.M1.ME.WEAK[i] <- sum(RANK.M1.ME$O.WEAK[LB:UB])/(UB-LB)
#  
#  LB <- floor(quantile(1:length(RANK.M1.ME$O.STRONG),0.01*(i-1)))
#  UB <- floor(quantile(1:length(RANK.M1.ME$O.STRONG),0.01*i))
#  PER.RATE.M1.ME.STRONG[i] <- sum(RANK.M1.ME$O.STRONG[LB:UB])/(UB-LB)
#  
#  LB <- floor(quantile(1:length(RANK.M1.ME$O.ACTIVATE),0.01*(i-1)))
#  UB <- floor(quantile(1:length(RANK.M1.ME$O.ACTIVATE),0.01*i))
#  PER.RATE.M1.ME.ACTIVATE[i] <- sum(RANK.M1.ME$O.ACTIVATE[LB:UB])/(UB-LB)
#  
#  LB <- floor(quantile(1:length(RANK.M1.ME$O.NSTRONG),0.01*(i-1)))
#  UB <- floor(quantile(1:length(RANK.M1.ME$O.NSTRONG),0.01*i))
#  PER.RATE.M1.ME.NSTRONG[i] <- sum(RANK.M1.ME$O.NSTRONG[LB:UB])/(UB-LB)
#}
#plot(  I,PER.RATE.M1.ME.NULL     /AVG.RATE.M1.ME.NULL     ,type="l",col="red",ylab="Lift",xlab="Percentile",main="M1.ME",ylim=c(0,80))
#points(I,PER.RATE.M1.ME.WEAK     /AVG.RATE.M1.ME.WEAK     ,type="l",col="blue")
#points(I,PER.RATE.M1.ME.STRONG   /AVG.RATE.M1.ME.STRONG   ,type="l",col="green")
#points(I,PER.RATE.M1.ME.ACTIVATE /AVG.RATE.M1.ME.ACTIVATE ,type="l",col="black")
#points(I,PER.RATE.M1.ME.NSTRONG  /AVG.RATE.M1.ME.NSTRONG  ,type="l",col="orange")
#legend(20,60,c("Null","Weak","Strong","Activator","Non-Strong"),fill=c("red","blue","green","black","orange"),cex=0.4)
#dev.off()

####PE model
##Use null, weak, strong classes
##Columns: Observed
##Rows: Predicted
#M1.PE <- matrix(0,nrow=3,ncol=3)
#M1.PE[1,3] <- sum(PREDICT.PE.CLASS$class == "strong" & DT.JOINT$pooled.class == "strong", na.rm=TRUE)
#M1.PE[2,3] <- sum(PREDICT.PE.CLASS$class == "weak"   & DT.JOINT$pooled.class == "strong", na.rm=TRUE)
#M1.PE[3,3] <- sum(PREDICT.PE.CLASS$class == "null"   & DT.JOINT$pooled.class == "strong", na.rm=TRUE)
#M1.PE[1,2] <- sum(PREDICT.PE.CLASS$class == "strong" & DT.JOINT$pooled.class == "weak",   na.rm=TRUE)
#M1.PE[2,2] <- sum(PREDICT.PE.CLASS$class == "weak"   & DT.JOINT$pooled.class == "weak",   na.rm=TRUE)
#M1.PE[3,2] <- sum(PREDICT.PE.CLASS$class == "null"   & DT.JOINT$pooled.class == "weak",   na.rm=TRUE)
#M1.PE[1,1] <- sum(PREDICT.PE.CLASS$class == "strong" & DT.JOINT$pooled.class == "null",   na.rm=TRUE)
#M1.PE[2,1] <- sum(PREDICT.PE.CLASS$class == "weak"   & DT.JOINT$pooled.class == "null",   na.rm=TRUE)
#M1.PE[3,1] <- sum(PREDICT.PE.CLASS$class == "null"   & DT.JOINT$pooled.class == "null",   na.rm=TRUE)
#colnames(M1.PE) <- c("O.Null","O.Weak","O.Strong")
#rownames(M1.PE) <- c("P.Strong","P.Weak","P.Null")

##True/False Positives/Negatives
##Double counting to account for trinary classification
#M1.PE.TP <- 2*M1.PE[1,3] +   M1.PE[2,2] 
#M1.PE.TN <-   M1.PE[2,2] + 2*M1.PE[3,1]
#M1.PE.FP <- 2*M1.PE[1,1] +   M1.PE[1,2] +   M1.PE[2,1]
#M1.PE.FN <-   M1.PE[2,3] +   M1.PE[3,2] + 2*M1.PE[3,3]
#M1.PE.CM <- matrix(c(M1.PE.TP,M1.PE.FP,M1.PE.FN,M1.PE.TN),nrow=2,byrow=TRUE);colnames(M1.PE.CM) <- c("O.TRUE","O.FALSE");rownames(M1.PE.CM) <- c("P.TRUE","P.FALSE")

##Sensitivity/Specificity and Precision/Recall for default cutoff
##Null
#M1.PE.NULL.TP   <- M1.PE[3,1]
#M1.PE.NULL.TN   <- M1.PE[1,2] +  M1.PE[1,3] + M1.PE[2,2] +  M1.PE[2,3] 
#M1.PE.NULL.FP   <- M1.PE[3,2] +  M1.PE[3,3]
#M1.PE.NULL.FN   <- M1.PE[1,1] +  M1.PE[2,1] 
#M1.PE.NULL.SENS <- M1.PE.NULL.TP/(M1.PE.NULL.TP + M1.PE.NULL.FN)
#M1.PE.NULL.SPEC <- M1.PE.NULL.TN/(M1.PE.NULL.TN + M1.PE.NULL.FP) 
#M1.PE.NULL.PR   <- M1.PE.NULL.TP/(M1.PE.NULL.TP + M1.PE.NULL.FP) 

##Weak
#M1.PE.WEAK.TP   <- M1.PE[2,2]
#M1.PE.WEAK.TN   <- M1.PE[1,1] +  M1.PE[1,3] + M1.PE[3,1] +  M1.PE[3,3] 
#M1.PE.WEAK.FP   <- M1.PE[2,1] +  M1.PE[2,3]
#M1.PE.WEAK.FN   <- M1.PE[1,2] +  M1.PE[3,2] 
#M1.PE.WEAK.SENS <- M1.PE.WEAK.TP/(M1.PE.WEAK.TP + M1.PE.WEAK.FN)
#M1.PE.WEAK.SPEC <- M1.PE.WEAK.TN/(M1.PE.WEAK.TN + M1.PE.WEAK.FP) 
#M1.PE.WEAK.PR   <- M1.PE.WEAK.TP/(M1.PE.WEAK.TP + M1.PE.WEAK.FP) 

##Strong
#M1.PE.STRONG.TP   <- M1.PE[1,3]
#M1.PE.STRONG.TN   <- M1.PE[2,1] +  M1.PE[2,2] + M1.PE[3,1] +  M1.PE[3,2] 
#M1.PE.STRONG.FP   <- M1.PE[1,1] +  M1.PE[1,2]
#M1.PE.STRONG.FN   <- M1.PE[2,3] +  M1.PE[3,3] 
#M1.PE.STRONG.SENS <- M1.PE.STRONG.TP/(M1.PE.STRONG.TP + M1.PE.STRONG.FN)
#M1.PE.STRONG.SPEC <- M1.PE.STRONG.TN/(M1.PE.STRONG.TN + M1.PE.STRONG.FP) 
#M1.PE.STRONG.PR   <- M1.PE.STRONG.TP/(M1.PE.STRONG.TP + M1.PE.STRONG.FP) 

##Activator (weak + strong)
#M1.PE.ACTIVATE.TP   <- M1.PE[1,2] +  M1.PE[1,3] + M1.PE[2,2] +  M1.PE[2,3] 
#M1.PE.ACTIVATE.TN   <- M1.PE[3,1] 
#M1.PE.ACTIVATE.FP   <- M1.PE[1,1] +  M1.PE[2,1]
#M1.PE.ACTIVATE.FN   <- M1.PE[3,2] +  M1.PE[3,3] 
#M1.PE.ACTIVATE.SENS <- M1.PE.ACTIVATE.TP/(M1.PE.ACTIVATE.TP + M1.PE.ACTIVATE.FN)
#M1.PE.ACTIVATE.SPEC <- M1.PE.ACTIVATE.TN/(M1.PE.ACTIVATE.TN + M1.PE.ACTIVATE.FP) 
#M1.PE.ACTIVATE.PR   <- M1.PE.ACTIVATE.TP/(M1.PE.ACTIVATE.TP + M1.PE.ACTIVATE.FP) 

##Not-Strong (null + weak)
#M1.PE.NSTRONG.TP   <- M1.PE[2,1] +  M1.PE[2,2] + M1.PE[3,1] +  M1.PE[3,2] 
#M1.PE.NSTRONG.TN   <- M1.PE[1,3] 
#M1.PE.NSTRONG.FP   <- M1.PE[2,3] +  M1.PE[3,3]
#M1.PE.NSTRONG.FN   <- M1.PE[1,1] +  M1.PE[1,2] 
#M1.PE.NSTRONG.SENS <- M1.PE.NSTRONG.TP/(M1.PE.NSTRONG.TP + M1.PE.NSTRONG.FN)
#M1.PE.NSTRONG.SPEC <- M1.PE.NSTRONG.TN/(M1.PE.NSTRONG.TN + M1.PE.NSTRONG.FP) 
#M1.PE.NSTRONG.PR   <- M1.PE.NSTRONG.TP/(M1.PE.NSTRONG.TP + M1.PE.NSTRONG.FP) 

##Rank data
#RANK.M1.PE <- data.frame(matrix(0,nrow=nrow(PREDICT.PE.CLASS$probs),ncol=10));colnames(RANK.M1.PE) = c("P.NULL","P.WEAK","P.STRONG","P.ACTIVATE","P.NSTRONG","O.NULL","O.WEAK","O.STRONG","O.ACTIVATE","O.NSTRONG")
#RANK.M1.PE$P.NULL     <- PREDICT.PE.CLASS$probs[,1,] 
#RANK.M1.PE$P.WEAK     <- PREDICT.PE.CLASS$probs[,2,] 
#RANK.M1.PE$P.STRONG   <- PREDICT.PE.CLASS$probs[,3,] 
#RANK.M1.PE$P.ACTIVATE <- PREDICT.PE.CLASS$probs[,2,] + PREDICT.PE.CLASS$probs[,3,] 
#RANK.M1.PE$P.NSTRONG  <- PREDICT.PE.CLASS$probs[,1,] + PREDICT.PE.CLASS$probs[,2,] 
#RANK.M1.PE[which(DT.JOINT$pooled.class %in% c("null")),6] <- 1
#RANK.M1.PE[which(DT.JOINT$pooled.class %in% c("weak")),7] <- 1
#RANK.M1.PE[which(DT.JOINT$pooled.class %in% c("strong")),8] <- 1
#RANK.M1.PE[which(DT.JOINT$pooled.class %in% c("weak","strong")),9] <- 1
#RANK.M1.PE[which(DT.JOINT$pooled.class %in% c("null","weak")),10] <- 1
#RANK.M1.PE[is.na(DT.JOINT$pooled.class),6:10] <- NA
#RANK.M1.PE <- RANK.M1.PE[order(RANK.M1.PE$P.STRONG),]
#RANK.M1.PE <- na.omit(RANK.M1.PE)

##ROC calculation
#ROC.M1.PE.NULL     <- roc(RANK.M1.PE$O.NULL,RANK.M1.PE$P.NULL,        ci=TRUE,plot=FALSE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,print.auc=TRUE,main="M1.PE",xlab="FP/(TN + FP)",ylab="TP/(TP + FN)")
#ROC.M1.PE.WEAK     <- roc(RANK.M1.PE$O.WEAK,RANK.M1.PE$P.WEAK,        ci=TRUE,plot=FALSE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,print.auc=TRUE,main="M1.PE",xlab="FP/(TN + FP)",ylab="TP/(TP + FN)")
#ROC.M1.PE.STRONG   <- roc(RANK.M1.PE$O.STRONG,RANK.M1.PE$P.STRONG,    ci=TRUE,plot=FALSE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,print.auc=TRUE,main="M1.PE",xlab="FP/(TN + FP)",ylab="TP/(TP + FN)")
#ROC.M1.PE.ACTIVATE <- roc(RANK.M1.PE$O.ACTIVATE,RANK.M1.PE$P.ACTIVATE,ci=TRUE,plot=FALSE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,print.auc=TRUE,main="M1.PE",xlab="FP/(TN + FP)",ylab="TP/(TP + FN)")
#ROC.M1.PE.NSTRONG  <- roc(RANK.M1.PE$O.NSTRONG,RANK.M1.PE$P.NSTRONG,  ci=TRUE,plot=FALSE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,print.auc=TRUE,main="M1.PE",xlab="FP/(TN + FP)",ylab="TP/(TP + FN)")

##Extract metrics
#METRICS.M1.PE.NULL     <- coords(ROC.M1.PE.NULL,     ret = "all", transpose = FALSE)
#METRICS.M1.PE.WEAK     <- coords(ROC.M1.PE.WEAK,     ret = "all", transpose = FALSE)
#METRICS.M1.PE.STRONG   <- coords(ROC.M1.PE.STRONG,   ret = "all", transpose = FALSE)
#METRICS.M1.PE.ACTIVATE <- coords(ROC.M1.PE.ACTIVATE, ret = "all", transpose = FALSE)
#METRICS.M1.PE.NSTRONG  <- coords(ROC.M1.PE.NSTRONG,  ret = "all", transpose = FALSE)

#pdf("Accuracy.PE.pdf")
#par(mfrow=c(2,2))
#par(pty="s")
##ROC
#plot(  1-rev(ROC.M1.PE.NULL$specificities),      rev(ROC.M1.PE.NULL$sensitivities),      type="l",col="red", xlab="FP/(TN + FP)",ylab="TP/(TP + FN)",main="M1.PE")
#points(1-rev(ROC.M1.PE.WEAK$specificities),      rev(ROC.M1.PE.WEAK$sensitivities),      type="l",col="blue")
#points(1-rev(ROC.M1.PE.STRONG$specificities),    rev(ROC.M1.PE.STRONG$sensitivities),    type="l",col="green")
#points(1-rev(ROC.M1.PE.ACTIVATE$specificities),  rev(ROC.M1.PE.ACTIVATE$sensitivities),  type="l",col="black")
#points(1-rev(ROC.M1.PE.NSTRONG$specificities),   rev(ROC.M1.PE.NSTRONG$sensitivities),   type="l",col="orange")
#points(1-M1.PE.NULL.SPEC,      M1.PE.NULL.SENS,      pch=19,cex=1,col="red")
#points(1-M1.PE.WEAK.SPEC,      M1.PE.WEAK.SENS,      pch=19,cex=1,col="blue")
#points(1-M1.PE.STRONG.SPEC,    M1.PE.STRONG.SENS,    pch=19,cex=1,col="green")
#points(1-M1.PE.ACTIVATE.SPEC,  M1.PE.ACTIVATE.SENS,  pch=19,cex=1,col="black")
#points(1-M1.PE.NSTRONG.SPEC,   M1.PE.NSTRONG.SENS,   pch=19,cex=1,col="orange")
#abline(a=0,b=1)
#legend(.7,.6,c("Null","Weak","Strong","Activator","Non-Strong"),fill=c("red","blue","green","black","orange"),cex=0.4)

##Precision/Recall
#plot(  METRICS.M1.PE.NULL$recall,     METRICS.M1.PE.NULL$precision,      type="l",col="red",xlab="TP/(TP+FN)",ylab="TP/(TP+FP)",main="M1.PE",ylim=c(0,1))
#points(METRICS.M1.PE.WEAK$recall,     METRICS.M1.PE.WEAK$precision,      type="l",col="blue")
#points(METRICS.M1.PE.STRONG$recall,   METRICS.M1.PE.STRONG$precision,    type="l",col="green")
#points(METRICS.M1.PE.ACTIVATE$recall, METRICS.M1.PE.ACTIVATE$precision,  type="l",col="black")
#points(METRICS.M1.PE.NSTRONG$recall,  METRICS.M1.PE.NSTRONG$precision,   type="l",col="orange")
#points(M1.PE.NULL.SENS,    M1.PE.NULL.PR,    pch=19,cex=0.6)
#points(M1.PE.WEAK.SENS,    M1.PE.WEAK.PR,    pch=19,cex=0.6)
#points(M1.PE.STRONG.SENS,  M1.PE.STRONG.PR,  pch=19,cex=0.6)
#points(M1.PE.ACTIVATE.SENS,M1.PE.ACTIVATE.PR,pch=19,cex=0.6)
#points(M1.PE.NSTRONG.SENS, M1.PE.NSTRONG.PR, pch=19,cex=0.6)
#legend(.3,.3,c("Null","Weak","Strong","Activator","Non-Strong"),fill=c("red","blue","green","black","orange"),cex=0.4)

##Lift
#AVG.RATE.M1.PE.NULL     <- sum(RANK.M1.PE$O.NULL    )/nrow(RANK.M1.PE)
#AVG.RATE.M1.PE.WEAK     <- sum(RANK.M1.PE$O.WEAK    )/nrow(RANK.M1.PE)
#AVG.RATE.M1.PE.STRONG   <- sum(RANK.M1.PE$O.STRONG  )/nrow(RANK.M1.PE)
#AVG.RATE.M1.PE.ACTIVATE <- sum(RANK.M1.PE$O.ACTIVATE)/nrow(RANK.M1.PE)
#AVG.RATE.M1.PE.NSTRONG  <- sum(RANK.M1.PE$O.NSTRONG )/nrow(RANK.M1.PE)

#I <- 1:100
#PER.RATE.M1.PE.NULL     <- numeric(length(I))
#PER.RATE.M1.PE.WEAK     <- numeric(length(I))
#PER.RATE.M1.PE.STRONG   <- numeric(length(I))
#PER.RATE.M1.PE.ACTIVATE <- numeric(length(I))
#PER.RATE.M1.PE.NSTRONG  <- numeric(length(I))
#for(i in I) {
#  LB <- floor(quantile(1:length(RANK.M1.PE$O.NULL),0.01*(i-1)))
#  UB <- floor(quantile(1:length(RANK.M1.PE$O.NULL),0.01*i))
#  PER.RATE.M1.PE.NULL[i] <- sum(RANK.M1.PE$O.NULL[LB:UB])/(UB-LB)
#  
#  LB <- floor(quantile(1:length(RANK.M1.PE$O.WEAK),0.01*(i-1)))
#  UB <- floor(quantile(1:length(RANK.M1.PE$O.WEAK),0.01*i))
#  PER.RATE.M1.PE.WEAK[i] <- sum(RANK.M1.PE$O.WEAK[LB:UB])/(UB-LB)
#  
#  LB <- floor(quantile(1:length(RANK.M1.PE$O.STRONG),0.01*(i-1)))
#  UB <- floor(quantile(1:length(RANK.M1.PE$O.STRONG),0.01*i))
#  PER.RATE.M1.PE.STRONG[i] <- sum(RANK.M1.PE$O.STRONG[LB:UB])/(UB-LB)
#  
#  LB <- floor(quantile(1:length(RANK.M1.PE$O.ACTIVATE),0.01*(i-1)))
#  UB <- floor(quantile(1:length(RANK.M1.PE$O.ACTIVATE),0.01*i))
#  PER.RATE.M1.PE.ACTIVATE[i] <- sum(RANK.M1.PE$O.ACTIVATE[LB:UB])/(UB-LB)
#  
#  LB <- floor(quantile(1:length(RANK.M1.PE$O.NSTRONG),0.01*(i-1)))
#  UB <- floor(quantile(1:length(RANK.M1.PE$O.NSTRONG),0.01*i))
#  PER.RATE.M1.PE.NSTRONG[i] <- sum(RANK.M1.PE$O.NSTRONG[LB:UB])/(UB-LB)
#}
#plot(  I,PER.RATE.M1.PE.NULL     /AVG.RATE.M1.PE.NULL     ,type="l",col="red",ylab="Lift",xlab="Percentile",main="M1.PE",ylim=c(0,100))
#points(I,PER.RATE.M1.PE.WEAK     /AVG.RATE.M1.PE.WEAK     ,type="l",col="blue")
#points(I,PER.RATE.M1.PE.STRONG   /AVG.RATE.M1.PE.STRONG   ,type="l",col="green")
#points(I,PER.RATE.M1.PE.ACTIVATE /AVG.RATE.M1.PE.ACTIVATE ,type="l",col="black")
#points(I,PER.RATE.M1.PE.NSTRONG  /AVG.RATE.M1.PE.NSTRONG  ,type="l",col="orange")
#legend(20,60,c("Null","Weak","Strong","Activator","Non-Strong"),fill=c("red","blue","green","black","orange"),cex=0.4)
#dev.off()

####TE model
##Use null, weak, strong classes
##Columns: Observed
##Rows: Predicted
M1.TE <- matrix(0,nrow=3,ncol=3)
M1.TE[1,3] <- sum(PREDICT.TE.CLASS$class == "strong" & DT.JOINT$pooled.class == "strong", na.rm=TRUE)
M1.TE[2,3] <- sum(PREDICT.TE.CLASS$class == "weak"   & DT.JOINT$pooled.class == "strong", na.rm=TRUE)
M1.TE[3,3] <- sum(PREDICT.TE.CLASS$class == "null"   & DT.JOINT$pooled.class == "strong", na.rm=TRUE)
M1.TE[1,2] <- sum(PREDICT.TE.CLASS$class == "strong" & DT.JOINT$pooled.class == "weak",   na.rm=TRUE)
M1.TE[2,2] <- sum(PREDICT.TE.CLASS$class == "weak"   & DT.JOINT$pooled.class == "weak",   na.rm=TRUE)
M1.TE[3,2] <- sum(PREDICT.TE.CLASS$class == "null"   & DT.JOINT$pooled.class == "weak",   na.rm=TRUE)
M1.TE[1,1] <- sum(PREDICT.TE.CLASS$class == "strong" & DT.JOINT$pooled.class == "null",   na.rm=TRUE)
M1.TE[2,1] <- sum(PREDICT.TE.CLASS$class == "weak"   & DT.JOINT$pooled.class == "null",   na.rm=TRUE)
M1.TE[3,1] <- sum(PREDICT.TE.CLASS$class == "null"   & DT.JOINT$pooled.class == "null",   na.rm=TRUE)
colnames(M1.TE) <- c("O.Null","O.Weak","O.Strong")
rownames(M1.TE) <- c("P.Strong","P.Weak","P.Null")

##True/False Positives/Negatives
#Double counting to account for trinary classification
M1.TE.TP <- 2*M1.TE[1,3] +   M1.TE[2,2] 
M1.TE.TN <-   M1.TE[2,2] + 2*M1.TE[3,1]
M1.TE.FP <- 2*M1.TE[1,1] +   M1.TE[1,2] +   M1.TE[2,1]
M1.TE.FN <-   M1.TE[2,3] +   M1.TE[3,2] + 2*M1.TE[3,3]
M1.TE.CM <- matrix(c(M1.TE.TP,M1.TE.FP,M1.TE.FN,M1.TE.TN),nrow=2,byrow=TRUE);colnames(M1.TE.CM) <- c("O.TRUE","O.FALSE");rownames(M1.TE.CM) <- c("P.TRUE","P.FALSE")

##Sensitivity/Specificity and Precision/Recall for default cutoff
##Null
M1.TE.NULL.TP   <- M1.TE[3,1]
M1.TE.NULL.TN   <- M1.TE[1,2] +  M1.TE[1,3] + M1.TE[2,2] +  M1.TE[2,3] 
M1.TE.NULL.FP   <- M1.TE[3,2] +  M1.TE[3,3]
M1.TE.NULL.FN   <- M1.TE[1,1] +  M1.TE[2,1] 
M1.TE.NULL.SENS <- M1.TE.NULL.TP/(M1.TE.NULL.TP + M1.TE.NULL.FN)
M1.TE.NULL.SPEC <- M1.TE.NULL.TN/(M1.TE.NULL.TN + M1.TE.NULL.FP) 
M1.TE.NULL.PR   <- M1.TE.NULL.TP/(M1.TE.NULL.TP + M1.TE.NULL.FP) 

##Weak
M1.TE.WEAK.TP   <- M1.TE[2,2]
M1.TE.WEAK.TN   <- M1.TE[1,1] +  M1.TE[1,3] + M1.TE[3,1] +  M1.TE[3,3] 
M1.TE.WEAK.FP   <- M1.TE[2,1] +  M1.TE[2,3]
M1.TE.WEAK.FN   <- M1.TE[1,2] +  M1.TE[3,2] 
M1.TE.WEAK.SENS <- M1.TE.WEAK.TP/(M1.TE.WEAK.TP + M1.TE.WEAK.FN)
M1.TE.WEAK.SPEC <- M1.TE.WEAK.TN/(M1.TE.WEAK.TN + M1.TE.WEAK.FP) 
M1.TE.WEAK.PR   <- M1.TE.WEAK.TP/(M1.TE.WEAK.TP + M1.TE.WEAK.FP) 

##Strong
M1.TE.STRONG.TP   <- M1.TE[1,3]
M1.TE.STRONG.TN   <- M1.TE[2,1] +  M1.TE[2,2] + M1.TE[3,1] +  M1.TE[3,2] 
M1.TE.STRONG.FP   <- M1.TE[1,1] +  M1.TE[1,2]
M1.TE.STRONG.FN   <- M1.TE[2,3] +  M1.TE[3,3] 
M1.TE.STRONG.SENS <- M1.TE.STRONG.TP/(M1.TE.STRONG.TP + M1.TE.STRONG.FN)
M1.TE.STRONG.SPEC <- M1.TE.STRONG.TN/(M1.TE.STRONG.TN + M1.TE.STRONG.FP) 
M1.TE.STRONG.PR   <- M1.TE.STRONG.TP/(M1.TE.STRONG.TP + M1.TE.STRONG.FP) 

##Activator
M1.TE.ACTIVATE.TP   <- M1.TE[1,2] +  M1.TE[1,3] + M1.TE[2,2] +  M1.TE[2,3] 
M1.TE.ACTIVATE.TN   <- M1.TE[3,1] 
M1.TE.ACTIVATE.FP   <- M1.TE[1,1] +  M1.TE[2,1]
M1.TE.ACTIVATE.FN   <- M1.TE[3,2] +  M1.TE[3,3] 
M1.TE.ACTIVATE.SENS <- M1.TE.ACTIVATE.TP/(M1.TE.ACTIVATE.TP + M1.TE.ACTIVATE.FN)
M1.TE.ACTIVATE.SPEC <- M1.TE.ACTIVATE.TN/(M1.TE.ACTIVATE.TN + M1.TE.ACTIVATE.FP) 
M1.TE.ACTIVATE.PR   <- M1.TE.ACTIVATE.TP/(M1.TE.ACTIVATE.TP + M1.TE.ACTIVATE.FP) 

##Not-Strong
M1.TE.NSTRONG.TP   <- M1.TE[2,1] +  M1.TE[2,2] + M1.TE[3,1] +  M1.TE[3,2] 
M1.TE.NSTRONG.TN   <- M1.TE[1,3] 
M1.TE.NSTRONG.FP   <- M1.TE[2,3] +  M1.TE[3,3]
M1.TE.NSTRONG.FN   <- M1.TE[1,1] +  M1.TE[1,2] 
M1.TE.NSTRONG.SENS <- M1.TE.NSTRONG.TP/(M1.TE.NSTRONG.TP + M1.TE.NSTRONG.FN)
M1.TE.NSTRONG.SPEC <- M1.TE.NSTRONG.TN/(M1.TE.NSTRONG.TN + M1.TE.NSTRONG.FP) 
M1.TE.NSTRONG.PR   <- M1.TE.NSTRONG.TP/(M1.TE.NSTRONG.TP + M1.TE.NSTRONG.FP) 

##Rank data
RANK.M1.TE <- data.frame(matrix(0,nrow=nrow(PREDICT.TE.CLASS$probs),ncol=10));colnames(RANK.M1.TE) = c("P.NULL","P.WEAK","P.STRONG","P.ACTIVATE","P.NSTRONG","O.NULL","O.WEAK","O.STRONG","O.ACTIVATE","O.NSTRONG")
RANK.M1.TE$P.NULL     <- PREDICT.TE.CLASS$probs[,1,] 
RANK.M1.TE$P.WEAK     <- PREDICT.TE.CLASS$probs[,2,] 
RANK.M1.TE$P.STRONG   <- PREDICT.TE.CLASS$probs[,3,] 
RANK.M1.TE$P.ACTIVATE <- PREDICT.TE.CLASS$probs[,2,] + PREDICT.TE.CLASS$probs[,3,] 
RANK.M1.TE$P.NSTRONG  <- PREDICT.TE.CLASS$probs[,1,] + PREDICT.TE.CLASS$probs[,2,] 
RANK.M1.TE[which(DT.JOINT$pooled.class %in% c("null")),6] <- 1
RANK.M1.TE[which(DT.JOINT$pooled.class %in% c("weak")),7] <- 1
RANK.M1.TE[which(DT.JOINT$pooled.class %in% c("strong")),8] <- 1
RANK.M1.TE[which(DT.JOINT$pooled.class %in% c("weak","strong")),9] <- 1
RANK.M1.TE[which(DT.JOINT$pooled.class %in% c("null","weak")),10] <- 1
RANK.M1.TE[is.na(DT.JOINT$pooled.class),6:10] <- NA
RANK.M1.TE <- RANK.M1.TE[order(RANK.M1.TE$P.STRONG),]
RANK.M1.TE <- na.omit(RANK.M1.TE)

##ROC calculation
ROC.M1.TE.NULL     <- roc(RANK.M1.TE$O.NULL,RANK.M1.TE$P.NULL,        ci=TRUE,plot=FALSE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,print.auc=TRUE,main="M1.TE",xlab="FP/(TN + FP)",ylab="TP/(TP + FN)")
ROC.M1.TE.WEAK     <- roc(RANK.M1.TE$O.WEAK,RANK.M1.TE$P.WEAK,        ci=TRUE,plot=FALSE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,print.auc=TRUE,main="M1.TE",xlab="FP/(TN + FP)",ylab="TP/(TP + FN)")
ROC.M1.TE.STRONG   <- roc(RANK.M1.TE$O.STRONG,RANK.M1.TE$P.STRONG,    ci=TRUE,plot=FALSE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,print.auc=TRUE,main="M1.TE",xlab="FP/(TN + FP)",ylab="TP/(TP + FN)")
ROC.M1.TE.ACTIVATE <- roc(RANK.M1.TE$O.ACTIVATE,RANK.M1.TE$P.ACTIVATE,ci=TRUE,plot=FALSE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,print.auc=TRUE,main="M1.TE",xlab="FP/(TN + FP)",ylab="TP/(TP + FN)")
ROC.M1.TE.NSTRONG  <- roc(RANK.M1.TE$O.NSTRONG,RANK.M1.TE$P.NSTRONG,  ci=TRUE,plot=FALSE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,print.auc=TRUE,main="M1.TE",xlab="FP/(TN + FP)",ylab="TP/(TP + FN)")

##Extract metrics
METRICS.M1.TE.NULL     <- coords(ROC.M1.TE.NULL,     ret = "all", transpose = FALSE)
METRICS.M1.TE.WEAK     <- coords(ROC.M1.TE.WEAK,     ret = "all", transpose = FALSE)
METRICS.M1.TE.STRONG   <- coords(ROC.M1.TE.STRONG,   ret = "all", transpose = FALSE)
METRICS.M1.TE.ACTIVATE <- coords(ROC.M1.TE.ACTIVATE, ret = "all", transpose = FALSE)
METRICS.M1.TE.NSTRONG  <- coords(ROC.M1.TE.NSTRONG,  ret = "all", transpose = FALSE)

#pdf("Accuracy.TE.pdf")
#par(mfrow=c(2,2))
#par(pty="s")
##ROC
#plot(  1-rev(ROC.M1.TE.NULL$specificities),      rev(ROC.M1.TE.NULL$sensitivities),      type="l",col="red", xlab="FP/(TN + FP)",ylab="TP/(TP + FN)",main="M1.TE")
#points(1-rev(ROC.M1.TE.WEAK$specificities),      rev(ROC.M1.TE.WEAK$sensitivities),      type="l",col="blue")
#points(1-rev(ROC.M1.TE.STRONG$specificities),    rev(ROC.M1.TE.STRONG$sensitivities),    type="l",col="green")
#points(1-rev(ROC.M1.TE.ACTIVATE$specificities),  rev(ROC.M1.TE.ACTIVATE$sensitivities),  type="l",col="black")
#points(1-rev(ROC.M1.TE.NSTRONG$specificities),   rev(ROC.M1.TE.NSTRONG$sensitivities),   type="l",col="orange")
#points(1-M1.TE.NULL.SPEC,      M1.TE.NULL.SENS,      pch=19,cex=1,col="red")
#points(1-M1.TE.WEAK.SPEC,      M1.TE.WEAK.SENS,      pch=19,cex=1,col="blue")
#points(1-M1.TE.STRONG.SPEC,    M1.TE.STRONG.SENS,    pch=19,cex=1,col="green")
#points(1-M1.TE.ACTIVATE.SPEC,  M1.TE.ACTIVATE.SENS,  pch=19,cex=1,col="black")
#points(1-M1.TE.NSTRONG.SPEC,   M1.TE.NSTRONG.SENS,   pch=19,cex=1,col="orange")
#abline(a=0,b=1)
#legend(.7,.6,c("Null","Weak","Strong","Activator","Non-Strong"),fill=c("red","blue","green","black","orange"),cex=0.4)

##ROC zoomed in on low FPR
#plot(  1-rev(ROC.M1.TE.NULL$specificities),      rev(ROC.M1.TE.NULL$sensitivities),      type="l",col="red", xlab="FP/(TN + FP)",ylab="TP/(TP + FN)",main="M1.TE",xlim=c(0,0.05))
#points(1-rev(ROC.M1.TE.WEAK$specificities),      rev(ROC.M1.TE.WEAK$sensitivities),      type="l",col="blue")
#points(1-rev(ROC.M1.TE.STRONG$specificities),    rev(ROC.M1.TE.STRONG$sensitivities),    type="l",col="green")
#points(1-rev(ROC.M1.TE.ACTIVATE$specificities),  rev(ROC.M1.TE.ACTIVATE$sensitivities),  type="l",col="black")
#points(1-rev(ROC.M1.TE.NSTRONG$specificities),   rev(ROC.M1.TE.NSTRONG$sensitivities),   type="l",col="orange")
#points(1-M1.TE.NULL.SPEC,      M1.TE.NULL.SENS,      pch=19,cex=1,col="red")
#points(1-M1.TE.WEAK.SPEC,      M1.TE.WEAK.SENS,      pch=19,cex=1,col="blue")
#points(1-M1.TE.STRONG.SPEC,    M1.TE.STRONG.SENS,    pch=19,cex=1,col="green")
#points(1-M1.TE.ACTIVATE.SPEC,  M1.TE.ACTIVATE.SENS,  pch=19,cex=1,col="black")
#points(1-M1.TE.NSTRONG.SPEC,   M1.TE.NSTRONG.SENS,   pch=19,cex=1,col="orange")
#abline(a=0,b=1)
#legend(.03,.5,c("Null","Weak","Strong","Activator","Non-Strong"),fill=c("red","blue","green","black","orange"),cex=0.4)

##Precision/Recall
#plot(  METRICS.M1.TE.NULL$recall,     METRICS.M1.TE.NULL$precision,      type="l",col="red",xlab="TP/(TP+FN)",ylab="TP/(TP+FP)",main="M1.TE",ylim=c(0,1))
#points(METRICS.M1.TE.WEAK$recall,     METRICS.M1.TE.WEAK$precision,      type="l",col="blue")
#points(METRICS.M1.TE.STRONG$recall,   METRICS.M1.TE.STRONG$precision,    type="l",col="green")
#points(METRICS.M1.TE.ACTIVATE$recall, METRICS.M1.TE.ACTIVATE$precision,  type="l",col="black")
#points(METRICS.M1.TE.NSTRONG$recall,  METRICS.M1.TE.NSTRONG$precision,   type="l",col="orange")
#points(M1.TE.NULL.SENS,    M1.TE.NULL.PR,    pch=19,cex=0.6)
#points(M1.TE.WEAK.SENS,    M1.TE.WEAK.PR,    pch=19,cex=0.6)
#points(M1.TE.STRONG.SENS,  M1.TE.STRONG.PR,  pch=19,cex=0.6)
#points(M1.TE.ACTIVATE.SENS,M1.TE.ACTIVATE.PR,pch=19,cex=0.6)
#points(M1.TE.NSTRONG.SENS, M1.TE.NSTRONG.PR, pch=19,cex=0.6)
#legend(.3,.3,c("Null","Weak","Strong","Activator","Non-Strong"),fill=c("red","blue","green","black","orange"),cex=0.4)

##Lift
AVG.RATE.M1.TE.NULL     <- sum(RANK.M1.TE$O.NULL    )/nrow(RANK.M1.TE)
AVG.RATE.M1.TE.WEAK     <- sum(RANK.M1.TE$O.WEAK    )/nrow(RANK.M1.TE)
AVG.RATE.M1.TE.STRONG   <- sum(RANK.M1.TE$O.STRONG  )/nrow(RANK.M1.TE)
AVG.RATE.M1.TE.ACTIVATE <- sum(RANK.M1.TE$O.ACTIVATE)/nrow(RANK.M1.TE)
AVG.RATE.M1.TE.NSTRONG  <- sum(RANK.M1.TE$O.NSTRONG )/nrow(RANK.M1.TE)

#I <- 1:100
#PER.RATE.M1.TE.NULL     <- numeric(length(I))
#PER.RATE.M1.TE.WEAK     <- numeric(length(I))
#PER.RATE.M1.TE.STRONG   <- numeric(length(I))
#PER.RATE.M1.TE.ACTIVATE <- numeric(length(I))
#PER.RATE.M1.TE.NSTRONG  <- numeric(length(I))
#for(i in I) {
#  LB <- floor(quantile(1:length(RANK.M1.TE$O.NULL),0.01*(i-1)))
#  UB <- floor(quantile(1:length(RANK.M1.TE$O.NULL),0.01*i))
#  PER.RATE.M1.TE.NULL[i] <- sum(RANK.M1.TE$O.NULL[LB:UB])/(UB-LB)
#  
#  LB <- floor(quantile(1:length(RANK.M1.TE$O.WEAK),0.01*(i-1)))
#  UB <- floor(quantile(1:length(RANK.M1.TE$O.WEAK),0.01*i))
#  PER.RATE.M1.TE.WEAK[i] <- sum(RANK.M1.TE$O.WEAK[LB:UB])/(UB-LB)
#  
#  LB <- floor(quantile(1:length(RANK.M1.TE$O.STRONG),0.01*(i-1)))
#  UB <- floor(quantile(1:length(RANK.M1.TE$O.STRONG),0.01*i))
#  PER.RATE.M1.TE.STRONG[i] <- sum(RANK.M1.TE$O.STRONG[LB:UB])/(UB-LB)
#  
#  LB <- floor(quantile(1:length(RANK.M1.TE$O.ACTIVATE),0.01*(i-1)))
#  UB <- floor(quantile(1:length(RANK.M1.TE$O.ACTIVATE),0.01*i))
#  PER.RATE.M1.TE.ACTIVATE[i] <- sum(RANK.M1.TE$O.ACTIVATE[LB:UB])/(UB-LB)
#  
#  LB <- floor(quantile(1:length(RANK.M1.TE$O.NSTRONG),0.01*(i-1)))
#  UB <- floor(quantile(1:length(RANK.M1.TE$O.NSTRONG),0.01*i))
#  PER.RATE.M1.TE.NSTRONG[i] <- sum(RANK.M1.TE$O.NSTRONG[LB:UB])/(UB-LB)
#}
#plot(  I,PER.RATE.M1.TE.NULL     /AVG.RATE.M1.TE.NULL     ,type="l",col="red",ylab="Lift",xlab="Percentile",main="M1.TE",ylim=c(0,100))
#points(I,PER.RATE.M1.TE.WEAK     /AVG.RATE.M1.TE.WEAK     ,type="l",col="blue")
#points(I,PER.RATE.M1.TE.STRONG   /AVG.RATE.M1.TE.STRONG   ,type="l",col="green")
#points(I,PER.RATE.M1.TE.ACTIVATE /AVG.RATE.M1.TE.ACTIVATE ,type="l",col="black")
#points(I,PER.RATE.M1.TE.NSTRONG  /AVG.RATE.M1.TE.NSTRONG  ,type="l",col="orange")
#legend(20,60,c("Null","Weak","Strong","Activator","Non-Strong"),fill=c("red","blue","green","black","orange"),cex=0.4)
#dev.off()

##################################
##11. Post-hoc coefficient check##
##################################
####Compare coefficient estimates with the reference free assumption that they are averages over sequences with particular amino acid states

####ME Model
##Collect sets of coefficients
ME.COEFS <- head(ME.MODEL$beta[,ME.BEST.LAMBDA.INDEX],-2)
ME.COEFS.INDEX.B1 <- grep("^X..$",     names(ME.COEFS))
ME.COEFS.INDEX.S1 <- grep("^X..:RE1$", names(ME.COEFS))
ME.COEFS.INDEX.S0 <- grep("^RE1$",     names(ME.COEFS))

##Get average values
#ME.ERE.MEAN <- mean(PREDICT.ME.LINK[1:160000] - ME.COEFS[ME.COEFS.INDEX.S0])
#ME.SRE.MEAN <- mean(PREDICT.ME.LINK[160001:320000] + ME.COEFS[ME.COEFS.INDEX.S0])
#ME.GLOBAL.MEAN <- mean(c(ME.ERE.MEAN,ME.SRE.MEAN))

#ME.CHECK <- matrix(0,nrow=length(ME.COEFS - 1)/2,ncol=8)
#colnames(ME.CHECK) <- c("Estimate.Bind","Estimate.Spec","Estimate.ERE","Estimate.SRE","Post.Bind","Post.Spec","Post.ERE","Post.SRE")
#rownames(ME.CHECK) <- names(ME.COEFS[ME.COEFS.INDEX.B1 ])
#ME.CHECK[,1] <- ME.COEFS[ME.COEFS.INDEX.B1]
#ME.CHECK[,2] <- ME.COEFS[ME.COEFS.INDEX.S1]
#ME.CHECK[,3] <- ME.CHECK[,1] + ME.CHECK[,2]
#ME.CHECK[,4] <- ME.CHECK[,1] - ME.CHECK[,2]

##Main Effects
#I <- 1:length(AAs)
#J <- c("AA1","AA2","AA3","AA4")
#for(j in J) {
#  for(i in I) {
#    INDEX <- (which(J==j) - 1)*20+i
#    SET.E <- which(DT.JOINT[,..j] == AAs[i] & DT.JOINT$RE == "E")
#    SET.S <- which(DT.JOINT[,..j] == AAs[i] & DT.JOINT$RE == "S")
#    SET.E.MEAN <- mean(PREDICT.ME.LINK[SET.E] - ME.COEFS[ME.COEFS.INDEX.S0]) - ME.ERE.MEAN 
#    SET.S.MEAN <- mean(PREDICT.ME.LINK[SET.S] + ME.COEFS[ME.COEFS.INDEX.S0]) - ME.SRE.MEAN 
#    ME.CHECK[INDEX,7] <- SET.E.MEAN 
#    ME.CHECK[INDEX,8] <- SET.S.MEAN 
#    ME.CHECK[INDEX,5] <- (ME.CHECK[INDEX,7] + ME.CHECK[INDEX,8])/2
#    ME.CHECK[INDEX,6] <- (ME.CHECK[INDEX,7] - ME.CHECK[INDEX,8])/2
#  }
#}
#pdf("Post.Hoc.ME.pdf")
#par(pty="s")
#par(mfrow=c(2,2))
#plot(ME.CHECK[1:80,1],ME.CHECK[1:80,5],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="ME Main Bind",col=c(rep("red",20),rep("blue",20),rep("green",20),rep("black",20)))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(ME.CHECK[1:80,2],ME.CHECK[1:80,6],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="ME Main Spec",col=c(rep("red",20),rep("blue",20),rep("green",20),rep("black",20)))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(ME.CHECK[1:80,3],ME.CHECK[1:80,7],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="ME Main ERE",col=c(rep("red",20),rep("blue",20),rep("green",20),rep("black",20)))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(ME.CHECK[1:80,4],ME.CHECK[1:80,8],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="ME Main SRE",col=c(rep("red",20),rep("blue",20),rep("green",20),rep("black",20)))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#dev.off()

##PE Model
PE.COEFS <- head(PE.MODEL$beta[,PE.BEST.LAMBDA.INDEX],-2)
PE.COEFS.INDEX.B1 <- grep("^X..$",         names(PE.COEFS))
PE.COEFS.INDEX.B2 <- grep("^X..:X..$",     names(PE.COEFS))
PE.COEFS.INDEX.S1 <- grep("^X..:RE1$",     names(PE.COEFS))
PE.COEFS.INDEX.S2 <- grep("^X..:X..:RE1$", names(PE.COEFS))
PE.COEFS.INDEX.S0 <- grep("^RE1$",         names(PE.COEFS))

#PE.ERE.MEAN <- mean(PREDICT.PE.LINK[1:160000] - PE.COEFS[PE.COEFS.INDEX.S0])
#PE.SRE.MEAN <- mean(PREDICT.PE.LINK[160001:320000] + PE.COEFS[PE.COEFS.INDEX.S0])
#PE.GLOBAL.MEAN <- mean(c(PE.ERE.MEAN,PE.SRE.MEAN))

#PE.CHECK <- matrix(0,nrow=(length(PE.COEFS)-1)/2,ncol=8)
#colnames(PE.CHECK) <- c("Estimate.Bind","Estimate.Spec","Estimate.ERE","Estimate.SRE","Post.Bind","Post.Spec","Post.ERE","Post.SRE")
#rownames(PE.CHECK) <- names(PE.COEFS[c(PE.COEFS.INDEX.B1,PE.COEFS.INDEX.B2)])
#PE.CHECK[,1] <- PE.COEFS[c(PE.COEFS.INDEX.B1,PE.COEFS.INDEX.B2)]
#PE.CHECK[,2] <- PE.COEFS[c(PE.COEFS.INDEX.S1,PE.COEFS.INDEX.S2)]
#PE.CHECK[,3] <- PE.CHECK[,1] + PE.CHECK[,2] 
#PE.CHECK[,4] <- PE.CHECK[,1] - PE.CHECK[,2] 

##Main Effects
#I <- 1:length(AAs)
#J <- c("AA1","AA2","AA3","AA4")
#for(j in J) {
#  for(i in I) {
#    INDEX <- (which(J==j) - 1)*20+i
#    SET.E <- which(DT.JOINT[,..j] == AAs[i] & DT.JOINT$RE == "E")
#    SET.S <- which(DT.JOINT[,..j] == AAs[i] & DT.JOINT$RE == "S")
#    SET.E.MEAN <- mean(PREDICT.PE.LINK[SET.E] - PE.COEFS[PE.COEFS.INDEX.S0]) - PE.ERE.MEAN
#    SET.S.MEAN <- mean(PREDICT.PE.LINK[SET.S] + PE.COEFS[PE.COEFS.INDEX.S0]) - PE.SRE.MEAN
#    PE.CHECK[INDEX,7] <- SET.E.MEAN 
#    PE.CHECK[INDEX,8] <- SET.S.MEAN 
#    PE.CHECK[INDEX,5] <- (PE.CHECK[INDEX,7] + PE.CHECK[INDEX,8])/2
#    PE.CHECK[INDEX,6] <- (PE.CHECK[INDEX,7] - PE.CHECK[INDEX,8])/2
#  }
#}

##Pairwise Epistasis
#SITES <- c("AA1","AA2","AA3","AA4")
#I <- 1:length(AAs)
#J <- 1:length(AAs)
#K <- 1:(length(SITES)-1)
#COUNT <- 0
#for(k in K) {
#  M <- (k+1):(length(K)+1)
#  for(m in M) {
#    COUNT <- COUNT + 1
#    for(j in J) {
#      for(i in I) {
#        INDEX <- 80 + (COUNT-1)*400+(j-1)*20+i
#        SITE1 <- SITES[k]
#        SITE2 <- SITES[m]
#        SET.E1  <- which(DT.JOINT[,..SITE1] == AAs[i] & DT.JOINT$RE == "E")
#        SET.E2  <- which(DT.JOINT[,..SITE2] == AAs[j] & DT.JOINT$RE == "E")
#        SET.E12 <- which(DT.JOINT[,..SITE1] == AAs[i] & DT.JOINT[,..SITE2] == AAs[j] & DT.JOINT$RE == "E")
#        
#        SET.S1  <- which(DT.JOINT[,..SITE1] == AAs[i] & DT.JOINT$RE == "S")
#        SET.S2  <- which(DT.JOINT[,..SITE2] == AAs[j] & DT.JOINT$RE == "S")
#        SET.S12 <- which(DT.JOINT[,..SITE1] == AAs[i] & DT.JOINT[,..SITE2] == AAs[j] & DT.JOINT$RE == "S")
#        
#        SET.E1.MEAN  <- mean(PREDICT.PE.LINK[SET.E1]  - PE.COEFS[PE.COEFS.INDEX.S0]) - PE.ERE.MEAN
#        SET.E2.MEAN  <- mean(PREDICT.PE.LINK[SET.E2]  - PE.COEFS[PE.COEFS.INDEX.S0]) - PE.ERE.MEAN
#        SET.E12.MEAN <- mean(PREDICT.PE.LINK[SET.E12] - PE.COEFS[PE.COEFS.INDEX.S0]) - PE.ERE.MEAN
#        SET.S1.MEAN  <- mean(PREDICT.PE.LINK[SET.S1]  + PE.COEFS[PE.COEFS.INDEX.S0]) - PE.SRE.MEAN
#        SET.S2.MEAN  <- mean(PREDICT.PE.LINK[SET.S2]  + PE.COEFS[PE.COEFS.INDEX.S0]) - PE.SRE.MEAN
#        SET.S12.MEAN <- mean(PREDICT.PE.LINK[SET.S12] + PE.COEFS[PE.COEFS.INDEX.S0]) - PE.SRE.MEAN
#
#        PE.CHECK[INDEX,7] <- SET.E12.MEAN - SET.E1.MEAN - SET.E2.MEAN
#        PE.CHECK[INDEX,8] <- SET.S12.MEAN - SET.S1.MEAN - SET.S2.MEAN
#        PE.CHECK[INDEX,5] <- (PE.CHECK[INDEX,7] + PE.CHECK[INDEX,8])/2
#        PE.CHECK[INDEX,6] <- (PE.CHECK[INDEX,7] - PE.CHECK[INDEX,8])/2
#      }
#    }
#  }
#}
#pdf("Post.Hoc.PE.pdf")
#par(pty="s")
#par(mfrow=c(2,2))
#plot(PE.CHECK[1:80,1],PE.CHECK[1:80,5],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="PE Main Bind",col=c(rep("red",20),rep("blue",20),rep("green",20),rep("black",20)))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(PE.CHECK[1:80,2],PE.CHECK[1:80,6],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="PE Main Spec",col=c(rep("red",20),rep("blue",20),rep("green",20),rep("black",20)))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(PE.CHECK[1:80,3],PE.CHECK[1:80,7],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="PE Main ERE",col=c(rep("red",20),rep("blue",20),rep("green",20),rep("black",20)))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(PE.CHECK[1:80,4],PE.CHECK[1:80,8],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="PE Main SRE",col=c(rep("red",20),rep("blue",20),rep("green",20),rep("black",20)))
#abline(a=0,b=1);abline(h=0);abline(v=0)

#par(mfrow=c(3,2))
#plot(  PE.CHECK[81:480,   1],PE.CHECK[81:480,   5],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="PE BIND 12",col=c("red"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  PE.CHECK[481:880,  1],PE.CHECK[481:880,  5],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="PE BIND 13",col=c("blue"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  PE.CHECK[881:1280, 1],PE.CHECK[881:1280, 5],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="PE BIND 14",col=c("green"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  PE.CHECK[1281:1680,1],PE.CHECK[1281:1680,5],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="PE BIND 23",col=c("black"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  PE.CHECK[1681:2080,1],PE.CHECK[1681:2080,5],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="PE BIND 24",col=c("purple"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  PE.CHECK[2081:2480,1],PE.CHECK[2081:2480,5],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="PE BIND 34",col=c("orange"))
#abline(a=0,b=1);abline(h=0);abline(v=0)

#plot(  PE.CHECK[81:480,   2],PE.CHECK[81:480,   6],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="PE SPEC 12",col=c("red"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  PE.CHECK[481:880,  2],PE.CHECK[481:880,  6],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="PE SPEC 13",col=c("blue"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  PE.CHECK[881:1280, 2],PE.CHECK[881:1280, 6],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="PE SPEC 14",col=c("green"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  PE.CHECK[1281:1680,2],PE.CHECK[1281:1680,6],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="PE SPEC 23",col=c("black"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  PE.CHECK[1681:2080,2],PE.CHECK[1681:2080,6],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="PE SPEC 24",col=c("purple"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  PE.CHECK[2081:2480,2],PE.CHECK[2081:2480,6],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="PE SPEC 34",col=c("orange"))
#abline(a=0,b=1);abline(h=0);abline(v=0)

#plot(  PE.CHECK[81:480,   3],PE.CHECK[81:480,   7],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="PE ERE 12",col=c("red"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  PE.CHECK[481:880,  3],PE.CHECK[481:880,  7],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="PE ERE 13",col=c("blue"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  PE.CHECK[881:1280, 3],PE.CHECK[881:1280, 7],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="PE ERE 14",col=c("green"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  PE.CHECK[1281:1680,3],PE.CHECK[1281:1680,7],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="PE ERE 23",col=c("black"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  PE.CHECK[1681:2080,3],PE.CHECK[1681:2080,7],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="PE ERE 24",col=c("purple"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  PE.CHECK[2081:2480,3],PE.CHECK[2081:2480,7],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="PE ERE 34",col=c("orange"))
#abline(a=0,b=1);abline(h=0);abline(v=0)

#plot(  PE.CHECK[81:480,   4],PE.CHECK[81:480,   8],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="PE SRE 12",col=c("red"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  PE.CHECK[481:880,  4],PE.CHECK[481:880,  8],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="PE SRE 13",col=c("blue"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  PE.CHECK[881:1280, 4],PE.CHECK[881:1280, 8],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="PE SRE 14",col=c("green"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  PE.CHECK[1281:1680,4],PE.CHECK[1281:1680,8],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="PE SRE 23",col=c("black"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  PE.CHECK[1681:2080,4],PE.CHECK[1681:2080,8],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="PE SRE 24",col=c("purple"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  PE.CHECK[2081:2480,4],PE.CHECK[2081:2480,8],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="PE SRE 34",col=c("orange"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#dev.off()

##TE Model
TE.COEFS <- head(TE.MODEL$beta[,TE.BEST.LAMBDA.INDEX],-2)
TE.COEFS.INDEX.B1 <- grep("^X..$",             names(TE.COEFS))
TE.COEFS.INDEX.B2 <- grep("^X..:X..$",         names(TE.COEFS))
TE.COEFS.INDEX.B3 <- grep("^X..:X..:X..$",     names(TE.COEFS))
TE.COEFS.INDEX.S1 <- grep("^X..:RE1$",         names(TE.COEFS))
TE.COEFS.INDEX.S2 <- grep("^X..:X..:RE1$",     names(TE.COEFS))
TE.COEFS.INDEX.S3 <- grep("^X..:X..:X..:RE1$", names(TE.COEFS))
TE.COEFS.INDEX.S0 <- grep("^RE1$",             names(TE.COEFS))

#TE.ERE.MEAN <- mean(PREDICT.TE.LINK[1:160000] - TE.COEFS[TE.COEFS.INDEX.S0])
#TE.SRE.MEAN <- mean(PREDICT.TE.LINK[160001:320000] + TE.COEFS[TE.COEFS.INDEX.S0])
#TE.GLOBAL.MEAN <- mean(c(TE.ERE.MEAN,TE.SRE.MEAN))

#TE.CHECK <- matrix(0,nrow=(length(TE.COEFS)-1)/2,ncol=8)
#colnames(TE.CHECK) <- c("Estimate.Bind","Estimate.Spec","Estimate.ERE","Estimate.SRE","Post.Bind","Post.Spec","Post.ERE","Post.SRE")
#rownames(TE.CHECK) <- names(TE.COEFS[c(TE.COEFS.INDEX.B1,TE.COEFS.INDEX.B2,TE.COEFS.INDEX.B3)])
#TE.CHECK[,1] <- TE.COEFS[c(TE.COEFS.INDEX.B1,TE.COEFS.INDEX.B2,TE.COEFS.INDEX.B3)]
#TE.CHECK[,2] <- TE.COEFS[c(TE.COEFS.INDEX.S1,TE.COEFS.INDEX.S2,TE.COEFS.INDEX.S3)]
#TE.CHECK[,3] <- TE.CHECK[,1] + TE.CHECK[,2] 
#TE.CHECK[,4] <- TE.CHECK[,1] - TE.CHECK[,2] 

##Main Effects
#I <- 1:length(AAs)
#J <- c("AA1","AA2","AA3","AA4")
#for(j in J) {
#  for(i in I) {
#    INDEX <- (which(J==j) - 1)*20+i
#    SET.E <- which(DT.JOINT[,..j] == AAs[i] & DT.JOINT$RE == "E")
#    SET.S <- which(DT.JOINT[,..j] == AAs[i] & DT.JOINT$RE == "S")
#    SET.E.MEAN <- mean(PREDICT.TE.LINK[SET.E] - TE.COEFS[TE.COEFS.INDEX.S0]) - TE.ERE.MEAN
#    SET.S.MEAN <- mean(PREDICT.TE.LINK[SET.S] + TE.COEFS[TE.COEFS.INDEX.S0]) - TE.SRE.MEAN
#    TE.CHECK[INDEX,7] <- SET.E.MEAN 
#    TE.CHECK[INDEX,8] <- SET.S.MEAN 
#    TE.CHECK[INDEX,5] <- (TE.CHECK[INDEX,7] + TE.CHECK[INDEX,8])/2
#    TE.CHECK[INDEX,6] <- (TE.CHECK[INDEX,7] - TE.CHECK[INDEX,8])/2
#  }
#}

##Pairwise Epistasis
#SITES <- c("AA1","AA2","AA3","AA4")
#I <- 1:length(AAs)
#J <- 1:length(AAs)
#K <- 1:(length(SITES)-1)
#COUNT <- 0
#for(k in K) {
#  M <- (k+1):(length(K)+1)
#  for(m in M) {
#    COUNT <- COUNT + 1
#    for(j in J) {
#      for(i in I) {
#        INDEX <- 80 + (COUNT-1)*400+(j-1)*20+i
#        SITE1 <- SITES[k]
#        SITE2 <- SITES[m]
#        SET.E1  <- which(DT.JOINT[,..SITE1] == AAs[i] & DT.JOINT$RE == "E")
#        SET.E2  <- which(DT.JOINT[,..SITE2] == AAs[j] & DT.JOINT$RE == "E")
#        SET.E12 <- which(DT.JOINT[,..SITE1] == AAs[i] & DT.JOINT[,..SITE2] == AAs[j] & DT.JOINT$RE == "E")
#        
#        SET.S1  <- which(DT.JOINT[,..SITE1] == AAs[i] & DT.JOINT$RE == "S")
#        SET.S2  <- which(DT.JOINT[,..SITE2] == AAs[j] & DT.JOINT$RE == "S")
#        SET.S12 <- which(DT.JOINT[,..SITE1] == AAs[i] & DT.JOINT[,..SITE2] == AAs[j] & DT.JOINT$RE == "S")
#        
#        SET.E1.MEAN  <- mean(PREDICT.TE.LINK[SET.E1]  - TE.COEFS[TE.COEFS.INDEX.S0]) - TE.ERE.MEAN
#        SET.E2.MEAN  <- mean(PREDICT.TE.LINK[SET.E2]  - TE.COEFS[TE.COEFS.INDEX.S0]) - TE.ERE.MEAN
#        SET.E12.MEAN <- mean(PREDICT.TE.LINK[SET.E12] - TE.COEFS[TE.COEFS.INDEX.S0]) - TE.ERE.MEAN
#        SET.S1.MEAN  <- mean(PREDICT.TE.LINK[SET.S1]  + TE.COEFS[TE.COEFS.INDEX.S0]) - TE.SRE.MEAN
#        SET.S2.MEAN  <- mean(PREDICT.TE.LINK[SET.S2]  + TE.COEFS[TE.COEFS.INDEX.S0]) - TE.SRE.MEAN
#        SET.S12.MEAN <- mean(PREDICT.TE.LINK[SET.S12] + TE.COEFS[TE.COEFS.INDEX.S0]) - TE.SRE.MEAN
#        
#        TE.CHECK[INDEX,7] <- SET.E12.MEAN - SET.E1.MEAN - SET.E2.MEAN
#        TE.CHECK[INDEX,8] <- SET.S12.MEAN - SET.S1.MEAN - SET.S2.MEAN
#        TE.CHECK[INDEX,5] <- (TE.CHECK[INDEX,7] + TE.CHECK[INDEX,8])/2
#        TE.CHECK[INDEX,6] <- (TE.CHECK[INDEX,7] - TE.CHECK[INDEX,8])/2
#      }
#    }
#  }
#}

##Third order epistasis
#SITES <- c("AA1","AA2","AA3","AA4")
#I <- 1:length(AAs)
#J <- 1:length(AAs)
#K <- 1:length(AAs)
#L <- 1:(length(SITES)-2)
#COUNT <- 0
#for(l in L) {
#  M <- (l+1):3
#  for(m in M) {
#    N <- (m+1):4
#    for(n in N) {
#      COUNT <- COUNT + 1
#      for(k in K) {
#        for(j in J) {
#          for(i in I) {
#            INDEX <- 2480 + (COUNT-1)*8000+(k-1)*400+(j-1)*20+i
#            SITE1 <- SITES[l]
#            SITE2 <- SITES[m]
#            SITE3 <- SITES[n]
#            
#            SET.E1   <- which(DT.JOINT[,..SITE1] == AAs[i] & DT.JOINT$RE == "E")
#            SET.E2   <- which(DT.JOINT[,..SITE2] == AAs[j] & DT.JOINT$RE == "E")
#            SET.E3   <- which(DT.JOINT[,..SITE3] == AAs[k] & DT.JOINT$RE == "E")
#            SET.E12  <- which(DT.JOINT[,..SITE1] == AAs[i] & DT.JOINT[,..SITE2] == AAs[j] & DT.JOINT$RE == "E")
#            SET.E13  <- which(DT.JOINT[,..SITE1] == AAs[i] & DT.JOINT[,..SITE3] == AAs[k] & DT.JOINT$RE == "E")
#            SET.E23  <- which(DT.JOINT[,..SITE2] == AAs[j] & DT.JOINT[,..SITE3] == AAs[k] & DT.JOINT$RE == "E")
#            SET.E123 <- which(DT.JOINT[,..SITE1] == AAs[i] & DT.JOINT[,..SITE2] == AAs[j] & DT.JOINT[,..SITE3] == AAs[k] & DT.JOINT$RE == "E") 
#            
#            SET.S1   <- which(DT.JOINT[,..SITE1] == AAs[i] & DT.JOINT$RE == "S")
#            SET.S2   <- which(DT.JOINT[,..SITE2] == AAs[j] & DT.JOINT$RE == "S")
#            SET.S3   <- which(DT.JOINT[,..SITE3] == AAs[k] & DT.JOINT$RE == "S")
#            SET.S12  <- which(DT.JOINT[,..SITE1] == AAs[i] & DT.JOINT[,..SITE2] == AAs[j] & DT.JOINT$RE == "S")
#            SET.S13  <- which(DT.JOINT[,..SITE1] == AAs[i] & DT.JOINT[,..SITE3] == AAs[k] & DT.JOINT$RE == "S")
#            SET.S23  <- which(DT.JOINT[,..SITE2] == AAs[j] & DT.JOINT[,..SITE3] == AAs[k] & DT.JOINT$RE == "S")
#            SET.S123 <- which(DT.JOINT[,..SITE1] == AAs[i] & DT.JOINT[,..SITE2] == AAs[j] & DT.JOINT[,..SITE3] == AAs[k] & DT.JOINT$RE == "S") 
#            
#            SET.E1.MEAN   <- mean(PREDICT.TE.LINK[SET.E1]   - TE.COEFS[TE.COEFS.INDEX.S0]) - TE.ERE.MEAN
#            SET.E2.MEAN   <- mean(PREDICT.TE.LINK[SET.E2]   - TE.COEFS[TE.COEFS.INDEX.S0]) - TE.ERE.MEAN
#            SET.E3.MEAN   <- mean(PREDICT.TE.LINK[SET.E3]   - TE.COEFS[TE.COEFS.INDEX.S0]) - TE.ERE.MEAN
#            SET.E12.MEAN  <- mean(PREDICT.TE.LINK[SET.E12]  - TE.COEFS[TE.COEFS.INDEX.S0]) - TE.ERE.MEAN
#            SET.E13.MEAN  <- mean(PREDICT.TE.LINK[SET.E13]  - TE.COEFS[TE.COEFS.INDEX.S0]) - TE.ERE.MEAN
#            SET.E23.MEAN  <- mean(PREDICT.TE.LINK[SET.E23]  - TE.COEFS[TE.COEFS.INDEX.S0]) - TE.ERE.MEAN
#            SET.E123.MEAN <- mean(PREDICT.TE.LINK[SET.E123] - TE.COEFS[TE.COEFS.INDEX.S0]) - TE.ERE.MEAN
#            
#            SET.S1.MEAN   <- mean(PREDICT.TE.LINK[SET.S1]   + TE.COEFS[TE.COEFS.INDEX.S0]) - TE.SRE.MEAN
#            SET.S2.MEAN   <- mean(PREDICT.TE.LINK[SET.S2]   + TE.COEFS[TE.COEFS.INDEX.S0]) - TE.SRE.MEAN
#            SET.S3.MEAN   <- mean(PREDICT.TE.LINK[SET.S3]   + TE.COEFS[TE.COEFS.INDEX.S0]) - TE.SRE.MEAN
#            SET.S12.MEAN  <- mean(PREDICT.TE.LINK[SET.S12]  + TE.COEFS[TE.COEFS.INDEX.S0]) - TE.SRE.MEAN
#            SET.S13.MEAN  <- mean(PREDICT.TE.LINK[SET.S13]  + TE.COEFS[TE.COEFS.INDEX.S0]) - TE.SRE.MEAN
#            SET.S23.MEAN  <- mean(PREDICT.TE.LINK[SET.S23]  + TE.COEFS[TE.COEFS.INDEX.S0]) - TE.SRE.MEAN
#            SET.S123.MEAN <- mean(PREDICT.TE.LINK[SET.S123] + TE.COEFS[TE.COEFS.INDEX.S0]) - TE.SRE.MEAN
#            
#            TE.CHECK[INDEX,7] <- SET.E123.MEAN - (SET.E12.MEAN + SET.E13.MEAN + SET.E23.MEAN - SET.E1.MEAN - SET.E2.MEAN - SET.E3.MEAN)
#            TE.CHECK[INDEX,8] <- SET.S123.MEAN - (SET.S12.MEAN + SET.S13.MEAN + SET.S23.MEAN - SET.S1.MEAN - SET.S2.MEAN - SET.S3.MEAN)
#            TE.CHECK[INDEX,5] <- (TE.CHECK[INDEX,7] + TE.CHECK[INDEX,8])/2
#            TE.CHECK[INDEX,6] <- (TE.CHECK[INDEX,7] - TE.CHECK[INDEX,8])/2
#          }
#        }
#      }
#    }
#  }
#}
#pdf("Post.Hoc.TE.pdf")
#par(pty="s")
#par(mfrow=c(2,2))
#plot(TE.CHECK[1:80,1],TE.CHECK[1:80,5],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE Main Bind",col=c(rep("red",20),rep("blue",20),rep("green",20),rep("black",20)))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(TE.CHECK[1:80,2],TE.CHECK[1:80,6],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE Main Spec",col=c(rep("red",20),rep("blue",20),rep("green",20),rep("black",20)))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(TE.CHECK[1:80,3],TE.CHECK[1:80,7],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE Main ERE",col=c(rep("red",20),rep("blue",20),rep("green",20),rep("black",20)))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(TE.CHECK[1:80,4],TE.CHECK[1:80,8],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE Main SRE",col=c(rep("red",20),rep("blue",20),rep("green",20),rep("black",20)))
#abline(a=0,b=1);abline(h=0);abline(v=0)

#par(mfrow=c(3,2))
#plot(  TE.CHECK[81:480,   1],TE.CHECK[81:480,   5],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE BIND 12",col=c("red"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  TE.CHECK[481:880,  1],TE.CHECK[481:880,  5],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE BIND 13",col=c("blue"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  TE.CHECK[881:1280, 1],TE.CHECK[881:1280, 5],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE BIND 14",col=c("green"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  TE.CHECK[1281:1680,1],TE.CHECK[1281:1680,5],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE BIND 23",col=c("black"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  TE.CHECK[1681:2080,1],TE.CHECK[1681:2080,5],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE BIND 24",col=c("purple"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  TE.CHECK[2081:2480,1],TE.CHECK[2081:2480,5],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE BIND 34",col=c("orange"))
#abline(a=0,b=1);abline(h=0);abline(v=0)

#plot(  TE.CHECK[81:480,   2],TE.CHECK[81:480,   6],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE SPEC 12",col=c("red"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  TE.CHECK[481:880,  2],TE.CHECK[481:880,  6],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE SPEC 13",col=c("blue"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  TE.CHECK[881:1280, 2],TE.CHECK[881:1280, 6],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE SPEC 14",col=c("green"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  TE.CHECK[1281:1680,2],TE.CHECK[1281:1680,6],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE SPEC 23",col=c("black"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  TE.CHECK[1681:2080,2],TE.CHECK[1681:2080,6],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE SPEC 24",col=c("purple"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  TE.CHECK[2081:2480,2],TE.CHECK[2081:2480,6],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE SPEC 34",col=c("orange"))
#abline(a=0,b=1);abline(h=0);abline(v=0)

#plot(  TE.CHECK[81:480,   3],TE.CHECK[81:480,   7],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE ERE 12",col=c("red"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  TE.CHECK[481:880,  3],TE.CHECK[481:880,  7],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE ERE 13",col=c("blue"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  TE.CHECK[881:1280, 3],TE.CHECK[881:1280, 7],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE ERE 14",col=c("green"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  TE.CHECK[1281:1680,3],TE.CHECK[1281:1680,7],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE ERE 23",col=c("black"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  TE.CHECK[1681:2080,3],TE.CHECK[1681:2080,7],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE ERE 24",col=c("purple"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  TE.CHECK[2081:2480,3],TE.CHECK[2081:2480,7],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE ERE 34",col=c("orange"))
#abline(a=0,b=1);abline(h=0);abline(v=0)

#plot(  TE.CHECK[81:480,   4],TE.CHECK[81:480,   8],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE SRE 12",col=c("red"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  TE.CHECK[481:880,  4],TE.CHECK[481:880,  8],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE SRE 13",col=c("blue"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  TE.CHECK[881:1280, 4],TE.CHECK[881:1280, 8],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE SRE 14",col=c("green"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  TE.CHECK[1281:1680,4],TE.CHECK[1281:1680,8],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE SRE 23",col=c("black"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  TE.CHECK[1681:2080,4],TE.CHECK[1681:2080,8],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE SRE 24",col=c("purple"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  TE.CHECK[2081:2480,4],TE.CHECK[2081:2480,8],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE SRE 34",col=c("orange"))
#abline(a=0,b=1);abline(h=0);abline(v=0)

#par(mfrow=c(2,2))
#plot(  TE.CHECK[2481:10480, 1],TE.CHECK[2481:10480, 5],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE BIND 123",col=c("red"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  TE.CHECK[10481:18480,1],TE.CHECK[10481:18480,5],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE BIND 124",col=c("blue"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  TE.CHECK[18481:26480,1],TE.CHECK[18481:26480,5],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE BIND 134",col=c("green"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  TE.CHECK[26481:34480,1],TE.CHECK[26481:34480,5],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE BIND 234",col=c("black"))
#abline(a=0,b=1);abline(h=0);abline(v=0)

#plot(  TE.CHECK[2481:10480, 2],TE.CHECK[2481:10480, 6],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE SPEC 123",col=c("red"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  TE.CHECK[10481:18480,2],TE.CHECK[10481:18480,6],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE SPEC 124",col=c("blue"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  TE.CHECK[18481:26480,2],TE.CHECK[18481:26480,6],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE SPEC 134",col=c("green"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  TE.CHECK[26481:34480,2],TE.CHECK[26481:34480,6],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE SPEC 234",col=c("black"))
#abline(a=0,b=1);abline(h=0);abline(v=0)

#plot(  TE.CHECK[2481:10480, 3],TE.CHECK[2481:10480, 7],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE ERE 123",col=c("red"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  TE.CHECK[10481:18480,3],TE.CHECK[10481:18480,7],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE ERE 124",col=c("blue"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  TE.CHECK[18481:26480,3],TE.CHECK[18481:26480,7],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE ERE 134",col=c("green"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  TE.CHECK[26481:34480,3],TE.CHECK[26481:34480,7],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE ERE 234",col=c("black"))
#abline(a=0,b=1);abline(h=0);abline(v=0)

#plot(  TE.CHECK[2481:10480, 4],TE.CHECK[2481:10480, 8],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE SRE 123",col=c("red"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  TE.CHECK[10481:18480,4],TE.CHECK[10481:18480,8],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE SRE 124",col=c("blue"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  TE.CHECK[18481:26480,4],TE.CHECK[18481:26480,8],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE SRE 134",col=c("green"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#plot(  TE.CHECK[26481:34480,4],TE.CHECK[26481:34480,8],pch=19,cex=0.6,xlim=c(-6,4),ylim=c(-6,4),xlab="Model Estimate",ylab="Post hoc check",main="TE SRE 234",col=c("black"))
#abline(a=0,b=1);abline(h=0);abline(v=0)
#dev.off()

###############################
##12. Coefficient Constraints##
###############################
####Adjust coefficients to fit the reference free assumption that the average of effects at a site/site combination is zero

####ME Model
##Index site specific coefficients
ME.COEFS.INDEX.B1.X1 <- grep("^X1.$",    names(ME.COEFS))
ME.COEFS.INDEX.B1.X2 <- grep("^X2.$",    names(ME.COEFS))
ME.COEFS.INDEX.B1.X3 <- grep("^X3.$",    names(ME.COEFS))
ME.COEFS.INDEX.B1.X4 <- grep("^X4.$",    names(ME.COEFS))
ME.COEFS.INDEX.S1.X1 <- grep("^X1.:RE1$",names(ME.COEFS))
ME.COEFS.INDEX.S1.X2 <- grep("^X2.:RE1$",names(ME.COEFS))
ME.COEFS.INDEX.S1.X3 <- grep("^X3.:RE1$",names(ME.COEFS))
ME.COEFS.INDEX.S1.X4 <- grep("^X4.:RE1$",names(ME.COEFS))

##Extract coefficient effects
ME.COEFS.EFFECT.B0  <- 0
ME.COEFS.EFFECT.S0  <- ME.COEFS[ME.COEFS.INDEX.S0]

ME.COEFS.EFFECT.B1.X1  <- ME.COEFS[ME.COEFS.INDEX.B1.X1]
ME.COEFS.EFFECT.B1.X2  <- ME.COEFS[ME.COEFS.INDEX.B1.X2]
ME.COEFS.EFFECT.B1.X3  <- ME.COEFS[ME.COEFS.INDEX.B1.X3]
ME.COEFS.EFFECT.B1.X4  <- ME.COEFS[ME.COEFS.INDEX.B1.X4]
ME.COEFS.EFFECT.S1.X1  <- ME.COEFS[ME.COEFS.INDEX.S1.X1]
ME.COEFS.EFFECT.S1.X2  <- ME.COEFS[ME.COEFS.INDEX.S1.X2]
ME.COEFS.EFFECT.S1.X3  <- ME.COEFS[ME.COEFS.INDEX.S1.X3]
ME.COEFS.EFFECT.S1.X4  <- ME.COEFS[ME.COEFS.INDEX.S1.X4]

##Impose constraints-Average of terms per site should equal zero
ME.COEFS.EFFECT.ADJ.B1.X1  <- ME.COEFS.EFFECT.B1.X1 - mean(ME.COEFS.EFFECT.B1.X1) 
ME.COEFS.EFFECT.ADJ.B1.X2  <- ME.COEFS.EFFECT.B1.X2 - mean(ME.COEFS.EFFECT.B1.X2) 
ME.COEFS.EFFECT.ADJ.B1.X3  <- ME.COEFS.EFFECT.B1.X3 - mean(ME.COEFS.EFFECT.B1.X3) 
ME.COEFS.EFFECT.ADJ.B1.X4  <- ME.COEFS.EFFECT.B1.X4 - mean(ME.COEFS.EFFECT.B1.X4)
ME.COEFS.EFFECT.ADJ.S1.X1  <- ME.COEFS.EFFECT.S1.X1 - mean(ME.COEFS.EFFECT.S1.X1) 
ME.COEFS.EFFECT.ADJ.S1.X2  <- ME.COEFS.EFFECT.S1.X2 - mean(ME.COEFS.EFFECT.S1.X2) 
ME.COEFS.EFFECT.ADJ.S1.X3  <- ME.COEFS.EFFECT.S1.X3 - mean(ME.COEFS.EFFECT.S1.X3) 
ME.COEFS.EFFECT.ADJ.S1.X4  <- ME.COEFS.EFFECT.S1.X4 - mean(ME.COEFS.EFFECT.S1.X4)

##Global effects are constrained so that genetic scores stay the same across all sequences after constraint
ME.COEFS.EFFECT.ADJ.B0     <- ME.COEFS.EFFECT.B0 + mean(ME.COEFS.EFFECT.B1.X1) + mean(ME.COEFS.EFFECT.B1.X2) + mean(ME.COEFS.EFFECT.B1.X3) + mean(ME.COEFS.EFFECT.B1.X4)
ME.COEFS.EFFECT.ADJ.S0     <- ME.COEFS.EFFECT.S0 + mean(ME.COEFS.EFFECT.S1.X1) + mean(ME.COEFS.EFFECT.S1.X2) + mean(ME.COEFS.EFFECT.S1.X3) + mean(ME.COEFS.EFFECT.S1.X4)

##Collect adjusted terms
ME.COEFS.ADJ <- c(ME.COEFS.EFFECT.ADJ.B1.X1,ME.COEFS.EFFECT.ADJ.B1.X2,ME.COEFS.EFFECT.ADJ.B1.X3,ME.COEFS.EFFECT.ADJ.B1.X4,
                  ME.COEFS.EFFECT.ADJ.S1.X1,ME.COEFS.EFFECT.ADJ.S1.X2,ME.COEFS.EFFECT.ADJ.S1.X3,ME.COEFS.EFFECT.ADJ.S1.X4);names(ME.COEFS.ADJ) <- names(ME.COEFS[-ME.COEFS.INDEX.S0])
#save(ME.COEFS.ADJ,file="ME.COEFS.ADJ.rda")
#save(ME.COEFS.EFFECT.ADJ.B0,file="ME.COEFS.EFFECT.ADJ.B0.rda")
#save(ME.COEFS.EFFECT.ADJ.S0,file="ME.COEFS.EFFECT.ADJ.S0.rda")

####PE Model
##Index site specific coefficients
PE.COEFS.INDEX.B1.X1   <- grep("^X1.$",        names(PE.COEFS))
PE.COEFS.INDEX.B1.X2   <- grep("^X2.$",        names(PE.COEFS))
PE.COEFS.INDEX.B1.X3   <- grep("^X3.$",        names(PE.COEFS))
PE.COEFS.INDEX.B1.X4   <- grep("^X4.$",        names(PE.COEFS))
PE.COEFS.INDEX.S1.X1   <- grep("^X1.:RE1$",    names(PE.COEFS))
PE.COEFS.INDEX.S1.X2   <- grep("^X2.:RE1$",    names(PE.COEFS))
PE.COEFS.INDEX.S1.X3   <- grep("^X3.:RE1$",    names(PE.COEFS))
PE.COEFS.INDEX.S1.X4   <- grep("^X4.:RE1$",    names(PE.COEFS))
PE.COEFS.INDEX.B2.X1X2 <- grep("^X1.:X2.$",    names(PE.COEFS))
PE.COEFS.INDEX.B2.X1X3 <- grep("^X1.:X3.$",    names(PE.COEFS))
PE.COEFS.INDEX.B2.X1X4 <- grep("^X1.:X4.$",    names(PE.COEFS))
PE.COEFS.INDEX.B2.X2X3 <- grep("^X2.:X3.$",    names(PE.COEFS))
PE.COEFS.INDEX.B2.X2X4 <- grep("^X2.:X4.$",    names(PE.COEFS))
PE.COEFS.INDEX.B2.X3X4 <- grep("^X3.:X4.$",    names(PE.COEFS))
PE.COEFS.INDEX.S2.X1X2 <- grep("^X1.:X2.:RE1$",names(PE.COEFS))
PE.COEFS.INDEX.S2.X1X3 <- grep("^X1.:X3.:RE1$",names(PE.COEFS))
PE.COEFS.INDEX.S2.X1X4 <- grep("^X1.:X4.:RE1$",names(PE.COEFS))
PE.COEFS.INDEX.S2.X2X3 <- grep("^X2.:X3.:RE1$",names(PE.COEFS))
PE.COEFS.INDEX.S2.X2X4 <- grep("^X2.:X4.:RE1$",names(PE.COEFS))
PE.COEFS.INDEX.S2.X3X4 <- grep("^X3.:X4.:RE1$",names(PE.COEFS))

##Extract coefficient effects
PE.COEFS.EFFECT.B0  <- 0
PE.COEFS.EFFECT.S0  <- PE.COEFS[PE.COEFS.INDEX.S0]
PE.COEFS.EFFECT.B1.X1  <- PE.COEFS[PE.COEFS.INDEX.B1.X1]
PE.COEFS.EFFECT.B1.X2  <- PE.COEFS[PE.COEFS.INDEX.B1.X2]
PE.COEFS.EFFECT.B1.X3  <- PE.COEFS[PE.COEFS.INDEX.B1.X3]
PE.COEFS.EFFECT.B1.X4  <- PE.COEFS[PE.COEFS.INDEX.B1.X4]
PE.COEFS.EFFECT.S1.X1  <- PE.COEFS[PE.COEFS.INDEX.S1.X1]
PE.COEFS.EFFECT.S1.X2  <- PE.COEFS[PE.COEFS.INDEX.S1.X2]
PE.COEFS.EFFECT.S1.X3  <- PE.COEFS[PE.COEFS.INDEX.S1.X3]
PE.COEFS.EFFECT.S1.X4  <- PE.COEFS[PE.COEFS.INDEX.S1.X4]
PE.COEFS.EFFECT.B2.X1X2  <- matrix(PE.COEFS[PE.COEFS.INDEX.B2.X1X2],nrow=20,byrow=TRUE)
PE.COEFS.EFFECT.B2.X1X3  <- matrix(PE.COEFS[PE.COEFS.INDEX.B2.X1X3],nrow=20,byrow=TRUE)
PE.COEFS.EFFECT.B2.X1X4  <- matrix(PE.COEFS[PE.COEFS.INDEX.B2.X1X4],nrow=20,byrow=TRUE)
PE.COEFS.EFFECT.B2.X2X3  <- matrix(PE.COEFS[PE.COEFS.INDEX.B2.X2X3],nrow=20,byrow=TRUE)
PE.COEFS.EFFECT.B2.X2X4  <- matrix(PE.COEFS[PE.COEFS.INDEX.B2.X2X4],nrow=20,byrow=TRUE)
PE.COEFS.EFFECT.B2.X3X4  <- matrix(PE.COEFS[PE.COEFS.INDEX.B2.X3X4],nrow=20,byrow=TRUE)
PE.COEFS.EFFECT.S2.X1X2  <- matrix(PE.COEFS[PE.COEFS.INDEX.S2.X1X2],nrow=20,byrow=TRUE)
PE.COEFS.EFFECT.S2.X1X3  <- matrix(PE.COEFS[PE.COEFS.INDEX.S2.X1X3],nrow=20,byrow=TRUE)
PE.COEFS.EFFECT.S2.X1X4  <- matrix(PE.COEFS[PE.COEFS.INDEX.S2.X1X4],nrow=20,byrow=TRUE)
PE.COEFS.EFFECT.S2.X2X3  <- matrix(PE.COEFS[PE.COEFS.INDEX.S2.X2X3],nrow=20,byrow=TRUE)
PE.COEFS.EFFECT.S2.X2X4  <- matrix(PE.COEFS[PE.COEFS.INDEX.S2.X2X4],nrow=20,byrow=TRUE)
PE.COEFS.EFFECT.S2.X3X4  <- matrix(PE.COEFS[PE.COEFS.INDEX.S2.X3X4],nrow=20,byrow=TRUE)

##Impose constraints-Average of terms per site should equal zero
#Starts with highest order interactions, constraining each dimension to sum to zero
#Lower order constraints account for constraints imposed on higher order terms
PE.COEFS.EFFECT.ADJ.B2.X1X2 <- PE.COEFS.EFFECT.B2.X1X2 - rowMeans(PE.COEFS.EFFECT.B2.X1X2) - rep(colMeans(PE.COEFS.EFFECT.B2.X1X2),each=20) + mean(PE.COEFS.EFFECT.B2.X1X2)
PE.COEFS.EFFECT.ADJ.B2.X1X3 <- PE.COEFS.EFFECT.B2.X1X3 - rowMeans(PE.COEFS.EFFECT.B2.X1X3) - rep(colMeans(PE.COEFS.EFFECT.B2.X1X3),each=20) + mean(PE.COEFS.EFFECT.B2.X1X3)
PE.COEFS.EFFECT.ADJ.B2.X1X4 <- PE.COEFS.EFFECT.B2.X1X4 - rowMeans(PE.COEFS.EFFECT.B2.X1X4) - rep(colMeans(PE.COEFS.EFFECT.B2.X1X4),each=20) + mean(PE.COEFS.EFFECT.B2.X1X4)
PE.COEFS.EFFECT.ADJ.B2.X2X3 <- PE.COEFS.EFFECT.B2.X2X3 - rowMeans(PE.COEFS.EFFECT.B2.X2X3) - rep(colMeans(PE.COEFS.EFFECT.B2.X2X3),each=20) + mean(PE.COEFS.EFFECT.B2.X2X3)
PE.COEFS.EFFECT.ADJ.B2.X2X4 <- PE.COEFS.EFFECT.B2.X2X4 - rowMeans(PE.COEFS.EFFECT.B2.X2X4) - rep(colMeans(PE.COEFS.EFFECT.B2.X2X4),each=20) + mean(PE.COEFS.EFFECT.B2.X2X4)
PE.COEFS.EFFECT.ADJ.B2.X3X4 <- PE.COEFS.EFFECT.B2.X3X4 - rowMeans(PE.COEFS.EFFECT.B2.X3X4) - rep(colMeans(PE.COEFS.EFFECT.B2.X3X4),each=20) + mean(PE.COEFS.EFFECT.B2.X3X4)
PE.COEFS.EFFECT.ADJ.B1.X1   <- PE.COEFS.EFFECT.B1.X1 - mean(PE.COEFS.EFFECT.B1.X1) + (colMeans(PE.COEFS.EFFECT.B2.X1X2) - mean(PE.COEFS.EFFECT.B2.X1X2)) + (colMeans(PE.COEFS.EFFECT.B2.X1X3) - mean(PE.COEFS.EFFECT.B2.X1X3)) + (colMeans(PE.COEFS.EFFECT.B2.X1X4) - mean(PE.COEFS.EFFECT.B2.X1X4))
PE.COEFS.EFFECT.ADJ.B1.X2   <- PE.COEFS.EFFECT.B1.X2 - mean(PE.COEFS.EFFECT.B1.X2) + (rowMeans(PE.COEFS.EFFECT.B2.X1X2) - mean(PE.COEFS.EFFECT.B2.X1X2)) + (colMeans(PE.COEFS.EFFECT.B2.X2X3) - mean(PE.COEFS.EFFECT.B2.X2X3)) + (colMeans(PE.COEFS.EFFECT.B2.X2X4) - mean(PE.COEFS.EFFECT.B2.X2X4))
PE.COEFS.EFFECT.ADJ.B1.X3   <- PE.COEFS.EFFECT.B1.X3 - mean(PE.COEFS.EFFECT.B1.X3) + (rowMeans(PE.COEFS.EFFECT.B2.X1X3) - mean(PE.COEFS.EFFECT.B2.X1X3)) + (rowMeans(PE.COEFS.EFFECT.B2.X2X3) - mean(PE.COEFS.EFFECT.B2.X2X3)) + (colMeans(PE.COEFS.EFFECT.B2.X3X4) - mean(PE.COEFS.EFFECT.B2.X3X4))
PE.COEFS.EFFECT.ADJ.B1.X4   <- PE.COEFS.EFFECT.B1.X4 - mean(PE.COEFS.EFFECT.B1.X4) + (rowMeans(PE.COEFS.EFFECT.B2.X1X4) - mean(PE.COEFS.EFFECT.B2.X1X4)) + (rowMeans(PE.COEFS.EFFECT.B2.X2X4) - mean(PE.COEFS.EFFECT.B2.X2X4)) + (rowMeans(PE.COEFS.EFFECT.B2.X3X4) - mean(PE.COEFS.EFFECT.B2.X3X4))
PE.COEFS.EFFECT.ADJ.B0      <- PE.COEFS.EFFECT.B0 + mean(PE.COEFS.EFFECT.B1.X1) + mean(PE.COEFS.EFFECT.B1.X2) + mean(PE.COEFS.EFFECT.B1.X3) + mean(PE.COEFS.EFFECT.B1.X4) + mean(PE.COEFS.EFFECT.B2.X1X2) + mean(PE.COEFS.EFFECT.B2.X1X3) + mean(PE.COEFS.EFFECT.B2.X1X4) + mean(PE.COEFS.EFFECT.B2.X2X3) + mean(PE.COEFS.EFFECT.B2.X2X4) + mean(PE.COEFS.EFFECT.B2.X3X4)

PE.COEFS.EFFECT.ADJ.S2.X1X2 <- PE.COEFS.EFFECT.S2.X1X2 - rowMeans(PE.COEFS.EFFECT.S2.X1X2) - rep(colMeans(PE.COEFS.EFFECT.S2.X1X2),each=20) + mean(PE.COEFS.EFFECT.S2.X1X2)
PE.COEFS.EFFECT.ADJ.S2.X1X3 <- PE.COEFS.EFFECT.S2.X1X3 - rowMeans(PE.COEFS.EFFECT.S2.X1X3) - rep(colMeans(PE.COEFS.EFFECT.S2.X1X3),each=20) + mean(PE.COEFS.EFFECT.S2.X1X3)
PE.COEFS.EFFECT.ADJ.S2.X1X4 <- PE.COEFS.EFFECT.S2.X1X4 - rowMeans(PE.COEFS.EFFECT.S2.X1X4) - rep(colMeans(PE.COEFS.EFFECT.S2.X1X4),each=20) + mean(PE.COEFS.EFFECT.S2.X1X4)
PE.COEFS.EFFECT.ADJ.S2.X2X3 <- PE.COEFS.EFFECT.S2.X2X3 - rowMeans(PE.COEFS.EFFECT.S2.X2X3) - rep(colMeans(PE.COEFS.EFFECT.S2.X2X3),each=20) + mean(PE.COEFS.EFFECT.S2.X2X3)
PE.COEFS.EFFECT.ADJ.S2.X2X4 <- PE.COEFS.EFFECT.S2.X2X4 - rowMeans(PE.COEFS.EFFECT.S2.X2X4) - rep(colMeans(PE.COEFS.EFFECT.S2.X2X4),each=20) + mean(PE.COEFS.EFFECT.S2.X2X4)
PE.COEFS.EFFECT.ADJ.S2.X3X4 <- PE.COEFS.EFFECT.S2.X3X4 - rowMeans(PE.COEFS.EFFECT.S2.X3X4) - rep(colMeans(PE.COEFS.EFFECT.S2.X3X4),each=20) + mean(PE.COEFS.EFFECT.S2.X3X4)
PE.COEFS.EFFECT.ADJ.S1.X1   <- PE.COEFS.EFFECT.S1.X1 - mean(PE.COEFS.EFFECT.S1.X1) + (colMeans(PE.COEFS.EFFECT.S2.X1X2) - mean(PE.COEFS.EFFECT.S2.X1X2)) + (colMeans(PE.COEFS.EFFECT.S2.X1X3) - mean(PE.COEFS.EFFECT.S2.X1X3)) + (colMeans(PE.COEFS.EFFECT.S2.X1X4) - mean(PE.COEFS.EFFECT.S2.X1X4))
PE.COEFS.EFFECT.ADJ.S1.X2   <- PE.COEFS.EFFECT.S1.X2 - mean(PE.COEFS.EFFECT.S1.X2) + (rowMeans(PE.COEFS.EFFECT.S2.X1X2) - mean(PE.COEFS.EFFECT.S2.X1X2)) + (colMeans(PE.COEFS.EFFECT.S2.X2X3) - mean(PE.COEFS.EFFECT.S2.X2X3)) + (colMeans(PE.COEFS.EFFECT.S2.X2X4) - mean(PE.COEFS.EFFECT.S2.X2X4))
PE.COEFS.EFFECT.ADJ.S1.X3   <- PE.COEFS.EFFECT.S1.X3 - mean(PE.COEFS.EFFECT.S1.X3) + (rowMeans(PE.COEFS.EFFECT.S2.X1X3) - mean(PE.COEFS.EFFECT.S2.X1X3)) + (rowMeans(PE.COEFS.EFFECT.S2.X2X3) - mean(PE.COEFS.EFFECT.S2.X2X3)) + (colMeans(PE.COEFS.EFFECT.S2.X3X4) - mean(PE.COEFS.EFFECT.S2.X3X4))
PE.COEFS.EFFECT.ADJ.S1.X4   <- PE.COEFS.EFFECT.S1.X4 - mean(PE.COEFS.EFFECT.S1.X4) + (rowMeans(PE.COEFS.EFFECT.S2.X1X4) - mean(PE.COEFS.EFFECT.S2.X1X4)) + (rowMeans(PE.COEFS.EFFECT.S2.X2X4) - mean(PE.COEFS.EFFECT.S2.X2X4)) + (rowMeans(PE.COEFS.EFFECT.S2.X3X4) - mean(PE.COEFS.EFFECT.S2.X3X4))
PE.COEFS.EFFECT.ADJ.S0      <- PE.COEFS.EFFECT.S0 + mean(PE.COEFS.EFFECT.S1.X1) + mean(PE.COEFS.EFFECT.S1.X2) + mean(PE.COEFS.EFFECT.S1.X3) + mean(PE.COEFS.EFFECT.S1.X4) + mean(PE.COEFS.EFFECT.S2.X1X2) + mean(PE.COEFS.EFFECT.S2.X1X3) + mean(PE.COEFS.EFFECT.S2.X1X4) + mean(PE.COEFS.EFFECT.S2.X2X3) + mean(PE.COEFS.EFFECT.S2.X2X4) + mean(PE.COEFS.EFFECT.S2.X3X4)

##Collect adjusted terms
PE.COEFS.ADJ <- c(PE.COEFS.EFFECT.ADJ.B1.X1,PE.COEFS.EFFECT.ADJ.B1.X2,PE.COEFS.EFFECT.ADJ.B1.X3,PE.COEFS.EFFECT.ADJ.B1.X4,
                  PE.COEFS.EFFECT.ADJ.S1.X1,PE.COEFS.EFFECT.ADJ.S1.X2,PE.COEFS.EFFECT.ADJ.S1.X3,PE.COEFS.EFFECT.ADJ.S1.X4,
                  c(t(PE.COEFS.EFFECT.ADJ.B2.X1X2)),c(t(PE.COEFS.EFFECT.ADJ.B2.X1X3)),c(t(PE.COEFS.EFFECT.ADJ.B2.X1X4)),c(t(PE.COEFS.EFFECT.ADJ.B2.X2X3)),c(t(PE.COEFS.EFFECT.ADJ.B2.X2X4)),c(t(PE.COEFS.EFFECT.ADJ.B2.X3X4)),
                  c(t(PE.COEFS.EFFECT.ADJ.S2.X1X2)),c(t(PE.COEFS.EFFECT.ADJ.S2.X1X3)),c(t(PE.COEFS.EFFECT.ADJ.S2.X1X4)),c(t(PE.COEFS.EFFECT.ADJ.S2.X2X3)),c(t(PE.COEFS.EFFECT.ADJ.S2.X2X4)),c(t(PE.COEFS.EFFECT.ADJ.S2.X3X4)));names(PE.COEFS.ADJ) <- names(PE.COEFS[-PE.COEFS.INDEX.S0])
#save(PE.COEFS.ADJ,file="PE.COEFS.ADJ.rda")
#save(PE.COEFS.EFFECT.ADJ.B0,file="PE.COEFS.EFFECT.ADJ.B0.rda")
#save(PE.COEFS.EFFECT.ADJ.S0,file="PE.COEFS.EFFECT.ADJ.S0.rda")

####TE Model
##Index site specific coefficients
TE.COEFS.INDEX.B1.X1     <- grep("^X1.$",            names(TE.COEFS))
TE.COEFS.INDEX.B1.X2     <- grep("^X2.$",            names(TE.COEFS))
TE.COEFS.INDEX.B1.X3     <- grep("^X3.$",            names(TE.COEFS))
TE.COEFS.INDEX.B1.X4     <- grep("^X4.$",            names(TE.COEFS))
TE.COEFS.INDEX.S1.X1     <- grep("^X1.:RE1$",        names(TE.COEFS))
TE.COEFS.INDEX.S1.X2     <- grep("^X2.:RE1$",        names(TE.COEFS))
TE.COEFS.INDEX.S1.X3     <- grep("^X3.:RE1$",        names(TE.COEFS))
TE.COEFS.INDEX.S1.X4     <- grep("^X4.:RE1$",        names(TE.COEFS))
TE.COEFS.INDEX.B2.X1X2   <- grep("^X1.:X2.$",        names(TE.COEFS))
TE.COEFS.INDEX.B2.X1X3   <- grep("^X1.:X3.$",        names(TE.COEFS))
TE.COEFS.INDEX.B2.X1X4   <- grep("^X1.:X4.$",        names(TE.COEFS))
TE.COEFS.INDEX.B2.X2X3   <- grep("^X2.:X3.$",        names(TE.COEFS))
TE.COEFS.INDEX.B2.X2X4   <- grep("^X2.:X4.$",        names(TE.COEFS))
TE.COEFS.INDEX.B2.X3X4   <- grep("^X3.:X4.$",        names(TE.COEFS))
TE.COEFS.INDEX.S2.X1X2   <- grep("^X1.:X2.:RE1$",    names(TE.COEFS))
TE.COEFS.INDEX.S2.X1X3   <- grep("^X1.:X3.:RE1$",    names(TE.COEFS))
TE.COEFS.INDEX.S2.X1X4   <- grep("^X1.:X4.:RE1$",    names(TE.COEFS))
TE.COEFS.INDEX.S2.X2X3   <- grep("^X2.:X3.:RE1$",    names(TE.COEFS))
TE.COEFS.INDEX.S2.X2X4   <- grep("^X2.:X4.:RE1$",    names(TE.COEFS))
TE.COEFS.INDEX.S2.X3X4   <- grep("^X3.:X4.:RE1$",    names(TE.COEFS))
TE.COEFS.INDEX.B3.X1X2X3 <- grep("^X1.:X2.:X3.$",    names(TE.COEFS))
TE.COEFS.INDEX.B3.X1X2X4 <- grep("^X1.:X2.:X4.$",    names(TE.COEFS))
TE.COEFS.INDEX.B3.X1X3X4 <- grep("^X1.:X3.:X4.$",    names(TE.COEFS))
TE.COEFS.INDEX.B3.X2X3X4 <- grep("^X2.:X3.:X4.$",    names(TE.COEFS))
TE.COEFS.INDEX.S3.X1X2X3 <- grep("^X1.:X2.:X3.:RE1$",names(TE.COEFS))
TE.COEFS.INDEX.S3.X1X2X4 <- grep("^X1.:X2.:X4.:RE1$",names(TE.COEFS))
TE.COEFS.INDEX.S3.X1X3X4 <- grep("^X1.:X3.:X4.:RE1$",names(TE.COEFS))
TE.COEFS.INDEX.S3.X2X3X4 <- grep("^X2.:X3.:X4.:RE1$",names(TE.COEFS))

##Index order specific coefficients
TE.COEFS.INDEX.O1  <- c(TE.COEFS.INDEX.B1,TE.COEFS.INDEX.S1)
TE.COEFS.INDEX.O2  <- c(TE.COEFS.INDEX.B2,TE.COEFS.INDEX.S2)
TE.COEFS.INDEX.O3  <- c(TE.COEFS.INDEX.B3,TE.COEFS.INDEX.S3)
TE.COEFS.INDEX.B   <- c(TE.COEFS.INDEX.B1,TE.COEFS.INDEX.B2,TE.COEFS.INDEX.B3)
TE.COEFS.INDEX.S   <- c(TE.COEFS.INDEX.S1,TE.COEFS.INDEX.S2,TE.COEFS.INDEX.S3)

##Extract coefficient effects
TE.COEFS.EFFECT.B0         <- 0
TE.COEFS.EFFECT.S0         <- TE.COEFS[TE.COEFS.INDEX.S0]
TE.COEFS.EFFECT.B1.X1      <- TE.COEFS[TE.COEFS.INDEX.B1.X1]
TE.COEFS.EFFECT.B1.X2      <- TE.COEFS[TE.COEFS.INDEX.B1.X2]
TE.COEFS.EFFECT.B1.X3      <- TE.COEFS[TE.COEFS.INDEX.B1.X3]
TE.COEFS.EFFECT.B1.X4      <- TE.COEFS[TE.COEFS.INDEX.B1.X4]
TE.COEFS.EFFECT.S1.X1      <- TE.COEFS[TE.COEFS.INDEX.S1.X1]
TE.COEFS.EFFECT.S1.X2      <- TE.COEFS[TE.COEFS.INDEX.S1.X2]
TE.COEFS.EFFECT.S1.X3      <- TE.COEFS[TE.COEFS.INDEX.S1.X3]
TE.COEFS.EFFECT.S1.X4      <- TE.COEFS[TE.COEFS.INDEX.S1.X4]
TE.COEFS.EFFECT.B2.X1X2    <- array(TE.COEFS[TE.COEFS.INDEX.B2.X1X2],  dim=c(20,20),   dimnames=list(names(TE.COEFS.EFFECT.B1.X1),names(TE.COEFS.EFFECT.B1.X2)))
TE.COEFS.EFFECT.B2.X1X3    <- array(TE.COEFS[TE.COEFS.INDEX.B2.X1X3],  dim=c(20,20),   dimnames=list(names(TE.COEFS.EFFECT.B1.X1),names(TE.COEFS.EFFECT.B1.X3)))
TE.COEFS.EFFECT.B2.X1X4    <- array(TE.COEFS[TE.COEFS.INDEX.B2.X1X4],  dim=c(20,20),   dimnames=list(names(TE.COEFS.EFFECT.B1.X1),names(TE.COEFS.EFFECT.B1.X4)))
TE.COEFS.EFFECT.B2.X2X3    <- array(TE.COEFS[TE.COEFS.INDEX.B2.X2X3],  dim=c(20,20),   dimnames=list(names(TE.COEFS.EFFECT.B1.X2),names(TE.COEFS.EFFECT.B1.X3)))
TE.COEFS.EFFECT.B2.X2X4    <- array(TE.COEFS[TE.COEFS.INDEX.B2.X2X4],  dim=c(20,20),   dimnames=list(names(TE.COEFS.EFFECT.B1.X2),names(TE.COEFS.EFFECT.B1.X4)))
TE.COEFS.EFFECT.B2.X3X4    <- array(TE.COEFS[TE.COEFS.INDEX.B2.X3X4],  dim=c(20,20),   dimnames=list(names(TE.COEFS.EFFECT.B1.X3),names(TE.COEFS.EFFECT.B1.X4)))
TE.COEFS.EFFECT.S2.X1X2    <- array(TE.COEFS[TE.COEFS.INDEX.S2.X1X2],  dim=c(20,20),   dimnames=list(names(TE.COEFS.EFFECT.S1.X1),names(TE.COEFS.EFFECT.S1.X2)))
TE.COEFS.EFFECT.S2.X1X3    <- array(TE.COEFS[TE.COEFS.INDEX.S2.X1X3],  dim=c(20,20),   dimnames=list(names(TE.COEFS.EFFECT.S1.X1),names(TE.COEFS.EFFECT.S1.X3)))
TE.COEFS.EFFECT.S2.X1X4    <- array(TE.COEFS[TE.COEFS.INDEX.S2.X1X4],  dim=c(20,20),   dimnames=list(names(TE.COEFS.EFFECT.S1.X1),names(TE.COEFS.EFFECT.S1.X4)))
TE.COEFS.EFFECT.S2.X2X3    <- array(TE.COEFS[TE.COEFS.INDEX.S2.X2X3],  dim=c(20,20),   dimnames=list(names(TE.COEFS.EFFECT.S1.X2),names(TE.COEFS.EFFECT.S1.X3)))
TE.COEFS.EFFECT.S2.X2X4    <- array(TE.COEFS[TE.COEFS.INDEX.S2.X2X4],  dim=c(20,20),   dimnames=list(names(TE.COEFS.EFFECT.S1.X2),names(TE.COEFS.EFFECT.S1.X4)))
TE.COEFS.EFFECT.S2.X3X4    <- array(TE.COEFS[TE.COEFS.INDEX.S2.X3X4],  dim=c(20,20),   dimnames=list(names(TE.COEFS.EFFECT.S1.X3),names(TE.COEFS.EFFECT.S1.X4)))
TE.COEFS.EFFECT.B3.X1X2X3  <- array(TE.COEFS[TE.COEFS.INDEX.B3.X1X2X3],dim=c(20,20,20),dimnames=list(names(TE.COEFS.EFFECT.B1.X1),names(TE.COEFS.EFFECT.B1.X2),names(TE.COEFS.EFFECT.B1.X3)))
TE.COEFS.EFFECT.B3.X1X2X4  <- array(TE.COEFS[TE.COEFS.INDEX.B3.X1X2X4],dim=c(20,20,20),dimnames=list(names(TE.COEFS.EFFECT.B1.X1),names(TE.COEFS.EFFECT.B1.X2),names(TE.COEFS.EFFECT.B1.X4)))
TE.COEFS.EFFECT.B3.X1X3X4  <- array(TE.COEFS[TE.COEFS.INDEX.B3.X1X3X4],dim=c(20,20,20),dimnames=list(names(TE.COEFS.EFFECT.B1.X1),names(TE.COEFS.EFFECT.B1.X3),names(TE.COEFS.EFFECT.B1.X4)))
TE.COEFS.EFFECT.B3.X2X3X4  <- array(TE.COEFS[TE.COEFS.INDEX.B3.X2X3X4],dim=c(20,20,20),dimnames=list(names(TE.COEFS.EFFECT.B1.X2),names(TE.COEFS.EFFECT.B1.X3),names(TE.COEFS.EFFECT.B1.X4)))
TE.COEFS.EFFECT.S3.X1X2X3  <- array(TE.COEFS[TE.COEFS.INDEX.S3.X1X2X3],dim=c(20,20,20),dimnames=list(names(TE.COEFS.EFFECT.S1.X1),names(TE.COEFS.EFFECT.S1.X2),names(TE.COEFS.EFFECT.S1.X3)))
TE.COEFS.EFFECT.S3.X1X2X4  <- array(TE.COEFS[TE.COEFS.INDEX.S3.X1X2X4],dim=c(20,20,20),dimnames=list(names(TE.COEFS.EFFECT.S1.X1),names(TE.COEFS.EFFECT.S1.X2),names(TE.COEFS.EFFECT.S1.X4)))
TE.COEFS.EFFECT.S3.X1X3X4  <- array(TE.COEFS[TE.COEFS.INDEX.S3.X1X3X4],dim=c(20,20,20),dimnames=list(names(TE.COEFS.EFFECT.S1.X1),names(TE.COEFS.EFFECT.S1.X3),names(TE.COEFS.EFFECT.S1.X4)))
TE.COEFS.EFFECT.S3.X2X3X4  <- array(TE.COEFS[TE.COEFS.INDEX.S3.X2X3X4],dim=c(20,20,20),dimnames=list(names(TE.COEFS.EFFECT.S1.X2),names(TE.COEFS.EFFECT.S1.X3),names(TE.COEFS.EFFECT.S1.X4)))

##Impose constraints-Average of terms per site should equal zero
#Starts with highest order interactions, constraining each dimension to sum to zero
#Lower order constraints account for constraints imposed on higher order terms
TE.COEFS.EFFECT.ADJ.B3.X1X2X3 <- TE.COEFS.EFFECT.B3.X1X2X3 - array(rep(apply(TE.COEFS.EFFECT.B3.X1X2X3,c(2,3),mean),each=20),dim=c(20,20,20)) - array(apply(apply(TE.COEFS.EFFECT.B3.X1X2X3,c(1,3),mean),2,rep,times=20),dim=c(20,20,20)) - array(rep(apply(TE.COEFS.EFFECT.B3.X1X2X3,c(1,2),mean),times=20),dim=c(20,20,20)) + array(rep(apply(TE.COEFS.EFFECT.B3.X1X2X3,1,mean),times=400),dim=c(20,20,20)) + array(rep(rep(apply(TE.COEFS.EFFECT.B3.X1X2X3,2,mean),20),each=20),dim=c(20,20,20)) + array(rep(apply(TE.COEFS.EFFECT.B3.X1X2X3,3,mean),each=400),dim=c(20,20,20)) - mean(TE.COEFS.EFFECT.B3.X1X2X3) 
TE.COEFS.EFFECT.ADJ.B3.X1X2X4 <- TE.COEFS.EFFECT.B3.X1X2X4 - array(rep(apply(TE.COEFS.EFFECT.B3.X1X2X4,c(2,3),mean),each=20),dim=c(20,20,20)) - array(apply(apply(TE.COEFS.EFFECT.B3.X1X2X4,c(1,3),mean),2,rep,times=20),dim=c(20,20,20)) - array(rep(apply(TE.COEFS.EFFECT.B3.X1X2X4,c(1,2),mean),times=20),dim=c(20,20,20)) + array(rep(apply(TE.COEFS.EFFECT.B3.X1X2X4,1,mean),times=400),dim=c(20,20,20)) + array(rep(rep(apply(TE.COEFS.EFFECT.B3.X1X2X4,2,mean),20),each=20),dim=c(20,20,20)) + array(rep(apply(TE.COEFS.EFFECT.B3.X1X2X4,3,mean),each=400),dim=c(20,20,20)) - mean(TE.COEFS.EFFECT.B3.X1X2X4) 
TE.COEFS.EFFECT.ADJ.B3.X1X3X4 <- TE.COEFS.EFFECT.B3.X1X3X4 - array(rep(apply(TE.COEFS.EFFECT.B3.X1X3X4,c(2,3),mean),each=20),dim=c(20,20,20)) - array(apply(apply(TE.COEFS.EFFECT.B3.X1X3X4,c(1,3),mean),2,rep,times=20),dim=c(20,20,20)) - array(rep(apply(TE.COEFS.EFFECT.B3.X1X3X4,c(1,2),mean),times=20),dim=c(20,20,20)) + array(rep(apply(TE.COEFS.EFFECT.B3.X1X3X4,1,mean),times=400),dim=c(20,20,20)) + array(rep(rep(apply(TE.COEFS.EFFECT.B3.X1X3X4,2,mean),20),each=20),dim=c(20,20,20)) + array(rep(apply(TE.COEFS.EFFECT.B3.X1X3X4,3,mean),each=400),dim=c(20,20,20)) - mean(TE.COEFS.EFFECT.B3.X1X3X4) 
TE.COEFS.EFFECT.ADJ.B3.X2X3X4 <- TE.COEFS.EFFECT.B3.X2X3X4 - array(rep(apply(TE.COEFS.EFFECT.B3.X2X3X4,c(2,3),mean),each=20),dim=c(20,20,20)) - array(apply(apply(TE.COEFS.EFFECT.B3.X2X3X4,c(1,3),mean),2,rep,times=20),dim=c(20,20,20)) - array(rep(apply(TE.COEFS.EFFECT.B3.X2X3X4,c(1,2),mean),times=20),dim=c(20,20,20)) + array(rep(apply(TE.COEFS.EFFECT.B3.X2X3X4,1,mean),times=400),dim=c(20,20,20)) + array(rep(rep(apply(TE.COEFS.EFFECT.B3.X2X3X4,2,mean),20),each=20),dim=c(20,20,20)) + array(rep(apply(TE.COEFS.EFFECT.B3.X2X3X4,3,mean),each=400),dim=c(20,20,20)) - mean(TE.COEFS.EFFECT.B3.X2X3X4) 
TE.COEFS.EFFECT.ADJ.B2.X1X2   <- TE.COEFS.EFFECT.B2.X1X2 - array(rep(apply(TE.COEFS.EFFECT.B2.X1X2,1,mean),times=20),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.B2.X1X2,2,mean),each=20),dim=c(20,20)) + mean(TE.COEFS.EFFECT.B2.X1X2) + array(apply(TE.COEFS.EFFECT.B3.X1X2X3,c(1,2),mean),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.B3.X1X2X3,1,mean),times=20),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.B3.X1X2X3,2,mean),each=20),dim=c(20,20)) + mean(TE.COEFS.EFFECT.B3.X1X2X3) + array(apply(TE.COEFS.EFFECT.B3.X1X2X4,c(1,2),mean),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.B3.X1X2X4,1,mean),times=20),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.B3.X1X2X4,2,mean),each=20),dim=c(20,20)) + mean(TE.COEFS.EFFECT.B3.X1X2X4)
TE.COEFS.EFFECT.ADJ.B2.X1X3   <- TE.COEFS.EFFECT.B2.X1X3 - array(rep(apply(TE.COEFS.EFFECT.B2.X1X3,1,mean),times=20),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.B2.X1X3,2,mean),each=20),dim=c(20,20)) + mean(TE.COEFS.EFFECT.B2.X1X3) + array(apply(TE.COEFS.EFFECT.B3.X1X2X3,c(1,3),mean),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.B3.X1X2X3,1,mean),times=20),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.B3.X1X2X3,3,mean),each=20),dim=c(20,20)) + mean(TE.COEFS.EFFECT.B3.X1X2X3) + array(apply(TE.COEFS.EFFECT.B3.X1X3X4,c(1,2),mean),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.B3.X1X3X4,1,mean),times=20),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.B3.X1X3X4,2,mean),each=20),dim=c(20,20)) + mean(TE.COEFS.EFFECT.B3.X1X3X4)
TE.COEFS.EFFECT.ADJ.B2.X1X4   <- TE.COEFS.EFFECT.B2.X1X4 - array(rep(apply(TE.COEFS.EFFECT.B2.X1X4,1,mean),times=20),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.B2.X1X4,2,mean),each=20),dim=c(20,20)) + mean(TE.COEFS.EFFECT.B2.X1X4) + array(apply(TE.COEFS.EFFECT.B3.X1X2X4,c(1,3),mean),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.B3.X1X2X4,1,mean),times=20),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.B3.X1X2X4,3,mean),each=20),dim=c(20,20)) + mean(TE.COEFS.EFFECT.B3.X1X2X4) + array(apply(TE.COEFS.EFFECT.B3.X1X3X4,c(1,3),mean),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.B3.X1X3X4,1,mean),times=20),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.B3.X1X3X4,3,mean),each=20),dim=c(20,20)) + mean(TE.COEFS.EFFECT.B3.X1X3X4)
TE.COEFS.EFFECT.ADJ.B2.X2X3   <- TE.COEFS.EFFECT.B2.X2X3 - array(rep(apply(TE.COEFS.EFFECT.B2.X2X3,1,mean),times=20),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.B2.X2X3,2,mean),each=20),dim=c(20,20)) + mean(TE.COEFS.EFFECT.B2.X2X3) + array(apply(TE.COEFS.EFFECT.B3.X1X2X3,c(2,3),mean),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.B3.X1X2X3,2,mean),times=20),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.B3.X1X2X3,3,mean),each=20),dim=c(20,20)) + mean(TE.COEFS.EFFECT.B3.X1X2X3) + array(apply(TE.COEFS.EFFECT.B3.X2X3X4,c(1,2),mean),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.B3.X2X3X4,1,mean),times=20),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.B3.X2X3X4,2,mean),each=20),dim=c(20,20)) + mean(TE.COEFS.EFFECT.B3.X2X3X4)
TE.COEFS.EFFECT.ADJ.B2.X2X4   <- TE.COEFS.EFFECT.B2.X2X4 - array(rep(apply(TE.COEFS.EFFECT.B2.X2X4,1,mean),times=20),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.B2.X2X4,2,mean),each=20),dim=c(20,20)) + mean(TE.COEFS.EFFECT.B2.X2X4) + array(apply(TE.COEFS.EFFECT.B3.X1X2X4,c(2,3),mean),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.B3.X1X2X4,2,mean),times=20),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.B3.X1X2X4,3,mean),each=20),dim=c(20,20)) + mean(TE.COEFS.EFFECT.B3.X1X2X4) + array(apply(TE.COEFS.EFFECT.B3.X2X3X4,c(1,3),mean),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.B3.X2X3X4,1,mean),times=20),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.B3.X2X3X4,3,mean),each=20),dim=c(20,20)) + mean(TE.COEFS.EFFECT.B3.X2X3X4)
TE.COEFS.EFFECT.ADJ.B2.X3X4   <- TE.COEFS.EFFECT.B2.X3X4 - array(rep(apply(TE.COEFS.EFFECT.B2.X3X4,1,mean),times=20),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.B2.X3X4,2,mean),each=20),dim=c(20,20)) + mean(TE.COEFS.EFFECT.B2.X3X4) + array(apply(TE.COEFS.EFFECT.B3.X1X3X4,c(2,3),mean),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.B3.X1X3X4,2,mean),times=20),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.B3.X1X3X4,3,mean),each=20),dim=c(20,20)) + mean(TE.COEFS.EFFECT.B3.X1X3X4) + array(apply(TE.COEFS.EFFECT.B3.X2X3X4,c(2,3),mean),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.B3.X2X3X4,2,mean),times=20),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.B3.X2X3X4,3,mean),each=20),dim=c(20,20)) + mean(TE.COEFS.EFFECT.B3.X2X3X4)
TE.COEFS.EFFECT.ADJ.B1.X1     <- TE.COEFS.EFFECT.B1.X1 - mean(TE.COEFS.EFFECT.B1.X1) + apply(TE.COEFS.EFFECT.B2.X1X2,1,mean) - mean(TE.COEFS.EFFECT.B2.X1X2) + apply(TE.COEFS.EFFECT.B2.X1X3,1,mean) - mean(TE.COEFS.EFFECT.B2.X1X3) + apply(TE.COEFS.EFFECT.B2.X1X4,1,mean) - mean(TE.COEFS.EFFECT.B2.X1X4) + apply(TE.COEFS.EFFECT.B3.X1X2X3,1,mean) - mean(TE.COEFS.EFFECT.B3.X1X2X3) + apply(TE.COEFS.EFFECT.B3.X1X2X4,1,mean) - mean(TE.COEFS.EFFECT.B3.X1X2X4) + apply(TE.COEFS.EFFECT.B3.X1X3X4,1,mean) - mean(TE.COEFS.EFFECT.B3.X1X3X4) 
TE.COEFS.EFFECT.ADJ.B1.X2     <- TE.COEFS.EFFECT.B1.X2 - mean(TE.COEFS.EFFECT.B1.X2) + apply(TE.COEFS.EFFECT.B2.X1X2,2,mean) - mean(TE.COEFS.EFFECT.B2.X1X2) + apply(TE.COEFS.EFFECT.B2.X2X3,1,mean) - mean(TE.COEFS.EFFECT.B2.X2X3) + apply(TE.COEFS.EFFECT.B2.X2X4,1,mean) - mean(TE.COEFS.EFFECT.B2.X2X4) + apply(TE.COEFS.EFFECT.B3.X1X2X3,2,mean) - mean(TE.COEFS.EFFECT.B3.X1X2X3) + apply(TE.COEFS.EFFECT.B3.X1X2X4,2,mean) - mean(TE.COEFS.EFFECT.B3.X1X2X4) + apply(TE.COEFS.EFFECT.B3.X2X3X4,1,mean) - mean(TE.COEFS.EFFECT.B3.X2X3X4) 
TE.COEFS.EFFECT.ADJ.B1.X3     <- TE.COEFS.EFFECT.B1.X3 - mean(TE.COEFS.EFFECT.B1.X3) + apply(TE.COEFS.EFFECT.B2.X1X3,2,mean) - mean(TE.COEFS.EFFECT.B2.X1X3) + apply(TE.COEFS.EFFECT.B2.X2X3,2,mean) - mean(TE.COEFS.EFFECT.B2.X2X3) + apply(TE.COEFS.EFFECT.B2.X3X4,1,mean) - mean(TE.COEFS.EFFECT.B2.X3X4) + apply(TE.COEFS.EFFECT.B3.X1X2X3,3,mean) - mean(TE.COEFS.EFFECT.B3.X1X2X3) + apply(TE.COEFS.EFFECT.B3.X1X3X4,2,mean) - mean(TE.COEFS.EFFECT.B3.X1X3X4) + apply(TE.COEFS.EFFECT.B3.X2X3X4,2,mean) - mean(TE.COEFS.EFFECT.B3.X2X3X4) 
TE.COEFS.EFFECT.ADJ.B1.X4     <- TE.COEFS.EFFECT.B1.X4 - mean(TE.COEFS.EFFECT.B1.X4) + apply(TE.COEFS.EFFECT.B2.X1X4,2,mean) - mean(TE.COEFS.EFFECT.B2.X1X4) + apply(TE.COEFS.EFFECT.B2.X2X4,2,mean) - mean(TE.COEFS.EFFECT.B2.X2X4) + apply(TE.COEFS.EFFECT.B2.X3X4,2,mean) - mean(TE.COEFS.EFFECT.B2.X3X4) + apply(TE.COEFS.EFFECT.B3.X1X2X4,3,mean) - mean(TE.COEFS.EFFECT.B3.X1X2X4) + apply(TE.COEFS.EFFECT.B3.X1X3X4,3,mean) - mean(TE.COEFS.EFFECT.B3.X1X3X4) + apply(TE.COEFS.EFFECT.B3.X2X3X4,3,mean) - mean(TE.COEFS.EFFECT.B3.X2X3X4) 
TE.COEFS.EFFECT.ADJ.B0        <- TE.COEFS.EFFECT.B0 + mean(TE.COEFS.EFFECT.B1.X1) + mean(TE.COEFS.EFFECT.B1.X2) + mean(TE.COEFS.EFFECT.B1.X3) + mean(TE.COEFS.EFFECT.B1.X4) + mean(TE.COEFS.EFFECT.B2.X1X2) + mean(TE.COEFS.EFFECT.B2.X1X3) + mean(TE.COEFS.EFFECT.B2.X1X4) + mean(TE.COEFS.EFFECT.B2.X2X3) + mean(TE.COEFS.EFFECT.B2.X2X4) + mean(TE.COEFS.EFFECT.B2.X3X4) + mean(TE.COEFS.EFFECT.B3.X1X2X3) + mean(TE.COEFS.EFFECT.B3.X1X2X4) + mean(TE.COEFS.EFFECT.B3.X1X3X4) + mean(TE.COEFS.EFFECT.B3.X2X3X4)

TE.COEFS.EFFECT.ADJ.S3.X1X2X3 <- TE.COEFS.EFFECT.S3.X1X2X3 - array(rep(apply(TE.COEFS.EFFECT.S3.X1X2X3,c(2,3),mean),each=20),dim=c(20,20,20)) - array(apply(apply(TE.COEFS.EFFECT.S3.X1X2X3,c(1,3),mean),2,rep,times=20),dim=c(20,20,20)) - array(rep(apply(TE.COEFS.EFFECT.S3.X1X2X3,c(1,2),mean),times=20),dim=c(20,20,20)) + array(rep(apply(TE.COEFS.EFFECT.S3.X1X2X3,1,mean),times=400),dim=c(20,20,20)) + array(rep(rep(apply(TE.COEFS.EFFECT.S3.X1X2X3,2,mean),20),each=20),dim=c(20,20,20)) + array(rep(apply(TE.COEFS.EFFECT.S3.X1X2X3,3,mean),each=400),dim=c(20,20,20)) - mean(TE.COEFS.EFFECT.S3.X1X2X3) 
TE.COEFS.EFFECT.ADJ.S3.X1X2X4 <- TE.COEFS.EFFECT.S3.X1X2X4 - array(rep(apply(TE.COEFS.EFFECT.S3.X1X2X4,c(2,3),mean),each=20),dim=c(20,20,20)) - array(apply(apply(TE.COEFS.EFFECT.S3.X1X2X4,c(1,3),mean),2,rep,times=20),dim=c(20,20,20)) - array(rep(apply(TE.COEFS.EFFECT.S3.X1X2X4,c(1,2),mean),times=20),dim=c(20,20,20)) + array(rep(apply(TE.COEFS.EFFECT.S3.X1X2X4,1,mean),times=400),dim=c(20,20,20)) + array(rep(rep(apply(TE.COEFS.EFFECT.S3.X1X2X4,2,mean),20),each=20),dim=c(20,20,20)) + array(rep(apply(TE.COEFS.EFFECT.S3.X1X2X4,3,mean),each=400),dim=c(20,20,20)) - mean(TE.COEFS.EFFECT.S3.X1X2X4) 
TE.COEFS.EFFECT.ADJ.S3.X1X3X4 <- TE.COEFS.EFFECT.S3.X1X3X4 - array(rep(apply(TE.COEFS.EFFECT.S3.X1X3X4,c(2,3),mean),each=20),dim=c(20,20,20)) - array(apply(apply(TE.COEFS.EFFECT.S3.X1X3X4,c(1,3),mean),2,rep,times=20),dim=c(20,20,20)) - array(rep(apply(TE.COEFS.EFFECT.S3.X1X3X4,c(1,2),mean),times=20),dim=c(20,20,20)) + array(rep(apply(TE.COEFS.EFFECT.S3.X1X3X4,1,mean),times=400),dim=c(20,20,20)) + array(rep(rep(apply(TE.COEFS.EFFECT.S3.X1X3X4,2,mean),20),each=20),dim=c(20,20,20)) + array(rep(apply(TE.COEFS.EFFECT.S3.X1X3X4,3,mean),each=400),dim=c(20,20,20)) - mean(TE.COEFS.EFFECT.S3.X1X3X4) 
TE.COEFS.EFFECT.ADJ.S3.X2X3X4 <- TE.COEFS.EFFECT.S3.X2X3X4 - array(rep(apply(TE.COEFS.EFFECT.S3.X2X3X4,c(2,3),mean),each=20),dim=c(20,20,20)) - array(apply(apply(TE.COEFS.EFFECT.S3.X2X3X4,c(1,3),mean),2,rep,times=20),dim=c(20,20,20)) - array(rep(apply(TE.COEFS.EFFECT.S3.X2X3X4,c(1,2),mean),times=20),dim=c(20,20,20)) + array(rep(apply(TE.COEFS.EFFECT.S3.X2X3X4,1,mean),times=400),dim=c(20,20,20)) + array(rep(rep(apply(TE.COEFS.EFFECT.S3.X2X3X4,2,mean),20),each=20),dim=c(20,20,20)) + array(rep(apply(TE.COEFS.EFFECT.S3.X2X3X4,3,mean),each=400),dim=c(20,20,20)) - mean(TE.COEFS.EFFECT.S3.X2X3X4) 
TE.COEFS.EFFECT.ADJ.S2.X1X2   <- TE.COEFS.EFFECT.S2.X1X2 - array(rep(apply(TE.COEFS.EFFECT.S2.X1X2,1,mean),times=20),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.S2.X1X2,2,mean),each=20),dim=c(20,20)) + mean(TE.COEFS.EFFECT.S2.X1X2) + array(apply(TE.COEFS.EFFECT.S3.X1X2X3,c(1,2),mean),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.S3.X1X2X3,1,mean),times=20),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.S3.X1X2X3,2,mean),each=20),dim=c(20,20)) + mean(TE.COEFS.EFFECT.S3.X1X2X3) + array(apply(TE.COEFS.EFFECT.S3.X1X2X4,c(1,2),mean),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.S3.X1X2X4,1,mean),times=20),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.S3.X1X2X4,2,mean),each=20),dim=c(20,20)) + mean(TE.COEFS.EFFECT.S3.X1X2X4)
TE.COEFS.EFFECT.ADJ.S2.X1X3   <- TE.COEFS.EFFECT.S2.X1X3 - array(rep(apply(TE.COEFS.EFFECT.S2.X1X3,1,mean),times=20),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.S2.X1X3,2,mean),each=20),dim=c(20,20)) + mean(TE.COEFS.EFFECT.S2.X1X3) + array(apply(TE.COEFS.EFFECT.S3.X1X2X3,c(1,3),mean),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.S3.X1X2X3,1,mean),times=20),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.S3.X1X2X3,3,mean),each=20),dim=c(20,20)) + mean(TE.COEFS.EFFECT.S3.X1X2X3) + array(apply(TE.COEFS.EFFECT.S3.X1X3X4,c(1,2),mean),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.S3.X1X3X4,1,mean),times=20),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.S3.X1X3X4,2,mean),each=20),dim=c(20,20)) + mean(TE.COEFS.EFFECT.S3.X1X3X4)
TE.COEFS.EFFECT.ADJ.S2.X1X4   <- TE.COEFS.EFFECT.S2.X1X4 - array(rep(apply(TE.COEFS.EFFECT.S2.X1X4,1,mean),times=20),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.S2.X1X4,2,mean),each=20),dim=c(20,20)) + mean(TE.COEFS.EFFECT.S2.X1X4) + array(apply(TE.COEFS.EFFECT.S3.X1X2X4,c(1,3),mean),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.S3.X1X2X4,1,mean),times=20),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.S3.X1X2X4,3,mean),each=20),dim=c(20,20)) + mean(TE.COEFS.EFFECT.S3.X1X2X4) + array(apply(TE.COEFS.EFFECT.S3.X1X3X4,c(1,3),mean),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.S3.X1X3X4,1,mean),times=20),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.S3.X1X3X4,3,mean),each=20),dim=c(20,20)) + mean(TE.COEFS.EFFECT.S3.X1X3X4)
TE.COEFS.EFFECT.ADJ.S2.X2X3   <- TE.COEFS.EFFECT.S2.X2X3 - array(rep(apply(TE.COEFS.EFFECT.S2.X2X3,1,mean),times=20),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.S2.X2X3,2,mean),each=20),dim=c(20,20)) + mean(TE.COEFS.EFFECT.S2.X2X3) + array(apply(TE.COEFS.EFFECT.S3.X1X2X3,c(2,3),mean),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.S3.X1X2X3,2,mean),times=20),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.S3.X1X2X3,3,mean),each=20),dim=c(20,20)) + mean(TE.COEFS.EFFECT.S3.X1X2X3) + array(apply(TE.COEFS.EFFECT.S3.X2X3X4,c(1,2),mean),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.S3.X2X3X4,1,mean),times=20),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.S3.X2X3X4,2,mean),each=20),dim=c(20,20)) + mean(TE.COEFS.EFFECT.S3.X2X3X4)
TE.COEFS.EFFECT.ADJ.S2.X2X4   <- TE.COEFS.EFFECT.S2.X2X4 - array(rep(apply(TE.COEFS.EFFECT.S2.X2X4,1,mean),times=20),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.S2.X2X4,2,mean),each=20),dim=c(20,20)) + mean(TE.COEFS.EFFECT.S2.X2X4) + array(apply(TE.COEFS.EFFECT.S3.X1X2X4,c(2,3),mean),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.S3.X1X2X4,2,mean),times=20),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.S3.X1X2X4,3,mean),each=20),dim=c(20,20)) + mean(TE.COEFS.EFFECT.S3.X1X2X4) + array(apply(TE.COEFS.EFFECT.S3.X2X3X4,c(1,3),mean),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.S3.X2X3X4,1,mean),times=20),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.S3.X2X3X4,3,mean),each=20),dim=c(20,20)) + mean(TE.COEFS.EFFECT.S3.X2X3X4)
TE.COEFS.EFFECT.ADJ.S2.X3X4   <- TE.COEFS.EFFECT.S2.X3X4 - array(rep(apply(TE.COEFS.EFFECT.S2.X3X4,1,mean),times=20),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.S2.X3X4,2,mean),each=20),dim=c(20,20)) + mean(TE.COEFS.EFFECT.S2.X3X4) + array(apply(TE.COEFS.EFFECT.S3.X1X3X4,c(2,3),mean),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.S3.X1X3X4,2,mean),times=20),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.S3.X1X3X4,3,mean),each=20),dim=c(20,20)) + mean(TE.COEFS.EFFECT.S3.X1X3X4) + array(apply(TE.COEFS.EFFECT.S3.X2X3X4,c(2,3),mean),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.S3.X2X3X4,2,mean),times=20),dim=c(20,20)) - array(rep(apply(TE.COEFS.EFFECT.S3.X2X3X4,3,mean),each=20),dim=c(20,20)) + mean(TE.COEFS.EFFECT.S3.X2X3X4)
TE.COEFS.EFFECT.ADJ.S1.X1     <- TE.COEFS.EFFECT.S1.X1 - mean(TE.COEFS.EFFECT.S1.X1) + apply(TE.COEFS.EFFECT.S2.X1X2,1,mean) - mean(TE.COEFS.EFFECT.S2.X1X2) + apply(TE.COEFS.EFFECT.S2.X1X3,1,mean) - mean(TE.COEFS.EFFECT.S2.X1X3) + apply(TE.COEFS.EFFECT.S2.X1X4,1,mean) - mean(TE.COEFS.EFFECT.S2.X1X4) + apply(TE.COEFS.EFFECT.S3.X1X2X3,1,mean) - mean(TE.COEFS.EFFECT.S3.X1X2X3) + apply(TE.COEFS.EFFECT.S3.X1X2X4,1,mean) - mean(TE.COEFS.EFFECT.S3.X1X2X4) + apply(TE.COEFS.EFFECT.S3.X1X3X4,1,mean) - mean(TE.COEFS.EFFECT.S3.X1X3X4) 
TE.COEFS.EFFECT.ADJ.S1.X2     <- TE.COEFS.EFFECT.S1.X2 - mean(TE.COEFS.EFFECT.S1.X2) + apply(TE.COEFS.EFFECT.S2.X1X2,2,mean) - mean(TE.COEFS.EFFECT.S2.X1X2) + apply(TE.COEFS.EFFECT.S2.X2X3,1,mean) - mean(TE.COEFS.EFFECT.S2.X2X3) + apply(TE.COEFS.EFFECT.S2.X2X4,1,mean) - mean(TE.COEFS.EFFECT.S2.X2X4) + apply(TE.COEFS.EFFECT.S3.X1X2X3,2,mean) - mean(TE.COEFS.EFFECT.S3.X1X2X3) + apply(TE.COEFS.EFFECT.S3.X1X2X4,2,mean) - mean(TE.COEFS.EFFECT.S3.X1X2X4) + apply(TE.COEFS.EFFECT.S3.X2X3X4,1,mean) - mean(TE.COEFS.EFFECT.S3.X2X3X4) 
TE.COEFS.EFFECT.ADJ.S1.X3     <- TE.COEFS.EFFECT.S1.X3 - mean(TE.COEFS.EFFECT.S1.X3) + apply(TE.COEFS.EFFECT.S2.X1X3,2,mean) - mean(TE.COEFS.EFFECT.S2.X1X3) + apply(TE.COEFS.EFFECT.S2.X2X3,2,mean) - mean(TE.COEFS.EFFECT.S2.X2X3) + apply(TE.COEFS.EFFECT.S2.X3X4,1,mean) - mean(TE.COEFS.EFFECT.S2.X3X4) + apply(TE.COEFS.EFFECT.S3.X1X2X3,3,mean) - mean(TE.COEFS.EFFECT.S3.X1X2X3) + apply(TE.COEFS.EFFECT.S3.X1X3X4,2,mean) - mean(TE.COEFS.EFFECT.S3.X1X3X4) + apply(TE.COEFS.EFFECT.S3.X2X3X4,2,mean) - mean(TE.COEFS.EFFECT.S3.X2X3X4) 
TE.COEFS.EFFECT.ADJ.S1.X4     <- TE.COEFS.EFFECT.S1.X4 - mean(TE.COEFS.EFFECT.S1.X4) + apply(TE.COEFS.EFFECT.S2.X1X4,2,mean) - mean(TE.COEFS.EFFECT.S2.X1X4) + apply(TE.COEFS.EFFECT.S2.X2X4,2,mean) - mean(TE.COEFS.EFFECT.S2.X2X4) + apply(TE.COEFS.EFFECT.S2.X3X4,2,mean) - mean(TE.COEFS.EFFECT.S2.X3X4) + apply(TE.COEFS.EFFECT.S3.X1X2X4,3,mean) - mean(TE.COEFS.EFFECT.S3.X1X2X4) + apply(TE.COEFS.EFFECT.S3.X1X3X4,3,mean) - mean(TE.COEFS.EFFECT.S3.X1X3X4) + apply(TE.COEFS.EFFECT.S3.X2X3X4,3,mean) - mean(TE.COEFS.EFFECT.S3.X2X3X4) 
TE.COEFS.EFFECT.ADJ.S0        <- TE.COEFS.EFFECT.S0 + mean(TE.COEFS.EFFECT.S1.X1) + mean(TE.COEFS.EFFECT.S1.X2) + mean(TE.COEFS.EFFECT.S1.X3) + mean(TE.COEFS.EFFECT.S1.X4) + mean(TE.COEFS.EFFECT.S2.X1X2) + mean(TE.COEFS.EFFECT.S2.X1X3) + mean(TE.COEFS.EFFECT.S2.X1X4) + mean(TE.COEFS.EFFECT.S2.X2X3) + mean(TE.COEFS.EFFECT.S2.X2X4) + mean(TE.COEFS.EFFECT.S2.X3X4) + mean(TE.COEFS.EFFECT.S3.X1X2X3) + mean(TE.COEFS.EFFECT.S3.X1X2X4) + mean(TE.COEFS.EFFECT.S3.X1X3X4) + mean(TE.COEFS.EFFECT.S3.X2X3X4)

##Collect adjusted terms
TE.COEFS.ADJ <- c(TE.COEFS.EFFECT.ADJ.B1.X1,TE.COEFS.EFFECT.ADJ.B1.X2,TE.COEFS.EFFECT.ADJ.B1.X3,TE.COEFS.EFFECT.ADJ.B1.X4,
                  TE.COEFS.EFFECT.ADJ.S1.X1,TE.COEFS.EFFECT.ADJ.S1.X2,TE.COEFS.EFFECT.ADJ.S1.X3,TE.COEFS.EFFECT.ADJ.S1.X4,
                  c(TE.COEFS.EFFECT.ADJ.B2.X1X2),c(TE.COEFS.EFFECT.ADJ.B2.X1X3),c(TE.COEFS.EFFECT.ADJ.B2.X1X4),c(TE.COEFS.EFFECT.ADJ.B2.X2X3),c(TE.COEFS.EFFECT.ADJ.B2.X2X4),c(TE.COEFS.EFFECT.ADJ.B2.X3X4),
                  c(TE.COEFS.EFFECT.ADJ.S2.X1X2),c(TE.COEFS.EFFECT.ADJ.S2.X1X3),c(TE.COEFS.EFFECT.ADJ.S2.X1X4),c(TE.COEFS.EFFECT.ADJ.S2.X2X3),c(TE.COEFS.EFFECT.ADJ.S2.X2X4),c(TE.COEFS.EFFECT.ADJ.S2.X3X4),
                  c(TE.COEFS.EFFECT.ADJ.B3.X1X2X3),c(TE.COEFS.EFFECT.ADJ.B3.X1X2X4),c(TE.COEFS.EFFECT.ADJ.B3.X1X3X4),c(TE.COEFS.EFFECT.ADJ.B3.X2X3X4),
                  c(TE.COEFS.EFFECT.ADJ.S3.X1X2X3),c(TE.COEFS.EFFECT.ADJ.S3.X1X2X4),c(TE.COEFS.EFFECT.ADJ.S3.X1X3X4),c(TE.COEFS.EFFECT.ADJ.S3.X2X3X4));names(TE.COEFS.ADJ) <- names(TE.COEFS[-TE.COEFS.INDEX.S0])
#save(TE.COEFS.ADJ,file="TE.COEFS.ADJ.rda")
#save(TE.COEFS.EFFECT.ADJ.B0,file="TE.COEFS.EFFECT.ADJ.B0.rda")
#save(TE.COEFS.EFFECT.ADJ.S0,file="TE.COEFS.EFFECT.ADJ.S0.rda")

##Compare pre and post constraint estimates
#pdf("Constraint.Cor.pdf")
#par(pty="s")

##ME Model
#par(mfrow=c(1,2))
#plot(c(ME.COEFS.EFFECT.B1.X1,    ME.COEFS.EFFECT.B1.X2,    ME.COEFS.EFFECT.B1.X3,    ME.COEFS.EFFECT.B1.X4),
#     c(ME.COEFS.EFFECT.ADJ.B1.X1,ME.COEFS.EFFECT.ADJ.B1.X2,ME.COEFS.EFFECT.ADJ.B1.X3,ME.COEFS.EFFECT.ADJ.B1.X4),
#     pch=19,cex=0.6,xlab="L2 fit Effect",ylab="Constrained Effect",main="ME B1",
#     col=c(rep("red",20),rep("blue",20),rep("green",20),rep("black",20)));abline(a=0,b=1)

#plot(c(ME.COEFS.EFFECT.S1.X1,    ME.COEFS.EFFECT.S1.X2,    ME.COEFS.EFFECT.S1.X3,    ME.COEFS.EFFECT.S1.X4),
#     c(ME.COEFS.EFFECT.ADJ.S1.X1,ME.COEFS.EFFECT.ADJ.S1.X2,ME.COEFS.EFFECT.ADJ.S1.X3,ME.COEFS.EFFECT.ADJ.S1.X4),
#     pch=19,cex=0.6,xlab="L2 fit Effect",ylab="Constrained Effect",main="ME S1",
#     col=c(rep("red",20),rep("blue",20),rep("green",20),rep("black",20)));abline(a=0,b=1)

###PE Model
#par(mfrow=c(2,2))
#plot(c(PE.COEFS.EFFECT.B1.X1,    PE.COEFS.EFFECT.B1.X2,    PE.COEFS.EFFECT.B1.X3,    PE.COEFS.EFFECT.B1.X4),
#     c(PE.COEFS.EFFECT.ADJ.B1.X1,PE.COEFS.EFFECT.ADJ.B1.X2,PE.COEFS.EFFECT.ADJ.B1.X3,PE.COEFS.EFFECT.ADJ.B1.X4),
#     pch=19,cex=0.6,xlab="L2 fit Effect",ylab="Constrained Effect",main="PE B1",
#     col=c(rep("red",20),rep("blue",20),rep("green",20),rep("black",20)));abline(a=0,b=1)

#plot(c(PE.COEFS.EFFECT.S1.X1,    PE.COEFS.EFFECT.S1.X2,    PE.COEFS.EFFECT.S1.X3,    PE.COEFS.EFFECT.S1.X4),
#     c(PE.COEFS.EFFECT.ADJ.S1.X1,PE.COEFS.EFFECT.ADJ.S1.X2,PE.COEFS.EFFECT.ADJ.S1.X3,PE.COEFS.EFFECT.ADJ.S1.X4),
#     pch=19,cex=0.6,xlab="L2 fit Effect",ylab="Constrained Effect",main="PE S1",
#     col=c(rep("red",20),rep("blue",20),rep("green",20),rep("black",20)));abline(a=0,b=1)

#plot(c(PE.COEFS.EFFECT.B2.X1X2,    PE.COEFS.EFFECT.B2.X1X3,    PE.COEFS.EFFECT.B2.X1X4,    PE.COEFS.EFFECT.B2.X2X3,    PE.COEFS.EFFECT.B2.X2X4,    PE.COEFS.EFFECT.B2.X3X4),
#     c(PE.COEFS.EFFECT.ADJ.B2.X1X2,PE.COEFS.EFFECT.ADJ.B2.X1X3,PE.COEFS.EFFECT.ADJ.B2.X1X4,PE.COEFS.EFFECT.ADJ.B2.X2X3,PE.COEFS.EFFECT.ADJ.B2.X2X4,PE.COEFS.EFFECT.ADJ.B2.X3X4),
#     pch=19,cex=0.6,xlab="L2 fit Effect",ylab="Constrained Effect",main="PE B2",
#     col=c(rep("red",20),rep("blue",20),rep("green",20),rep("black",20)));abline(a=0,b=1)

#plot(c(PE.COEFS.EFFECT.S2.X1X2,    PE.COEFS.EFFECT.S2.X1X3,    PE.COEFS.EFFECT.S2.X1X4,    PE.COEFS.EFFECT.S2.X2X3,    PE.COEFS.EFFECT.S2.X2X4,    PE.COEFS.EFFECT.S2.X3X4),
#     c(PE.COEFS.EFFECT.ADJ.S2.X1X2,PE.COEFS.EFFECT.ADJ.S2.X1X3,PE.COEFS.EFFECT.ADJ.S2.X1X4,PE.COEFS.EFFECT.ADJ.S2.X2X3,PE.COEFS.EFFECT.ADJ.S2.X2X4,PE.COEFS.EFFECT.ADJ.S2.X3X4),
#     pch=19,cex=0.6,xlab="L2 fit Effect",ylab="Constrained Effect",main="PE S2",
#     col=c(rep("red",20),rep("blue",20),rep("green",20),rep("black",20)));abline(a=0,b=1)

###TE Model
#par(mfrow=c(3,2))
pdf("CONSTRAINT.B1.pdf")
par(pty="s")
plot(c(TE.COEFS.EFFECT.B1.X1,    TE.COEFS.EFFECT.B1.X2,    TE.COEFS.EFFECT.B1.X3,    TE.COEFS.EFFECT.B1.X4),
    c(TE.COEFS.EFFECT.ADJ.B1.X1,TE.COEFS.EFFECT.ADJ.B1.X2,TE.COEFS.EFFECT.ADJ.B1.X3,TE.COEFS.EFFECT.ADJ.B1.X4),
    pch=19,cex=0.6,xlab="L2 fit Effect",ylab="Constrained Effect",main="TE B1",
    col=c(rep("red",20),rep("blue",20),rep("green",20),rep("black",20)),xlim=c(-4,3),ylim=c(-4,3));abline(a=0,b=1)
dev.off()
pdf("CONSTRAINT.S1.pdf")
par(pty="s")
plot(c(TE.COEFS.EFFECT.S1.X1,    TE.COEFS.EFFECT.S1.X2,    TE.COEFS.EFFECT.S1.X3,    TE.COEFS.EFFECT.S1.X4),
     c(TE.COEFS.EFFECT.ADJ.S1.X1,TE.COEFS.EFFECT.ADJ.S1.X2,TE.COEFS.EFFECT.ADJ.S1.X3,TE.COEFS.EFFECT.ADJ.S1.X4),
     pch=19,cex=0.6,xlab="L2 fit Effect",ylab="Constrained Effect",main="TE S1",
     col=c(rep("red",20),rep("blue",20),rep("green",20),rep("black",20)),xlim=c(-4,3),ylim=c(-4,3));abline(a=0,b=1)
dev.off()
pdf("CONSTRAINT.B2.pdf")
par(pty="s")
plot(c(TE.COEFS.EFFECT.B2.X1X2,    TE.COEFS.EFFECT.B2.X1X3,    TE.COEFS.EFFECT.B2.X1X4,    TE.COEFS.EFFECT.B2.X2X3,    TE.COEFS.EFFECT.B2.X2X4,    TE.COEFS.EFFECT.B2.X3X4),
     c(TE.COEFS.EFFECT.ADJ.B2.X1X2,TE.COEFS.EFFECT.ADJ.B2.X1X3,TE.COEFS.EFFECT.ADJ.B2.X1X4,TE.COEFS.EFFECT.ADJ.B2.X2X3,TE.COEFS.EFFECT.ADJ.B2.X2X4,TE.COEFS.EFFECT.ADJ.B2.X3X4),
     pch=19,cex=0.6,xlab="L2 fit Effect",ylab="Constrained Effect",main="TE B2",
     col=c(rep("red",400),rep("blue",400),rep("green",400),rep("black",400),rep("purple",400),rep("orange",400)),xlim=c(-4,3),ylim=c(-4,3));abline(a=0,b=1)
dev.off()
pdf("CONSTRAINT.S2.pdf")
par(pty="s")
plot(c(TE.COEFS.EFFECT.S2.X1X2,    TE.COEFS.EFFECT.S2.X1X3,    TE.COEFS.EFFECT.S2.X1X4,    TE.COEFS.EFFECT.S2.X2X3,    TE.COEFS.EFFECT.S2.X2X4,    TE.COEFS.EFFECT.S2.X3X4),
     c(TE.COEFS.EFFECT.ADJ.S2.X1X2,TE.COEFS.EFFECT.ADJ.S2.X1X3,TE.COEFS.EFFECT.ADJ.S2.X1X4,TE.COEFS.EFFECT.ADJ.S2.X2X3,TE.COEFS.EFFECT.ADJ.S2.X2X4,TE.COEFS.EFFECT.ADJ.S2.X3X4),
     pch=19,cex=0.6,xlab="L2 fit Effect",ylab="Constrained Effect",main="TE S2",
     col=c(rep("red",400),rep("blue",400),rep("green",400),rep("black",400),rep("purple",400),rep("orange",400)),xlim=c(-4,3),ylim=c(-4,3));abline(a=0,b=1)
dev.off()
pdf("CONSTRAINT.B3.pdf")
par(pty="s")
plot(c(TE.COEFS.EFFECT.B3.X1X2X3,    TE.COEFS.EFFECT.B3.X1X2X4,    TE.COEFS.EFFECT.B3.X1X3X4,    TE.COEFS.EFFECT.B3.X2X3X4),
     c(TE.COEFS.EFFECT.ADJ.B3.X1X2X3,TE.COEFS.EFFECT.ADJ.B3.X1X2X4,TE.COEFS.EFFECT.ADJ.B3.X1X3X4,TE.COEFS.EFFECT.ADJ.B3.X2X3X4),
     pch=19,cex=0.6,xlab="L2 fit Effect",ylab="Constrained Effect",main="TE B3",
     col=c(rep("red",20),rep("blue",20),rep("green",20),rep("black",20)),xlim=c(-4,3),ylim=c(-4,3));abline(a=0,b=1)
dev.off()
pdf("CONSTRAINT.S3.pdf")
par(pty="s")
plot(c(TE.COEFS.EFFECT.S3.X1X2X3,    TE.COEFS.EFFECT.S3.X1X2X4,    TE.COEFS.EFFECT.S3.X1X3X4,    TE.COEFS.EFFECT.S3.X2X3X4),
     c(TE.COEFS.EFFECT.ADJ.S3.X1X2X3,TE.COEFS.EFFECT.ADJ.S3.X1X2X4,TE.COEFS.EFFECT.ADJ.S3.X1X3X4,TE.COEFS.EFFECT.ADJ.S3.X2X3X4),
     pch=19,cex=0.6,xlab="L2 fit Effect",ylab="Constrained Effect",main="TE S3",
     col=c(rep("red",20),rep("blue",20),rep("green",20),rep("black",20)),xlim=c(-4,3),ylim=c(-4,3));abline(a=0,b=1)
dev.off()
#dev.off()

###########################
##13. Genetic score vs dG##
###########################
####Compare genetic score with previously measured dG values
##Ka values from Anderson and Mckeown
dGs <- data.frame("AAseq"=c("EGKA","GGKA","ESKA","EGKV","GSKA","ESKV","GGKV","GSKV"))
dGs$ERE.ka <- c(1.74E+16,8.87E+16,1.16E+13,8.85E+12,8.20E+13,5.32E+11,1.40E+14,1.84E+13)
dGs$SRE.ka <- c(1.58E+14,3.19E+16,4.85E+12,2.02E+12,2.14E+14,1.90E+12,4.90E+15,5.06E+14)

##Order sequences
dGs <- dGs[c(1,4,3,6,2,7,5,8),]

##Convert to kcal
RT <- 0.592

dGs$ERE.dG <- -RT*log(dGs$ERE.ka)
dGs$SRE.dG <- -RT*log(dGs$SRE.ka)
dGs$ERE.meanF <- DT.JOINT$pooled.meanF[which(DT.JOINT$AAseq %in% c("EEGKA","EGGKA","EESKA","EEGKV","EGSKA","EESKV","EGGKV","EGSKV"))]
dGs$SRE.meanF <- DT.JOINT$pooled.meanF[which(DT.JOINT$AAseq %in% c("SEGKA","SGGKA","SESKA","SEGKV","SGSKA","SESKV","SGGKV","SGSKV"))]
dGs$ERE.TE.link  <- -1*PREDICT.TE.LINK[which(DT.JOINT$AAseq %in% c("EEGKA","EGGKA","EESKA","EEGKV","EGSKA","EESKV","EGGKV","EGSKV"))]
dGs$SRE.TE.link  <- -1*PREDICT.TE.LINK[which(DT.JOINT$AAseq %in% c("SEGKA","SGGKA","SESKA","SEGKV","SGSKA","SESKV","SGGKV","SGSKV"))]
dGs$ERE.TE.link <- dGs$ERE.TE.link/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)
dGs$SRE.TE.link <- dGs$SRE.TE.link/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)

##Remove non-binders
dGs$ERE.TE.link[which(dGs$ERE.TE.link < -1*TE.THRESH.NULL/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0))] <- NA
dGs$SRE.TE.link[which(dGs$SRE.TE.link < -1*TE.THRESH.NULL/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0))] <- NA

##pdf("Latent.KMAC.meanF.Cor.pdf")
#par(mfrow=c(1,2))
#par(pty="s")
#plot(c(dGs$ERE.TE.link,dGs$SRE.TE.link),c(-1*dGs$ERE.dG,-1*dGs$SRE.dG),ylab="dG",xlab="TE Genetic Score",pch=19)
#abline(lm(c(-1*dGs$ERE.dG,-1*dGs$SRE.dG) ~ c(dGs$ERE.TE.link,dGs$SRE.TE.link)))

#plot(c(dGs$ERE.meanF,dGs$SRE.meanF),c(-1*dGs$ERE.dG,-1*dGs$SRE.dG),ylab="dG",xlab="meanF",pch=19)
#abline(lm(c(-1*dGs$ERE.dG,-1*dGs$SRE.dG) ~ c(dGs$ERE.meanF,dGs$SRE.meanF)))
##dev.off()

##Simple model to convert change in fluorescence or genetic score to ddG
SCORE.TO.dG <- lm(c(-1*dGs$ERE.dG,-1*dGs$SRE.dG) ~ c(dGs$ERE.TE.link,dGs$SRE.TE.link))
MEANF.TO.dG <- lm(c(-1*dGs$ERE.dG,-1*dGs$SRE.dG) ~ c(dGs$ERE.meanF,dGs$SRE.meanF))

cor.test(c(-1*dGs$ERE.dG,-1*dGs$SRE.dG),c(dGs$ERE.meanF,dGs$SRE.meanF))
cor.test(c(-1*dGs$ERE.dG,-1*dGs$SRE.dG),c(dGs$ERE.TE.link,dGs$SRE.TE.link))

##Comparison to model thresholds
RELATIVE.TE.LINK <- PREDICT.TE.LINK/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)
ddG <- -1*RELATIVE.TE.LINK*SCORE.TO.dG$coefficients[2]

##Plot distribution of genetic scores and ddG values
##pdf("DISTRIBUTION.pdf")
#par(pty="s")
#X <- seq(-19,10,0.01)
#plot(  -1*X/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)*SCORE.TO.dG$coefficients[2],PREDICT.TE.PROBS.NULL,  type="l",ylim=c(0,1),xlim=c(-5,15),ylab="Probability",xlab="ddG")
#points(-1*X/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)*SCORE.TO.dG$coefficients[2],PREDICT.TE.PROBS.WEAK,  type="l",col="blue")
#points(-1*X/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)*SCORE.TO.dG$coefficients[2],PREDICT.TE.PROBS.STRONG,type="l",col="red")
#hist(ddG,probability = TRUE,breaks=50,add=TRUE)
#abline(v=-1*TE.THRESH.NULL/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)*SCORE.TO.dG$coefficients[2])
#abline(v=-1*TE.THRESH.WEAK/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)*SCORE.TO.dG$coefficients[2])

#X <- seq(-19,10,0.01)
#plot(  -1*X,PREDICT.TE.PROBS.NULL,  type="l",ylim=c(0,1),xlim=c(-10,20),ylab="Probability",xlab="Genetic Score")
#points(-1*X,PREDICT.TE.PROBS.WEAK,  type="l",col="blue")
#points(-1*X,PREDICT.TE.PROBS.STRONG,type="l",col="red")
#hist(-1*PREDICT.TE.LINK,probability = FALSE,breaks=50,xlim=c(-10,20),ylim=c(0,40000),main="",xlab="Genetic Score",ylab="Count")
#abline(v=-1*TE.THRESH.NULL)
#abline(v=-1*TE.THRESH.WEAK)
##dev.off()

##########################
##14. Variance Explained##
##########################
####Variance explained by a term is proportional to its squared value, divided by the number of terms of that order of interaction per site/site combination

####ME Model
##Total phenotypic variance
ME.TOTAL.VAR <- mean((PREDICT.ME.LINK - mean(PREDICT.ME.LINK))^2)
##Number of model terms at a given site
ME.N.VAR <- rep(20,160)
##Calculate individual contributions to total variance
ME.ADJ.VAR <- 1/c(ME.N.VAR,1)*c(ME.COEFS.EFFECT.ADJ.B1.X1^2,ME.COEFS.EFFECT.ADJ.B1.X2^2,ME.COEFS.EFFECT.ADJ.B1.X3^2,ME.COEFS.EFFECT.ADJ.B1.X4^2,
                                ME.COEFS.EFFECT.ADJ.S1.X1^2,ME.COEFS.EFFECT.ADJ.S1.X2^2,ME.COEFS.EFFECT.ADJ.S1.X3^2,ME.COEFS.EFFECT.ADJ.S1.X4^2,ME.COEFS.EFFECT.ADJ.S0^2)
##Scale contributions by total variance
ME.REL.VAR <- ME.ADJ.VAR/sum(ME.ADJ.VAR)

####PE Model
##Total phenotypic variance
PE.TOTAL.VAR <- mean((PREDICT.PE.LINK - mean(PREDICT.PE.LINK))^2)
##Number of model terms at a given site
PE.N.VAR <- c(rep(20,160),rep(400,4800))
##Calculate individual contributions to total variance
PE.ADJ.VAR <- 1/c(PE.N.VAR,1)*c(PE.COEFS.EFFECT.ADJ.B1.X1^2,  PE.COEFS.EFFECT.ADJ.B1.X2^2,  PE.COEFS.EFFECT.ADJ.B1.X3^2,  PE.COEFS.EFFECT.ADJ.B1.X4^2,
                                PE.COEFS.EFFECT.ADJ.S1.X1^2,  PE.COEFS.EFFECT.ADJ.S1.X2^2,  PE.COEFS.EFFECT.ADJ.S1.X3^2,  PE.COEFS.EFFECT.ADJ.S1.X4^2,
                                PE.COEFS.EFFECT.ADJ.B2.X1X2^2,PE.COEFS.EFFECT.ADJ.B2.X1X3^2,PE.COEFS.EFFECT.ADJ.B2.X1X4^2,PE.COEFS.EFFECT.ADJ.B2.X2X3^2,PE.COEFS.EFFECT.ADJ.B2.X2X4^2,PE.COEFS.EFFECT.ADJ.B2.X3X4^2,
                                PE.COEFS.EFFECT.ADJ.S2.X1X2^2,PE.COEFS.EFFECT.ADJ.S2.X1X3^2,PE.COEFS.EFFECT.ADJ.S2.X1X4^2,PE.COEFS.EFFECT.ADJ.S2.X2X3^2,PE.COEFS.EFFECT.ADJ.S2.X2X4^2,PE.COEFS.EFFECT.ADJ.S2.X3X4^2,PE.COEFS.EFFECT.ADJ.S0^2)
names(PE.ADJ.VAR) <- names(PE.COEFS)
##Scale contributions by total variance
PE.REL.VAR <- PE.ADJ.VAR/sum(PE.ADJ.VAR)

####TE Model
##Total phenotypic variance
TE.TOTAL.VAR <- mean((PREDICT.TE.LINK - mean(PREDICT.TE.LINK))^2)
##Number of model terms at a given site
TE.N.VAR <- c(rep(20,160),rep(400,4800),rep(8000,64000))
##Calculate individual contributions to total variance
TE.ADJ.VAR <- 1/c(TE.N.VAR,1)*c(TE.COEFS.EFFECT.ADJ.B1.X1^2,    TE.COEFS.EFFECT.ADJ.B1.X2^2,    TE.COEFS.EFFECT.ADJ.B1.X3^2,    TE.COEFS.EFFECT.ADJ.B1.X4^2,
                                TE.COEFS.EFFECT.ADJ.S1.X1^2,    TE.COEFS.EFFECT.ADJ.S1.X2^2,    TE.COEFS.EFFECT.ADJ.S1.X3^2,    TE.COEFS.EFFECT.ADJ.S1.X4^2,
                                TE.COEFS.EFFECT.ADJ.B2.X1X2^2,  TE.COEFS.EFFECT.ADJ.B2.X1X3^2,  TE.COEFS.EFFECT.ADJ.B2.X1X4^2,  TE.COEFS.EFFECT.ADJ.B2.X2X3^2,TE.COEFS.EFFECT.ADJ.B2.X2X4^2,TE.COEFS.EFFECT.ADJ.B2.X3X4^2,
                                TE.COEFS.EFFECT.ADJ.S2.X1X2^2,  TE.COEFS.EFFECT.ADJ.S2.X1X3^2,  TE.COEFS.EFFECT.ADJ.S2.X1X4^2,  TE.COEFS.EFFECT.ADJ.S2.X2X3^2,TE.COEFS.EFFECT.ADJ.S2.X2X4^2,TE.COEFS.EFFECT.ADJ.S2.X3X4^2,
                                TE.COEFS.EFFECT.ADJ.B3.X1X2X3^2,TE.COEFS.EFFECT.ADJ.B3.X1X2X4^2,TE.COEFS.EFFECT.ADJ.B3.X1X3X4^2,TE.COEFS.EFFECT.ADJ.B3.X2X3X4^2,
                                TE.COEFS.EFFECT.ADJ.S3.X1X2X3^2,TE.COEFS.EFFECT.ADJ.S3.X1X2X4^2,TE.COEFS.EFFECT.ADJ.S3.X1X3X4^2,TE.COEFS.EFFECT.ADJ.S3.X2X3X4^2,TE.COEFS.EFFECT.ADJ.S0^2)
names(TE.ADJ.VAR) <- names(TE.COEFS)
##Scale contributions by total variance
TE.REL.VAR <- TE.ADJ.VAR/sum(TE.ADJ.VAR)

####Calculate % of variance explained by different categories and groupings
##Remove S0 term and renormalize
TE.REL.VAR.NOS0 <- TE.REL.VAR[-TE.COEFS.INDEX.S0]/sum(TE.REL.VAR[-TE.COEFS.INDEX.S0])

##Order terms by % of variance explained
REL.VAR.ORDER <- TE.REL.VAR.NOS0[order(TE.REL.VAR.NOS0,decreasing = TRUE)]
REL.VAR.ORDER.SUM  <- cumsum(REL.VAR.ORDER)

##Index ordered effects 
TE.COEFS.INDEX.ORDER.B1        <- grep("^X..$",            names(REL.VAR.ORDER))
TE.COEFS.INDEX.ORDER.B2        <- grep("^X..:X..$",        names(REL.VAR.ORDER))
TE.COEFS.INDEX.ORDER.B3        <- grep("^X..:X..:X..$",    names(REL.VAR.ORDER))
TE.COEFS.INDEX.ORDER.S1        <- grep("^X..:RE1$",        names(REL.VAR.ORDER))
TE.COEFS.INDEX.ORDER.S2        <- grep("^X..:X..:RE1$",    names(REL.VAR.ORDER))
TE.COEFS.INDEX.ORDER.S3        <- grep("^X..:X..:X..:RE1$",names(REL.VAR.ORDER))
TE.COEFS.INDEX.ORDER.B1.X1     <- grep("^X1.$",            names(REL.VAR.ORDER))
TE.COEFS.INDEX.ORDER.B1.X2     <- grep("^X2.$",            names(REL.VAR.ORDER))
TE.COEFS.INDEX.ORDER.B1.X3     <- grep("^X3.$",            names(REL.VAR.ORDER))
TE.COEFS.INDEX.ORDER.B1.X4     <- grep("^X4.$",            names(REL.VAR.ORDER))
TE.COEFS.INDEX.ORDER.S1.X1     <- grep("^X1.:RE1$",        names(REL.VAR.ORDER))
TE.COEFS.INDEX.ORDER.S1.X2     <- grep("^X2.:RE1$",        names(REL.VAR.ORDER))
TE.COEFS.INDEX.ORDER.S1.X3     <- grep("^X3.:RE1$",        names(REL.VAR.ORDER))
TE.COEFS.INDEX.ORDER.S1.X4     <- grep("^X4.:RE1$",        names(REL.VAR.ORDER))
TE.COEFS.INDEX.ORDER.B2.X1X2   <- grep("^X1.:X2.$",        names(REL.VAR.ORDER))
TE.COEFS.INDEX.ORDER.B2.X1X3   <- grep("^X1.:X3.$",        names(REL.VAR.ORDER))
TE.COEFS.INDEX.ORDER.B2.X1X4   <- grep("^X1.:X4.$",        names(REL.VAR.ORDER))
TE.COEFS.INDEX.ORDER.B2.X2X3   <- grep("^X2.:X3.$",        names(REL.VAR.ORDER))
TE.COEFS.INDEX.ORDER.B2.X2X4   <- grep("^X2.:X4.$",        names(REL.VAR.ORDER))
TE.COEFS.INDEX.ORDER.B2.X3X4   <- grep("^X3.:X4.$",        names(REL.VAR.ORDER))
TE.COEFS.INDEX.ORDER.S2.X1X2   <- grep("^X1.:X2.:RE1$",    names(REL.VAR.ORDER))
TE.COEFS.INDEX.ORDER.S2.X1X3   <- grep("^X1.:X3.:RE1$",    names(REL.VAR.ORDER))
TE.COEFS.INDEX.ORDER.S2.X1X4   <- grep("^X1.:X4.:RE1$",    names(REL.VAR.ORDER))
TE.COEFS.INDEX.ORDER.S2.X2X3   <- grep("^X2.:X3.:RE1$",    names(REL.VAR.ORDER))
TE.COEFS.INDEX.ORDER.S2.X2X4   <- grep("^X2.:X4.:RE1$",    names(REL.VAR.ORDER))
TE.COEFS.INDEX.ORDER.S2.X3X4   <- grep("^X3.:X4.:RE1$",    names(REL.VAR.ORDER))
TE.COEFS.INDEX.ORDER.B3.X1X2X3 <- grep("^X1.:X2.:X3.$",    names(REL.VAR.ORDER))
TE.COEFS.INDEX.ORDER.B3.X1X2X4 <- grep("^X1.:X2.:X4.$",    names(REL.VAR.ORDER))
TE.COEFS.INDEX.ORDER.B3.X1X3X4 <- grep("^X1.:X3.:X4.$",    names(REL.VAR.ORDER))
TE.COEFS.INDEX.ORDER.B3.X2X3X4 <- grep("^X2.:X3.:X4.$",    names(REL.VAR.ORDER))
TE.COEFS.INDEX.ORDER.S3.X1X2X3 <- grep("^X1.:X2.:X3.:RE1$",names(REL.VAR.ORDER))
TE.COEFS.INDEX.ORDER.S3.X1X2X4 <- grep("^X1.:X2.:X4.:RE1$",names(REL.VAR.ORDER))
TE.COEFS.INDEX.ORDER.S3.X1X3X4 <- grep("^X1.:X3.:X4.:RE1$",names(REL.VAR.ORDER))
TE.COEFS.INDEX.ORDER.S3.X2X3X4 <- grep("^X2.:X3.:X4.:RE1$",names(REL.VAR.ORDER))

#pdf("Explained.percent.category.pdf")
#par(pty="s")
##Percent explained as the number of terms included increases
##Main and pairwise terms only
#REL.VAR.ORDER.PAIRWISE <- REL.VAR.ORDER[names(REL.VAR.ORDER) %in% names(TE.COEFS[c(TE.COEFS.INDEX.B1,TE.COEFS.INDEX.B2,TE.COEFS.INDEX.S1,TE.COEFS.INDEX.S2)])]
#REL.VAR.ORDER.SUM.PAIRWISE  <- cumsum(REL.VAR.ORDER.PAIRWISE)

#plot(0:length(REL.VAR.ORDER.SUM), c(0,REL.VAR.ORDER.SUM),xlab="# of Coefficients included",ylab="% variance explained",ylim=c(0,1),pch=19,cex=0.4,xlim=c(0,10000))
#points(0:length(REL.VAR.ORDER.SUM.PAIRWISE), c(0,REL.VAR.ORDER.SUM.PAIRWISE),xlab="# of Coefficients included (Main + Pairwise only)",ylab="% variance explained",ylim=c(0,1),type="l",lwd=3,col="gray")
#abline(h=0.9)
#abline(v=min(which(REL.VAR.ORDER.SUM > .9)))
#abline(h=0.99,lty=2)
#abline(v=min(which(REL.VAR.ORDER.SUM > .99)),lty=2)

####Compare variation due to different types of effects
##AAs (binding terms), DNA (S0), and AA-DNA interaction (specificity terms)
#par(pty="s")
#barplot(c(sum(TE.REL.VAR[TE.COEFS.INDEX.B]),sum(TE.REL.VAR[TE.COEFS.INDEX.S0]),sum(TE.REL.VAR[TE.COEFS.INDEX.S])),col=c("firebrick","black","deepskyblue1"),ylim=c(0,1),ylab="% Explained")
#axis(1,at=seq(0.7,3.1,1.2),labels=c("B1-B3","S0","S1-S3"))

##Binding (B) vs specificity (S)
#Percent explained
#barplot(cbind(c(sum(TE.REL.VAR[TE.COEFS.INDEX.B1]),sum(TE.REL.VAR[TE.COEFS.INDEX.B2]),sum(TE.REL.VAR[TE.COEFS.INDEX.B3])),c(sum(TE.REL.VAR[TE.COEFS.INDEX.S1]),sum(TE.REL.VAR[TE.COEFS.INDEX.S2]),sum(TE.REL.VAR[TE.COEFS.INDEX.S3]))),ylim=c(0,1),col=c("firebrick","deepskyblue1","brown"),ylab="% Explained")
#axis(1,at=seq(0.7,1.9,1.2),labels=c("B","S"))

##Binding and specificity split by order and site/site combination
#par(mfrow=c(2,1))
##Percent explained
#barplot(cbind(c(sum(TE.REL.VAR[TE.COEFS.INDEX.B1.X1]),sum(TE.REL.VAR[TE.COEFS.INDEX.B1.X2]),
#          sum(TE.REL.VAR[TE.COEFS.INDEX.B1.X3]),sum(TE.REL.VAR[TE.COEFS.INDEX.B1.X4]),0,0),
#        c(sum(TE.REL.VAR[TE.COEFS.INDEX.B2.X1X2]),sum(TE.REL.VAR[TE.COEFS.INDEX.B2.X1X3]),
#          sum(TE.REL.VAR[TE.COEFS.INDEX.B2.X1X4]),sum(TE.REL.VAR[TE.COEFS.INDEX.B2.X2X3]),
#          sum(TE.REL.VAR[TE.COEFS.INDEX.B2.X2X4]),sum(TE.REL.VAR[TE.COEFS.INDEX.B2.X3X4])),
#        c(sum(TE.REL.VAR[TE.COEFS.INDEX.B3.X1X2X3]),sum(TE.REL.VAR[TE.COEFS.INDEX.B3.X1X2X4]),
#          sum(TE.REL.VAR[TE.COEFS.INDEX.B3.X1X3X4]),sum(TE.REL.VAR[TE.COEFS.INDEX.B3.X2X3X4]),0,0))  
#,ylim=c(0,0.5),col=c("black","firebrick","deepskyblue1","orange","forestgreen","brown"),ylab="% Explained", main="Binding")
#axis(1,at=seq(0.7,3.1,1.2),labels=c("B1","B2","B3"),cex.axis=0.6)

#barplot(cbind(c(sum(TE.REL.VAR[TE.COEFS.INDEX.S1.X1]),    sum(TE.REL.VAR[TE.COEFS.INDEX.S1.X2]),
#                sum(TE.REL.VAR[TE.COEFS.INDEX.S1.X3]),    sum(TE.REL.VAR[TE.COEFS.INDEX.S1.X4]),0,0),
#              c(sum(TE.REL.VAR[TE.COEFS.INDEX.S2.X1X2]),  sum(TE.REL.VAR[TE.COEFS.INDEX.S2.X1X3]),
#                sum(TE.REL.VAR[TE.COEFS.INDEX.S2.X1X4]),  sum(TE.REL.VAR[TE.COEFS.INDEX.S2.X2X3]),
#                sum(TE.REL.VAR[TE.COEFS.INDEX.S2.X2X4]),  sum(TE.REL.VAR[TE.COEFS.INDEX.S2.X3X4])),
#              c(sum(TE.REL.VAR[TE.COEFS.INDEX.S3.X1X2X3]),sum(TE.REL.VAR[TE.COEFS.INDEX.S3.X1X2X4]),
#                sum(TE.REL.VAR[TE.COEFS.INDEX.S3.X1X3X4]),sum(TE.REL.VAR[TE.COEFS.INDEX.S3.X2X3X4]),0,0))  
#        ,ylim=c(0,0.2),col=c("black","firebrick","deepskyblue1","orange","forestgreen","brown"),ylab="% Explained", main="Specificity")
#axis(1,at=seq(0.7,3.1,1.2),labels=c("S1","S2","S3"),cex.axis=0.6)

##Compare main effects vs. variance in epistatic effects
#par(mfrow=c(2,2))
#plot(0,0,type="n",pch=19,cex=0.6,xlim=c(-0.15,0.4),ylim=c(0,.07),xlab="Effect on Binding",ylab="2º variance explained")
#text(paste0(AAs,"1"),x=-1*TE.COEFS.ADJ[TE.COEFS.INDEX.B1.X1]/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0),y=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2.X1X2],dim=c(20,20)),1,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2.X1X3],dim=c(20,20)),1,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2.X1X4],dim=c(20,20)),1,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2])),cex=0.6,col="firebrick")
#text(paste0(AAs,"2"),x=-1*TE.COEFS.ADJ[TE.COEFS.INDEX.B1.X2]/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0),y=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2.X1X2],dim=c(20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2.X2X3],dim=c(20,20)),1,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2.X2X4],dim=c(20,20)),1,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2])),cex=0.6,col="deepskyblue1")
#text(paste0(AAs,"3"),x=-1*TE.COEFS.ADJ[TE.COEFS.INDEX.B1.X3]/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0),y=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2.X1X3],dim=c(20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2.X2X3],dim=c(20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2.X3X4],dim=c(20,20)),1,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2])),cex=0.6,col="forestgreen")
#text(paste0(AAs,"4"),x=-1*TE.COEFS.ADJ[TE.COEFS.INDEX.B1.X4]/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0),y=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2.X1X4],dim=c(20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2.X2X4],dim=c(20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2.X3X4],dim=c(20,20)),2,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2])),cex=0.6,col="black")

#plot(0,0,type="n",pch=19,cex=0.6,xlim=c(-0.15,0.4),ylim=c(0,.07),xlab="Effect on Binding",ylab="3º variance explained")
#text(paste0(AAs,"1"),x=-1*TE.COEFS.ADJ[TE.COEFS.INDEX.B1.X1]/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0),y=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3.X1X2X3],dim=c(20,20,20)),1,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3.X1X2X4],dim=c(20,20,20)),1,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3.X1X3X4],dim=c(20,20,20)),1,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3])),cex=0.6,col="firebrick")
#text(paste0(AAs,"2"),x=-1*TE.COEFS.ADJ[TE.COEFS.INDEX.B1.X2]/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0),y=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3.X1X2X3],dim=c(20,20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3.X1X2X4],dim=c(20,20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3.X2X3X4],dim=c(20,20,20)),1,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3])),cex=0.6,col="deepskyblue1")
#text(paste0(AAs,"3"),x=-1*TE.COEFS.ADJ[TE.COEFS.INDEX.B1.X3]/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0),y=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3.X1X2X3],dim=c(20,20,20)),3,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3.X1X3X4],dim=c(20,20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3.X2X3X4],dim=c(20,20,20)),2,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3])),cex=0.6,col="forestgreen")
#text(paste0(AAs,"4"),x=-1*TE.COEFS.ADJ[TE.COEFS.INDEX.B1.X4]/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0),y=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3.X1X2X4],dim=c(20,20,20)),3,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3.X1X3X4],dim=c(20,20,20)),3,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3.X2X3X4],dim=c(20,20,20)),3,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3])),cex=0.6,col="black")

#plot(0,0,type="n",pch=19,cex=0.6,xlim=c(-0.1,0.1),ylim=c(0,.07),xlab="Effect on Specificity",ylab="2º variance explained")
#text(paste0(AAs,"1"),x=-1*TE.COEFS.ADJ[TE.COEFS.INDEX.S1.X1]/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0),y=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2.X1X2],dim=c(20,20)),1,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2.X1X3],dim=c(20,20)),1,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2.X1X4],dim=c(20,20)),1,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2])),cex=0.6,col="firebrick")
#text(paste0(AAs,"2"),x=-1*TE.COEFS.ADJ[TE.COEFS.INDEX.S1.X2]/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0),y=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2.X1X2],dim=c(20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2.X2X3],dim=c(20,20)),1,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2.X2X4],dim=c(20,20)),1,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2])),cex=0.6,col="deepskyblue1")
#text(paste0(AAs,"3"),x=-1*TE.COEFS.ADJ[TE.COEFS.INDEX.S1.X3]/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0),y=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2.X1X3],dim=c(20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2.X2X3],dim=c(20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2.X3X4],dim=c(20,20)),1,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2])),cex=0.6,col="forestgreen")
#text(paste0(AAs,"4"),x=-1*TE.COEFS.ADJ[TE.COEFS.INDEX.S1.X4]/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0),y=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2.X1X4],dim=c(20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2.X2X4],dim=c(20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2.X3X4],dim=c(20,20)),2,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2])),cex=0.6,col="black")

#plot(0,0,type="n",pch=19,cex=0.6,xlim=c(-0.1,0.1),ylim=c(0,.07),xlab="Effect on Specificity",ylab="3º variance explained")
#text(paste0(AAs,"1"),x=-1*TE.COEFS.ADJ[TE.COEFS.INDEX.S1.X1]/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0),y=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3.X1X2X3],dim=c(20,20,20)),1,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3.X1X2X4],dim=c(20,20,20)),1,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3.X1X3X4],dim=c(20,20,20)),1,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3])),cex=0.6,col="firebrick")
#text(paste0(AAs,"2"),x=-1*TE.COEFS.ADJ[TE.COEFS.INDEX.S1.X2]/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0),y=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3.X1X2X3],dim=c(20,20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3.X1X2X4],dim=c(20,20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3.X2X3X4],dim=c(20,20,20)),1,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3])),cex=0.6,col="deepskyblue1")
#text(paste0(AAs,"3"),x=-1*TE.COEFS.ADJ[TE.COEFS.INDEX.S1.X3]/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0),y=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3.X1X2X3],dim=c(20,20,20)),3,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3.X1X3X4],dim=c(20,20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3.X2X3X4],dim=c(20,20,20)),2,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3])),cex=0.6,col="forestgreen")
#text(paste0(AAs,"4"),x=-1*TE.COEFS.ADJ[TE.COEFS.INDEX.S1.X4]/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0),y=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3.X1X2X4],dim=c(20,20,20)),3,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3.X1X3X4],dim=c(20,20,20)),3,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3.X2X3X4],dim=c(20,20,20)),3,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3])),cex=0.6,col="black")

##Compare variance in 2º and 3º effects
#par(mfrow=c(1,2))
#plot(0,0,type="n",pch=19,cex=0.6,xlim=c(0,.07),ylim=c(0,.07),xlab="2º binding variance explained",ylab="3º binding variance explained")
#text(paste0(AAs,"1"),x=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2.X1X2],dim=c(20,20)),1,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2.X1X3],dim=c(20,20)),1,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2.X1X4],dim=c(20,20)),1,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2])),y=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3.X1X2X3],dim=c(20,20,20)),1,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3.X1X2X4],dim=c(20,20,20)),1,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3.X1X3X4],dim=c(20,20,20)),1,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3])),cex=0.6,col="firebrick")
#text(paste0(AAs,"2"),x=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2.X1X2],dim=c(20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2.X2X3],dim=c(20,20)),1,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2.X2X4],dim=c(20,20)),1,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2])),y=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3.X1X2X3],dim=c(20,20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3.X1X2X4],dim=c(20,20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3.X2X3X4],dim=c(20,20,20)),1,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3])),cex=0.6,col="deepskyblue1")
#text(paste0(AAs,"3"),x=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2.X1X3],dim=c(20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2.X2X3],dim=c(20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2.X3X4],dim=c(20,20)),1,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2])),y=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3.X1X2X3],dim=c(20,20,20)),3,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3.X1X3X4],dim=c(20,20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3.X2X3X4],dim=c(20,20,20)),2,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3])),cex=0.6,col="forestgreen")
#text(paste0(AAs,"4"),x=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2.X1X4],dim=c(20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2.X2X4],dim=c(20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2.X3X4],dim=c(20,20)),2,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2])),y=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3.X1X2X4],dim=c(20,20,20)),3,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3.X1X3X4],dim=c(20,20,20)),3,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3.X2X3X4],dim=c(20,20,20)),3,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3])),cex=0.6,col="black")

#plot(0,0,type="n",pch=19,cex=0.6,xlim=c(0,.07),ylim=c(0,.07),xlab="2º specificity variance explained",ylab="3º specificity variance explained")
#text(paste0(AAs,"1"),x=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2.X1X2],dim=c(20,20)),1,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2.X1X3],dim=c(20,20)),1,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2.X1X4],dim=c(20,20)),1,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2])),y=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3.X1X2X3],dim=c(20,20,20)),1,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3.X1X2X4],dim=c(20,20,20)),1,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3.X1X3X4],dim=c(20,20,20)),1,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3])),cex=0.6,col="firebrick")
#text(paste0(AAs,"2"),x=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2.X1X2],dim=c(20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2.X2X3],dim=c(20,20)),1,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2.X2X4],dim=c(20,20)),1,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2])),y=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3.X1X2X3],dim=c(20,20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3.X1X2X4],dim=c(20,20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3.X2X3X4],dim=c(20,20,20)),1,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3])),cex=0.6,col="deepskyblue1")
#text(paste0(AAs,"3"),x=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2.X1X3],dim=c(20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2.X2X3],dim=c(20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2.X3X4],dim=c(20,20)),1,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2])),y=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3.X1X2X3],dim=c(20,20,20)),3,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3.X1X3X4],dim=c(20,20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3.X2X3X4],dim=c(20,20,20)),2,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3])),cex=0.6,col="forestgreen")
#text(paste0(AAs,"4"),x=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2.X1X4],dim=c(20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2.X2X4],dim=c(20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2.X3X4],dim=c(20,20)),2,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2])),y=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3.X1X2X4],dim=c(20,20,20)),3,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3.X1X3X4],dim=c(20,20,20)),3,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3.X2X3X4],dim=c(20,20,20)),3,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3])),cex=0.6,col="black")

##Compare variance in  binding and specificity epistatic effects
#par(mfrow=c(1,2))
#plot(0,0,type="n",pch=19,cex=0.6,xlim=c(0,.07),ylim=c(0,.07),xlab="2º binding variance explained",ylab="2º specificity variance explained")
#text(paste0(AAs,"1"),x=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2.X1X2],dim=c(20,20)),1,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2.X1X3],dim=c(20,20)),1,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2.X1X4],dim=c(20,20)),1,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2])),y=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2.X1X2],dim=c(20,20)),1,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2.X1X3],dim=c(20,20)),1,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2.X1X4],dim=c(20,20)),1,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2])),cex=0.6,col="firebrick")
#text(paste0(AAs,"2"),x=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2.X1X2],dim=c(20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2.X2X3],dim=c(20,20)),1,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2.X2X4],dim=c(20,20)),1,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2])),y=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2.X1X2],dim=c(20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2.X2X3],dim=c(20,20)),1,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2.X2X4],dim=c(20,20)),1,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2])),cex=0.6,col="deepskyblue1")
#text(paste0(AAs,"3"),x=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2.X1X3],dim=c(20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2.X2X3],dim=c(20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2.X3X4],dim=c(20,20)),1,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2])),y=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2.X1X3],dim=c(20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2.X2X3],dim=c(20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2.X3X4],dim=c(20,20)),1,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2])),cex=0.6,col="forestgreen")
#text(paste0(AAs,"4"),x=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2.X1X4],dim=c(20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2.X2X4],dim=c(20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2.X3X4],dim=c(20,20)),2,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B2])),y=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2.X1X4],dim=c(20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2.X2X4],dim=c(20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2.X3X4],dim=c(20,20)),2,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S2])),cex=0.6,col="black")

#plot(0,0,type="n",pch=19,cex=0.6,xlim=c(0,.07),ylim=c(0,.07),xlab="3º binding variance explained",ylab="3º specificity variance explained")
#text(paste0(AAs,"1"),x=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3.X1X2X3],dim=c(20,20,20)),1,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3.X1X2X4],dim=c(20,20,20)),1,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3.X1X3X4],dim=c(20,20,20)),1,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3])),y=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3.X1X2X3],dim=c(20,20,20)),1,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3.X1X2X4],dim=c(20,20,20)),1,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3.X1X3X4],dim=c(20,20,20)),1,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3])),cex=0.6,col="firebrick")
#text(paste0(AAs,"2"),x=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3.X1X2X3],dim=c(20,20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3.X1X2X4],dim=c(20,20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3.X2X3X4],dim=c(20,20,20)),1,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3])),y=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3.X1X2X3],dim=c(20,20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3.X1X2X4],dim=c(20,20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3.X2X3X4],dim=c(20,20,20)),1,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3])),cex=0.6,col="deepskyblue1")
#text(paste0(AAs,"3"),x=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3.X1X2X3],dim=c(20,20,20)),3,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3.X1X3X4],dim=c(20,20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3.X2X3X4],dim=c(20,20,20)),2,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3])),y=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3.X1X2X3],dim=c(20,20,20)),3,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3.X1X3X4],dim=c(20,20,20)),2,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3.X2X3X4],dim=c(20,20,20)),2,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3])),cex=0.6,col="forestgreen")
#text(paste0(AAs,"4"),x=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3.X1X2X4],dim=c(20,20,20)),3,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3.X1X3X4],dim=c(20,20,20)),3,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3.X2X3X4],dim=c(20,20,20)),3,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.B3])),y=(apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3.X1X2X4],dim=c(20,20,20)),3,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3.X1X3X4],dim=c(20,20,20)),3,sum) + apply(array(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3.X2X3X4],dim=c(20,20,20)),3,sum))/(3*sum(TE.REL.VAR.NOS0[TE.COEFS.INDEX.S3])),cex=0.6,col="black")

####Identity of important effects
IMPORTANT.EFFECTS.INDEX <- 1:min(which(REL.VAR.ORDER.SUM > .99))
IMPORTANT.EFFECTS <- REL.VAR.ORDER[IMPORTANT.EFFECTS.INDEX]

##Main effects
IMPORTANT.EFFECTS.B1 <- IMPORTANT.EFFECTS[TE.COEFS.INDEX.ORDER.B1]; IMPORTANT.EFFECTS.B1 <- IMPORTANT.EFFECTS.B1[!is.na(IMPORTANT.EFFECTS.B1)]
IMPORTANT.EFFECTS.S1 <- IMPORTANT.EFFECTS[TE.COEFS.INDEX.ORDER.S1]; IMPORTANT.EFFECTS.S1 <- IMPORTANT.EFFECTS.B1[!is.na(IMPORTANT.EFFECTS.S1)]

##Epistatic effects
IMPORTANT.EFFECTS.B2 <- IMPORTANT.EFFECTS[TE.COEFS.INDEX.ORDER.B2]; IMPORTANT.EFFECTS.B2 <- IMPORTANT.EFFECTS.B2[!is.na(IMPORTANT.EFFECTS.B2)]
IMPORTANT.EFFECTS.S2 <- IMPORTANT.EFFECTS[TE.COEFS.INDEX.ORDER.S2]; IMPORTANT.EFFECTS.S2 <- IMPORTANT.EFFECTS.B2[!is.na(IMPORTANT.EFFECTS.S2)]
IMPORTANT.EFFECTS.B3 <- IMPORTANT.EFFECTS[TE.COEFS.INDEX.ORDER.B3]; IMPORTANT.EFFECTS.B3 <- IMPORTANT.EFFECTS.B3[!is.na(IMPORTANT.EFFECTS.B3)]
IMPORTANT.EFFECTS.S3 <- IMPORTANT.EFFECTS[TE.COEFS.INDEX.ORDER.S3]; IMPORTANT.EFFECTS.S3 <- IMPORTANT.EFFECTS.B3[!is.na(IMPORTANT.EFFECTS.S3)]

I <- 1:4
J <- 1:length(AAs)
IMPORTANT.EFFECTS.B2.TABLE.POS <- matrix(0,nrow=length(J),ncol=length(I)); rownames(IMPORTANT.EFFECTS.B2.TABLE.POS) <- AAs; colnames(IMPORTANT.EFFECTS.B2.TABLE.POS) <- c("X1","X2","X3","X4")
IMPORTANT.EFFECTS.B2.TABLE.NEG <- matrix(0,nrow=length(J),ncol=length(I)); rownames(IMPORTANT.EFFECTS.B2.TABLE.NEG) <- AAs; colnames(IMPORTANT.EFFECTS.B2.TABLE.NEG) <- c("X1","X2","X3","X4")
IMPORTANT.EFFECTS.S2.TABLE.POS <- matrix(0,nrow=length(J),ncol=length(I)); rownames(IMPORTANT.EFFECTS.S2.TABLE.POS) <- AAs; colnames(IMPORTANT.EFFECTS.S2.TABLE.POS) <- c("X1","X2","X3","X4")
IMPORTANT.EFFECTS.S2.TABLE.NEG <- matrix(0,nrow=length(J),ncol=length(I)); rownames(IMPORTANT.EFFECTS.S2.TABLE.NEG) <- AAs; colnames(IMPORTANT.EFFECTS.S2.TABLE.NEG) <- c("X1","X2","X3","X4")
IMPORTANT.EFFECTS.B3.TABLE.POS <- matrix(0,nrow=length(J),ncol=length(I)); rownames(IMPORTANT.EFFECTS.B3.TABLE.POS) <- AAs; colnames(IMPORTANT.EFFECTS.B3.TABLE.POS) <- c("X1","X2","X3","X4")
IMPORTANT.EFFECTS.B3.TABLE.NEG <- matrix(0,nrow=length(J),ncol=length(I)); rownames(IMPORTANT.EFFECTS.B3.TABLE.NEG) <- AAs; colnames(IMPORTANT.EFFECTS.B3.TABLE.NEG) <- c("X1","X2","X3","X4")
IMPORTANT.EFFECTS.S3.TABLE.POS <- matrix(0,nrow=length(J),ncol=length(I)); rownames(IMPORTANT.EFFECTS.S3.TABLE.POS) <- AAs; colnames(IMPORTANT.EFFECTS.S3.TABLE.POS) <- c("X1","X2","X3","X4")
IMPORTANT.EFFECTS.S3.TABLE.NEG <- matrix(0,nrow=length(J),ncol=length(I)); rownames(IMPORTANT.EFFECTS.S3.TABLE.NEG) <- AAs; colnames(IMPORTANT.EFFECTS.S3.TABLE.NEG) <- c("X1","X2","X3","X4")
for(i in I) {
  for(j in J) {
    INDEX <- grep(paste("X",i,AAs[j],sep="")  ,names(IMPORTANT.EFFECTS.B2))
    IMPORTANT.EFFECTS.B2.TABLE.POS[j,i] <- sum(TE.COEFS.ADJ[TE.COEFS.INDEX.ORDER.B2][INDEX] < 0)
    IMPORTANT.EFFECTS.B2.TABLE.NEG[j,i] <- sum(TE.COEFS.ADJ[TE.COEFS.INDEX.ORDER.B2][INDEX] > 0)
    INDEX <- grep(paste("X",i,AAs[j],sep="")  ,names(IMPORTANT.EFFECTS.S2))
    IMPORTANT.EFFECTS.S2.TABLE.POS[j,i] <- sum(TE.COEFS.ADJ[TE.COEFS.INDEX.ORDER.S2][INDEX] < 0)
    IMPORTANT.EFFECTS.S2.TABLE.NEG[j,i] <- sum(TE.COEFS.ADJ[TE.COEFS.INDEX.ORDER.S2][INDEX] > 0)
    INDEX <- grep(paste("X",i,AAs[j],sep="")  ,names(IMPORTANT.EFFECTS.B3))
    IMPORTANT.EFFECTS.B3.TABLE.POS[j,i] <- sum(TE.COEFS.ADJ[TE.COEFS.INDEX.ORDER.B3][INDEX] < 0)
    IMPORTANT.EFFECTS.B3.TABLE.NEG[j,i] <- sum(TE.COEFS.ADJ[TE.COEFS.INDEX.ORDER.B3][INDEX] > 0)
    INDEX <- grep(paste("X",i,AAs[j],sep="")  ,names(IMPORTANT.EFFECTS.S3))
    IMPORTANT.EFFECTS.S3.TABLE.POS[j,i] <- sum(TE.COEFS.ADJ[TE.COEFS.INDEX.ORDER.S3][INDEX] < 0)
    IMPORTANT.EFFECTS.S3.TABLE.NEG[j,i] <- sum(TE.COEFS.ADJ[TE.COEFS.INDEX.ORDER.S3][INDEX] > 0)
  }
}

####################
##15. Effect sizes##
####################
####Median effect relative to threshold to be a strong activator
median(abs(TE.COEFS.ADJ[TE.COEFS.INDEX.O1]))/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)
median(abs(TE.COEFS.ADJ[TE.COEFS.INDEX.O2]))/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)
median(abs(TE.COEFS.ADJ[TE.COEFS.INDEX.O3]))/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)

median(abs(TE.COEFS.ADJ[TE.COEFS.INDEX.B1]))/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)
median(abs(TE.COEFS.ADJ[TE.COEFS.INDEX.B2]))/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)
median(abs(TE.COEFS.ADJ[TE.COEFS.INDEX.B3]))/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)

median(abs(TE.COEFS.ADJ[TE.COEFS.INDEX.S1]))/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)
median(abs(TE.COEFS.ADJ[TE.COEFS.INDEX.S2]))/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)
median(abs(TE.COEFS.ADJ[TE.COEFS.INDEX.S3]))/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)

####Table of effects per genotype
EFFECT.TABLE.TE <- data.frame(matrix(0,nrow=320000,ncol=32))
colnames(EFFECT.TABLE.TE) <- c("SEQ","LINK","B0","S0","BX1","BX2","BX3","BX4","SX1","SX2","SX3","SX4","BX1X2","BX1X3","BX1X4","BX2X3","BX2X4","BX3X4","SX1X2","SX1X3","SX1X4","SX2X3","SX2X4","SX3X4","BX1X2X3","BX1X2X4","BX1X3X4","BX2X3X4","SX1X2X3","SX1X2X4","SX1X3X4","SX2X3X4")
EFFECT.TABLE.TE$SEQ <- DT.JOINT$AAseq
EFFECT.TABLE.TE$LINK <- PREDICT.TE.LINK
EFFECT.TABLE.TE$B0 <- TE.COEFS.EFFECT.ADJ.B0
EFFECT.TABLE.TE$S0 <- c(rep(TE.COEFS.EFFECT.ADJ.S0,160000),rep(-1*TE.COEFS.EFFECT.ADJ.S0,160000))

TE.COEF.MATRIX.ADJ    <- t(t(TE.MATRIX[,-TE.COEFS.INDEX.S0])*TE.COEFS.ADJ)
#save(TE.COEF.MATRIX.ADJ, file = "TE.COEF.MATRIX.ADJ.rda")

#I <- 1:320000
#for(i in I) {
#  SET <- which(TE.COEF.MATRIX.ADJ[i,] != 0)
#  EFFECT.TABLE.TE[i,5:32] <- TE.COEF.MATRIX.ADJ[i,SET]
#}
#EFFECT.TABLE.TE[,2:32] <- -1*EFFECT.TABLE.TE[,3:32]
#save(EFFECT.TABLE.TE,file="EFFECT.TABLE.TE.rda")
load("EFFECT.TABLE.TE.rda")


##Single sites sufficient
X1 <- EFFECT.TABLE.TE$B0 + EFFECT.TABLE.TE$BX1 + EFFECT.TABLE.TE$S0 + EFFECT.TABLE.TE$SX1; X1 <- X1/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)
X2 <- EFFECT.TABLE.TE$B0 + EFFECT.TABLE.TE$BX2 + EFFECT.TABLE.TE$S0 + EFFECT.TABLE.TE$SX2; X2 <- X2/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)
X3 <- EFFECT.TABLE.TE$B0 + EFFECT.TABLE.TE$BX3 + EFFECT.TABLE.TE$S0 + EFFECT.TABLE.TE$SX3; X3 <- X3/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)
X4 <- EFFECT.TABLE.TE$B0 + EFFECT.TABLE.TE$BX4 + EFFECT.TABLE.TE$S0 + EFFECT.TABLE.TE$SX4; X4 <- X4/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)

##Pairs of sites sufficient
X1X2 <- EFFECT.TABLE.TE$B0 + EFFECT.TABLE.TE$BX1 + EFFECT.TABLE.TE$BX2 + EFFECT.TABLE.TE$BX1X2 + EFFECT.TABLE.TE$S0 + EFFECT.TABLE.TE$SX1 + EFFECT.TABLE.TE$SX2 + EFFECT.TABLE.TE$SX1X2; X1X2 <- X1X2/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)
X1X3 <- EFFECT.TABLE.TE$B0 + EFFECT.TABLE.TE$BX1 + EFFECT.TABLE.TE$BX3 + EFFECT.TABLE.TE$BX1X3 + EFFECT.TABLE.TE$S0 + EFFECT.TABLE.TE$SX1 + EFFECT.TABLE.TE$SX3 + EFFECT.TABLE.TE$SX1X3; X1X3 <- X1X3/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)
X1X4 <- EFFECT.TABLE.TE$B0 + EFFECT.TABLE.TE$BX1 + EFFECT.TABLE.TE$BX4 + EFFECT.TABLE.TE$BX1X4 + EFFECT.TABLE.TE$S0 + EFFECT.TABLE.TE$SX1 + EFFECT.TABLE.TE$SX4 + EFFECT.TABLE.TE$SX1X4; X1X4 <- X1X4/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)
X2X3 <- EFFECT.TABLE.TE$B0 + EFFECT.TABLE.TE$BX2 + EFFECT.TABLE.TE$BX4 + EFFECT.TABLE.TE$BX2X3 + EFFECT.TABLE.TE$S0 + EFFECT.TABLE.TE$SX2 + EFFECT.TABLE.TE$SX3 + EFFECT.TABLE.TE$SX2X3; X2X3 <- X2X3/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)
X2X4 <- EFFECT.TABLE.TE$B0 + EFFECT.TABLE.TE$BX2 + EFFECT.TABLE.TE$BX4 + EFFECT.TABLE.TE$BX2X4 + EFFECT.TABLE.TE$S0 + EFFECT.TABLE.TE$SX2 + EFFECT.TABLE.TE$SX4 + EFFECT.TABLE.TE$SX2X4; X2X4 <- X2X4/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)
X3X4 <- EFFECT.TABLE.TE$B0 + EFFECT.TABLE.TE$BX3 + EFFECT.TABLE.TE$BX4 + EFFECT.TABLE.TE$BX3X4 + EFFECT.TABLE.TE$S0 + EFFECT.TABLE.TE$SX3 + EFFECT.TABLE.TE$SX4 + EFFECT.TABLE.TE$SX3X4; X3X4 <- X3X4/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)

##Triples of sites sufficient
X1X2X3 <- EFFECT.TABLE.TE$B0 + EFFECT.TABLE.TE$BX1 + EFFECT.TABLE.TE$BX2 + EFFECT.TABLE.TE$BX3 + EFFECT.TABLE.TE$BX1X2 + EFFECT.TABLE.TE$BX1X3 + EFFECT.TABLE.TE$BX2X3 + EFFECT.TABLE.TE$BX1X2X3 + EFFECT.TABLE.TE$S0 + EFFECT.TABLE.TE$SX1 + EFFECT.TABLE.TE$SX2 + EFFECT.TABLE.TE$SX3 + EFFECT.TABLE.TE$SX1X2 + EFFECT.TABLE.TE$SX1X3 + EFFECT.TABLE.TE$SX2X3 + EFFECT.TABLE.TE$SX1X2X3; X1X2X3 <- X1X2X3/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)
X1X2X4 <- EFFECT.TABLE.TE$B0 + EFFECT.TABLE.TE$BX1 + EFFECT.TABLE.TE$BX2 + EFFECT.TABLE.TE$BX4 + EFFECT.TABLE.TE$BX1X2 + EFFECT.TABLE.TE$BX1X4 + EFFECT.TABLE.TE$BX2X4 + EFFECT.TABLE.TE$BX1X2X4 + EFFECT.TABLE.TE$S0 + EFFECT.TABLE.TE$SX1 + EFFECT.TABLE.TE$SX2 + EFFECT.TABLE.TE$SX4 + EFFECT.TABLE.TE$SX1X2 + EFFECT.TABLE.TE$SX1X4 + EFFECT.TABLE.TE$SX2X4 + EFFECT.TABLE.TE$SX1X2X4; X1X2X4 <- X1X2X4/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)
X1X3X4 <- EFFECT.TABLE.TE$B0 + EFFECT.TABLE.TE$BX1 + EFFECT.TABLE.TE$BX3 + EFFECT.TABLE.TE$BX4 + EFFECT.TABLE.TE$BX1X3 + EFFECT.TABLE.TE$BX1X4 + EFFECT.TABLE.TE$BX3X4 + EFFECT.TABLE.TE$BX1X3X4 + EFFECT.TABLE.TE$S0 + EFFECT.TABLE.TE$SX1 + EFFECT.TABLE.TE$SX3 + EFFECT.TABLE.TE$SX4 + EFFECT.TABLE.TE$SX1X3 + EFFECT.TABLE.TE$SX1X4 + EFFECT.TABLE.TE$SX3X4 + EFFECT.TABLE.TE$SX1X3X4; X1X3X4 <- X1X3X4/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)
X2X3X4 <- EFFECT.TABLE.TE$B0 + EFFECT.TABLE.TE$BX2 + EFFECT.TABLE.TE$BX3 + EFFECT.TABLE.TE$BX4 + EFFECT.TABLE.TE$BX2X3 + EFFECT.TABLE.TE$BX2X4 + EFFECT.TABLE.TE$BX3X4 + EFFECT.TABLE.TE$BX2X3X4 + EFFECT.TABLE.TE$S0 + EFFECT.TABLE.TE$SX2 + EFFECT.TABLE.TE$SX3 + EFFECT.TABLE.TE$SX4 + EFFECT.TABLE.TE$SX2X3 + EFFECT.TABLE.TE$SX2X4 + EFFECT.TABLE.TE$SX3X4 + EFFECT.TABLE.TE$SX2X3X4; X2X3X4 <- X2X3X4/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)


###Effect correlations
#pdf("Cor.plots.pdf")
X1 <- -1*TE.COEFS.ADJ[TE.COEFS.INDEX.B1.X1]/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)
X2 <- -1*TE.COEFS.ADJ[TE.COEFS.INDEX.B1.X2]/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)
X3 <- -1*TE.COEFS.ADJ[TE.COEFS.INDEX.B1.X3]/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)
X4 <- -1*TE.COEFS.ADJ[TE.COEFS.INDEX.B1.X4]/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)
X12 <- matrix(-1*TE.COEFS.ADJ[TE.COEFS.INDEX.B2.X1X2]/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0),nrow=20)
X13 <- matrix(-1*TE.COEFS.ADJ[TE.COEFS.INDEX.B2.X1X3]/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0),nrow=20)
X14 <- matrix(-1*TE.COEFS.ADJ[TE.COEFS.INDEX.B2.X1X4]/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0),nrow=20)
X23 <- matrix(-1*TE.COEFS.ADJ[TE.COEFS.INDEX.B2.X2X3]/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0),nrow=20)
X24 <- matrix(-1*TE.COEFS.ADJ[TE.COEFS.INDEX.B2.X2X4]/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0),nrow=20)
X34 <- matrix(-1*TE.COEFS.ADJ[TE.COEFS.INDEX.B2.X3X4]/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0),nrow=20)
COUNT <- 0
CORS <- numeric(length=20*3*4)
I <- 1:20
#for(i in I) {
#  COUNT <- COUNT + 1
#  plot(X1,X12[,i],cex=0.0001,main=AAs[i],xlab="Main effect 1",ylab="Interaction 2")
#  text(X1,X12[,i],label=AAs)
#  CORS[COUNT] <- cor(X1,X12[,i])
#}
#for(i in I) {
#  COUNT <- COUNT + 1
#  plot(X1,X13[,i],cex=0.0001,main=AAs[i],xlab="Main effect 1",ylab="Interaction 3")
#  text(X1,X13[,i],label=AAs)
#  CORS[COUNT] <- cor(X1,X13[,i])
#}
#for(i in I) {
#  COUNT <- COUNT + 1
#  plot(X1,X14[,i],cex=0.0001,main=AAs[i],xlab="Main effect 1",ylab="Interaction 4")
#  text(X1,X14[,i],label=AAs)
#  CORS[COUNT] <- cor(X1,X14[,i])
#}
#for(i in I) {
#  COUNT <- COUNT + 1
#  plot(X2,X12[i,],cex=0.0001,main=AAs[i],xlab="Main effect 2",ylab="Interaction 1")
#  text(X2,X12[i,],label=AAs)
#  CORS[COUNT] <- cor(X2,X12[i,])
#}
#for(i in I) {
#  COUNT <- COUNT + 1
#  plot(X2,X23[,i],cex=0.0001,main=AAs[i],xlab="Main effect 2",ylab="Interaction 3")
#  text(X2,X23[,i],label=AAs)
#  CORS[COUNT] <- cor(X2,X23[,i])
#}
#for(i in I) {
#  COUNT <- COUNT + 1
#  plot(X2,X24[,i],cex=0.0001,main=AAs[i],xlab="Main effect 2",ylab="Interaction 4")
#  text(X2,X24[,i],label=AAs)
#  CORS[COUNT] <- cor(X2,X24[,i])
#}
#for(i in I) {
#  COUNT <- COUNT + 1
#  plot(X3,X13[i,],cex=0.0001,main=AAs[i],xlab="Main effect 3",ylab="Interaction 1")
#  text(X3,X13[i,],label=AAs)
#  CORS[COUNT] <- cor(X3,X13[i,])
#}
#for(i in I) {
#  COUNT <- COUNT + 1
#  plot(X3,X23[i,],cex=0.0001,main=AAs[i],xlab="Main effect 3",ylab="Interaction 2")
#  text(X3,X23[i,],label=AAs)
#  CORS[COUNT] <- cor(X3,X23[i,])
#}
#for(i in I) {
#  COUNT <- COUNT + 1
#  plot(X3,X34[,i],cex=0.0001,main=AAs[i],xlab="Main effect 3",ylab="Interaction 4")
#  text(X3,X34[,i],label=AAs)
#  CORS[COUNT] <- cor(X3,X34[,i])
#}
#for(i in I) {
#  COUNT <- COUNT + 1
#  plot(X4,X14[i,],cex=0.0001,main=AAs[i],xlab="Main effect 4",ylab="Interaction 1")
#  text(X4,X14[i,],label=AAs)
#  CORS[COUNT] <- cor(X4,X14[i,])
#}
#for(i in I) {
#  COUNT <- COUNT + 1
#  plot(X4,X24[i,],cex=0.0001,main=AAs[i],xlab="Main effect 4",ylab="Interaction 2")
#  text(X4,X24[i,],label=AAs)
#  CORS[COUNT] <- cor(X4,X24[i,])
#}
#for(i in I) {
#  COUNT <- COUNT + 1
#  plot(X4,X34[i,],cex=0.0001,main=AAs[i],xlab="Main effect 4",ylab="Interaction 3")
#  text(X4,X34[i,],label=AAs)
#  CORS[COUNT] <- cor(X4,X34[i,])
#}
#dev.off()


############################
##16. Sequence Differences##
############################
####Sequence logos of weak and strong activators
ERE.STRONG.SEQS <- subset(DT.JOINT,DT.JOINT$RE=="E" & as.vector(PREDICT.TE.CLASS$class) == "strong")$AAseq; ERE.STRONG.SEQS <- substr(ERE.STRONG.SEQS,2,5)
ERE.WEAK.SEQS   <- subset(DT.JOINT,DT.JOINT$RE=="E" & as.vector(PREDICT.TE.CLASS$class) == "weak"  )$AAseq; ERE.WEAK.SEQS   <- substr(ERE.WEAK.SEQS,2,5)
SRE.STRONG.SEQS <- subset(DT.JOINT,DT.JOINT$RE=="S" & as.vector(PREDICT.TE.CLASS$class) == "strong")$AAseq; SRE.STRONG.SEQS <- substr(SRE.STRONG.SEQS,2,5)
SRE.WEAK.SEQS   <- subset(DT.JOINT,DT.JOINT$RE=="S" & as.vector(PREDICT.TE.CLASS$class) == "weak"  )$AAseq; SRE.WEAK.SEQS   <- substr(SRE.WEAK.SEQS,2,5)

ERE.STRONG.LOGO <- ggplot() + geom_logo(ERE.STRONG.SEQS,method="p") + theme(legend.position="none")
ERE.WEAK.LOGO   <- ggplot() + geom_logo(ERE.WEAK.SEQS,method="p")   + theme(legend.position="none")
SRE.STRONG.LOGO <- ggplot() + geom_logo(SRE.STRONG.SEQS,method="p") + theme(legend.position="none")
SRE.WEAK.LOGO   <- ggplot() + geom_logo(SRE.WEAK.SEQS,method="p")   + theme(legend.position="none")

#pdf("SEQ.LOGOS.pdf")
#print(ERE.STRONG.LOGO) 
#print(ERE.WEAK.LOGO)   
#print(SRE.STRONG.LOGO) 
#print(SRE.WEAK.LOGO)   
#dev.off()

###Distribution of sequences differences for strong binders
I <- 1:(length(ERE.STRONG.SEQS)-1)
ERE.SEQ.DIFF <- NULL
for(i in I) {
  J <- (i+1):length(ERE.STRONG.SEQS)
  for(j in J) {
    ERE.SEQ.DIFF <- c(ERE.SEQ.DIFF,get.HAM.distance(ERE.STRONG.SEQS[i],ERE.STRONG.SEQS[j]))
  }
}
ERE.SEQ.DIFF <- (4 - ERE.SEQ.DIFF)

I <- 1:(length(SRE.STRONG.SEQS)-1)
SRE.SEQ.DIFF <- NULL
for(i in I) {
  J <- (i+1):length(SRE.STRONG.SEQS)
  for(j in J) {
    SRE.SEQ.DIFF <- c(SRE.SEQ.DIFF,get.HAM.distance(SRE.STRONG.SEQS[i],SRE.STRONG.SEQS[j]))
  }
}
SRE.SEQ.DIFF <- (4 - SRE.SEQ.DIFF)

#pdf("Seq.Diffs.pdf")
#par(mfrow=c(2,1))
#par(pty="s")
#hist(ERE.SEQ.DIFF,xlab="Sequence Similarity",main="ERE")
#hist(SRE.SEQ.DIFF,xlab="Sequence Similarity",main="SRE")
#dev.off()

##########################
##17. Types of Epistasis##
###########################
####Classify, count, and determine percent of variance explained by different forms of epistasis
##Second order binding interactions
TE.B2.EPISTASIS <- data.frame(matrix(0,nrow=2400,ncol=9))
colnames(TE.B2.EPISTASIS) <- c("A","B","AB","Expected","Observed","Sign.A","Sign.B","Sign.A.AB","Sign.B.AB")
rownames(TE.B2.EPISTASIS) <- names(TE.COEFS.ADJ[TE.COEFS.INDEX.B2])
##Get first amino acid state effect by itself
TE.B2.EPISTASIS$A <- c(rep(TE.COEFS.ADJ[TE.COEFS.INDEX.B1.X1],60),rep(TE.COEFS.ADJ[TE.COEFS.INDEX.B1.X2],40),rep(TE.COEFS.ADJ[TE.COEFS.INDEX.B1.X3],20))
##Get second amino acid state effect by itself
TE.B2.EPISTASIS$B <- c(rep(TE.COEFS.ADJ[TE.COEFS.INDEX.B1.X2],each=20),rep(TE.COEFS.ADJ[TE.COEFS.INDEX.B1.X3],each=20),rep(TE.COEFS.ADJ[TE.COEFS.INDEX.B1.X4],each=20),rep(TE.COEFS.ADJ[TE.COEFS.INDEX.B1.X3],each=20),rep(TE.COEFS.ADJ[TE.COEFS.INDEX.B1.X4],each=20),rep(TE.COEFS.ADJ[TE.COEFS.INDEX.B1.X4],each=20))
##Get pairwise epistastic effect between amino acid states
TE.B2.EPISTASIS$AB <- TE.COEFS.ADJ[TE.COEFS.INDEX.B2]
##Expected phenotype in absence of epistasis (A+B)
TE.B2.EPISTASIS$Expected <- TE.B2.EPISTASIS[,1] + TE.B2.EPISTASIS[,2]
##Observed phenotype in presence of epistasis (A+B+AB)
TE.B2.EPISTASIS$Observed <- TE.B2.EPISTASIS[,1] + TE.B2.EPISTASIS[,2] + TE.B2.EPISTASIS[,3]
##Direction of amino acid A by itself
TE.B2.EPISTASIS$Sign.A <- sign(TE.B2.EPISTASIS[,1])
##Direction of amino acid B by itself
TE.B2.EPISTASIS$Sign.B <- sign(TE.B2.EPISTASIS[,2])
##Direction of amino acid A in presence of B
TE.B2.EPISTASIS$Sign.A.AB <- sign(TE.B2.EPISTASIS[,1] + TE.B2.EPISTASIS[,3])
##Direction of amino acid B in presence of A
TE.B2.EPISTASIS$Sign.B.AB <- sign(TE.B2.EPISTASIS[,2] + TE.B2.EPISTASIS[,3])
##Flip all signs so that an increase in function is positive
TE.B2.EPISTASIS <- -1*TE.B2.EPISTASIS

####Classify epistasis
TE.B2.EPISTASIS.LIST <- list()
TE.B2.EPISTASIS.LIST$PPP  <- which(  TE.B2.EPISTASIS$Sign.A > 0 & TE.B2.EPISTASIS$Sign.B > 0 & TE.B2.EPISTASIS$Observed > TE.B2.EPISTASIS$Expected)
TE.B2.EPISTASIS.LIST$PPN  <- which(  TE.B2.EPISTASIS$Sign.A > 0 & TE.B2.EPISTASIS$Sign.B > 0 & TE.B2.EPISTASIS$Observed < TE.B2.EPISTASIS$Expected)
TE.B2.EPISTASIS.LIST$PNP  <- which(((TE.B2.EPISTASIS$Sign.A > 0 & TE.B2.EPISTASIS$Sign.B < 0) | (TE.B2.EPISTASIS$Sign.A < 0 & TE.B2.EPISTASIS$Sign.B > 0)) & TE.B2.EPISTASIS$Observed > TE.B2.EPISTASIS$Expected)
TE.B2.EPISTASIS.LIST$PNN  <- which(((TE.B2.EPISTASIS$Sign.A > 0 & TE.B2.EPISTASIS$Sign.B < 0) | (TE.B2.EPISTASIS$Sign.A < 0 & TE.B2.EPISTASIS$Sign.B > 0)) & TE.B2.EPISTASIS$Observed < TE.B2.EPISTASIS$Expected)
TE.B2.EPISTASIS.LIST$NNP  <- which(  TE.B2.EPISTASIS$Sign.A < 0 & TE.B2.EPISTASIS$Sign.B < 0 & TE.B2.EPISTASIS$Observed > TE.B2.EPISTASIS$Expected)
TE.B2.EPISTASIS.LIST$NNN  <- which(  TE.B2.EPISTASIS$Sign.A < 0 & TE.B2.EPISTASIS$Sign.B < 0 & TE.B2.EPISTASIS$Observed < TE.B2.EPISTASIS$Expected)

TE.B2.EPISTASIS.REL.VAR <- matrix(0,nrow=3,ncol=2);rownames(TE.B2.EPISTASIS.REL.VAR) <- c("PP","PN","NN"); colnames(TE.B2.EPISTASIS.REL.VAR) <- c("+","-")
TE.B2.EPISTASIS.REL.VAR[1,1] <- sum(TE.REL.VAR[TE.COEFS.INDEX.B2][TE.B2.EPISTASIS.LIST$PPP])
TE.B2.EPISTASIS.REL.VAR[1,2] <- sum(TE.REL.VAR[TE.COEFS.INDEX.B2][TE.B2.EPISTASIS.LIST$PPN])
TE.B2.EPISTASIS.REL.VAR[2,1] <- sum(TE.REL.VAR[TE.COEFS.INDEX.B2][TE.B2.EPISTASIS.LIST$PNP])
TE.B2.EPISTASIS.REL.VAR[2,2] <- sum(TE.REL.VAR[TE.COEFS.INDEX.B2][TE.B2.EPISTASIS.LIST$PNN])
TE.B2.EPISTASIS.REL.VAR[3,1] <- sum(TE.REL.VAR[TE.COEFS.INDEX.B2][TE.B2.EPISTASIS.LIST$NNP])
TE.B2.EPISTASIS.REL.VAR[3,2] <- sum(TE.REL.VAR[TE.COEFS.INDEX.B2][TE.B2.EPISTASIS.LIST$NNN])

###Second order specificity interactions
TE.S2.EPISTASIS <- data.frame(matrix(0,nrow=2400,ncol=9))
colnames(TE.S2.EPISTASIS) <- c("A","B","AB","Expected","Observed","Sign.A","Sign.B","Sign.A.AB","Sign.B.AB")
rownames(TE.S2.EPISTASIS) <- names(TE.COEFS.ADJ[TE.COEFS.INDEX.S2])
##Get first amino acid state effect by itself
TE.S2.EPISTASIS$A <- c(rep(TE.COEFS.ADJ[TE.COEFS.INDEX.S1.X1],60),rep(TE.COEFS.ADJ[TE.COEFS.INDEX.S1.X2],40),rep(TE.COEFS.ADJ[TE.COEFS.INDEX.S1.X3],20))
##Get second amino acid state effect by itself
TE.S2.EPISTASIS$B <- c(rep(TE.COEFS.ADJ[TE.COEFS.INDEX.S1.X2],each=20),rep(TE.COEFS.ADJ[TE.COEFS.INDEX.S1.X3],each=20),rep(TE.COEFS.ADJ[TE.COEFS.INDEX.S1.X4],each=20),rep(TE.COEFS.ADJ[TE.COEFS.INDEX.S1.X3],each=20),rep(TE.COEFS.ADJ[TE.COEFS.INDEX.S1.X4],each=20),rep(TE.COEFS.ADJ[TE.COEFS.INDEX.S1.X4],each=20))
##Get pairwise epistastic effect between amino acid states
TE.S2.EPISTASIS$AB <- TE.COEFS.ADJ[TE.COEFS.INDEX.S2]
##Expected phenotype in absence of epistasis (A+B)
TE.S2.EPISTASIS$Expected <- TE.S2.EPISTASIS[,1] + TE.S2.EPISTASIS[,2]
##Observed phenotype in presence of epistasis (A+B+AB)
TE.S2.EPISTASIS$Observed <- TE.S2.EPISTASIS[,1] + TE.S2.EPISTASIS[,2] + TE.S2.EPISTASIS[,3]
##Direction of amino acid A by itself
TE.S2.EPISTASIS$Sign.A <- sign(TE.S2.EPISTASIS[,1])
##Direction of amino acid B by itself
TE.S2.EPISTASIS$Sign.B <- sign(TE.S2.EPISTASIS[,2])
##Direction of amino acid A in presence of B
TE.S2.EPISTASIS$Sign.A.AB <- sign(TE.S2.EPISTASIS[,1] + TE.S2.EPISTASIS[,3])
##Direction of amino acid B in presence of A
TE.S2.EPISTASIS$Sign.B.AB <- sign(TE.S2.EPISTASIS[,2] + TE.S2.EPISTASIS[,3])
##Flip all signs so that an increase in function is positive
TE.S2.EPISTASIS <- -1*TE.S2.EPISTASIS

###Third order binding interactions
TE.B3.EPISTASIS <- data.frame(matrix(0,nrow=3*length(TE.COEFS.ADJ[TE.COEFS.INDEX.B3]),ncol=9))
colnames(TE.B3.EPISTASIS) <- c("A","B","AB","Expected","Observed","Sign.A","Sign.B","Sign.A.AB","Sign.B.AB")
rownames(TE.B3.EPISTASIS) <- paste0(rep(names(TE.COEFS.ADJ[TE.COEFS.INDEX.B3]),each=3),rep(c(".E1",".E2",".E3"),length(TE.COEFS.ADJ[TE.COEFS.INDEX.B3])))
##Get first amino acid state effect by itself
TE.B3.EPISTASIS[seq(1,96000,3),1] <- c(rep(TE.COEFS.ADJ[TE.COEFS.INDEX.B1.X1],1200),rep(TE.COEFS.ADJ[TE.COEFS.INDEX.B1.X2],400))
##Get second amino acid state effect by itself
TE.B3.EPISTASIS[seq(2,96000,3),1] <- c(rep(rep(TE.COEFS.ADJ[TE.COEFS.INDEX.B1.X2],each=20),40),rep(rep(TE.COEFS.ADJ[TE.COEFS.INDEX.B1.X3],each=20),40))
##Get third amino acid state effect by itself
TE.B3.EPISTASIS[seq(3,96000,3),1] <- c(rep(TE.COEFS.ADJ[TE.COEFS.INDEX.B1.X3],each=400),rep(rep(TE.COEFS.ADJ[TE.COEFS.INDEX.B1.X4],each=400),3))
##Get pairwise epistastic effect between first and second amino acid states
TE.B3.EPISTASIS[seq(1,96000,3),2] <- c(rep(TE.COEFS.ADJ[TE.COEFS.INDEX.B2.X2X3],each=20),rep(TE.COEFS.ADJ[TE.COEFS.INDEX.B2.X2X4],each=20),rep(TE.COEFS.ADJ[TE.COEFS.INDEX.B2.X3X4],each=20),rep(TE.COEFS.ADJ[TE.COEFS.INDEX.B2.X3X4],each=20))
##Get pairwise epistastic effect between first and third amino acid states
TE.B3.EPISTASIS[seq(2,96000,3),2] <- c(rep(TE.COEFS.ADJ[TE.COEFS.INDEX.B2.X1X3],20),rep(TE.COEFS.ADJ[TE.COEFS.INDEX.B2.X1X4],40),rep(TE.COEFS.ADJ[TE.COEFS.INDEX.B2.X2X4],20))
##Get pairwise epistastic effect between second and third amino acid states
TE.B3.EPISTASIS[seq(3,96000,3),2] <- c(rep(TE.COEFS.ADJ[TE.COEFS.INDEX.B2.X1X2],40),rep(TE.COEFS.ADJ[TE.COEFS.INDEX.B2.X1X3],20),rep(TE.COEFS.ADJ[TE.COEFS.INDEX.B2.X2X3],20))
##Get third order epistatic effect between amino acids
TE.B3.EPISTASIS[,3] <- rep(TE.COEFS.ADJ[TE.COEFS.INDEX.B3],each=3)
##Expected phenotype in absence of epistasis (A+B)
TE.B3.EPISTASIS[,4] <- TE.B3.EPISTASIS[,1] + TE.B3.EPISTASIS[,2]
##Observed phenotype in presence of epistasis (A+B+AB)
TE.B3.EPISTASIS[,5] <- TE.B3.EPISTASIS[,1] + TE.B3.EPISTASIS[,2] + TE.B3.EPISTASIS[,3]
##Direction of amino acid A by itself
TE.B3.EPISTASIS[,6] <- sign(TE.B3.EPISTASIS[,1])
##Direction of amino acid B by itself
TE.B3.EPISTASIS[,7] <- sign(TE.B3.EPISTASIS[,2])
##Direction of amino acid A in presence of B
TE.B3.EPISTASIS[,8] <- sign(TE.B3.EPISTASIS[,1] + TE.B3.EPISTASIS[,3])
##Direction of amino acid B in presence of A
TE.B3.EPISTASIS[,9] <- sign(TE.B3.EPISTASIS[,2] + TE.B3.EPISTASIS[,3])
##Flip all signs so that an increase in function is positive
TE.B3.EPISTASIS <- -1*TE.B3.EPISTASIS

###Classify
####Classify epistasis
TE.B3.EPISTASIS.LIST <- list()
TE.B3.EPISTASIS.LIST$PPP  <- which(  TE.B3.EPISTASIS$Sign.A > 0 & TE.B3.EPISTASIS$Sign.B > 0 & TE.B3.EPISTASIS$Observed > TE.B3.EPISTASIS$Expected)
TE.B3.EPISTASIS.LIST$PPN  <- which(  TE.B3.EPISTASIS$Sign.A > 0 & TE.B3.EPISTASIS$Sign.B > 0 & TE.B3.EPISTASIS$Observed < TE.B3.EPISTASIS$Expected)
TE.B3.EPISTASIS.LIST$PNP  <- which(((TE.B3.EPISTASIS$Sign.A > 0 & TE.B3.EPISTASIS$Sign.B < 0) | (TE.B3.EPISTASIS$Sign.A < 0 & TE.B3.EPISTASIS$Sign.B > 0)) & TE.B3.EPISTASIS$Observed > TE.B3.EPISTASIS$Expected)
TE.B3.EPISTASIS.LIST$PNN  <- which(((TE.B3.EPISTASIS$Sign.A > 0 & TE.B3.EPISTASIS$Sign.B < 0) | (TE.B3.EPISTASIS$Sign.A < 0 & TE.B3.EPISTASIS$Sign.B > 0)) & TE.B3.EPISTASIS$Observed < TE.B3.EPISTASIS$Expected)
TE.B3.EPISTASIS.LIST$NNP  <- which(  TE.B3.EPISTASIS$Sign.A < 0 & TE.B3.EPISTASIS$Sign.B < 0 & TE.B3.EPISTASIS$Observed > TE.B3.EPISTASIS$Expected)
TE.B3.EPISTASIS.LIST$NNN  <- which(  TE.B3.EPISTASIS$Sign.A < 0 & TE.B3.EPISTASIS$Sign.B < 0 & TE.B3.EPISTASIS$Observed < TE.B3.EPISTASIS$Expected)

TE.B3.EPISTASIS.REL.VAR <- matrix(0,nrow=3,ncol=2);rownames(TE.B3.EPISTASIS.REL.VAR) <- c("PP","PN","NN"); colnames(TE.B3.EPISTASIS.REL.VAR) <- c("+","-")
TE.B3.EPISTASIS.REL.VAR[1,1] <- sum(rep(TE.REL.VAR[TE.COEFS.INDEX.B3],each=3)[TE.B3.EPISTASIS.LIST$PPP])/3
TE.B3.EPISTASIS.REL.VAR[1,2] <- sum(rep(TE.REL.VAR[TE.COEFS.INDEX.B3],each=3)[TE.B3.EPISTASIS.LIST$PPN])/3
TE.B3.EPISTASIS.REL.VAR[2,1] <- sum(rep(TE.REL.VAR[TE.COEFS.INDEX.B3],each=3)[TE.B3.EPISTASIS.LIST$PNP])/3
TE.B3.EPISTASIS.REL.VAR[2,2] <- sum(rep(TE.REL.VAR[TE.COEFS.INDEX.B3],each=3)[TE.B3.EPISTASIS.LIST$PNN])/3
TE.B3.EPISTASIS.REL.VAR[3,1] <- sum(rep(TE.REL.VAR[TE.COEFS.INDEX.B3],each=3)[TE.B3.EPISTASIS.LIST$NNP])/3
TE.B3.EPISTASIS.REL.VAR[3,2] <- sum(rep(TE.REL.VAR[TE.COEFS.INDEX.B3],each=3)[TE.B3.EPISTASIS.LIST$NNN])/3

###Third order specificity interactions
TE.S3.EPISTASIS <- data.frame(matrix(0,nrow=3*length(TE.COEFS.ADJ[TE.COEFS.INDEX.S3]),ncol=9))
colnames(TE.S3.EPISTASIS) <- c("A","B","AB","Expected","Observed","Sign.A","Sign.B","Sign.A.AB","Sign.B.AB")
rownames(TE.S3.EPISTASIS) <- paste0(rep(names(TE.COEFS.ADJ[TE.COEFS.INDEX.S3]),each=3),rep(c(".E1",".E2",".E3"),length(TE.COEFS.ADJ[TE.COEFS.INDEX.S3])))
##Get first amino acid state effect by itself
TE.S3.EPISTASIS[seq(1,96000,3),1] <- c(rep(TE.COEFS.ADJ[TE.COEFS.INDEX.S1.X1],1200),rep(TE.COEFS.ADJ[TE.COEFS.INDEX.S1.X2],400))
##Get second amino acid state effect by itself
TE.S3.EPISTASIS[seq(2,96000,3),1] <- c(rep(rep(TE.COEFS.ADJ[TE.COEFS.INDEX.S1.X2],each=20),40),rep(rep(TE.COEFS.ADJ[TE.COEFS.INDEX.S1.X3],each=20),40))
##Get third amino acid state effect by itself
TE.S3.EPISTASIS[seq(3,96000,3),1] <- c(rep(TE.COEFS.ADJ[TE.COEFS.INDEX.S1.X3],each=400),rep(rep(TE.COEFS.ADJ[TE.COEFS.INDEX.S1.X4],each=400),3))
##Get pairwise epistastic effect between first and second amino acid states
TE.S3.EPISTASIS[seq(1,96000,3),2] <- c(rep(TE.COEFS.ADJ[TE.COEFS.INDEX.S2.X2X3],each=20),rep(TE.COEFS.ADJ[TE.COEFS.INDEX.S2.X2X4],each=20),rep(TE.COEFS.ADJ[TE.COEFS.INDEX.S2.X3X4],each=20),rep(TE.COEFS.ADJ[TE.COEFS.INDEX.S2.X3X4],each=20))
##Get pairwise epistastic effect between first and third amino acid states
TE.S3.EPISTASIS[seq(2,96000,3),2] <- c(rep(TE.COEFS.ADJ[TE.COEFS.INDEX.S2.X1X3],20),rep(TE.COEFS.ADJ[TE.COEFS.INDEX.S2.X1X4],40),rep(TE.COEFS.ADJ[TE.COEFS.INDEX.S2.X2X4],20))
##Get pairwise epistastic effect between second and third amino acid states
TE.S3.EPISTASIS[seq(3,96000,3),2] <- c(rep(TE.COEFS.ADJ[TE.COEFS.INDEX.S2.X1X2],40),rep(TE.COEFS.ADJ[TE.COEFS.INDEX.S2.X1X3],20),rep(TE.COEFS.ADJ[TE.COEFS.INDEX.S2.X2X3],20))
##Get third order epistatic effect between amino acids
TE.S3.EPISTASIS[,3] <- rep(TE.COEFS.ADJ[TE.COEFS.INDEX.S3],each=3)
##Expected phenotype in absence of epistasis (A+B)
TE.S3.EPISTASIS[,4] <- TE.S3.EPISTASIS[,1] + TE.S3.EPISTASIS[,2]
##Observed phenotype in presence of epistasis (A+B+AB)
TE.S3.EPISTASIS[,5] <- TE.S3.EPISTASIS[,1] + TE.S3.EPISTASIS[,2] + TE.S3.EPISTASIS[,3]
##Direction of amino acid A by itself
TE.S3.EPISTASIS[,6] <- sign(TE.S3.EPISTASIS[,1])
##Direction of amino acid B by itself
TE.S3.EPISTASIS[,7] <- sign(TE.S3.EPISTASIS[,2])
##Direction of amino acid A in presence of B
TE.S3.EPISTASIS[,8] <- sign(TE.S3.EPISTASIS[,1] + TE.S3.EPISTASIS[,3])
##Direction of amino acid B in presence of A
TE.S3.EPISTASIS[,9] <- sign(TE.S3.EPISTASIS[,2] + TE.S3.EPISTASIS[,3])
##Flip all signs so that an increase in function is positive
TE.S3.EPISTASIS <- -1*TE.S3.EPISTASIS

##pdf("EPISTASIS.TYPES.pdf")
par(pty="s")
par(mfrow=c(1,2))
barplot(t(TE.B2.EPISTASIS.REL.VAR))
barplot(t(TE.B3.EPISTASIS.REL.VAR))
##dev.off()

################################
##18. Heatmaps of coefficients##
################################
####ME Model
#Rescale to be in units towards the weak/strong intercept
ME.COEFS.ADJ.SCALED <- ME.COEFS.ADJ/(-1*ME.THRESH.WEAK + ME.COEFS.EFFECT.ADJ.B0)

#Censor values
ME.COEFS.ADJ.SCALED[which(ME.COEFS.ADJ.SCALED >  0.000 & ME.COEFS.ADJ.SCALED <  0.025)] <-  0.000
ME.COEFS.ADJ.SCALED[which(ME.COEFS.ADJ.SCALED >  0.025 & ME.COEFS.ADJ.SCALED <  0.050)] <-  0.025
ME.COEFS.ADJ.SCALED[which(ME.COEFS.ADJ.SCALED >  0.050 & ME.COEFS.ADJ.SCALED <  0.075)] <-  0.050
ME.COEFS.ADJ.SCALED[which(ME.COEFS.ADJ.SCALED >  0.075 & ME.COEFS.ADJ.SCALED <  0.100)] <-  0.075
ME.COEFS.ADJ.SCALED[which(ME.COEFS.ADJ.SCALED >  0.100 & ME.COEFS.ADJ.SCALED <  0.125)] <-  0.100
ME.COEFS.ADJ.SCALED[which(ME.COEFS.ADJ.SCALED >  0.125 & ME.COEFS.ADJ.SCALED <  0.150)] <-  0.125
ME.COEFS.ADJ.SCALED[which(ME.COEFS.ADJ.SCALED >  0.150 & ME.COEFS.ADJ.SCALED <  0.175)] <-  0.150
ME.COEFS.ADJ.SCALED[which(ME.COEFS.ADJ.SCALED >  0.175 & ME.COEFS.ADJ.SCALED <  0.200)] <-  0.175
ME.COEFS.ADJ.SCALED[which(ME.COEFS.ADJ.SCALED >  0.200 & ME.COEFS.ADJ.SCALED <  0.225)] <-  0.200
ME.COEFS.ADJ.SCALED[which(ME.COEFS.ADJ.SCALED >  0.225 & ME.COEFS.ADJ.SCALED <  0.250)] <-  0.225
ME.COEFS.ADJ.SCALED[which(ME.COEFS.ADJ.SCALED >  0.250)] <- 0.250
ME.COEFS.ADJ.SCALED[which(ME.COEFS.ADJ.SCALED <  0.000 & ME.COEFS.ADJ.SCALED > -0.025)] <-  0.000
ME.COEFS.ADJ.SCALED[which(ME.COEFS.ADJ.SCALED < -0.025 & ME.COEFS.ADJ.SCALED > -0.050)] <- -0.025
ME.COEFS.ADJ.SCALED[which(ME.COEFS.ADJ.SCALED < -0.050 & ME.COEFS.ADJ.SCALED > -0.075)] <- -0.050
ME.COEFS.ADJ.SCALED[which(ME.COEFS.ADJ.SCALED < -0.075 & ME.COEFS.ADJ.SCALED > -0.100)] <- -0.075
ME.COEFS.ADJ.SCALED[which(ME.COEFS.ADJ.SCALED < -0.100 & ME.COEFS.ADJ.SCALED > -0.125)] <- -0.100
ME.COEFS.ADJ.SCALED[which(ME.COEFS.ADJ.SCALED < -0.125 & ME.COEFS.ADJ.SCALED > -0.150)] <- -0.125
ME.COEFS.ADJ.SCALED[which(ME.COEFS.ADJ.SCALED < -0.150 & ME.COEFS.ADJ.SCALED > -0.175)] <- -0.150
ME.COEFS.ADJ.SCALED[which(ME.COEFS.ADJ.SCALED < -0.175 & ME.COEFS.ADJ.SCALED > -0.200)] <- -0.175
ME.COEFS.ADJ.SCALED[which(ME.COEFS.ADJ.SCALED < -0.200 & ME.COEFS.ADJ.SCALED > -0.225)] <- -0.200
ME.COEFS.ADJ.SCALED[which(ME.COEFS.ADJ.SCALED < -0.225 & ME.COEFS.ADJ.SCALED > -0.250)] <- -0.225
ME.COEFS.ADJ.SCALED[which(ME.COEFS.ADJ.SCALED < -0.250)] <- -0.250

##Binding effects
ME.MAIN.BIND.COEFS.df <- data.frame(state=as.character(c(AAs,AAs,AAs,AAs)),site=factor(c(rep(1,length(AAs)),rep(2,length(AAs)),rep(3,length(AAs)),rep(4,length(AAs))),levels=c("4","3","2","1")),value=-1*(ME.COEFS.ADJ.SCALED[ME.COEFS.INDEX.B1]),row.names=NULL)

#Main binding effects for position 1
ME.BIND.M1 <- ggplot(ME.MAIN.BIND.COEFS.df[ME.MAIN.BIND.COEFS.df$site==1,], aes(state, site))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+coord_equal()+scale_x_discrete(position = "top") + theme(legend.position = "none",panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
#Main binding effects for position 2
ME.BIND.M2 <- ggplot(ME.MAIN.BIND.COEFS.df[ME.MAIN.BIND.COEFS.df$site==2,], aes(state, site))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+coord_equal()+scale_x_discrete(position = "top") + theme(legend.position = "none",panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
#Main binding effects for position 3
ME.BIND.M3 <- ggplot(ME.MAIN.BIND.COEFS.df[ME.MAIN.BIND.COEFS.df$site==3,], aes(state, site))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+coord_equal()+scale_x_discrete(position = "top") + theme(legend.position = "none",panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
#Main binding effects for position 4
ME.BIND.M4 <- ggplot(ME.MAIN.BIND.COEFS.df[ME.MAIN.BIND.COEFS.df$site==4,], aes(state, site))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+coord_equal()+scale_x_discrete(position = "top") + theme(legend.position = "none",panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

#Specificity effects
ME.MAIN.SPEC.COEFS.df <- data.frame(state=as.character(c(AAs,AAs,AAs,AAs)),site=factor(c(rep(1,length(AAs)),rep(2,length(AAs)),rep(3,length(AAs)),rep(4,length(AAs))),levels=c("4","3","2","1")),value=-1*(ME.COEFS.ADJ.SCALED[ME.COEFS.INDEX.S1]),row.names=NULL)

#Main specificity effects for position 1
ME.SPEC.M1 <- ggplot(ME.MAIN.SPEC.COEFS.df[ME.MAIN.SPEC.COEFS.df$site==1,], aes(state, site))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#00B300","#96FF63","#FFFFFF","#DB94FF","#7A00CC"),limits=c(-0.25,0.25)) +labs(x="",y="")+
  theme_classic()+ theme(legend.position = "none", panel.grid = element_blank(), axis.text.x = element_blank()) +
  coord_flip() + theme(aspect.ratio=18) + scale_x_discrete(limits = rev(levels(ME.MAIN.SPEC.COEFS.df$state)))
#Main specificity effects for position 2
ME.SPEC.M2 <- ggplot(ME.MAIN.SPEC.COEFS.df[ME.MAIN.SPEC.COEFS.df$site==2,], aes(state, site))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#00B300","#96FF63","#FFFFFF","#DB94FF","#7A00CC"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+ theme(legend.position = "none", panel.grid = element_blank(), axis.text.x = element_blank()) +
  coord_flip() + theme(aspect.ratio=18) + scale_x_discrete(limits = rev(levels(ME.MAIN.SPEC.COEFS.df$state)))
#Main specificity effects for position 3
ME.SPEC.M3 <- ggplot(ME.MAIN.SPEC.COEFS.df[ME.MAIN.SPEC.COEFS.df$site==3,], aes(state, site))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#00B300","#96FF63","#FFFFFF","#DB94FF","#7A00CC"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+ theme(legend.position = "none", panel.grid = element_blank(), axis.text.x = element_blank()) +
  coord_flip() + theme(aspect.ratio=18) + scale_x_discrete(limits = rev(levels(ME.MAIN.SPEC.COEFS.df$state)))
#Main specificity effects for position 4
ME.SPEC.M4 <- ggplot(ME.MAIN.SPEC.COEFS.df[ME.MAIN.SPEC.COEFS.df$site==4,], aes(state, site))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#00B300","#96FF63","#FFFFFF","#DB94FF","#7A00CC"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+ theme(legend.position = "none", panel.grid = element_blank(), axis.text.x = element_blank()) +
  coord_flip() + theme(aspect.ratio=18) + scale_x_discrete(limits = rev(levels(ME.MAIN.SPEC.COEFS.df$state)))

##Main binding effects across positions
##pdf("ME.BIND.MAIN.pdf",5,3)
#margin = theme(plot.margin = unit(c(0,0,0,0), "cm"))
#pl <- list(ME.BIND.M1,ME.BIND.M2,ME.BIND.M3,ME.BIND.M4)
#grid.arrange(grobs = lapply(pl, "+", margin),layout_matrix = rbind(c(1,2,3,4)))
##dev.off()

##Main specificty effects across positions
##pdf("ME.SPEC.MAIN.pdf",5,3)
#margin = theme(plot.margin = unit(c(0,0,0,0), "cm"))
#pl <- list(ME.SPEC.M1,ME.SPEC.M2,ME.SPEC.M3,ME.SPEC.M4)
#grid.arrange(grobs = lapply(pl, "+", margin),layout_matrix = cbind(c(1,2,3,4)))
##dev.off()

###PE Model
##Rescale to be in units towards the weak/strong intercept
PE.COEFS.ADJ.SCALED <- PE.COEFS.ADJ/(-1*PE.THRESH.WEAK+PE.COEFS.EFFECT.ADJ.B0)

#Censor values
PE.COEFS.ADJ.SCALED[which(PE.COEFS.ADJ.SCALED >  0.000 & PE.COEFS.ADJ.SCALED <  0.025)] <-  0.000
PE.COEFS.ADJ.SCALED[which(PE.COEFS.ADJ.SCALED >  0.025 & PE.COEFS.ADJ.SCALED <  0.050)] <-  0.025
PE.COEFS.ADJ.SCALED[which(PE.COEFS.ADJ.SCALED >  0.050 & PE.COEFS.ADJ.SCALED <  0.075)] <-  0.050
PE.COEFS.ADJ.SCALED[which(PE.COEFS.ADJ.SCALED >  0.075 & PE.COEFS.ADJ.SCALED <  0.100)] <-  0.075
PE.COEFS.ADJ.SCALED[which(PE.COEFS.ADJ.SCALED >  0.100 & PE.COEFS.ADJ.SCALED <  0.125)] <-  0.100
PE.COEFS.ADJ.SCALED[which(PE.COEFS.ADJ.SCALED >  0.125 & PE.COEFS.ADJ.SCALED <  0.150)] <-  0.125
PE.COEFS.ADJ.SCALED[which(PE.COEFS.ADJ.SCALED >  0.150 & PE.COEFS.ADJ.SCALED <  0.175)] <-  0.150
PE.COEFS.ADJ.SCALED[which(PE.COEFS.ADJ.SCALED >  0.175 & PE.COEFS.ADJ.SCALED <  0.200)] <-  0.175
PE.COEFS.ADJ.SCALED[which(PE.COEFS.ADJ.SCALED >  0.200 & PE.COEFS.ADJ.SCALED <  0.225)] <-  0.200
PE.COEFS.ADJ.SCALED[which(PE.COEFS.ADJ.SCALED >  0.225 & PE.COEFS.ADJ.SCALED <  0.250)] <-  0.225
PE.COEFS.ADJ.SCALED[which(PE.COEFS.ADJ.SCALED >  0.250)] <- 0.250
PE.COEFS.ADJ.SCALED[which(PE.COEFS.ADJ.SCALED <  0.000 & PE.COEFS.ADJ.SCALED > -0.025)] <-  0.000
PE.COEFS.ADJ.SCALED[which(PE.COEFS.ADJ.SCALED < -0.025 & PE.COEFS.ADJ.SCALED > -0.050)] <- -0.025
PE.COEFS.ADJ.SCALED[which(PE.COEFS.ADJ.SCALED < -0.050 & PE.COEFS.ADJ.SCALED > -0.075)] <- -0.050
PE.COEFS.ADJ.SCALED[which(PE.COEFS.ADJ.SCALED < -0.075 & PE.COEFS.ADJ.SCALED > -0.100)] <- -0.075
PE.COEFS.ADJ.SCALED[which(PE.COEFS.ADJ.SCALED < -0.100 & PE.COEFS.ADJ.SCALED > -0.125)] <- -0.100
PE.COEFS.ADJ.SCALED[which(PE.COEFS.ADJ.SCALED < -0.125 & PE.COEFS.ADJ.SCALED > -0.150)] <- -0.125
PE.COEFS.ADJ.SCALED[which(PE.COEFS.ADJ.SCALED < -0.150 & PE.COEFS.ADJ.SCALED > -0.175)] <- -0.150
PE.COEFS.ADJ.SCALED[which(PE.COEFS.ADJ.SCALED < -0.175 & PE.COEFS.ADJ.SCALED > -0.200)] <- -0.175
PE.COEFS.ADJ.SCALED[which(PE.COEFS.ADJ.SCALED < -0.200 & PE.COEFS.ADJ.SCALED > -0.225)] <- -0.200
PE.COEFS.ADJ.SCALED[which(PE.COEFS.ADJ.SCALED < -0.225 & PE.COEFS.ADJ.SCALED > -0.250)] <- -0.225
PE.COEFS.ADJ.SCALED[which(PE.COEFS.ADJ.SCALED < -0.250)] <- -0.250

##Binding effects
PE.MAIN.BIND.COEFS.df <- data.frame(state=as.character(c(AAs,AAs,AAs,AAs)),site=factor(c(rep(1,length(AAs)),rep(2,length(AAs)),rep(3,length(AAs)),rep(4,length(AAs))),levels=c("4","3","2","1")),value=-1*(PE.COEFS.ADJ.SCALED[PE.COEFS.INDEX.B1]),row.names=NULL)
PE.EP12.BIND.COEFS.df <- data.frame(state1=factor(rep(AAs,length(AAs)),levels=rev(AAs)), state2=factor(as.vector(sapply(1:length(AAs), function(x) rep(AAs[x],length(AAs))))),value= -1*(PE.COEFS.ADJ.SCALED[PE.COEFS.INDEX.B2.X1X2]))
PE.EP13.BIND.COEFS.df <- data.frame(state1=factor(rep(AAs,length(AAs)),levels=rev(AAs)), state2=factor(as.vector(sapply(1:length(AAs), function(x) rep(AAs[x],length(AAs))))),value= -1*(PE.COEFS.ADJ.SCALED[PE.COEFS.INDEX.B2.X1X3]))
PE.EP14.BIND.COEFS.df <- data.frame(state1=factor(rep(AAs,length(AAs)),levels=rev(AAs)), state2=factor(as.vector(sapply(1:length(AAs), function(x) rep(AAs[x],length(AAs))))),value= -1*(PE.COEFS.ADJ.SCALED[PE.COEFS.INDEX.B2.X1X4]))
PE.EP23.BIND.COEFS.df <- data.frame(state1=factor(rep(AAs,length(AAs)),levels=rev(AAs)), state2=factor(as.vector(sapply(1:length(AAs), function(x) rep(AAs[x],length(AAs))))),value= -1*(PE.COEFS.ADJ.SCALED[PE.COEFS.INDEX.B2.X2X3]))
PE.EP24.BIND.COEFS.df <- data.frame(state1=factor(rep(AAs,length(AAs)),levels=rev(AAs)), state2=factor(as.vector(sapply(1:length(AAs), function(x) rep(AAs[x],length(AAs))))),value= -1*(PE.COEFS.ADJ.SCALED[PE.COEFS.INDEX.B2.X2X4]))
PE.EP34.BIND.COEFS.df <- data.frame(state1=factor(rep(AAs,length(AAs)),levels=rev(AAs)), state2=factor(as.vector(sapply(1:length(AAs), function(x) rep(AAs[x],length(AAs))))),value= -1*(PE.COEFS.ADJ.SCALED[PE.COEFS.INDEX.B2.X3X4]))

##Plots of binding effects
#Main binding effects for position 1
PE.BIND.M1 <- ggplot(PE.MAIN.BIND.COEFS.df[PE.MAIN.BIND.COEFS.df$site==1,], aes(state, site))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+coord_equal()+scale_x_discrete(position = "top") + theme(legend.position = "none",panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
#Main binding effects for position 2
PE.BIND.M2 <- ggplot(PE.MAIN.BIND.COEFS.df[PE.MAIN.BIND.COEFS.df$site==2,], aes(state, site))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+coord_equal()+scale_x_discrete(position = "top") + theme(legend.position = "none",panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
#Main binding effects for position 3
PE.BIND.M3 <- ggplot(PE.MAIN.BIND.COEFS.df[PE.MAIN.BIND.COEFS.df$site==3,], aes(state, site))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+coord_equal()+scale_x_discrete(position = "top") + theme(legend.position = "none",panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
#Main binding effects for position 4
PE.BIND.M4 <- ggplot(PE.MAIN.BIND.COEFS.df[PE.MAIN.BIND.COEFS.df$site==4,], aes(state, site))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+coord_equal()+scale_x_discrete(position = "top") + theme(legend.position = "none",panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
#Pairwise binding interaction effects between positions 1 and 2
PE.BIND.P12 <- ggplot(PE.EP12.BIND.COEFS.df ,aes(state2,state1))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+coord_equal()+ theme(legend.position = "none",panel.grid = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank()) + 
  xlim(rev(levels(PE.EP12.BIND.COEFS.df$state1))) + ylim(rev(levels(PE.EP12.BIND.COEFS.df$state2)))
#Pairwise binding interaction effects between positions 1 and 3
PE.BIND.P13 <- ggplot(PE.EP13.BIND.COEFS.df ,aes(state2,state1))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+coord_equal()+ theme(legend.position = "none",panel.grid = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank()) + 
  xlim(rev(levels(PE.EP13.BIND.COEFS.df$state1))) + ylim(rev(levels(PE.EP13.BIND.COEFS.df$state2)))
#Pairwise binding interaction effects between positions 1 and 4
PE.BIND.P14 <- ggplot(PE.EP14.BIND.COEFS.df ,aes(state2,state1))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+coord_equal()+ theme(legend.position = "none",panel.grid = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank()) + 
  xlim(rev(levels(PE.EP14.BIND.COEFS.df$state1))) + ylim(rev(levels(PE.EP14.BIND.COEFS.df$state2)))
#Pairwise binding interaction effects between positions 2 and 3
PE.BIND.P23 <- ggplot(PE.EP23.BIND.COEFS.df ,aes(state2,state1))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+coord_equal()+ theme(legend.position = "none",panel.grid = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank()) + 
  xlim(rev(levels(PE.EP24.BIND.COEFS.df$state1))) + ylim(rev(levels(PE.EP24.BIND.COEFS.df$state2)))
#Pairwise binding interaction effects between positions 2 and 4
PE.BIND.P24 <- ggplot(PE.EP24.BIND.COEFS.df ,aes(state2,state1))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+coord_equal()+ theme(legend.position = "none",panel.grid = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank()) + 
  xlim(rev(levels(PE.EP24.BIND.COEFS.df$state1))) + ylim(rev(levels(PE.EP24.BIND.COEFS.df$state2)))
#Pairwise binding interaction effects between positions 3 and 4
PE.BIND.P34 <- ggplot(PE.EP34.BIND.COEFS.df ,aes(state2,state1))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+coord_equal()+ theme(legend.position = "none",panel.grid = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank()) + 
  xlim(rev(levels(PE.EP34.BIND.COEFS.df$state1))) + ylim(rev(levels(PE.EP34.BIND.COEFS.df$state2)))

##Specificity effects
PE.MAIN.SPEC.COEFS.df <- data.frame(state=as.character(c(AAs,AAs,AAs,AAs)),site=factor(c(rep(1,length(AAs)),rep(2,length(AAs)),rep(3,length(AAs)),rep(4,length(AAs))),levels=c("4","3","2","1")),value=-1*(PE.COEFS.ADJ.SCALED[PE.COEFS.INDEX.S1]),row.names=NULL)
PE.EP12.SPEC.COEFS.df <- data.frame(state1=factor(rep(AAs,length(AAs)),levels=rev(AAs)), state2=factor(as.vector(sapply(1:length(AAs), function(x) rep(AAs[x],length(AAs))))),value= -1*(PE.COEFS.ADJ.SCALED[PE.COEFS.INDEX.S2.X1X2]))
PE.EP13.SPEC.COEFS.df <- data.frame(state1=factor(rep(AAs,length(AAs)),levels=rev(AAs)), state2=factor(as.vector(sapply(1:length(AAs), function(x) rep(AAs[x],length(AAs))))),value= -1*(PE.COEFS.ADJ.SCALED[PE.COEFS.INDEX.S2.X1X3]))
PE.EP14.SPEC.COEFS.df <- data.frame(state1=factor(rep(AAs,length(AAs)),levels=rev(AAs)), state2=factor(as.vector(sapply(1:length(AAs), function(x) rep(AAs[x],length(AAs))))),value= -1*(PE.COEFS.ADJ.SCALED[PE.COEFS.INDEX.S2.X1X4]))
PE.EP23.SPEC.COEFS.df <- data.frame(state1=factor(rep(AAs,length(AAs)),levels=rev(AAs)), state2=factor(as.vector(sapply(1:length(AAs), function(x) rep(AAs[x],length(AAs))))),value= -1*(PE.COEFS.ADJ.SCALED[PE.COEFS.INDEX.S2.X2X3]))
PE.EP24.SPEC.COEFS.df <- data.frame(state1=factor(rep(AAs,length(AAs)),levels=rev(AAs)), state2=factor(as.vector(sapply(1:length(AAs), function(x) rep(AAs[x],length(AAs))))),value= -1*(PE.COEFS.ADJ.SCALED[PE.COEFS.INDEX.S2.X2X4]))
PE.EP34.SPEC.COEFS.df <- data.frame(state1=factor(rep(AAs,length(AAs)),levels=rev(AAs)), state2=factor(as.vector(sapply(1:length(AAs), function(x) rep(AAs[x],length(AAs))))),value= -1*(PE.COEFS.ADJ.SCALED[PE.COEFS.INDEX.S2.X3X4]))

##Plots of specificty effects
#Main specificity effects for position 1
PE.SPEC.M1 <- ggplot(PE.MAIN.SPEC.COEFS.df[PE.MAIN.SPEC.COEFS.df$site==1,], aes(state, site))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#00B300","#96FF63","#FFFFFF","#DB94FF","#7A00CC"),limits=c(-0.25,0.25)) +labs(x="",y="")+
  theme_classic()+ theme(legend.position = "none", panel.grid = element_blank(), axis.text.x = element_blank()) +
  coord_flip() + theme(aspect.ratio=18) + scale_x_discrete(limits = rev(levels(PE.MAIN.SPEC.COEFS.df$state)))
#Main specificity effects for position 2
PE.SPEC.M2 <- ggplot(PE.MAIN.SPEC.COEFS.df[PE.MAIN.SPEC.COEFS.df$site==2,], aes(state, site))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#00B300","#96FF63","#FFFFFF","#DB94FF","#7A00CC"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+ theme(legend.position = "none", panel.grid = element_blank(), axis.text.x = element_blank()) +
  coord_flip() + theme(aspect.ratio=18) + scale_x_discrete(limits = rev(levels(PE.MAIN.SPEC.COEFS.df$state)))
#Main specificity effects for position 3
PE.SPEC.M3 <- ggplot(PE.MAIN.SPEC.COEFS.df[PE.MAIN.SPEC.COEFS.df$site==3,], aes(state, site))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#00B300","#96FF63","#FFFFFF","#DB94FF","#7A00CC"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+ theme(legend.position = "none", panel.grid = element_blank(), axis.text.x = element_blank()) +
  coord_flip() + theme(aspect.ratio=18) + scale_x_discrete(limits = rev(levels(PE.MAIN.SPEC.COEFS.df$state)))
#Main specificity effects for position 4
PE.SPEC.M4 <- ggplot(PE.MAIN.SPEC.COEFS.df[PE.MAIN.SPEC.COEFS.df$site==4,], aes(state, site))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#00B300","#96FF63","#FFFFFF","#DB94FF","#7A00CC"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+ theme(legend.position = "none", panel.grid = element_blank(), axis.text.x = element_blank()) +
  coord_flip() + theme(aspect.ratio=18) + scale_x_discrete(limits = rev(levels(PE.MAIN.SPEC.COEFS.df$state)))
#Pairwise specificity interaction effects between positions 1 and 2
PE.SPEC.P12 <- ggplot(PE.EP12.SPEC.COEFS.df ,aes(state1,state2))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#00B300","#96FF63","#FFFFFF","#DB94FF","#7A00CC"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+coord_equal() + theme(legend.position = "none",panel.grid = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank()) + 
  xlim(rev(levels(PE.EP12.SPEC.COEFS.df$state1))) + ylim(rev(levels(PE.EP12.SPEC.COEFS.df$state2)))
#Pairwise specificity interaction effects between positions 1 and 3
PE.SPEC.P13 <- ggplot(PE.EP13.SPEC.COEFS.df ,aes(state1,state2))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#00B300","#96FF63","#FFFFFF","#DB94FF","#7A00CC"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+coord_equal() + theme(legend.position = "none",panel.grid = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank()) + 
  xlim(rev(levels(PE.EP13.SPEC.COEFS.df$state1))) + ylim(rev(levels(PE.EP13.SPEC.COEFS.df$state2)))
#Pairwise specificity interaction effects between positions 1 and 4
PE.SPEC.P14 <- ggplot(PE.EP14.SPEC.COEFS.df ,aes(state1,state2))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#00B300","#96FF63","#FFFFFF","#DB94FF","#7A00CC"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+coord_equal() + theme(legend.position = "none",panel.grid = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank()) + 
  xlim(rev(levels(PE.EP14.SPEC.COEFS.df$state1))) + ylim(rev(levels(PE.EP14.SPEC.COEFS.df$state2)))
#Pairwise specificity interaction effects between positions 2 and 3
PE.SPEC.P23 <- ggplot(PE.EP23.SPEC.COEFS.df ,aes(state1,state2))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#00B300","#96FF63","#FFFFFF","#DB94FF","#7A00CC"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+coord_equal() + theme(legend.position = "none",panel.grid = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank()) + 
  xlim(rev(levels(PE.EP23.SPEC.COEFS.df$state1))) + ylim(rev(levels(PE.EP23.SPEC.COEFS.df$state2)))
#Pairwise specificity interaction effects between positions 2 and 4
PE.SPEC.P24 <- ggplot(PE.EP24.SPEC.COEFS.df ,aes(state1,state2))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#00B300","#96FF63","#FFFFFF","#DB94FF","#7A00CC"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+coord_equal() + theme(legend.position = "none",panel.grid = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank()) + 
  xlim(rev(levels(PE.EP24.SPEC.COEFS.df$state1))) + ylim(rev(levels(PE.EP24.SPEC.COEFS.df$state2)))
#Pairwise specificity interaction effects between positions 3 and 4
PE.SPEC.P34 <- ggplot(PE.EP34.SPEC.COEFS.df ,aes(state1,state2))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#00B300","#96FF63","#FFFFFF","#DB94FF","#7A00CC"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+coord_equal() + theme(legend.position = "none",panel.grid = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank()) + 
  xlim(rev(levels(PE.EP34.SPEC.COEFS.df$state1))) + ylim(rev(levels(PE.EP34.SPEC.COEFS.df$state2)))

##Combine binding and specificty interactions across positions into a single plot 
##pdf("PE.BIND.SPEC.PAIRWISE.pdf",5,3)
#margin = theme(plot.margin = unit(c(0,0,0,0), "cm"))
#pl <- list(PE.BIND.P12,PE.BIND.P13,PE.BIND.P14,PE.BIND.P23,PE.BIND.P24,PE.BIND.P34,
#           PE.SPEC.P12,PE.SPEC.P13,PE.SPEC.P14,PE.SPEC.P23,PE.SPEC.P24,PE.SPEC.P34)
#grid.arrange(grobs = lapply(pl, "+", margin),layout_matrix = rbind(c(NA,1,2,3),c(7,NA,4,5),c(8,10,NA,6),c(9,11,12,NA)))
##dev.off()  

##Main binding effects across positions
##pdf("PE.BIND.MAIN.pdf",5,3)
#margin = theme(plot.margin = unit(c(0,0,0,0), "cm"))
#pl <- list(PE.BIND.M1,PE.BIND.M2,PE.BIND.M3,PE.BIND.M4)
#grid.arrange(grobs = lapply(pl, "+", margin),layout_matrix = rbind(c(1,2,3,4)))
##dev.off()

##Main specificty effects across positions
##pdf("PE.SPEC.MAIN.pdf",5,3)
#margin = theme(plot.margin = unit(c(0,0,0,0), "cm"))
#pl <- list(PE.SPEC.M1,PE.SPEC.M2,PE.SPEC.M3,PE.SPEC.M4)
#grid.arrange(grobs = lapply(pl, "+", margin),layout_matrix = cbind(c(1,2,3,4)))
##dev.off()

###TE Model
##Amplify specificity effects for visualization
TE.COEFS.ADJ.AMP <- TE.COEFS.ADJ
TE.COEFS.ADJ.AMP[TE.COEFS.INDEX.S1] <- 2*TE.COEFS.ADJ.AMP[TE.COEFS.INDEX.S1] 
TE.COEFS.ADJ.AMP[TE.COEFS.INDEX.S2] <- 2*TE.COEFS.ADJ.AMP[TE.COEFS.INDEX.S2] 
TE.COEFS.ADJ.AMP[TE.COEFS.INDEX.S3] <- 2*TE.COEFS.ADJ.AMP[TE.COEFS.INDEX.S3]
#Rescale to be in units towards the weak/strong intercept
TE.COEFS.ADJ.SCALED <- TE.COEFS.ADJ.AMP/(-1*TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0)

#Censor values
TE.COEFS.ADJ.SCALED[which(TE.COEFS.ADJ.SCALED >  0.000 & TE.COEFS.ADJ.SCALED <  0.025)] <-  0.000
TE.COEFS.ADJ.SCALED[which(TE.COEFS.ADJ.SCALED >  0.025 & TE.COEFS.ADJ.SCALED <  0.050)] <-  0.025
TE.COEFS.ADJ.SCALED[which(TE.COEFS.ADJ.SCALED >  0.050 & TE.COEFS.ADJ.SCALED <  0.075)] <-  0.050
TE.COEFS.ADJ.SCALED[which(TE.COEFS.ADJ.SCALED >  0.075 & TE.COEFS.ADJ.SCALED <  0.100)] <-  0.075
TE.COEFS.ADJ.SCALED[which(TE.COEFS.ADJ.SCALED >  0.100 & TE.COEFS.ADJ.SCALED <  0.125)] <-  0.100
TE.COEFS.ADJ.SCALED[which(TE.COEFS.ADJ.SCALED >  0.125 & TE.COEFS.ADJ.SCALED <  0.150)] <-  0.125
TE.COEFS.ADJ.SCALED[which(TE.COEFS.ADJ.SCALED >  0.150 & TE.COEFS.ADJ.SCALED <  0.175)] <-  0.150
TE.COEFS.ADJ.SCALED[which(TE.COEFS.ADJ.SCALED >  0.175 & TE.COEFS.ADJ.SCALED <  0.200)] <-  0.175
TE.COEFS.ADJ.SCALED[which(TE.COEFS.ADJ.SCALED >  0.200 & TE.COEFS.ADJ.SCALED <  0.225)] <-  0.200
TE.COEFS.ADJ.SCALED[which(TE.COEFS.ADJ.SCALED >  0.225 & TE.COEFS.ADJ.SCALED <  0.250)] <-  0.225
TE.COEFS.ADJ.SCALED[which(TE.COEFS.ADJ.SCALED >  0.250)] <- 0.250
TE.COEFS.ADJ.SCALED[which(TE.COEFS.ADJ.SCALED <  0.000 & TE.COEFS.ADJ.SCALED > -0.025)] <-  0.000
TE.COEFS.ADJ.SCALED[which(TE.COEFS.ADJ.SCALED < -0.025 & TE.COEFS.ADJ.SCALED > -0.050)] <- -0.025
TE.COEFS.ADJ.SCALED[which(TE.COEFS.ADJ.SCALED < -0.050 & TE.COEFS.ADJ.SCALED > -0.075)] <- -0.050
TE.COEFS.ADJ.SCALED[which(TE.COEFS.ADJ.SCALED < -0.075 & TE.COEFS.ADJ.SCALED > -0.100)] <- -0.075
TE.COEFS.ADJ.SCALED[which(TE.COEFS.ADJ.SCALED < -0.100 & TE.COEFS.ADJ.SCALED > -0.125)] <- -0.100
TE.COEFS.ADJ.SCALED[which(TE.COEFS.ADJ.SCALED < -0.125 & TE.COEFS.ADJ.SCALED > -0.150)] <- -0.125
TE.COEFS.ADJ.SCALED[which(TE.COEFS.ADJ.SCALED < -0.150 & TE.COEFS.ADJ.SCALED > -0.175)] <- -0.150
TE.COEFS.ADJ.SCALED[which(TE.COEFS.ADJ.SCALED < -0.175 & TE.COEFS.ADJ.SCALED > -0.200)] <- -0.175
TE.COEFS.ADJ.SCALED[which(TE.COEFS.ADJ.SCALED < -0.200 & TE.COEFS.ADJ.SCALED > -0.225)] <- -0.200
TE.COEFS.ADJ.SCALED[which(TE.COEFS.ADJ.SCALED < -0.225 & TE.COEFS.ADJ.SCALED > -0.250)] <- -0.225
TE.COEFS.ADJ.SCALED[which(TE.COEFS.ADJ.SCALED < -0.250)] <- -0.250

##Main and pairwise effects
TE.MAIN.BIND.COEFS.df <- data.frame(state=as.character(c(AAs,AAs,AAs,AAs)),site=factor(c(rep(1,length(AAs)),rep(2,length(AAs)),rep(3,length(AAs)),rep(4,length(AAs))),levels=c("4","3","2","1")),value=-1*(TE.COEFS.ADJ.SCALED[TE.COEFS.INDEX.B1]),row.names=NULL)
TE.EP12.BIND.COEFS.df <- data.frame(state1=factor(rep(AAs,length(AAs)),levels=rev(AAs)), state2=factor(as.vector(sapply(1:length(AAs), function(x) rep(AAs[x],length(AAs))))),value= -1*(TE.COEFS.ADJ.SCALED[TE.COEFS.INDEX.B2.X1X2]))
TE.EP13.BIND.COEFS.df <- data.frame(state1=factor(rep(AAs,length(AAs)),levels=rev(AAs)), state2=factor(as.vector(sapply(1:length(AAs), function(x) rep(AAs[x],length(AAs))))),value= -1*(TE.COEFS.ADJ.SCALED[TE.COEFS.INDEX.B2.X1X3]))
TE.EP14.BIND.COEFS.df <- data.frame(state1=factor(rep(AAs,length(AAs)),levels=rev(AAs)), state2=factor(as.vector(sapply(1:length(AAs), function(x) rep(AAs[x],length(AAs))))),value= -1*(TE.COEFS.ADJ.SCALED[TE.COEFS.INDEX.B2.X1X4]))
TE.EP23.BIND.COEFS.df <- data.frame(state1=factor(rep(AAs,length(AAs)),levels=rev(AAs)), state2=factor(as.vector(sapply(1:length(AAs), function(x) rep(AAs[x],length(AAs))))),value= -1*(TE.COEFS.ADJ.SCALED[TE.COEFS.INDEX.B2.X2X3]))
TE.EP24.BIND.COEFS.df <- data.frame(state1=factor(rep(AAs,length(AAs)),levels=rev(AAs)), state2=factor(as.vector(sapply(1:length(AAs), function(x) rep(AAs[x],length(AAs))))),value= -1*(TE.COEFS.ADJ.SCALED[TE.COEFS.INDEX.B2.X2X4]))
TE.EP34.BIND.COEFS.df <- data.frame(state1=factor(rep(AAs,length(AAs)),levels=rev(AAs)), state2=factor(as.vector(sapply(1:length(AAs), function(x) rep(AAs[x],length(AAs))))),value= -1*(TE.COEFS.ADJ.SCALED[TE.COEFS.INDEX.B2.X3X4]))

##Plots of binding effects
#Main binding effects for position 1
TE.BIND.M1 <- ggplot(TE.MAIN.BIND.COEFS.df[TE.MAIN.BIND.COEFS.df$site==1,], aes(state, site))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+coord_equal()+scale_x_discrete(position = "top") + theme(legend.position = "none",panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
#Main binding effects for position 2
TE.BIND.M2 <- ggplot(TE.MAIN.BIND.COEFS.df[TE.MAIN.BIND.COEFS.df$site==2,], aes(state, site))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+coord_equal()+scale_x_discrete(position = "top") + theme(legend.position = "none",panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
#Main binding effects for position 3
TE.BIND.M3 <- ggplot(TE.MAIN.BIND.COEFS.df[TE.MAIN.BIND.COEFS.df$site==3,], aes(state, site))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+coord_equal()+scale_x_discrete(position = "top") + theme(legend.position = "none",panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
#Main binding effects for position 4
TE.BIND.M4 <- ggplot(TE.MAIN.BIND.COEFS.df[TE.MAIN.BIND.COEFS.df$site==4,], aes(state, site))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+coord_equal()+scale_x_discrete(position = "top") + theme(legend.position = "none",panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
#Pairwise binding interaction effects between positions 1 and 2
TE.BIND.P12 <- ggplot(TE.EP12.BIND.COEFS.df ,aes(state2,state1))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+coord_equal()+ theme(legend.position = "none",panel.grid = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank()) + 
  xlim(rev(levels(TE.EP12.BIND.COEFS.df$state1))) + ylim(rev(levels(TE.EP12.BIND.COEFS.df$state2)))
#Pairwise binding interaction effects between positions 1 and 3
TE.BIND.P13 <- ggplot(TE.EP13.BIND.COEFS.df ,aes(state2,state1))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+coord_equal()+ theme(legend.position = "none",panel.grid = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank()) + 
  xlim(rev(levels(TE.EP13.BIND.COEFS.df$state1))) + ylim(rev(levels(TE.EP13.BIND.COEFS.df$state2)))
#Pairwise binding interaction effects between positions 1 and 4
TE.BIND.P14 <- ggplot(TE.EP14.BIND.COEFS.df ,aes(state2,state1))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+coord_equal()+ theme(legend.position = "none",panel.grid = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank()) + 
  xlim(rev(levels(TE.EP14.BIND.COEFS.df$state1))) + ylim(rev(levels(TE.EP14.BIND.COEFS.df$state2)))
#Pairwise binding interaction effects between positions 2 and 3
TE.BIND.P23 <- ggplot(TE.EP23.BIND.COEFS.df ,aes(state2,state1))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+coord_equal()+ theme(legend.position = "none",panel.grid = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank()) + 
  xlim(rev(levels(TE.EP24.BIND.COEFS.df$state1))) + ylim(rev(levels(TE.EP24.BIND.COEFS.df$state2)))
#Pairwise binding interaction effects between positions 2 and 4
TE.BIND.P24 <- ggplot(TE.EP24.BIND.COEFS.df ,aes(state2,state1))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+coord_equal()+ theme(legend.position = "none",panel.grid = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank()) + 
  xlim(rev(levels(TE.EP24.BIND.COEFS.df$state1))) + ylim(rev(levels(TE.EP24.BIND.COEFS.df$state2)))
#Pairwise binding interaction effects between positions 3 and 4
TE.BIND.P34 <- ggplot(TE.EP34.BIND.COEFS.df ,aes(state2,state1))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+coord_equal()+ theme(legend.position = "none",panel.grid = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank()) + 
  xlim(rev(levels(TE.EP34.BIND.COEFS.df$state1))) + ylim(rev(levels(TE.EP34.BIND.COEFS.df$state2)))

##Specificity effects
TE.MAIN.SPEC.COEFS.df <- data.frame(state=as.character(c(AAs,AAs,AAs,AAs)),site=factor(c(rep(1,length(AAs)),rep(2,length(AAs)),rep(3,length(AAs)),rep(4,length(AAs))),levels=c("4","3","2","1")),value=-1*(TE.COEFS.ADJ.SCALED[TE.COEFS.INDEX.S1]),row.names=NULL)
TE.EP12.SPEC.COEFS.df <- data.frame(state1=factor(rep(AAs,length(AAs)),levels=rev(AAs)), state2=factor(as.vector(sapply(1:length(AAs), function(x) rep(AAs[x],length(AAs))))),value= -1*(TE.COEFS.ADJ.SCALED[TE.COEFS.INDEX.S2.X1X2]))
TE.EP13.SPEC.COEFS.df <- data.frame(state1=factor(rep(AAs,length(AAs)),levels=rev(AAs)), state2=factor(as.vector(sapply(1:length(AAs), function(x) rep(AAs[x],length(AAs))))),value= -1*(TE.COEFS.ADJ.SCALED[TE.COEFS.INDEX.S2.X1X3]))
TE.EP14.SPEC.COEFS.df <- data.frame(state1=factor(rep(AAs,length(AAs)),levels=rev(AAs)), state2=factor(as.vector(sapply(1:length(AAs), function(x) rep(AAs[x],length(AAs))))),value= -1*(TE.COEFS.ADJ.SCALED[TE.COEFS.INDEX.S2.X1X4]))
TE.EP23.SPEC.COEFS.df <- data.frame(state1=factor(rep(AAs,length(AAs)),levels=rev(AAs)), state2=factor(as.vector(sapply(1:length(AAs), function(x) rep(AAs[x],length(AAs))))),value= -1*(TE.COEFS.ADJ.SCALED[TE.COEFS.INDEX.S2.X2X3]))
TE.EP24.SPEC.COEFS.df <- data.frame(state1=factor(rep(AAs,length(AAs)),levels=rev(AAs)), state2=factor(as.vector(sapply(1:length(AAs), function(x) rep(AAs[x],length(AAs))))),value= -1*(TE.COEFS.ADJ.SCALED[TE.COEFS.INDEX.S2.X2X4]))
TE.EP34.SPEC.COEFS.df <- data.frame(state1=factor(rep(AAs,length(AAs)),levels=rev(AAs)), state2=factor(as.vector(sapply(1:length(AAs), function(x) rep(AAs[x],length(AAs))))),value= -1*(TE.COEFS.ADJ.SCALED[TE.COEFS.INDEX.S2.X3X4]))

##Plots of specificty effects
#Main specificty effects for position 1
TE.SPEC.M1 <- ggplot(TE.MAIN.SPEC.COEFS.df[TE.MAIN.SPEC.COEFS.df$site==1,], aes(state, site))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#00B300","#96FF63","#FFFFFF","#DB94FF","#7A00CC"),limits=c(-0.25,0.25)) +labs(x="",y="")+
  theme_classic()+ theme(legend.position = "none", panel.grid = element_blank(), axis.text.x = element_blank()) +
  coord_flip() + theme(aspect.ratio=18) + scale_x_discrete(limits = rev(levels(TE.MAIN.SPEC.COEFS.df$state)))
#Main binding effects for position 2
TE.SPEC.M2 <- ggplot(TE.MAIN.SPEC.COEFS.df[TE.MAIN.SPEC.COEFS.df$site==2,], aes(state, site))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#00B300","#96FF63","#FFFFFF","#DB94FF","#7A00CC"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+ theme(legend.position = "none", panel.grid = element_blank(), axis.text.x = element_blank()) +
  coord_flip() + theme(aspect.ratio=18) + scale_x_discrete(limits = rev(levels(TE.MAIN.SPEC.COEFS.df$state)))
#Main specificty effects for position 3
TE.SPEC.M3 <- ggplot(TE.MAIN.SPEC.COEFS.df[TE.MAIN.SPEC.COEFS.df$site==3,], aes(state, site))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#00B300","#96FF63","#FFFFFF","#DB94FF","#7A00CC"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+ theme(legend.position = "none", panel.grid = element_blank(), axis.text.x = element_blank()) +
  coord_flip() + theme(aspect.ratio=18) + scale_x_discrete(limits = rev(levels(TE.MAIN.SPEC.COEFS.df$state)))
#Main specificty effects for position 4
TE.SPEC.M4 <- ggplot(TE.MAIN.SPEC.COEFS.df[TE.MAIN.SPEC.COEFS.df$site==4,], aes(state, site))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#00B300","#96FF63","#FFFFFF","#DB94FF","#7A00CC"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+ theme(legend.position = "none", panel.grid = element_blank(), axis.text.x = element_blank()) +
  coord_flip() + theme(aspect.ratio=18) + scale_x_discrete(limits = rev(levels(TE.MAIN.SPEC.COEFS.df$state)))
#Pairwise specificty interaction effects between positions 1 and 2
TE.SPEC.P12 <- ggplot(TE.EP12.SPEC.COEFS.df ,aes(state1,state2))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#00B300","#96FF63","#FFFFFF","#DB94FF","#7A00CC"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+coord_equal() + theme(legend.position = "none",panel.grid = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank()) + 
  xlim(rev(levels(TE.EP12.SPEC.COEFS.df$state1))) + ylim(rev(levels(TE.EP12.SPEC.COEFS.df$state2)))
#Pairwise specificty interaction effects between positions 1 and 3
TE.SPEC.P13 <- ggplot(TE.EP13.SPEC.COEFS.df ,aes(state1,state2))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#00B300","#96FF63","#FFFFFF","#DB94FF","#7A00CC"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+coord_equal() + theme(legend.position = "none",panel.grid = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank()) + 
  xlim(rev(levels(TE.EP13.SPEC.COEFS.df$state1))) + ylim(rev(levels(TE.EP13.SPEC.COEFS.df$state2)))
#Pairwise specificty interaction effects between positions 1 and 4
TE.SPEC.P14 <- ggplot(TE.EP14.SPEC.COEFS.df ,aes(state1,state2))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#00B300","#96FF63","#FFFFFF","#DB94FF","#7A00CC"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+coord_equal() + theme(legend.position = "none",panel.grid = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank()) + 
  xlim(rev(levels(TE.EP14.SPEC.COEFS.df$state1))) + ylim(rev(levels(TE.EP14.SPEC.COEFS.df$state2)))
#Pairwise specificty interaction effects between positions 2 and 3
TE.SPEC.P23 <- ggplot(TE.EP23.SPEC.COEFS.df ,aes(state1,state2))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#00B300","#96FF63","#FFFFFF","#DB94FF","#7A00CC"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+coord_equal() + theme(legend.position = "none",panel.grid = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank()) + 
  xlim(rev(levels(TE.EP23.SPEC.COEFS.df$state1))) + ylim(rev(levels(TE.EP23.SPEC.COEFS.df$state2)))
#Pairwise specificty interaction effects between positions 2 and 4
TE.SPEC.P24 <- ggplot(TE.EP24.SPEC.COEFS.df ,aes(state1,state2))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#00B300","#96FF63","#FFFFFF","#DB94FF","#7A00CC"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+coord_equal() + theme(legend.position = "none",panel.grid = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank()) + 
  xlim(rev(levels(TE.EP24.SPEC.COEFS.df$state1))) + ylim(rev(levels(TE.EP24.SPEC.COEFS.df$state2)))
#Pairwise Specificity interaction effects between positions 3 and 4
TE.SPEC.P34 <- ggplot(TE.EP34.SPEC.COEFS.df ,aes(state1,state2))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#00B300","#96FF63","#FFFFFF","#DB94FF","#7A00CC"),limits=c(-0.25,0.25))+labs(x="",y="")+
  theme_classic()+coord_equal() + theme(legend.position = "none",panel.grid = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank()) + 
  xlim(rev(levels(TE.EP34.SPEC.COEFS.df$state1))) + ylim(rev(levels(TE.EP34.SPEC.COEFS.df$state2)))

###Combine binding and specificty interactions across positions into a single plot 
##pdf("TE.BIND.SPEC.PAIRWISE.pdf",5,3)
#margin = theme(plot.margin = unit(c(0,0,0,0), "cm"))
#pl <- list(TE.BIND.P12,TE.BIND.P13,TE.BIND.P14,TE.BIND.P23,TE.BIND.P24,TE.BIND.P34,
#           TE.SPEC.P12,TE.SPEC.P13,TE.SPEC.P14,TE.SPEC.P23,TE.SPEC.P24,TE.SPEC.P34)
#grid.arrange(grobs = lapply(pl, "+", margin),layout_matrix = rbind(c(NA,1,2,3),c(7,NA,4,5),c(8,10,NA,6),c(9,11,12,NA)))
##dev.off()  

###Main binding effects across positions
##pdf("TE.BIND.MAIN.pdf",5,3)
#margin = theme(plot.margin = unit(c(0,0,0,0), "cm"))
#pl <- list(TE.BIND.M1,TE.BIND.M2,TE.BIND.M3,TE.BIND.M4)
#grid.arrange(grobs = lapply(pl, "+", margin),layout_matrix = rbind(c(1,2,3,4)))
##dev.off()

###Main specificty effects across positions
##pdf("TE.SPEC.MAIN.pdf",5,3)
#margin = theme(plot.margin = unit(c(0,0,0,0), "cm"))
#pl <- list(TE.SPEC.M1,TE.SPEC.M2,TE.SPEC.M3,TE.SPEC.M4)
#grid.arrange(grobs = lapply(pl, "+", margin),layout_matrix = cbind(c(1,2,3,4)))
##dev.off()

##Number of + and - epistatic effects per amino acid
BIND.COUNT <- data.frame(matrix(0,nrow=80,ncol=2)); colnames(BIND.COUNT ) <- c("POS","NEG")
SPEC.COUNT <- data.frame(matrix(0,nrow=80,ncol=2)); colnames(SPEC.COUNT ) <- c("POS","NEG")
I <- 1:4
for(i in I) {
  J <- 1:20
  for(j in J) {
    TARGET <- paste0("X",i,AAs[j])
    SET <- TE.COEFS.ADJ.SCALED[TE.COEFS.INDEX.B2]
    BAIT <- SET[grep(TARGET,names(SET))]
    rownames(BIND.COUNT )[(i-1)*20+j] <- TARGET
    BIND.COUNT [(i-1)*20+j,1] <- sum(BAIT < 0)
    BIND.COUNT [(i-1)*20+j,2] <- sum(BAIT > 0)
    
    SET <- TE.COEFS.ADJ.SCALED[TE.COEFS.INDEX.S2]
    BAIT <- SET[grep(TARGET,names(SET))]
    rownames(SPEC.COUNT )[(i-1)*20+j] <- TARGET
    SPEC.COUNT [(i-1)*20+j,1] <- sum(BAIT < 0)
    SPEC.COUNT [(i-1)*20+j,2] <- sum(BAIT > 0)
  }
}

BIND.TOTAL <- apply(BIND.COUNT,1,sum)
SPEC.TOTAL <- apply(SPEC.COUNT,1,sum)

###3-way interaction effects
TE.COEF.3X <- data.frame(value = -1*TE.COEFS.ADJ.SCALED[c(TE.COEFS.INDEX.B3,TE.COEFS.INDEX.S3)])

TE.COEF.3X$X1 <- NA
TE.COEF.3X$X2 <- NA
TE.COEF.3X$X3 <- NA
TE.COEF.3X$X4 <- NA
TE.COEF.3X$RE <- rep(c(rep("B",8000),rep("S",8000)),4)
TE.COEF.3X$X12 <- NA
TE.COEF.3X$X13 <- NA
TE.COEF.3X$X14 <- NA
TE.COEF.3X$X23 <- NA
TE.COEF.3X$X24 <- NA
TE.COEF.3X$X34 <- NA

INDEX <- grep("X1",rownames(TE.COEF.3X)) 
TE.COEF.3X$X1[INDEX] <- substr(rownames(TE.COEF.3X)[INDEX],3,3)
INDEX <- grep("X1.:X2",rownames(TE.COEF.3X)) 
TE.COEF.3X$X2[INDEX] <- substr(rownames(TE.COEF.3X)[INDEX],7,7)
INDEX <- grep("X2.:X3.:X4",rownames(TE.COEF.3X)) 
TE.COEF.3X$X2[INDEX] <- substr(rownames(TE.COEF.3X)[INDEX],3,3)
INDEX <- grep("X1.:X2.:X3",rownames(TE.COEF.3X)) 
TE.COEF.3X$X3[INDEX] <- substr(rownames(TE.COEF.3X)[INDEX],11,11)
INDEX <- grep("X3.:X4",rownames(TE.COEF.3X)) 
TE.COEF.3X$X3[INDEX] <- substr(rownames(TE.COEF.3X)[INDEX],7,7)
INDEX <- grep("X4",rownames(TE.COEF.3X)) 
TE.COEF.3X$X4[INDEX] <- substr(rownames(TE.COEF.3X)[INDEX],11,11)

TE.COEF.3X$X12 <- paste0(TE.COEF.3X$X1,TE.COEF.3X$X2)
TE.COEF.3X$X13 <- paste0(TE.COEF.3X$X1,TE.COEF.3X$X3)
TE.COEF.3X$X14 <- paste0(TE.COEF.3X$X1,TE.COEF.3X$X4)
TE.COEF.3X$X23 <- paste0(TE.COEF.3X$X2,TE.COEF.3X$X3)
TE.COEF.3X$X24 <- paste0(TE.COEF.3X$X2,TE.COEF.3X$X4)
TE.COEF.3X$X34 <- paste0(TE.COEF.3X$X3,TE.COEF.3X$X4)

TE.COEF.3X$X1  <- factor(TE.COEF.3X$X1,levels=AAs)
TE.COEF.3X$X2  <- factor(TE.COEF.3X$X2,levels=AAs)
TE.COEF.3X$X3  <- factor(TE.COEF.3X$X3,levels=AAs)
TE.COEF.3X$X4  <- factor(TE.COEF.3X$X4,levels=AAs)
TE.COEF.3X$X12 <- factor(TE.COEF.3X$X12,levels=paste0(rep(AAs,each=20),AAs))
TE.COEF.3X$X13 <- factor(TE.COEF.3X$X13,levels=paste0(rep(AAs,each=20),AAs))
TE.COEF.3X$X14 <- factor(TE.COEF.3X$X14,levels=paste0(rep(AAs,each=20),AAs))
TE.COEF.3X$X23 <- factor(TE.COEF.3X$X23,levels=paste0(rep(AAs,each=20),AAs))
TE.COEF.3X$X24 <- factor(TE.COEF.3X$X24,levels=paste0(rep(AAs,each=20),AAs))
TE.COEF.3X$X34 <- factor(TE.COEF.3X$X34,levels=paste0(rep(AAs,each=20),AAs))

TE.COEF.3X.SUB <- TE.COEF.3X[TE.COEFS.INDEX.B3.X1X2X3-4960,]
TE.BIND.23.1 <- ggplot(TE.COEF.3X.SUB, aes(X1,X23))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-0.25,0.25))+
  labs(x="",y="")+theme_classic()+coord_equal()+scale_x_discrete(position = "top",drop=FALSE) + 
  scale_y_discrete(limits = rev(levels(TE.COEF.3X.SUB$X23))) +
  theme(legend.position = "none",panel.grid = element_blank(), text = element_text(size=5),axis.text = element_text(size = 1),axis.text.x = element_text(angle=90, hjust=1))
TE.BIND.23.1.PLOT <- set_panel_size(TE.BIND.23.1, height=unit(length(levels(TE.COEF.3X.SUB$X23))/2.5, "mm"),width=unit(length(levels(TE.COEF.3X.SUB$X1))/2.5, "mm") )

TE.COEF.3X.SUB <- TE.COEF.3X[TE.COEFS.INDEX.B3.X1X2X4-4960,]
TE.BIND.24.1 <- ggplot(TE.COEF.3X.SUB, aes(X1,X24))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-0.25,0.25))+
  labs(x="",y="")+theme_classic()+coord_equal()+scale_x_discrete(position = "top",drop=FALSE) + 
  scale_y_discrete(limits = rev(levels(TE.COEF.3X.SUB$X24))) +
  theme(legend.position = "none",panel.grid = element_blank(), text = element_text(size=5),axis.text = element_text(size = 1),axis.text.x = element_text(angle=90, hjust=1))
TE.BIND.24.1.PLOT <- set_panel_size( TE.BIND.24.1, height=unit(length(levels(TE.COEF.3X.SUB$X24))/2.5, "mm"),width=unit(length(levels(TE.COEF.3X.SUB$X1))/2.5, "mm") )

TE.COEF.3X.SUB <- TE.COEF.3X[TE.COEFS.INDEX.B3.X1X3X4-4960,]
TE.BIND.34.1 <- ggplot(TE.COEF.3X.SUB, aes(X1,X34))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-0.25,0.25))+
  labs(x="",y="")+theme_classic()+coord_equal()+scale_x_discrete(position = "top",drop=FALSE) + 
  scale_y_discrete(limits = rev(levels(TE.COEF.3X.SUB$X34))) +
  theme(legend.position = "none",panel.grid = element_blank(), text = element_text(size=5),axis.text = element_text(size = 1),axis.text.x = element_text(angle=90, hjust=1))
TE.BIND.34.1.PLOT <- set_panel_size( TE.BIND.34.1, height=unit(length(levels(TE.COEF.3X.SUB$X34))/2.5, "mm"),width=unit(length(levels(TE.COEF.3X.SUB$X1))/2.5, "mm") )

TE.COEF.3X.SUB <- TE.COEF.3X[TE.COEFS.INDEX.B3.X2X3X4-4960,]
TE.BIND.23.4 <- ggplot(TE.COEF.3X.SUB, aes(X4,X23))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-0.25,0.25))+
  labs(x="",y="")+theme_classic()+coord_equal()+scale_x_discrete(position = "top",drop=FALSE) + 
  scale_y_discrete(limits = rev(levels(TE.COEF.3X.SUB$X23))) +
  theme(legend.position = "none",panel.grid = element_blank(), text = element_text(size=5),axis.text = element_text(size = 1),axis.text.x = element_text(angle=90, hjust=1))
TE.BIND.23.4.PLOT <- set_panel_size( TE.BIND.23.4, height=unit(length(levels(TE.COEF.3X.SUB$X23))/2.5, "mm"),width=unit(length(levels(TE.COEF.3X.SUB$X4))/2.5, "mm") )

TE.COEF.3X.SUB <- TE.COEF.3X[TE.COEFS.INDEX.S3.X1X2X3-4960,]
TE.SPEC.23.1 <- ggplot(TE.COEF.3X.SUB, aes(X1,X23))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#00B300","#96FF63","#FFFFFF","#DB94FF","#7A00CC"),limits=c(-0.25,0.25))+
  labs(x="",y="")+theme_classic()+coord_equal()+scale_x_discrete(position = "top",drop=FALSE) + 
  scale_y_discrete(limits = rev(levels(TE.COEF.3X.SUB$X23))) +
  theme(legend.position = "none",panel.grid = element_blank(), text = element_text(size=5),axis.text = element_text(size = 1),axis.text.x = element_text(angle=90, hjust=1))
TE.SPEC.23.1.PLOT <- set_panel_size( TE.SPEC.23.1, height=unit(length(levels(TE.COEF.3X.SUB$X23))/2.5, "mm"),width=unit(length(levels(TE.COEF.3X.SUB$X1))/2.5, "mm") )

TE.COEF.3X.SUB <- TE.COEF.3X[TE.COEFS.INDEX.S3.X1X2X4-4960,]
TE.SPEC.24.1 <- ggplot(TE.COEF.3X.SUB, aes(X1,X24))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#00B300","#96FF63","#FFFFFF","#DB94FF","#7A00CC"),limits=c(-0.25,0.25))+
  labs(x="",y="")+theme_classic()+coord_equal()+scale_x_discrete(position = "top",drop=FALSE) + 
  scale_y_discrete(limits = rev(levels(TE.COEF.3X.SUB$X24))) +
  theme(legend.position = "none",panel.grid = element_blank(), text = element_text(size=5),axis.text = element_text(size = 1),axis.text.x = element_text(angle=90, hjust=1))
TE.SPEC.24.1.PLOT <- set_panel_size( TE.SPEC.24.1, height=unit(length(levels(TE.COEF.3X.SUB$X24))/2.5, "mm"),width=unit(length(levels(TE.COEF.3X.SUB$X1))/2.5, "mm") )

TE.COEF.3X.SUB <- TE.COEF.3X[TE.COEFS.INDEX.S3.X1X3X4-4960,]
TE.SPEC.34.1 <- ggplot(TE.COEF.3X.SUB, aes(X1,X34))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#00B300","#96FF63","#FFFFFF","#DB94FF","#7A00CC"),limits=c(-0.25,0.25))+
  labs(x="",y="")+theme_classic()+coord_equal()+scale_x_discrete(position = "top",drop=FALSE) + 
  scale_y_discrete(limits = rev(levels(TE.COEF.3X.SUB$X34))) +
  theme(legend.position = "none",panel.grid = element_blank(), text = element_text(size=5),axis.text = element_text(size = 1),axis.text.x = element_text(angle=90, hjust=1))
TE.SPEC.34.1.PLOT <- set_panel_size( TE.SPEC.34.1, height=unit(length(levels(TE.COEF.3X.SUB$X34))/2.5, "mm"),width=unit(length(levels(TE.COEF.3X.SUB$X1))/2.5, "mm") )

TE.COEF.3X.SUB <- TE.COEF.3X[TE.COEFS.INDEX.S3.X2X3X4-4960,]
TE.SPEC.23.4<- ggplot(TE.COEF.3X.SUB, aes(X4,X23))+geom_tile(aes(fill=value))+
  scale_fill_gradientn(colours=c("#00B300","#96FF63","#FFFFFF","#DB94FF","#7A00CC"),limits=c(-0.25,0.25))+
  labs(x="",y="")+theme_classic()+coord_equal()+scale_x_discrete(position = "top",drop=FALSE) + 
  scale_y_discrete(limits = rev(levels(TE.COEF.3X.SUB$X23))) +
  theme(legend.position = "none",panel.grid = element_blank(), text = element_text(size=5),axis.text = element_text(size = 1),axis.text.x = element_text(angle=90, hjust=1))
TE.SPEC.23.4.PLOT <- set_panel_size( TE.SPEC.23.4, height=unit(length(levels(TE.COEF.3X.SUB$X23))/2.5, "mm"),width=unit(length(levels(TE.COEF.3X.SUB$X4))/2.5, "mm") )

#pdf("TE.3D.BIND.pdf")
#cowplot::plot_grid(TE.BIND.23.1.PLOT,TE.BIND.24.1.PLOT,TE.BIND.34.1.PLOT,TE.BIND.23.4.PLOT,nrow=1)
#dev.off()
#pdf("TE.3D.SPEC.pdf")
#cowplot::plot_grid(TE.SPEC.23.1.PLOT,TE.SPEC.24.1.PLOT,TE.SPEC.34.1.PLOT,TE.SPEC.23.4.PLOT,nrow=1)
#dev.off()

#####################
##19. Net epistasis##
#####################
PRED.NULL <- rownames(TE.MATRIX)[PREDICT.TE.CLASS$class == "null"]
PRED.ACTIVATE <- rownames(TE.MATRIX)[PREDICT.TE.CLASS$class != "null"]
PRED.WEAK <- rownames(TE.MATRIX)[PREDICT.TE.CLASS$class == "weak"]
PRED.STRONG <- rownames(TE.MATRIX)[PREDICT.TE.CLASS$class == "strong"]

load("TE.COEF.MATRIX.ADJ.rda")
TE.COEF.MATRIX.ADJ    <- t(t(TE.MATRIX[,-TE.COEFS.INDEX.S0])*TE.COEFS.ADJ)

##Net First Order (i.e. net main effect)
NET.EPISTASIS.1.PHENOTYPE.WEAK    <- -1*rowSums(TE.COEF.MATRIX.ADJ[PRED.WEAK,   TE.COEFS.INDEX.O1])/(-1*TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0)
NET.EPISTASIS.1.PHENOTYPE.STRONG  <- -1*rowSums(TE.COEF.MATRIX.ADJ[PRED.STRONG, TE.COEFS.INDEX.O1])/(-1*TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0)

##Net Second Order Epistasis
NET.EPISTASIS.2.PHENOTYPE.WEAK    <- -1*rowSums(TE.COEF.MATRIX.ADJ[PRED.WEAK,   TE.COEFS.INDEX.O2])/(-1*TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0)
NET.EPISTASIS.2.PHENOTYPE.STRONG  <- -1*rowSums(TE.COEF.MATRIX.ADJ[PRED.STRONG, TE.COEFS.INDEX.O2])/(-1*TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0)

##Net Third order Epistasis
NET.EPISTASIS.3.PHENOTYPE.WEAK    <- -1*rowSums(TE.COEF.MATRIX.ADJ[PRED.WEAK,   TE.COEFS.INDEX.O3])/(-1*TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0)
NET.EPISTASIS.3.PHENOTYPE.STRONG  <- -1*rowSums(TE.COEF.MATRIX.ADJ[PRED.STRONG, TE.COEFS.INDEX.O3])/(-1*TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0)

##Net Second and Third order (i.e net epistasis)
NET.EPISTASIS.23.PHENOTYPE.WEAK    <- -1*rowSums(TE.COEF.MATRIX.ADJ[PRED.WEAK,    c(TE.COEFS.INDEX.O2,TE.COEFS.INDEX.O3)])/(-1*TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0)
NET.EPISTASIS.23.PHENOTYPE.STRONG  <- -1*rowSums(TE.COEF.MATRIX.ADJ[PRED.STRONG,  c(TE.COEFS.INDEX.O2,TE.COEFS.INDEX.O3)])/(-1*TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0)

##Net first and second order
NET.EPISTASIS.12.PHENOTYPE.WEAK    <- -1*rowSums(TE.COEF.MATRIX.ADJ[PRED.WEAK,    c(TE.COEFS.INDEX.O1,TE.COEFS.INDEX.O2)])/(-1*TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0)
NET.EPISTASIS.12.PHENOTYPE.STRONG  <- -1*rowSums(TE.COEF.MATRIX.ADJ[PRED.STRONG,  c(TE.COEFS.INDEX.O1,TE.COEFS.INDEX.O2)])/(-1*TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0)

##Net first and third order
NET.EPISTASIS.13.PHENOTYPE.WEAK    <- -1*rowSums(TE.COEF.MATRIX.ADJ[PRED.WEAK,    c(TE.COEFS.INDEX.O1,TE.COEFS.INDEX.O3)])/(-1*TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0)
NET.EPISTASIS.13.PHENOTYPE.STRONG  <- -1*rowSums(TE.COEF.MATRIX.ADJ[PRED.STRONG,  c(TE.COEFS.INDEX.O1,TE.COEFS.INDEX.O3)])/(-1*TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0)


##pdf("NET.EPISTASIS.pdf")
#par(mfrow=c(2,2))
#hist(NET.EPISTASIS.PHENOTYPE.STRONG.S,    breaks=seq(-1.5,1.6,.05),xlim=c(-0.5,1.6),xlab="Net Epistasis",ylim=c(0,200),  main="Strong 2+3 rel S")
#abline(v=mean(NET.EPISTASIS.PHENOTYPE.STRONG.S),col="red"); abline(v=1,col="orange");abline(v=(TE.THRESH.NULL+TE.COEFS.EFFECT.ADJ.B0)/(TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0),col="brown")
#hist(NET.EPISTASIS.PHENOTYPE.WEAK.S,      breaks=seq(-1.5,1.6,.05),xlim=c(-0.5,1.6),xlab="Net Epistasis",ylim=c(0,400),  main="Weak 2+3 rel S"  )
#abline(v=mean(NET.EPISTASIS.PHENOTYPE.WEAK.S),  col="red"); abline(v=1,col="orange");abline(v=(TE.THRESH.NULL+TE.COEFS.EFFECT.ADJ.B0)/(TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0),col="brown")
#hist(NET.EPISTASIS.PHENOTYPE.STRONG.W,    breaks=seq(-1.5,1.6,.05),xlim=c(-0.5,1.6),xlab="Net Epistasis",ylim=c(0,200),  main="Strong 2+3 rel W")
#abline(v=mean(NET.EPISTASIS.PHENOTYPE.STRONG.W),col="red"); abline(v=1,col="orange");abline(v=(TE.THRESH.NULL+TE.COEFS.EFFECT.ADJ.B0)/(TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0),col="brown")
#hist(NET.EPISTASIS.PHENOTYPE.WEAK.W,      breaks=seq(-1.5,1.6,.05),xlim=c(-0.5,1.6),xlab="Net Epistasis",ylim=c(0,400),  main="Weak 2+3 rel W"  )
#abline(v=mean(NET.EPISTASIS.PHENOTYPE.WEAK.W),  col="red"); abline(v=1,col="orange");abline(v=(TE.THRESH.NULL+TE.COEFS.EFFECT.ADJ.B0)/(TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0),col="brown")

#hist(NET.EPISTASIS.3.PHENOTYPE.STRONG.S,    breaks=seq(-1.5,1.6,.05),xlim=c(-0.5,1.6),xlab="Net Epistasis",ylim=c(0,400),  main="Strong 3 rel S")
#abline(v=mean(NET.EPISTASIS.3.PHENOTYPE.STRONG.S),col="red"); abline(v=1,col="orange");abline(v=(TE.THRESH.NULL+TE.COEFS.EFFECT.ADJ.B0)/(TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0),col="brown")
#hist(NET.EPISTASIS.3.PHENOTYPE.WEAK.S,      breaks=seq(-1.5,1.6,.05),xlim=c(-0.5,1.6),xlab="Net Epistasis",ylim=c(0,700),  main="Weak 3 rel S"  )
#abline(v=mean(NET.EPISTASIS.3.PHENOTYPE.WEAK.S),  col="red"); abline(v=1,col="orange");abline(v=(TE.THRESH.NULL+TE.COEFS.EFFECT.ADJ.B0)/(TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0),col="brown")
#hist(NET.EPISTASIS.3.PHENOTYPE.STRONG.W,    breaks=seq(-1.5,1.6,.05),xlim=c(-0.5,1.6),xlab="Net Epistasis",ylim=c(0,400),  main="Strong 3 rel W")
#abline(v=mean(NET.EPISTASIS.3.PHENOTYPE.STRONG.W),col="red"); abline(v=1,col="orange");abline(v=(TE.THRESH.NULL+TE.COEFS.EFFECT.ADJ.B0)/(TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0),col="brown")
#hist(NET.EPISTASIS.3.PHENOTYPE.WEAK.W,      breaks=seq(-1.5,1.6,.05),xlim=c(-0.5,1.6),xlab="Net Epistasis",ylim=c(0,700),  main="Weak 3 rel W"  )
#abline(v=mean(NET.EPISTASIS.3.PHENOTYPE.WEAK.W),  col="red"); abline(v=1,col="orange");abline(v=(TE.THRESH.NULL+TE.COEFS.EFFECT.ADJ.B0)/(TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0),col="brown")
##dev.off()

################################################
##20. Number, Identity, Function of Activators##
################################################
##Number
sum(PREDICT.TE.CLASS$class == "strong")
sum(PREDICT.PE.CLASS$class == "strong")
sum(PREDICT.ME.CLASS$class == "strong")

sum(PREDICT.TE.CLASS$class == "weak")
sum(PREDICT.PE.CLASS$class == "weak")
sum(PREDICT.ME.CLASS$class == "weak")

##Identity
sum(PREDICT.ME.CLASS$class[PREDICT.TE.CLASS$class == "strong"] == "strong")
sum(PREDICT.TE.CLASS$class[PREDICT.ME.CLASS$class == "strong"] == "strong")

#############################################
##20. Direction of effect for substitutions##
#############################################
###Determine if each amino acid has a benefit/detriment on binding in all genetic contexts

#K <- 1:4
#J <- 1:20
#I <- 1:20
#AA.EFFECT.START.TG <- matrix(0,nrow=1600,ncol=16000)
#AA.EFFECT.START.PG <- matrix(0,nrow=1600,ncol=16000)
#AA.EFFECT.START.MG <- matrix(0,nrow=1600,ncol=16000)
#AA.EFFECT.STOP.TG  <- matrix(0,nrow=1600,ncol=16000)
#AA.EFFECT.STOP.PG  <- matrix(0,nrow=1600,ncol=16000)
#AA.EFFECT.STOP.MG  <- matrix(0,nrow=1600,ncol=16000)
#rownames(AA.EFFECT.START.TG) <- c(paste0("X1:",rep(AAs,each=20),AAs),paste0("X2:",rep(AAs,each=20),AAs),paste0("X3:",rep(AAs,each=20),AAs),paste0("X4:",rep(AAs,each=20),AAs))
#rownames(AA.EFFECT.START.PG) <- c(paste0("X1:",rep(AAs,each=20),AAs),paste0("X2:",rep(AAs,each=20),AAs),paste0("X3:",rep(AAs,each=20),AAs),paste0("X4:",rep(AAs,each=20),AAs))
#rownames(AA.EFFECT.START.MG) <- c(paste0("X1:",rep(AAs,each=20),AAs),paste0("X2:",rep(AAs,each=20),AAs),paste0("X3:",rep(AAs,each=20),AAs),paste0("X4:",rep(AAs,each=20),AAs))
#rownames(AA.EFFECT.STOP.TG)  <- c(paste0("X1:",rep(AAs,each=20),AAs),paste0("X2:",rep(AAs,each=20),AAs),paste0("X3:",rep(AAs,each=20),AAs),paste0("X4:",rep(AAs,each=20),AAs))
#rownames(AA.EFFECT.STOP.PG)  <- c(paste0("X1:",rep(AAs,each=20),AAs),paste0("X2:",rep(AAs,each=20),AAs),paste0("X3:",rep(AAs,each=20),AAs),paste0("X4:",rep(AAs,each=20),AAs))
#rownames(AA.EFFECT.STOP.MG)  <- c(paste0("X1:",rep(AAs,each=20),AAs),paste0("X2:",rep(AAs,each=20),AAs),paste0("X3:",rep(AAs,each=20),AAs),paste0("X4:",rep(AAs,each=20),AAs))
#for(k in K) {
#  for(j in J) {
#    for(i in I) {
#      if(k == 1) {
#        S1 <- DT.JOINT[1:320000,][AA1==AAs[j]]
#        S2 <- DT.JOINT[1:320000,][AA1==AAs[i]]
#      }else if(k == 2) {
#        S1 <- DT.JOINT[1:320000,][AA2==AAs[j]]
#        S2 <- DT.JOINT[1:320000,][AA2==AAs[i]]
#      }else if(k == 3) {
#        S1 <- DT.JOINT[1:320000,][AA3==AAs[j]]
#        S2 <- DT.JOINT[1:320000,][AA3==AAs[i]]
#      }else if(k == 4) {
#        S1 <- DT.JOINT[1:320000,][AA4==AAs[j]]
#        S2 <- DT.JOINT[1:320000,][AA4==AAs[i]]
#      }
#      AA.EFFECT.START.TG[(k-1)*400+(j-1)*20+i,] <- S1$PRED.TE.LINK
#      AA.EFFECT.START.PG[(k-1)*400+(j-1)*20+i,] <- S1$PRED.PE.LINK
#      AA.EFFECT.START.MG[(k-1)*400+(j-1)*20+i,] <- S1$PRED.ME.LINK
#      AA.EFFECT.STOP.TG[ (k-1)*400+(j-1)*20+i,] <- S2$PRED.TE.LINK
#      AA.EFFECT.STOP.PG[ (k-1)*400+(j-1)*20+i,] <- S2$PRED.PE.LINK
#      AA.EFFECT.STOP.MG[ (k-1)*400+(j-1)*20+i,] <- S2$PRED.ME.LINK
#    }
#  }
#}
###Remove instances where i==j
#AA.EFFECT.START.TG <- AA.EFFECT.START.TG[-c((I-1)*20+I,400+(I-1)*20+I,800+(I-1)*20+I,1200+(I-1)*20+I),]
#AA.EFFECT.START.PG <- AA.EFFECT.START.PG[-c((I-1)*20+I,400+(I-1)*20+I,800+(I-1)*20+I,1200+(I-1)*20+I),]
#AA.EFFECT.START.MG <- AA.EFFECT.START.MG[-c((I-1)*20+I,400+(I-1)*20+I,800+(I-1)*20+I,1200+(I-1)*20+I),]
#AA.EFFECT.STOP.TG  <- AA.EFFECT.STOP.TG[ -c((I-1)*20+I,400+(I-1)*20+I,800+(I-1)*20+I,1200+(I-1)*20+I),]
#AA.EFFECT.STOP.PG  <- AA.EFFECT.STOP.PG[ -c((I-1)*20+I,400+(I-1)*20+I,800+(I-1)*20+I,1200+(I-1)*20+I),]
#AA.EFFECT.STOP.MG  <- AA.EFFECT.STOP.MG[ -c((I-1)*20+I,400+(I-1)*20+I,800+(I-1)*20+I,1200+(I-1)*20+I),]
#
#save(AA.EFFECT.START.TG,file="AA.EFFECT.START.TG.rda")
#save(AA.EFFECT.START.PG,file="AA.EFFECT.START.PG.rda")
#save(AA.EFFECT.START.MG,file="AA.EFFECT.START.MG.rda")
#save(AA.EFFECT.STOP.TG, file="AA.EFFECT.STOP.TG.rda")
#save(AA.EFFECT.STOP.PG, file="AA.EFFECT.STOP.PG.rda")
#save(AA.EFFECT.STOP.MG, file="AA.EFFECT.STOP.MG.rda")
#
load("AA.EFFECT.START.TG.rda")
load("AA.EFFECT.START.PG.rda")
load("AA.EFFECT.START.MG.rda")
load("AA.EFFECT.STOP.TG.rda")
load("AA.EFFECT.STOP.PG.rda")
load("AA.EFFECT.STOP.MG.rda")

AA.EFFECT.TG <- AA.EFFECT.START.TG - AA.EFFECT.STOP.TG
AA.EFFECT.PG <- AA.EFFECT.START.PG - AA.EFFECT.STOP.PG
AA.EFFECT.MG <- AA.EFFECT.START.MG - AA.EFFECT.STOP.MG

###Plot results for TE model
#hist(c(AA.EFFECT.TG),xlim=c(-25,25),xlab="Mutation Effect",main="TE")
#hist(apply(AA.EFFECT.TG,1,min),breaks=seq(-25,25,.5),xlab="Worst Mutation Effect",main="TE")
#hist(apply(AA.EFFECT.TG,1,max),breaks=seq(-25,25,.5),xlab="Best Mutation Effect",main="TE")
#hist(apply(AA.EFFECT.TG < 0,1,sum),breaks=seq(0,16000,800),xlab="# backgrounds reduced binding",main="TE")
#hist(apply(AA.EFFECT.TG > 0,1,sum),breaks=seq(0,16000,800),xlab="# backgrounds increased binding",main="TE")

#Never Decrease
sum(apply(AA.EFFECT.TG < 0,1,sum)==0)

#Never Increase
sum(apply(AA.EFFECT.TG > 0,1,sum)==0)

#plot(-AA.EFFECT.START.TG[1,],AA.EFFECT.TG[1,],pch=19,cex=0.4,col="#00000022",xlab="Starting Phenotype",ylab="Change in Phenotype")
#abline(h=0,col="red")
#abline(v=-TE.THRESH.NULL,col="red",lty=2)
#abline(v=-TE.THRESH.WEAK,col="red",lty=2)

##Look at substitutions among things that are activators or become activators
#How many backgrounds did each sub have a positive effect on
ACT.INCREASE.EFFECT <- sapply(1:1520,FUN = function(x) {
  sum(AA.EFFECT.TG[x,which(AA.EFFECT.START.TG[x,] < TE.THRESH.NULL | AA.EFFECT.STOP.TG[x,] < TE.THRESH.NULL)] > 0)
})
#How many backgrounds did each sub have a negative effect on
ACT.DECREASE.EFFECT <- sapply(1:1520,FUN = function(x) {
  sum(AA.EFFECT.TG[x,which(AA.EFFECT.START.TG[x,] < TE.THRESH.NULL | AA.EFFECT.STOP.TG[x,] < TE.THRESH.NULL)] < 0)
})
#How many background was each sub evaluated on
ACT.NUMBER <- sapply(1:1520,FUN = function(x) {
  length(which(AA.EFFECT.START.TG[x,] < TE.THRESH.NULL | AA.EFFECT.STOP.TG[x,] < TE.THRESH.NULL))
})

#Never increase
sum(ACT.INCREASE.EFFECT==0)
#Never decrease
sum(ACT.DECREASE.EFFECT==0)
#Never tested 
sum(ACT.NUMBER < 2)

#Fraction of tested on activator subs that:
#Never increase
(sum(ACT.INCREASE.EFFECT==0) - sum(ACT.NUMBER < 2))/(1520 - sum(ACT.NUMBER < 2))
#Never decrease
(sum(ACT.DECREASE.EFFECT==0) - sum(ACT.NUMBER < 2))/(1520 - sum(ACT.NUMBER < 2))
#Have + and - effects
1-(sum(ACT.INCREASE.EFFECT==0) - sum(ACT.NUMBER < 2) + sum(ACT.DECREASE.EFFECT==0) - sum(ACT.NUMBER < 2))/(1520 - sum(ACT.NUMBER < 2))

#####Check ERE and SRE binding separately
###Sign epistasis with small mutation threshold
#Get relative effect of each mutation
REL.EFFECT.TG <- AA.EFFECT.TG/(-1*TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)
#Set cutoff value for 'small' effects
CUTOFF <- .05
#Determine increases, decreases, and no effect of each mutation
REL.EFFECT.TG[which(REL.EFFECT.TG >= -CUTOFF & REL.EFFECT.TG <= CUTOFF)] <- 0 
REL.EFFECT.TG[which(REL.EFFECT.TG > CUTOFF)] <- 1  
REL.EFFECT.TG[which(REL.EFFECT.TG < -CUTOFF)] <- -1  

#Mutations that never/only increase/decrease on ERE
sum(apply(REL.EFFECT.TG[,1:8000] > 0,1,sum) == 0)
sum(apply(REL.EFFECT.TG[,1:8000] > 0,1,sum) == 8000)
sum(apply(REL.EFFECT.TG[,1:8000] < 0,1,sum) == 0)
sum(apply(REL.EFFECT.TG[,1:8000] < 0,1,sum) == 8000)

#Number of backgrounds increase/decrease on ERE
plot(apply(REL.EFFECT.TG[,1:8000] > 0,1,sum),apply(REL.EFFECT.TG[,1:8000] < 0,1,sum),pch=19,cex=0.6,col="#00000044",xlab="# of backgrounds increase",ylab="# of backgrounds decrease")
hist(apply(REL.EFFECT.TG[,1:8000] > 0,1,sum),breaks=seq(0,8000,200),xlab="# backgrounds increased binding on ERE",main="TE ERE")
hist(apply(REL.EFFECT.TG[,1:8000] < 0,1,sum),breaks=seq(0,8000,200),xlab="# backgrounds decreased binding on ERE",main="TE ERE")

#Number of mutation types that both increase and decrease on more than 10% of backgrounds for ERE
sum(apply(REL.EFFECT.TG[,1:8000] > 0,1,sum) > 400 & apply(REL.EFFECT.TG[,1:8000] > 0,1,sum) < 7600 & apply(REL.EFFECT.TG[,1:8000] < 0,1,sum) > 400 & apply(REL.EFFECT.TG[,1:8000] < 0,1,sum) < 7600)

#Mutations that never/only increase/decrease on SRE
sum(apply(REL.EFFECT.TG[,8001:16000] > 0,1,sum) == 0)
sum(apply(REL.EFFECT.TG[,8001:16000] > 0,1,sum) == 8000)
sum(apply(REL.EFFECT.TG[,8001:16000] < 0,1,sum) == 0)
sum(apply(REL.EFFECT.TG[,8001:16000] < 0,1,sum) == 8000)

#Number of backgrounds increase/decrease on SRE
plot(apply(REL.EFFECT.TG[,8001:16000] > 0,1,sum),apply(REL.EFFECT.TG[,8001:16000] < 0,1,sum),pch=19,cex=0.6,col="#00000044",xlab="# of backgrounds increase",ylab="# of backgrounds decrease")
hist(apply(REL.EFFECT.TG[,8001:16000] > 0,1,sum),breaks=seq(0,8000,200),xlab="# backgrounds increased binding on SRE",main="TE SRE")
hist(apply(REL.EFFECT.TG[,8001:16000] < 0,1,sum),breaks=seq(0,8000,200),xlab="# backgrounds decreased binding on SRE",main="TE SRE")

#Number of mutation types that both increase and decrease on more than 10% of backgrounds for SRE
sum(apply(REL.EFFECT.TG[,8001:16000] > 0,1,sum) > 400 & apply(REL.EFFECT.TG[,8001:16000] > 0,1,sum) < 7600 & apply(REL.EFFECT.TG[,8001:16000] < 0,1,sum) > 400 & apply(REL.EFFECT.TG[,8001:16000] < 0,1,sum) < 7600)

#Mutations that never/only increase/decrease regardless of RE
sum(apply(REL.EFFECT.TG > 0,1,sum) == 0)
sum(apply(REL.EFFECT.TG > 0,1,sum) == 16000)
sum(apply(REL.EFFECT.TG < 0,1,sum) == 0)
sum(apply(REL.EFFECT.TG < 0,1,sum) == 16000)

#Number of backgrounds increase/decrease across both REs
plot(apply(REL.EFFECT.TG > 0,1,sum),apply(REL.EFFECT.TG < 0,1,sum),pch=19,cex=0.6,col="#00000044",xlab="# of backgrounds increase",ylab="# of backgrounds decrease")
hist(apply(REL.EFFECT.TG > 0,1,sum),breaks=seq(0,16000,400),xlab="# backgrounds increased binding",main="TE")
hist(apply(REL.EFFECT.TG < 0,1,sum),breaks=seq(0,16000,400),xlab="# backgrounds decreased binding",main="TE")

#Number of mutation types that both increase and decrease on more than 10% of backgrounds
sum(apply(REL.EFFECT.TG > 0,1,sum) > 800 & apply(REL.EFFECT.TG > 0,1,sum) < 15200 & apply(REL.EFFECT.TG < 0,1,sum) > 800 & apply(REL.EFFECT.TG < 0,1,sum) < 15200)


###ERE
#How many backgrounds did each sub have a positive effect on for ERE
ACT.ERE.INCREASE.EFFECT <- sapply(1:1520,FUN = function(x) {
  sum(AA.EFFECT.TG[x,which(AA.EFFECT.START.TG[x,1:8000] < TE.THRESH.NULL | AA.EFFECT.STOP.TG[x,1:8000] < TE.THRESH.NULL)] > 0)
})
#How many backgrounds did each sub have a negative effect on for ERE
ACT.ERE.DECREASE.EFFECT <- sapply(1:1520,FUN = function(x) {
  sum(AA.EFFECT.TG[x,which(AA.EFFECT.START.TG[x,1:8000] < TE.THRESH.NULL | AA.EFFECT.STOP.TG[x,1:8000] < TE.THRESH.NULL)] < 0)
})
#How many background was each sub evaluated on
ACT.ERE.NUMBER <- sapply(1:1520,FUN = function(x) {
  length(which(AA.EFFECT.START.TG[x,1:8000] < TE.THRESH.NULL | AA.EFFECT.STOP.TG[x,1:8000] < TE.THRESH.NULL))
})
#Never increase
sum(ACT.ERE.INCREASE.EFFECT==0)
#Never decrease
sum(ACT.ERE.DECREASE.EFFECT==0)
#Never tested 
sum(ACT.ERE.NUMBER < 2)

#Fraction of tested on activator subs that:
#Never increase
(sum(ACT.ERE.INCREASE.EFFECT==0) - sum(ACT.ERE.NUMBER < 2))/(1520 - sum(ACT.ERE.NUMBER < 2))
#Never decrease
(sum(ACT.ERE.DECREASE.EFFECT==0) - sum(ACT.ERE.NUMBER < 2))/(1520 - sum(ACT.ERE.NUMBER < 2))
#Have + and - effects
1-(sum(ACT.ERE.INCREASE.EFFECT==0) - sum(ACT.ERE.NUMBER < 2) + sum(ACT.ERE.DECREASE.EFFECT==0) - sum(ACT.ERE.NUMBER < 2))/(1520 - sum(ACT.ERE.NUMBER < 2))

###SRE
#How many backgrounds did each sub have a positive effect on for SRE
ACT.SRE.INCREASE.EFFECT <- sapply(1:1520,FUN = function(x) {
  sum(AA.EFFECT.TG[x,which(AA.EFFECT.START.TG[x,8001:16000] < TE.THRESH.NULL | AA.EFFECT.STOP.TG[x,8001:16000] < TE.THRESH.NULL)] > 0)
})
#How many backgrounds did each sub have a negative effect on for SRE
ACT.SRE.DECREASE.EFFECT <- sapply(1:1520,FUN = function(x) {
  sum(AA.EFFECT.TG[x,which(AA.EFFECT.START.TG[x,8001:16000] < TE.THRESH.NULL | AA.EFFECT.STOP.TG[x,8001:16000] < TE.THRESH.NULL)] < 0)
})
#How many background was each sub evaluated on
ACT.SRE.NUMBER <- sapply(1:1520,FUN = function(x) {
  length(which(AA.EFFECT.START.TG[x,8001:16000] < TE.THRESH.NULL | AA.EFFECT.STOP.TG[x,8001:16000] < TE.THRESH.NULL))
})
#Never increase
sum(ACT.SRE.INCREASE.EFFECT==0)
#Never decrease
sum(ACT.SRE.DECREASE.EFFECT==0)
#Never tested 
sum(ACT.SRE.NUMBER < 2)

#Fraction of tested on activator subs that:
#Never increase
(sum(ACT.SRE.INCREASE.EFFECT==0) - sum(ACT.SRE.NUMBER < 2))/(1520 - sum(ACT.SRE.NUMBER < 2))
#Never decrease
(sum(ACT.SRE.DECREASE.EFFECT==0) - sum(ACT.SRE.NUMBER < 2))/(1520 - sum(ACT.SRE.NUMBER < 2))
#Have + and - effects
1-(sum(ACT.SRE.INCREASE.EFFECT==0) - sum(ACT.SRE.NUMBER < 2) + sum(ACT.SRE.DECREASE.EFFECT==0) - sum(ACT.SRE.NUMBER < 2))/(1520 - sum(ACT.SRE.NUMBER < 2))


###Collect data
SUB.FREQ <- matrix(0,nrow=6,ncol=3)
SUB.FREQ[,1] <- c(1520,1520,1520,1520-sum(ACT.NUMBER < 2),1520-sum(ACT.ERE.NUMBER < 2),1520-sum(ACT.SRE.NUMBER < 2))
SUB.FREQ[,2] <- c(sum(apply(AA.EFFECT.TG < 0,1,sum)==0),sum(apply(AA.EFFECT.TG[,1:8000] > 0,1,sum)==0),sum(apply(AA.EFFECT.TG[,8001:16000] > 0,1,sum)==0),2*(sum(ACT.INCREASE.EFFECT==0) - sum(ACT.NUMBER < 2)),2*(sum(ACT.ERE.INCREASE.EFFECT==0) - sum(ACT.ERE.NUMBER < 2)),2*(sum(ACT.SRE.INCREASE.EFFECT==0) - sum(ACT.SRE.NUMBER < 2)))
SUB.FREQ[,3] <- SUB.FREQ[,1] - SUB.FREQ[,2]
rownames(SUB.FREQ) <- c("ALL.BOTH","ALL.ERE","ALL.SRE","ACT.BOTH","ACT.ERE","ACT.SRE")
colnames(SUB.FREQ) <- c("AVAILABLE","UNIDIRECTIONAL","MULTIDIMENSIONAL")

##pdf("SUB.DIRECTION.pdf")
#par(pty="s")
#barplot(t(SUB.FREQ[,c(1,3)]),beside = TRUE)
#
#par(mfrow=c(3,2))
#hist(apply(AA.EFFECT.TG < 0,1,sum)/16000,breaks=seq(0,1,0.05),xlab="",main="TE")
#hist(ACT.INCREASE.EFFECT/(ACT.INCREASE.EFFECT + ACT.DECREASE.EFFECT),breaks=seq(0,1,0.05),xlab="",main="ACT TE")
#hist(apply(AA.EFFECT.TG[,1:8000] > 0,1,sum)/8000,breaks=seq(0,1,0.05),xlab="",main="TE ERE")
#hist(ACT.ERE.INCREASE.EFFECT/(ACT.ERE.INCREASE.EFFECT + ACT.ERE.DECREASE.EFFECT),breaks=seq(0,1,0.05),xlab="",main="ACT TE ERE")
#hist(apply(AA.EFFECT.TG[,8001:16000] > 0,1,sum)/8000,breaks=seq(0,1,0.05),xlab="",main="TE SRE")
#hist(ACT.SRE.INCREASE.EFFECT/(ACT.SRE.INCREASE.EFFECT + ACT.SRE.DECREASE.EFFECT),breaks=seq(0,1,0.05),xlab="",main="ACT TE SRE")
##dev.off()

###Alternative visualization
#Mutation effects across all 8000 backgrounds
ERE.ALL.DIR.MEAN <- apply(AA.EFFECT.TG[,1:8000],1,mean)
ERE.ALL.DIR.POS  <- apply(AA.EFFECT.TG[,1:8000],1,FUN=function(x) {sum(x > 0)})
ERE.ALL.DIR.NEG  <- apply(AA.EFFECT.TG[,1:8000],1,FUN=function(x) {sum(x < 0)})

SRE.ALL.DIR.MEAN <- apply(AA.EFFECT.TG[,8001:16000],1,mean)
SRE.ALL.DIR.POS  <- apply(AA.EFFECT.TG[,8001:16000],1,FUN=function(x) {sum(x > 0)})
SRE.ALL.DIR.NEG  <- apply(AA.EFFECT.TG[,8001:16000],1,FUN=function(x) {sum(x < 0)})

#Mutation effect across activator/near activator backgrounds
ERE.ACT.DIR.MEAN <- sapply(1:1520,FUN = function(x) {
  mean(AA.EFFECT.TG[x,which(AA.EFFECT.START.TG[x,1:8000] < TE.THRESH.NULL | AA.EFFECT.STOP.TG[x,1:8000] < TE.THRESH.NULL)])})
ERE.ACT.DIR.POS <- sapply(1:1520,FUN = function(x) {
  sum(AA.EFFECT.TG[x,which(AA.EFFECT.START.TG[x,1:8000] < TE.THRESH.NULL | AA.EFFECT.STOP.TG[x,1:8000] < TE.THRESH.NULL)] > 0) })
ERE.ACT.DIR.NEG <- sapply(1:1520,FUN = function(x) {
  sum(AA.EFFECT.TG[x,which(AA.EFFECT.START.TG[x,1:8000] < TE.THRESH.NULL | AA.EFFECT.STOP.TG[x,1:8000] < TE.THRESH.NULL)] < 0) })

SRE.ACT.DIR.MEAN <- sapply(1:1520,FUN = function(x) {
  mean(AA.EFFECT.TG[x,which(AA.EFFECT.START.TG[x,8001:16000] < TE.THRESH.NULL | AA.EFFECT.STOP.TG[x,8001:16000] < TE.THRESH.NULL)])})
SRE.ACT.DIR.POS <- sapply(1:1520,FUN = function(x) {
  sum(AA.EFFECT.TG[x,which(AA.EFFECT.START.TG[x,8001:16000] < TE.THRESH.NULL | AA.EFFECT.STOP.TG[x,8001:16000] < TE.THRESH.NULL)] > 0) })
SRE.ACT.DIR.NEG <- sapply(1:1520,FUN = function(x) {
  sum(AA.EFFECT.TG[x,which(AA.EFFECT.START.TG[x,8001:16000] < TE.THRESH.NULL | AA.EFFECT.STOP.TG[x,8001:16000] < TE.THRESH.NULL)] < 0) })

#Collect data
AA.EFFECT.LONG.TG <- data.frame(
  "STATE1"       = rep(factor(substr(rownames(AA.EFFECT.TG),4,4),levels = AAs),2),
  "STATE2"       = rep(factor(substr(rownames(AA.EFFECT.TG),5,5),levels = AAs),2),
  "SITE"         = rep(as.numeric(substr(rownames(AA.EFFECT.TG),2,2)),2),
  "RE"           = c(rep("S",1520),rep("E",1520)),
  "ALL.DIR.MEAN" = -1*c(SRE.ALL.DIR.MEAN,ERE.ALL.DIR.MEAN)/TE.THRESH.NULL,
  "ALL.DIR.POS"  = c(SRE.ALL.DIR.POS,ERE.ALL.DIR.POS),
  "ALL.DIR.NEG"  = c(SRE.ALL.DIR.NEG,ERE.ALL.DIR.POS),
  "ALL.DIR.SIZE" = c(SRE.ALL.DIR.POS + SRE.ALL.DIR.NEG,ERE.ALL.DIR.POS + ERE.ALL.DIR.NEG),
  "ALL.DIR.PROP" = c(SRE.ALL.DIR.POS/(SRE.ALL.DIR.POS + SRE.ALL.DIR.NEG),ERE.ALL.DIR.POS/(ERE.ALL.DIR.POS + ERE.ALL.DIR.NEG)),
  "ACT.DIR.MEAN" = -1*c(SRE.ACT.DIR.MEAN,ERE.ACT.DIR.MEAN)/TE.THRESH.NULL,
  "ACT.DIR.POS"  = c(SRE.ACT.DIR.POS,ERE.ACT.DIR.POS),
  "ACT.DIR.NEG"  = c(SRE.ACT.DIR.NEG,ERE.ACT.DIR.POS),
  "ACT.DIR.SIZE" = c(SRE.ACT.DIR.POS + SRE.ACT.DIR.NEG,ERE.ACT.DIR.POS + ERE.ACT.DIR.NEG),
  "ACT.DIR.PROP" = c(SRE.ACT.DIR.POS/(SRE.ACT.DIR.POS + SRE.ACT.DIR.NEG),ERE.ACT.DIR.POS/(ERE.ACT.DIR.POS + ERE.ACT.DIR.NEG))
)

#hist(AA.EFFECT.LONG.TG$ALL.DIR.PROP[AA.EFFECT.LONG.TG])

##pdf("SUB.DIRECTION.2.pdf")
#par(pty="s")
#theme_set(theme_bw(16))
#p1 <- AA.EFFECT.LONG.TG %>%
#  ggplot(aes(x=ALL.DIR.MEAN,ALL.DIR.PROP,color=RE)) +
#  geom_point(size=log10(ALL.DIR.SIZE)/4,fill="black") +
#  theme(legend.position = "none") +
#  scale_x_continuous(name=paste("Average Effect of","\n","Mutation on Binding"), breaks=seq(-0.1,0.1,0.02)) +
#  scale_y_continuous(name=paste("% of Backgrounds","\n","Mutation Increases Binding"), breaks=seq(0,1,0.1)) + 
#  scale_color_manual(values=c("darkorchid4","forestgreen"))
#ggMarginal(p1, type="histogram", size = 5, groupColour = TRUE, groupFill = TRUE, 
#           xparams = list(breaks=seq(-0.1,0.1,0.01)), 
#           yparams = list(breaks=c(-.1,.0001,.1,.2,.3,.4,.5,.6,.7,.8,.9,.9999,1.1)))
##dev.off()
#
##pdf("SUB.DIRECTION.3.pdf")
#par(pty="s")
#p2 <- AA.EFFECT.LONG.TG %>%
#  ggplot(aes(x=ACT.DIR.MEAN,ACT.DIR.PROP,color=RE)) +
#  geom_point(size=log10(ACT.DIR.SIZE)/2,fill="black") +
#  theme(legend.position = "none") +
#  scale_x_continuous(name=paste("Average Effect of","\n","Mutation on Binding")) +
#  scale_y_continuous(name=paste("% of Backgrounds","\n","Mutation Increases Binding"), breaks=seq(0,1,0.1)) + 
#  scale_color_manual(values=c("darkorchid4","forestgreen"))
#ggMarginal(p2, type="histogram", size = 5, groupColour = TRUE, groupFill = TRUE, 
#           xparams = list(breaks=seq(-0.15,0.15,0.025)), 
#           yparams = list(breaks=c(-.1,.0001,.1,.2,.3,.4,.5,.6,.7,.8,.9,.9999,1.1)))
##dev.off()

##pdf("Random.plots.pdf")
##Binding 3º model
par(pty="s")
EFFECT <- -1*TE.COEFS.ADJ[1:80]/(-1*TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0)
vioplot(EFFECT[1:20],EFFECT[21:40],EFFECT[41:60],EFFECT[61:80],ylim=c(-0.2,0.4),horizontal = TRUE,col="white",main="Main Binding")
text(-1*TE.COEFS.ADJ[1:80]/(-1*TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0),jitter(c(rep(1.25,20),rep(2.25,20),rep(3.25,20),rep(4.25,20))),labels=AAs,cex=0.4)
EFFECT <- -1*TE.COEFS.ADJ[161:2560]/(-1*TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0)
SHOW <- which(abs(EFFECT) > 0.05)
vioplot(EFFECT[1:400],EFFECT[401:800],EFFECT[801:1200],EFFECT[1201:1600],EFFECT[1601:2000],EFFECT[2001:2400],ylim=c(-0.2,0.4),horizontal = TRUE,col="white",main="Pairwise Binding")
text(-1*TE.COEFS.ADJ[161:2560][SHOW]/(-1*TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0),jitter(c(rep(1.25,400),rep(2.25,400),rep(3.25,400),rep(4.25,400),rep(5.25,400),rep(6.25,400))[SHOW]),labels=rep(paste0(AAs,rep(AAs,each=20)),6)[SHOW],cex=0.4)
EFFECT <- -1*TE.COEFS.ADJ[4961:36960]/(-1*TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0)
SHOW <- which(abs(EFFECT) > 0.05)
vioplot(EFFECT[1:8000],EFFECT[8001:16000],EFFECT[16001:24000],EFFECT[24001:32000],ylim=c(-0.2,0.4),horizontal = TRUE,col="white",main="Third order Binding")
text(-1*TE.COEFS.ADJ[4961:36960][SHOW]/(-1*TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0),jitter(c(rep(1.25,8000),rep(2.25,8000),rep(3.25,8000),rep(4.25,8000))[SHOW]),labels=rep(paste0(AAs,rep(AAs,each=20),rep(AAs,each=400)),4)[SHOW],cex=0.4)

EFFECT <- -1*TE.COEFS.ADJ[c(1:80,161:2560,4961:36960)]/(-1*TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0)
vioplot(EFFECT[1:80],EFFECT[81:2480],EFFECT[2481:34480],ylim=c(-0.2,0.4),col="white",main="Binding Effects")
EFFECT <- -1*TE.COEFS.ADJ[1:80]/(-1*TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0); SHOW <- which(EFFECT >= min(tail(sort(EFFECT),20)) | EFFECT <= max(head(sort(EFFECT),20)) )
text(jitter(rep(1,20),10),EFFECT[SHOW],  labels=c(paste0(AAs,"xxx"),paste0("x",AAs,"xx"),paste0("xx",AAs,"x"),paste0("xxx",AAs))[SHOW],cex=0.8)
EFFECT <- -1*TE.COEFS.ADJ[161:2460]/(-1*TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0); SHOW <- which(EFFECT >= min(tail(sort(EFFECT),20)) | EFFECT <= max(head(sort(EFFECT),20)) )
text(jitter(rep(2,20),6),EFFECT[SHOW],  labels=c(paste0(AAs,rep(AAs,each=20),"xx"),paste0(AAs,"x",rep(AAs,each=20),"x"),paste0(AAs,"xx",rep(AAs,each=20)),paste0("x",AAs,rep(AAs,each=20),"x"),paste0("x",AAs,"x",rep(AAs,each=20)),paste0("xx",AAs,rep(AAs,each=20)))[SHOW],cex=0.8)
EFFECT <- -1*TE.COEFS.ADJ[4961:36960]/(-1*TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0); SHOW <- which(EFFECT >= min(tail(sort(EFFECT),20)) | EFFECT <= max(head(sort(EFFECT),20)) )
text(jitter(rep(3,20),6),EFFECT[SHOW],  labels=c(paste0(AAs,rep(AAs,each=20),rep(AAs,each=400),"x"),paste0(AAs,rep(AAs,each=20),"x",rep(AAs,each=400)),paste0(AAs,"x",rep(AAs,each=20),rep(AAs,each=400)),paste0("x",AAs,rep(AAs,each=20),rep(AAs,each=400)))[SHOW],cex=0.8)

#Specificity 3º model
EFFECT <- -1*TE.COEFS.ADJ[81:160]/(-1*TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0)
vioplot(EFFECT[1:20],EFFECT[21:40],EFFECT[41:60],EFFECT[61:80],ylim=c(-0.2,0.2),horizontal = TRUE,col="white",main="Main Specificity")
text(-1*TE.COEFS.ADJ[81:160]/(-1*TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0),jitter(c(rep(1.25,20),rep(2.25,20),rep(3.25,20),rep(4.25,20))),labels=AAs,cex=0.4)
EFFECT <- -1*TE.COEFS.ADJ[2561:4960]/(-1*TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0)
SHOW <- which(abs(EFFECT) > 0.045)
vioplot(EFFECT[1:400],EFFECT[401:800],EFFECT[801:1200],EFFECT[1201:1600],EFFECT[1601:2000],EFFECT[2001:2400],ylim=c(-0.2,0.2),horizontal = TRUE,col="white",main="Pairwise Specificity")
text(-1*TE.COEFS.ADJ[2561:4960][SHOW]/(-1*TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0),jitter(c(rep(1.25,400),rep(2.25,400),rep(3.25,400),rep(4.25,400),rep(5.25,400),rep(6.25,400))[SHOW]),labels=rep(paste0(AAs,rep(AAs,each=20)),6)[SHOW],cex=0.4)
EFFECT <- -1*TE.COEFS.ADJ[36961:68960]/(-1*TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0)
SHOW <- which(abs(EFFECT) > .045)
vioplot(EFFECT[1:8000],EFFECT[8001:16000],EFFECT[16001:24000],EFFECT[24001:32000],ylim=c(-0.2,0.2),horizontal = TRUE,col="white",main="Third order Specificity")
text(-1*TE.COEFS.ADJ[36961:68960][SHOW]/(-1*TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0),jitter(c(rep(1.25,8000),rep(2.25,8000),rep(3.25,8000),rep(4.25,8000))[SHOW]),labels=rep(paste0(AAs,rep(AAs,each=20),rep(AAs,each=400)),4)[SHOW],cex=0.4)

EFFECT <- -1*TE.COEFS.ADJ[c(81:160,2561:4960,36961:68960)]/(-1*TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0)
vioplot(EFFECT[1:80],EFFECT[81:2480],EFFECT[2481:34480],ylim=c(-0.2,0.2),col="white",main="Specificity Effects")
EFFECT <- -1*TE.COEFS.ADJ[81:160]/(-1*TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0); SHOW <- which(EFFECT >= min(tail(sort(EFFECT),20)) | EFFECT <= max(head(sort(EFFECT),20)) )
text(jitter(rep(1,20),10),EFFECT[SHOW],  labels=c(paste0(AAs,"xxx"),paste0("x",AAs,"xx"),paste0("xx",AAs,"x"),paste0("xxx",AAs))[SHOW],cex=0.8)
EFFECT <- -1*TE.COEFS.ADJ[2561:4960]/(-1*TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0); SHOW <- which(EFFECT >= min(tail(sort(EFFECT),20)) | EFFECT <= max(head(sort(EFFECT),20)) )
text(jitter(rep(2,20),6),EFFECT[SHOW],  labels=c(paste0(AAs,rep(AAs,each=20),"xx"),paste0(AAs,"x",rep(AAs,each=20),"x"),paste0(AAs,"xx",rep(AAs,each=20)),paste0("x",AAs,rep(AAs,each=20),"x"),paste0("x",AAs,"x",rep(AAs,each=20)),paste0("xx",AAs,rep(AAs,each=20)))[SHOW],cex=0.8)
EFFECT <- -1*TE.COEFS.ADJ[36961:68960]/(-1*TE.THRESH.WEAK+TE.COEFS.EFFECT.ADJ.B0); SHOW <- which(EFFECT >= min(tail(sort(EFFECT),20)) | EFFECT <= max(head(sort(EFFECT),20)) )
text(jitter(rep(3,20),6),EFFECT[SHOW],  labels=c(paste0(AAs,rep(AAs,each=20),rep(AAs,each=400),"x"),paste0(AAs,rep(AAs,each=20),"x",rep(AAs,each=400)),paste0(AAs,"x",rep(AAs,each=20),rep(AAs,each=400)),paste0("x",AAs,rep(AAs,each=20),rep(AAs,each=400)))[SHOW],cex=0.8)

#Binding 1º Model
EFFECT <- -1*ME.COEFS.ADJ[1:80]/(-1*ME.THRESH.WEAK+ME.COEFS.EFFECT.ADJ.B0)
vioplot(EFFECT[1:20],EFFECT[21:40],EFFECT[41:60],EFFECT[61:80],ylim=c(-0.25,0.55),horizontal = TRUE,col="white",main="Main Binding")
text(-1*ME.COEFS.ADJ[1:80]/(-1*ME.THRESH.WEAK+ME.COEFS.EFFECT.ADJ.B0),jitter(c(rep(1.25,20),rep(2.25,20),rep(3.25,20),rep(4.25,20))),labels=AAs,cex=0.4)

vioplot(EFFECT,ylim=c(-0.3,0.6),col="white",main="1º Binding Effects")
SHOW <- which(EFFECT >= min(tail(sort(EFFECT),20)) | EFFECT <= max(head(sort(EFFECT),20)) )
text(jitter(rep(1,20),10),EFFECT[SHOW],  labels=c(paste0(AAs,"xxx"),paste0("x",AAs,"xx"),paste0("xx",AAs,"x"),paste0("xxx",AAs))[SHOW],cex=0.8)

#Specificity 1º model
EFFECT <- -1*ME.COEFS.ADJ[81:160]/(-1*ME.THRESH.WEAK+ME.COEFS.EFFECT.ADJ.B0)
vioplot(EFFECT[1:20],EFFECT[21:40],EFFECT[41:60],EFFECT[61:80],ylim=c(-0.2,0.2),horizontal = TRUE,col="white",main="Main Specificity")
text(-1*ME.COEFS.ADJ[81:160]/(-1*ME.THRESH.WEAK+ME.COEFS.EFFECT.ADJ.B0),jitter(c(rep(1.25,20),rep(2.25,20),rep(3.25,20),rep(4.25,20))),labels=AAs,cex=0.4)

vioplot(EFFECT,ylim=c(-0.2,0.2),col="white",main="1º Specificity Effects")
SHOW <- which(EFFECT >= min(tail(sort(EFFECT),20)) | EFFECT <= max(head(sort(EFFECT),20)) )
text(jitter(rep(1,20),10),EFFECT[SHOW],  labels=c(paste0(AAs,"xxx"),paste0("x",AAs,"xx"),paste0("xx",AAs,"x"),paste0("xxx",AAs))[SHOW],cex=0.8)
##dev.off()

#######################################################
##21. Effects over historical evolutionary trajectory##
#######################################################
HISTORICAL.GENOS <- c("EEGKA","EGGKA","EESKA","EEGKV","EGSKA","EGGKV","EESKV","EGSKV","SEGKA","SGGKA","SESKA","SEGKV","SGSKA","SGGKV","SESKV","SGSKV")
HISTORICAL.EFFECTS <- EFFECT.TABLE.TE[EFFECT.TABLE.TE$SEQ %in% HISTORICAL.GENOS,]
HISTORICAL.EFFECTS[,2:32] <- -1*HISTORICAL.EFFECTS[,2:32]

HISTORICAL.B0 <- mean(HISTORICAL.EFFECTS$LINK)
HISTORICAL.S0 <- mean(HISTORICAL.EFFECTS$LINK[9:16]) - HISTORICAL.B0
HISTORICAL.B1.X1 <- mean(HISTORICAL.EFFECTS$LINK[c(5:8,13:16)]) - HISTORICAL.B0
HISTORICAL.B1.X2 <- mean(HISTORICAL.EFFECTS$LINK[c(3:4,7:8,11:12,15:16)]) - HISTORICAL.B0
HISTORICAL.B1.X4 <- mean(HISTORICAL.EFFECTS$LINK[c(2,4,6,8,10,12,14,16)]) - HISTORICAL.B0
HISTORICAL.S1.X1 <- mean(HISTORICAL.EFFECTS$LINK[c(13:16)]) - mean(HISTORICAL.EFFECTS$LINK[c(5:8)]) - HISTORICAL.S0
HISTORICAL.S1.X2 <- mean(HISTORICAL.EFFECTS$LINK[c(11:12,15:16)]) - mean(HISTORICAL.EFFECTS$LINK[c(3:4,7:8)]) - HISTORICAL.S0
HISTORICAL.S1.X4 <- mean(HISTORICAL.EFFECTS$LINK[c(10,12,14,16)]) - mean(HISTORICAL.EFFECTS$LINK[c(2,4,6,8)]) - HISTORICAL.S0
HISTORICAL.B2.X1X2 <- mean(HISTORICAL.EFFECTS$LINK[c(7:8,15:16)]) - HISTORICAL.B1.X1 - HISTORICAL.B1.X2 - HISTORICAL.B0
HISTORICAL.B2.X1X4 <- mean(HISTORICAL.EFFECTS$LINK[c(6,8,14,16)]) - HISTORICAL.B1.X1 - HISTORICAL.B1.X4 - HISTORICAL.B0
HISTORICAL.B2.X2X4 <- mean(HISTORICAL.EFFECTS$LINK[c(4,8,12,16)]) - HISTORICAL.B1.X2 - HISTORICAL.B1.X4 - HISTORICAL.B0
HISTORICAL.S2.X1X2 <- mean(HISTORICAL.EFFECTS$LINK[c(15:16)]) - mean(HISTORICAL.EFFECTS$LINK[c(7:8)]) - HISTORICAL.S1.X1 - HISTORICAL.S1.X2 - HISTORICAL.S0
HISTORICAL.S2.X1X4 <- mean(HISTORICAL.EFFECTS$LINK[c(14,16)]) - mean(HISTORICAL.EFFECTS$LINK[c(6,8)]) - HISTORICAL.S1.X1 - HISTORICAL.S1.X4 - HISTORICAL.S0
HISTORICAL.S2.X2X4 <- mean(HISTORICAL.EFFECTS$LINK[c(12,16)]) - mean(HISTORICAL.EFFECTS$LINK[c(4,8)]) - HISTORICAL.S1.X2 - HISTORICAL.S1.X4 - HISTORICAL.S0
HISTORICAL.B3.X1X2X4 <- mean(HISTORICAL.EFFECTS$LINK[c(8,16)]) - HISTORICAL.B1.X1 - HISTORICAL.B1.X2 - HISTORICAL.B1.X4 - HISTORICAL.B2.X1X2 - HISTORICAL.B2.X1X4 - HISTORICAL.B2.X2X4- HISTORICAL.B0
HISTORICAL.S3.X1X2X4 <- mean(HISTORICAL.EFFECTS$LINK[16]) - mean(HISTORICAL.EFFECTS$LINK[8]) - HISTORICAL.S1.X1 - HISTORICAL.S1.X2 - HISTORICAL.S1.X4 - HISTORICAL.S2.X1X2 - HISTORICAL.S2.X1X4 - HISTORICAL.S2.X2X4- HISTORICAL.S0

HISTORICAL.EFFECTS[,2:32] <- HISTORICAL.EFFECTS[,2:32]/(TE.THRESH.WEAK + TE.COEFS.EFFECT.ADJ.B0)

##Historical route changes
HISTORICAL.EFFECTS[HISTORICAL.EFFECTS$SEQ=="EGGKA",3:32] - HISTORICAL.EFFECTS[HISTORICAL.EFFECTS$SEQ=="EEGKA",3:32]
HISTORICAL.EFFECTS[HISTORICAL.EFFECTS$SEQ=="SGGKA",3:32] - HISTORICAL.EFFECTS[HISTORICAL.EFFECTS$SEQ=="SEGKA",3:32]

HISTORICAL.EFFECTS[HISTORICAL.EFFECTS$SEQ=="EGGKV",3:32] - HISTORICAL.EFFECTS[HISTORICAL.EFFECTS$SEQ=="EGGKA",3:32]
HISTORICAL.EFFECTS[HISTORICAL.EFFECTS$SEQ=="SGGKV",3:32] - HISTORICAL.EFFECTS[HISTORICAL.EFFECTS$SEQ=="SGGKA",3:32]

HISTORICAL.EFFECTS[HISTORICAL.EFFECTS$SEQ=="EGSKV",3:32] - HISTORICAL.EFFECTS[HISTORICAL.EFFECTS$SEQ=="EGGKV",3:32]
HISTORICAL.EFFECTS[HISTORICAL.EFFECTS$SEQ=="SGSKV",3:32] - HISTORICAL.EFFECTS[HISTORICAL.EFFECTS$SEQ=="SGGKV",3:32]


##Epistasis necessary (check by removing epi)
HISTORY.ACTIVATORS <- rowSums(HISTORICAL.EFFECTS[,c(3:32)]) > 1; names(HISTORY.ACTIVATORS) <- HISTORICAL.EFFECTS$SEQ
HISTORY.EPI.NEC  <- rowSums(HISTORICAL.EFFECTS[HISTORY.ACTIVATORS,c(3:12)]) < 1; names(HISTORY.EPI.NEC) <- names(HISTORY.ACTIVATORS[which(HISTORY.ACTIVATORS == TRUE)])
HISTORY.2EPI.NEC <- rowSums(HISTORICAL.EFFECTS[HISTORY.ACTIVATORS,c(3:12,25:32)]) < 1; names(HISTORY.2EPI.NEC) <- names(HISTORY.ACTIVATORS[which(HISTORY.ACTIVATORS == TRUE)])
HISTORY.3EPI.NEC <- rowSums(HISTORICAL.EFFECTS[HISTORY.ACTIVATORS,c(3:12,13:24)]) < 1; names(HISTORY.3EPI.NEC) <- names(HISTORY.ACTIVATORS[which(HISTORY.ACTIVATORS == TRUE)])

#Epistasis necessary (check by lower order model)
DT.JOINT[which(AAseq %in% HISTORICAL.GENOS),list(AAseq,PRED.TE.CLASS,PRED.PE.CLASS,PRED.ME.CLASS)]
