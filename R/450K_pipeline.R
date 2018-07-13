####
# Describe: use appropriate packages to preprocess 450K idat
# Environment: R 3.2.2
# Author: Wang Kangli (SKLMG)
# Date: 2016-06-28
####

###Reading 450K idat files with the minfi package (idat path and basename needed)
library(minfi)
path <- "GSE66210_RAW"
list.files(path)
targets <- read.csv("GSE66210_tec_trait.csv")
names(targets)
rownames(targets) <- targets[,1]
targets$Basename <- file.path(path, targets$Basename)
targets$Basename
targets$batch <- substring(targets$Sample_name,1,10)
targets$batch
targets$position <- substring(targets$Sample_name,12)
targets$position

RGsetEx <- read.450k(targets$Basename, extended = TRUE, verbose = TRUE) # extended needed for get NBeads
pData(RGsetEx) <- targets
dim(getRed(RGsetEx))
dim(getGreen(RGsetEx))
manifest <- getManifest(RGsetEx)
manifest
head(getProbeInfo(manifest))
head(getRed(RGsetEx))

###sample filter by detectionP
detP <- detectionP(RGsetEx)
Pfailed <- detP>0.01
colMeans(Pfailed) # Fraction of failed positions per sample
sum(colMeans(Pfailed)>0.01)
cutSample <- colMeans(Pfailed)>0.01
RGsetEx_cutSam <- RGsetEx[,!cutSample]
###find failed probe by detectionP
detP <- detectionP(RGsetEx_cutSam)
Pfailed <- detP>0.01
sum(rowMeans(Pfailed)>0.1) # How many positions failed in >10% of samples?
failedProbesP <-rownames(Pfailed)[rowMeans(Pfailed)>0.1]
###find failed probe by Nbeads use wateRmelon packages
library(wateRmelon)
beadcount <- beadcount(RGsetEx_cutSam)
NBfailed <- is.na(beadcount)
sum(rowMeans(NBfailed)>0.05)
failedProbesNB <-rownames(NBfailed)[rowMeans(NBfailed)>0.05]
###get the list of failed probes
failedProbes <- union(failedProbesP, failedProbesNB) 
###one-step filter: you also can use wateRmelon to filter based on bead count and detection p-values, but the old version is not avaiable for RGsetEx，you can not get MSet after filter.
#RGsetEx.pf <- pfilter(RGsetEx, perCount=5, pnthresh = 0.01, perc = 1, pthresh = 1) 

###preprocess and get beta value
#MSet <- preprocessRaw(RGsetEx_cutSam)
Mset <-preprocessIllumina(RGsetEx_cutSam, bg.correct=TRUE, normalize = c("controls"), reference =1) 
MSet
###remove failed probes
Mset_cutProbe <- Mset[!rownames(Mset)%in% failedProbes,] 
#Mset_cutProbe2 <- Mset[!(rowMeans(Pfailed)>0.1),]
head(getMeth(Mset_cutProbe))
###get beta value
RSet <- ratioConvert(Mset_cutProbe, what = "both", keepCN = TRUE)
RSet
beta <- getBeta(RSet)

###BMIQ
#typeI <- getProbeInfo(Mset_cutProbe, type = "I")
#typeII <- getProbeInfo(Mset_cutProbe, type = "II")
#typeI$type <- 1
#typeII$type <- 2
#design.typeI <- typeI[,c(1,9)]
#design.typeII <- typeII[,c(1,5)]
#design.type <- rbind(design.typeI,design.typeII)
# #design.v <- subset(design.type,Name==rownames(beta),type)
#lookup <- rownames(beta)
#design.v <- design.type[match(lookup,design.type[,1]),2]
##another method
design.type <- got(Mset_cutProbe)
design.v <- sub("II",2,design.type)
design.v <- sub("I",1,design.v)
design.v1 <- as.numeric(design.v1)

BMIQ1 <- apply(beta, 2, function(x) BMIQ(beta.v=x, design.v=design.v, plots=FALSE))
BMIQ2 <- lapply(BMIQ1, function(y) y$nbeta)
beta_BMIQ <- as.matrix(data.frame(BMIQ2))
##BMIQ can be used for a MethyLumiSet or MethylSet, but it got error when I use MethySet
#MethySet_BMIQ <- BMIQ(MethySet) # got error, else it will be more useful
#MethylumiSet_BMIQ <- BMIQ(MethylumiSet)

###remove cross and snp probes
filter.probes<- read.csv(file="list167085_for45Kfilter.csv", header=TRUE)
dim(filter.probes)
filter.p <- filter.probes[,1]
length(filter.p) # 167085
beta_BMIQ_filter <- beta_BMIQ[rownames(beta_BMIQ)%in% filter.p,] 
beta_BMIQ_nonfilter <- beta_BMIQ[!rownames(beta_BMIQ)%in% filter.p,] 
dim(beta_BMIQ_filter) #features.460350 samples.146 

###PCA and PVCA plot
source("pca_for450K_pipeline.R")
source("pvca1_for450K_pipeline.R")
source("pvca2_for450K_pipeline.R")

###QC plot
qc <- getQC(Mset_cutProbe)
head(qc)
pdf("QC.pdf")
plotQC(qc)
densityPlot(Mset_cutProbe, sampGroups = pData(Mset_cutProbe)$rep)
densityBeanPlot(Mset_cutProbe, sampGroups = pData(Mset_cutProbe)$rep)
controlStripPlot(RGsetEx_cutSam, controls="BISULFITE CONVERSION II", sampNames=pData(RGsetEx_cutSam)$description)
#qcReport(RGsetEx_cutSam, sampNames=colnames(RGsetEx_cutSam), sampGroups=pData(RGsetEx_cutSam)$rep, pdf= "qcReport.pdf")
### Type plot and mds plot
plotBetasByType(Mset_cutProbe[,1], main=colnames(Mset_cutProbe)[1])
mdsPlot(Mset_cutProbe, numPositions = 1000, sampNames=colnames(Mset_cutProbe), sampGroups=pData(Mset_cutProbe)$rep)
dev.off()
