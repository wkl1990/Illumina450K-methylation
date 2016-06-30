#Programmer:  Chao Chen
#Location: UCHICAGO
#email  chenchaor@gmail.com
#Code: R
#Program name: pvca.R
#Date: May 26, 2010

########## load libraries ##########
library(lme4)
####### Edit these variables according to user defined parameters and the path to your data and data files names ##########

#myPath <- "C:/Users/WangKangli/Desktop/methylationdata"
#theGene_expression_file <- "RAW_62REPLICATES.TXT"
#theMethylation_file <- "methylation.csv"
#theExperiment_data_file <- "expinfo_tab_delimited.TXT"
#theMethylation_data_file <- "sample.csv"
pct_threshold = .5876 # Amount of variability desired to be explained by the principal components.  Set to match the results in book chapter and SAS code.  User can adjust this to a higher (>= 0.8) number but < 1.0

### In addition, be sure to modify the mixed linear model by adding the appropriate random effects terms in the model

################################################

#theGEDFilePath = paste(myPath,theMethylation_file, sep="/")
#theExpDataFilePath = paste(myPath,theMethylation_data_file, sep="/")

########## Load data ##########

theDataMatrix <- beta_BMIQ_filter
#theDataMatrix <- read.delim(theGEDFilePath, row.names = 1, header = TRUE, sep = ",")
dataRowN <- nrow(theDataMatrix)
dataColN <- ncol(theDataMatrix)

########## Center the data (center rows) ##########
theDataMatrixCentered <- matrix(data = 0, nrow = dataRowN, ncol = dataColN)
theDataMatrixCentered_transposed = apply(theDataMatrix, 1, scale, center = TRUE, scale = FALSE)
theDataMatrixCentered = t(theDataMatrixCentered_transposed)

#exp_design <- read.delim(theExpDataFilePath, sep = ",", header = TRUE, row.names = 1)
#exp_design <- exp_design[,-c(3,7)]
exp_design <- targets[,c("individual", "tissue", "pregnancy", "batch", "position")]
rownames(exp_design) <- rownames(targets)
expDesignRowN <- nrow(exp_design)
expDesignColN <- ncol(exp_design)
myColNames <- names(exp_design)


########## Compute correlation matrix ##########

theDataCor <- cor(theDataMatrixCentered)

########## Obtain eigenvalues ##########

eigenData <- eigen(theDataCor)
eigenValues = eigenData$values
ev_n <- length(eigenValues)
eigenVectorsMatrix = eigenData$vectors
eigenValuesSum = sum(eigenValues)
percents_PCs = eigenValues /eigenValuesSum 

########## Merge experimental file and eigenvectors for n components ##########

my_counter_2 = 0
my_sum_2 = 1
for (i in ev_n:1){
my_sum_2  = my_sum_2 - percents_PCs[i]
	if ((my_sum_2) <= pct_threshold ){
		my_counter_2 = my_counter_2 + 1
	}

}
if (my_counter_2 < 3){
	pc_n  = 3

}else {
	pc_n = my_counter_2 
}

# pc_n is the number of principal components to model

pc_data_matrix <- matrix(data = 0, nrow = (expDesignRowN*pc_n), ncol = 1)
mycounter = 0
for (i in 1:pc_n){
	for (j in 1:expDesignRowN){
	mycounter <- mycounter + 1
		pc_data_matrix[mycounter,1] = eigenVectorsMatrix[j,i]

	}
}

AAA <- exp_design[rep(1:expDesignRowN,pc_n),]

Data <- cbind(AAA,pc_data_matrix)
####### Edit these variables according to your factors #######

#Data$sex <- as.factor(Data$sex)
#Data$age <- as.factor(Data$age)
#Data$Index <- as.factor(Data$Index)

variables <- c(colnames(exp_design))
    for (i in 1:length(variables)) {
        Data$variables[i] <- as.factor(Data$variables[i])
    }

########## Mixed linear model ##########
op <- options(warn = (-1))
    effects_n = expDesignColN + 1
    #effects_n = expDesignColN + choose(expDesignColN, 2) + 1
    randomEffectsMatrix <- matrix(data = 0, nrow = pc_n, ncol = effects_n)

model.func <- c()
    index <- 1
    for (i in 1:length(variables)) {
        mod = paste("(1|", variables[i], ")", sep = "")
        model.func[index] = mod
        index = index + 1
    }
    #for (i in 1:(length(variables) - 1)) {
    #    for (j in (i + 1):length(variables)) {
    #        mod = paste("(1|", variables[i], ":", variables[j], 
    #            ")", sep = "")
    #        model.func[index] = mod
    #        index = index + 1
    #    }
    #}
    function.mods <- paste(model.func, collapse = " + ")
    
    for (i in 1:pc_n) {
        y = (((i - 1) * expDesignRowN) + 1)
        funct <- paste("pc_data_matrix", function.mods, sep = " ~ ")
        Rm1ML <- lmer(funct, Data[y:(((i - 1) * expDesignRowN) + 
            expDesignRowN), ], REML = TRUE, control=lmerControl(check.nobs.vs.nlev = "ignore",check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"),verbose = FALSE, 
            na.action = na.omit)
        randomEffects <- Rm1ML
        randomEffectsMatrix[i, ] <- c(unlist(VarCorr(Rm1ML)), 
            resid = sigma(Rm1ML)^2)
    }
    effectsNames <- c(names(getME(Rm1ML, "cnms")), "resid")
########## Standardize Variance ##########

randomEffectsMatrixStdze <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
for (i in 1:pc_n){
	mySum = sum(randomEffectsMatrix[i,])
	for (j in 1:effects_n){
		randomEffectsMatrixStdze[i,j] = randomEffectsMatrix[i,j]/mySum	
	}
}

########## Compute Weighted Proportions ##########

randomEffectsMatrixWtProp <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
for (i in 1:pc_n){
	weight = eigenValues[i]/eigenValuesSum
	for (j in 1:effects_n){
		randomEffectsMatrixWtProp[i,j] = randomEffectsMatrixStdze[i,j]*weight
	}
}

########## Compute Weighted Ave Proportions ##########

randomEffectsSums <- matrix(data = 0, nrow = 1, ncol = effects_n)
randomEffectsSums <-colSums(randomEffectsMatrixWtProp)
totalSum = sum(randomEffectsSums)
randomEffectsMatrixWtAveProp <- matrix(data = 0, nrow = 1, ncol = effects_n)

for (j in 1:effects_n){
	randomEffectsMatrixWtAveProp[j] = randomEffectsSums[j]/totalSum 	
	
}

pdf("pvca1.pdf")
bp <- barplot(randomEffectsMatrixWtAveProp,  xlab = "Effects", ylab = "Weighted average proportion variance", ylim= c(0,1.1),col = c("blue"), las=2)

axis(1, at = bp, labels = effectsNames, xlab = "Effects", cex.axis = 0.5, las=2)

## replace the above code of "axis(1, at = bp, labels = effectsNames, xlab = "Effects", cex.axis = 0.5, las=2)" if you want rotate the x axis labels. 
## text(bp, par("usr")[3]-0.02, srt = 45, adj = 1,labels = effectsNames, xpd = TRUE,cex=0.8)

values = randomEffectsMatrixWtAveProp
new_values = round(values , 3)
text(bp,randomEffectsMatrixWtAveProp,labels = new_values, pos=3, cex = 0.8) # place numbers on top of bars 
dev.off()



