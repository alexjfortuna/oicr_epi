library(data.table)

TissSourceSite.table <- read.delim("data/tissueSourceSite.tsv", header = TRUE)

GB.df <- TissSourceSite.table[grep("Glioblastoma",TissSourceSite.table$Study.Name),]

GB.TSS.list <- paste(GB.df[,1], collapse ="|")

BLGG.df <- TissSourceSite.table[grep("Glioma",TissSourceSite.table$Study.Name),]

BLGG.TSS.list <- paste(BLGG.df[,1], collapse ="|")

#script to take text file and process into R list
#read in the data
PatientIDs <- scan("data/PATIENTID.txt", what ="character")

#seperate elements by one or more whitespace
y <- strsplit(PatientIDs,"[[:space:]]+")

#split list at "-" to make list of lists
temp <- strsplit(as.character(y),"-")

#take list of lists and make matrix
mat  <- matrix(unlist(temp), ncol=7, byrow=TRUE)

#take matrix and make data frame
df  <- as.data.frame(mat)

colnames(df) <- c("TCGA","TSS","Participant","Sample","Portion","Plate","Center")

#create new data frame for Glioblasta TSS samples
GB.TCGA.df <- df[grep(GB.TSS.list,df$TSS),]


#create new data frame for Brain lower Grade Glioma
BLGG.TCGA.df <- df[grep(BLGG.TSS.list,df$TSS),]

library(data.table)
TCGA.clinical <- read.csv("data/clinical_PANCAN_patient_with_followup.tsv", na.strings=c("[Not Applicable]","[Not Available]"), sep = "\t", stringsAsFactors=FALSE, header = TRUE)

TCGA.clinical2 <- subset(TCGA.clinical, select = c(2,5,7,8))

Cat.BLGG.clinical.IDs <- paste(BLGG.TCGA.df$TCGA, BLGG.TCGA.df$TSS, BLGG.TCGA.df$Participant, sep = "-")

#Create Clinical Data Frame for Brain Lower Grade Glioma Patients
TCGA_BLGG_clinical <- subset(TCGA.clinical2, bcr_patient_barcode %in% Cat.BLGG.clinical.IDs)
#513 Patients
#!/usr/bin/Rscript
print(getwd())

#Survival_GeneMeth########################################
##########################################################

#Author: Alex Fortuna
#Date_started: January 26 2017

#Description:
#Investigate if there are any predictors of survival based on
#a median expression cutoff for Ion Channel Genes and Brain Cancer

#Preamble#################################################

###Libraries##############################################
library(data.table)
library(survival)
library(ggplot2)
library(ggfortify)
library(cluster)
source("bin/220217.r")
source("bin/200317_clinicalBLGG.r")
###Functions##############################################
insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  existingDF
}

###Data###################################################
clin <- as.data.frame(TCGA_BLGG_clinical)
rownames(clin) <- clin$bcr_patient_barcode
#Read in expresion data for BLGG Patients for ICGs
meth.dt <- fread("data/jhu_BLGG_icg.tsv", header = TRUE, sep = "\t")
#Merge first and second columns, GeneName and ProbesetID
names(meth.dt)[1] <- "GeneName"
names(meth.dt)[2] <- "ProbesetID"
Cat.BLGG.Pat.IDs <- paste(BLGG.TCGA.df$TCGA, BLGG.TCGA.df$TSS, BLGG.TCGA.df$Participant, sep="-")
colnames(meth.dt) <- c("GeneName", "ProbesetID", Cat.BLGG.Pat.IDs)
#select only columns which we have Clinical Info for
meth.clin <- subset(meth.dt, select = c("GeneName", "ProbesetID", clin$bcr_patient_barcode))
#change gene names to rows so have just a matrix
#of expression values
meth.clin$feature <- paste(meth.clin$GeneName, meth.clin$ProbesetID, sep=":")
#remove rows with NA values
clin.no.na <- na.omit(meth.clin)
rownames(meth.clin) <- meth.clin$feature
#After adding rownames as Gene and Probeset name of interest remove those columns
clin.no.na$GeneName <- NULL ; clin.no.na$ProbesetID <- NULL
#transform to matrix, do not include last column which contains meth.dt$feature
#matrix with rows as patients with relative clinical information.
#event and time columns must be marked
#conver data table to matrix
meth <- as.matrix(clin.no.na[,1:513])
rownames(meth) <- clin.no.na$feature
#> dim(meth) [1] 1665 513
"gene" <- rownames(meth)


##############################################################################
###Survival Analysis Matrix###################################################
##############################################################################

#storing data
data_results <- as.data.frame(matrix(ncol=5))
colnames(data_results) <- c("gene", "coxphHR", "coxPval", "logrankChisq", "logrankPval")


for(i in 1:nrow(meth)){

	#set up gene data matrix with expression and survival data
	###########################################################

	geneid <- rownames(meth)[i]
	data <- as.data.frame(matrix(ncol=8, nrow=513))
	colnames(data) <- c("patient", geneid, "vital_status", "event", "time", "cancer_type", "tag", "tag2")
	data[,2] <- meth[i,]
	data$patient <- colnames(meth)
  data$vital_status <- clin$vital_status
  data$cancer_type <- "Brain Lower Grade Glioma"
  for(i in 1:nrow(data)){
    pat <- data$patient[i]
    z <- which(clin$bcr_patient_barcode %in% pat)
    t <- clin$days_to_death[z]
    data$time[i] <- t
    if(is.na(t)){
      ot <- clin$days_to_last_followup[z]
      data$time[i] <- ot
    }
    data$event[i] <- clin$vital_status[z]
		if(data$event[i]=="Alive"){data$event[i] <- 0}
		if(data$event[i]=="Dead"){data$event[i] <- 1}
	}

	##add tag based on high or low methylation at probe via median
	##############################################################
	median <- median(data[,2])

	if(!(median == 0)){

	for(i in 1:nrow(data)){
		med <- median(data[i,2])
		if(med <= median){
			data$tag[i] <- "Low"
			data$tag2[i] <- 0
		}
		if(med > median){
			data$tag[i] <- "high"
			data$tag2[i] <- 1
		}
		}

	##add boxplot showing differences in methylation between groups
	###############################################################
	name_plot <- paste(geneid,".pdf", sep="")
    pdf(name_plot)
    plot <- ggplot(data, aes(data$tag, data[,2])) + ggtitle(name_plot)
    print(plot + geom_boxplot())
    dev.off()

    ###survival analsys############################################
    ###############################################################
    data$event <- as.numeric(data$event)
    data$time <- as.numeric(data$time)
    S <- Surv(data$time, data$event)
    model <- survfit(Surv(data$time, data$event) ~ data$tag2)

    ###cox proportional hazard test##
    res_cox <- coxph(formula = Surv(data$time, data$event) ~ as.factor(data$tag2))
    ###hazard ratio
    s <- summary(res_cox)
    cox_pvalue <- s$coefficients[5]
    hr <- s$coefficients[1]

    ###Log-rank proportional hazard test##
    res_logrank <- survdiff(formula = Surv(data$time, data$event) ~ data$tag2)
    log_pvalue <- 1 - pchisq(res_logrank$chisq, length(res_logrank$n) - 1)
    chisq <- res_logrank$chisq

    ###store results
	row <- c(as.character(geneid), as.character(hr), as.character(cox_pvalue), as.character(chisq), as.character(log_pvalue))
	data_results <- insertRow(data_results, row, 1)

	##plot survival plots
	###############################################################

	name <- paste(geneid, "surv.pdf", sep="_")
	pdf(name)
	fit <- survfit(Surv(time, event) ~ as.factor(data$tag2), data = data)
	print(autoplot(fit,conf.int = FALSE, censor.size = 0.5))
	dev.off()

	}}

	### WRITING AND ANALYSIS OF RESULTS ##################################################
	data_results <- as.data.frame(data_results)
	data_results <- data_results[-which(is.na(data_results[,1])),]

	data_results$coxphHR <- as.numeric(data_results$coxphHR)
	data_results$coxPval <- as.numeric(data_results$coxPval)
	data_results$logrankChisq <- as.numeric(data_results$logrankChisq)
	data_results$logrankPval <- as.numeric(data_results$logrankPval)

    ## multiple testing correction BH FDR##
    pvals <- data_results$coxPval
    pval_adj <- p.adjust(pvals, method = 'BH')
    data_results$coxPval_adj <- pval_adj

    ## subset to those with adj cox pvalues < 0.1 ##
    data_results_symbols_sig <- subset(data_results, data_results$coxPval_adj <= 0.1)


    ## What are signficant coxphHR genes?
    write.csv(data_results, file= "methylation_Surv_0317_BLGG.csv", quote=F)

##########################################################################33
drugtargets <- read.csv("data/targets_and_families.csv")
drugtargets.icg <- subset(drugtargets, Type %in% c("vgic", "lgic", "other_ic"))

library(FDb.InfiniumMethylation.hg19)

ProbenameIDs <- scan("data/PANCAN_HumanMethylationProbesetIDs.tsv", what ="character")
hm450 <- get450k()
ProbenameIDs <- ProbenameIDs[-1]
probes <- hm450[ProbenameIDs]
NearestGeneprobes <- getNearestGene(probes)
