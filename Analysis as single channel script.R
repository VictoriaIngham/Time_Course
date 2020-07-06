workingdirectory <-'/Users/vickyingham/Dropbox (LSTM)/Vicky Documents/PhD Year 1/GaGa'
workingdirectory2 <-'/Users/vickyingham/Dropbox (LSTM)/Vicky Documents/Post Doc/Exposure Time Course, Transcriptome/Aging Experiment/Aging and PANG'
targetfile <- 'targets.txt'
gene.name.file = '/AGAP_correct2.txt'

suppressMessages(library(limma))

##SET WORKING DIRECTORY AS AGILENT FILE LOCATION##
setwd(workingdirectory2)


##READ IN ALL DATA AND CONVERT PROBES TO AGAP##
targets<- readTargets(targetfile)
gene_names_file<- read.table(paste(workingdirectory2,gene.name.file,sep=""), header = TRUE, fill = TRUE, sep = "\t")
gene_names_file<-as.matrix(gene_names_file)
red_green_data<- read.maimages(targets, source = "agilent.median", annotation = c("Row", "Col", "FeatureNum", "ControlType", "ProbeName", "SystematicName")) 
positions<-grep("DETOX",red_green_data$genes$ProbeName,value=FALSE)
source(paste(workingdirectory,"/Functions/correct_gene_names3.R",sep=""))
new_new_red_green<- correct_gene_names3(red_green_data, gene_names_file,positions)
new_red_green<- new_new_red_green[which(new_new_red_green$genes$ControlType == 0),]


convert_to_MA<- normalizeWithinArrays(new_red_green, method = "loess")
within_norm_RG<- RG.MA(convert_to_MA)
within_norm_RG<-backgroundCorrect(within_norm_RG, method="normexp",normexp.method='mle',offset=50, verbose=TRUE)
between_norm<- normalizeBetweenArrays(within_norm_RG, method = "Aquantile")

MA = between_norm

data_m = MA$M
rownames(data_m) = new_red_green$genes$SystematicName
plotMA3by2(MA)
plotDensities(MA)

targets2 <- targetsA2C(targets)
u <- unique(targets2$Target)
f <- factor(targets2$Target, levels=u)
design <- model.matrix(~0+f)
colnames(design) <- u

MA$M[which(is.na(MA$M)==T)] = 0
MA$A[which(is.na(MA$A)==T)] = 0

corfit <- intraspotCorrelation(MA, design)
fit <- lmscFit(MA, design, correlation=corfit$consensus)

cont.matrix <- makeContrasts("seventytwo-Unexp",levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
topTable(fit2, adjust="BH")


tops <- topTable(fit2,adjust='BH', number=Inf)
write.table(tops,file = paste(workingdirectory2,'/72h - Unexposed.txt',sep=""),sep='\t',row.names=F)
