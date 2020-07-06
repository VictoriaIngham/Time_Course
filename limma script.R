workingdirectory <-'D:\\Dropbox Data\\Dropbox (LSTM)\\Vicky Documents\\PhD Year 1\\GaGa'
workingdirectory2 <-'D:\\Dropbox Data\\Dropbox (LSTM)\\Vicky Documents\\Post Doc\\Exposure Time Course, Transcriptome\\Aging Experiment\\Aging and PANG'
targetfile <- 'targets.txt'
filename<- '12h vs 8h'
reference = '8h'
gene.name.file = '\\AGAP_correct2.txt'

suppressMessages(library(limma))
suppressMessages(library(marray))
suppressMessages(library(stats))
suppressMessages(library(affycoretools))

##SET WORKING DIRECTORY AS AGILENT FILE LOCATION##
setwd(workingdirectory2)


##READ IN ALL DATA AND CONVERT PROBES TO AGAP##
targets<- readTargets(targetfile)
gene_names_file<- read.table(paste(workingdirectory2,gene.name.file,sep=""), header = TRUE, fill = TRUE, sep = "\t")
gene_names_file<-as.matrix(gene_names_file)
red_green_data<- read.maimages(targets, source = "agilent.median", annotation = c("Row", "Col", "FeatureNum", "ControlType", "ProbeName", "SystematicName")) 
positions<-grep("DETOX",red_green_data$genes$ProbeName,value=FALSE)
source(paste(workingdirectory,"\\Functions\\correct_gene_names3.R",sep=""))
new_new_red_green<- correct_gene_names3(red_green_data, gene_names_file,positions)
new_red_green<- new_new_red_green[which(new_new_red_green$genes$ControlType == 0),]


##NORMALIZE ALL DATA##
source(paste(workingdirectory,"\\Functions\\MyNormalizeWithinArrays.R",sep=""))
source(paste(workingdirectory,"\\Functions\\myMA.RG.R",sep=""))
source(paste(workingdirectory,"\\Functions\\myMA2.RG.R",sep=""))
source(paste(workingdirectory,"\\Functions\\MyNormalizeBetweenArrays.R",sep=""))
convert_to_MA<- mynormalizeWithinArrays(new_red_green, method = "loess")
within_norm_RG<- RG.MA(convert_to_MA)
within_norm_RG<-backgroundCorrect(within_norm_RG, method="normexp",normexp.method='mle',offset=50, verbose=TRUE)
between_norm<- mynormalizeBetweenArrays(within_norm_RG, method = "Aquantile")

boxplot(between_norm$M~col(between_norm$M),names=colnames(between_norm$M))
design.vector <- modelMatrix(targets, ref=reference)

out_RG<-RG.MA(between_norm)

##OUTPUT ALL NORMALISED FILES##

dir.create(paste(workingdirectory2,"\\",filename,sep=""))

for(i in 1:nrow(targets))
{
  output<-between_norm[,i]
  output2<-out_RG[,i]
  nameoffile<-output$targets$FileName
  nameoffile2<- paste(workingdirectory2,"\\",filename,"\\",nameoffile,sep="")
  
  if(file.exists(nameoffile2)==FALSE)
  {
    output2<-cbind(output$genes$ProbeName,output$genes$SystematicName,output$M[,1],output$A[,1],output2$R[,1],output2$G[,1])
    colnames(output2)<-c("Probe","ID_REF","M","A","R","G")
    write.table(output2, file = paste(workingdirectory2,"/",filename,"/",nameoffile,sep=""),sep='\t',row.names=F)
  }
}


###Normal design
fit<-lmFit(between_norm,design = design.vector)
limma::plotMA(fit)
abline(0,0,col="blue")
fit<-eBayes(fit)
qqt(fit$t,df=fit$df.prior+fit$df.residual,pch=16,cex=0.2)
abline(0,1)
topTable(fit,adjust.method ='BH')

##Write out top table and produce graphs##
tops <- topTable(fit,adjust='BH', number=Inf)
write.table(tops,file = paste(workingdirectory2,'\\',filename,'\\genesofinterest.txt',sep=""),sep='\t',row.names=F)

limma::plotMA(fit)
top30<-order(fit$lods,decreasing=TRUE)[1:30]
text(fit$Amean[top30],fit$coef[top30],labels=fit$genes[top30,"SystematicName"],cex=0.6,col='blue')



