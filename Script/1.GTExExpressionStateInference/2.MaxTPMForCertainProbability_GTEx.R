library(ggplot2)
library(mclust)

ExpressionProFilePath = 'D:/Project/TumorAntigen/TestData/Results/GTEx_ExpressPosterior'
TPMFilePath = 'D:/Project/TumorAntigen/TestData/GTExTissueData'

Pathout = 'D:/Project/TumorAntigen/TestData/Results'

MaxTPMForCertainProbability_GTEx_Path = paste0(Pathout,'/','MaxTPMForCertainProbability_GTEx')
if (!file.exists(MaxTPMForCertainProbability_GTEx_Path)){
  dir.create(MaxTPMForCertainProbability_GTEx_Path)
}

ThresholdSelectionPics_GTEx = MaxTPMForCertainProbability_GTEx_Path

MaxTPM = c()
TissueType = c()
SampleID = c()
FileNames = list.files(ExpressionProFilePath)

Warning_SampleID = c()
Warning_TissueID = c()
Warning_TPM = c()
Warning_Clust = c()
numrlt = 0
for (i in 1:length(FileNames)){
#for (i in 1:1){
  FileName = FileNames[i]
  Tissue = gsub(".ExpressPosterior.txt", "", FileName)
  FilePath = paste(ExpressionProFilePath,'/',FileName,sep = '')
  NonExpressPro = 1-read.delim(FilePath, header = TRUE,row.names = 1,sep = "\t")
  FilePath = paste(TPMFilePath,'/',Tissue,'.gene_tpm.gct',sep = '')
  TPM = read.delim(FilePath, header = TRUE, row.names = 1,sep = "\t")
  TPM = TPM[,-1]
  TPM = TPM[rownames(NonExpressPro),colnames(NonExpressPro)]
  
  Count = 0
  for (j in 1:ncol(NonExpressPro)){
    tmp = TPM[,j]
    tmp2 = log(tmp[tmp>0])
    m_clust = Mclust(as.matrix(tmp2),G=1:2)
    
    if (max(tmp[NonExpressPro[,j]>=0.5])>10 || m_clust$G == 1){
      numrlt = numrlt + 1
      Warning_SampleID = c(Warning_SampleID,colnames(NonExpressPro)[j])
      Warning_TissueID = c(Warning_TissueID,Tissue)
      Warning_TPM = c(Warning_TPM,max(tmp[NonExpressPro[,j]>=0.5]))
      Warning_Clust = c(Warning_Clust,m_clust$G)
      next
    }
    
    MaxTPM = c(MaxTPM,max(tmp[NonExpressPro[,j]>=0.5]))
    SampleID = c(SampleID,colnames(NonExpressPro)[j])
    Count = Count + 1
  }
  TissueType = c(TissueType,rep(Tissue,Count))
  
  print(Tissue)
}

Warning = data.frame(SampleID=Warning_SampleID,MaxTPM=Warning_TPM,TissueType=Warning_TissueID,Clust=Warning_Clust)

Data = data.frame(SampleID = SampleID, MaxTPM = MaxTPM, TissueType = TissueType)

saveRDS(Data, paste0(MaxTPMForCertainProbability_GTEx_Path,'/','MaxTPMForCertainProbability_GTEx.rds'))


Data <- readRDS(paste0(MaxTPMForCertainProbability_GTEx_Path,'/','MaxTPMForCertainProbability_GTEx.rds'))

pathout = ThresholdSelectionPics_GTEx
TissueType2 = unique(Data$TissueType)
for (i in 1:length(TissueType2)){
  TissueType3 = TissueType2[i]
  Data2 = Data[Data$TissueType==TissueType3,]
  p = ggplot(Data2, aes(x=rownames(Data2), y=MaxTPM)) + 
    geom_point() + 
    geom_hline(yintercept=median(Data2$MaxTPM), linetype="dashed", color = "blue", size=2) +
    geom_hline(yintercept=2, linetype="dashed", color = "red", size=2) + 
    geom_hline(yintercept=1, linetype="dashed", color = "red", size=2) +
    geom_hline(yintercept=0, linetype="dashed", color = "red", size=2) + 
    xlab('Sample ID') + 
    ylab('Max TPM of not express probability > 0.5')
  
  ggsave(paste(pathout,'/',TissueType3,'.png',sep = ''), width = 10, height = 5)
}

library(ggplot2)
p <- ggplot(Data, aes(x=TissueType, y=MaxTPM)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE)
p + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(paste(pathout,'/','AllTissue_BoxPlot','.png',sep = ''), width = 20, height = 5)



