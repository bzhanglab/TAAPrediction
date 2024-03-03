library(ggplot2)
library(mclust)

ExpressionProFilePath = 'D:/Project/TumorAntigen/TestData/Results/TCGA_ExpressPosterior'
TPMFilePath = 'D:/Project/TumorAntigen/TestData/TCGATumorTissueData'
Pathout = 'D:/Project/TumorAntigen/TestData/Results'

MaxTPMForCertainProbability_TCGA_Path = paste0(Pathout,'/','MaxTPMForCertainProbability_TCGA')
if (!file.exists(MaxTPMForCertainProbability_TCGA_Path)){
  dir.create(MaxTPMForCertainProbability_TCGA_Path)
}

ThresholdSelectionPics_TCGA = MaxTPMForCertainProbability_TCGA_Path

MaxTPM = c()
CancerType = c()
SampleID = c()
FileNames = list.files(ExpressionProFilePath)

Warning_SampleID = c()
Warning_CancerID = c()
Warning_TPM = c()
Warning_Clust = c()
numrlt = 0
for (i in 1:length(FileNames)){
#for (i in 1:1){
  FileName = FileNames[i]
  
  if (!grepl(".Primary Tumor.", FileName, fixed=TRUE) & !grepl('SKCM_M.Metastatic.ExpressPosterior.txt',FileName,fixed=TRUE) & !grepl('LAML.Primary Blood Derived Cancer - Peripheral Blood.ExpressPosterior.txt',FileName,fixed=TRUE)){
    next
  }
  
  Cancer = gsub(".ExpressPosterior.txt", "", FileName)
  FilePath = paste(ExpressionProFilePath,'/',FileName,sep = '')
  NonExpressPro = 1-read.delim(FilePath, header = TRUE,row.names = 1,sep = "\t")
  FilePath = paste(TPMFilePath,'/',Cancer,'.bam.gene_tpm.gct',sep = '')
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
      Warning_CancerID = c(Warning_CancerID,Cancer)
      Warning_TPM = c(Warning_TPM,max(tmp[NonExpressPro[,j]>=0.5]))
      Warning_Clust = c(Warning_Clust,m_clust$G)
      next
    }
    
    MaxTPM = c(MaxTPM,max(tmp[NonExpressPro[,j]>=0.5]))
    SampleID = c(SampleID,colnames(NonExpressPro)[j])
    Count = Count + 1
  }
  CancerType = c(CancerType,rep(Cancer,Count))
  
  print(i)
}

Warning = data.frame(SampleID=Warning_SampleID,MaxTPM=Warning_TPM,CancerType=Warning_CancerID,Clust=Warning_Clust)

Data = data.frame(SampleID = SampleID, MaxTPM = MaxTPM, CancerType = CancerType)

saveRDS(Data, paste0(MaxTPMForCertainProbability_TCGA_Path,'/MaxTPMForCertainProbability_TCGA.rds'))

Data <- readRDS(paste0(MaxTPMForCertainProbability_TCGA_Path,'/MaxTPMForCertainProbability_TCGA.rds'))

pathout = ThresholdSelectionPics_TCGA
CancerType2 = unique(Data$CancerType)
for (i in 1:length(CancerType2)){
  CancerType3 = CancerType2[i]
  Data2 = Data[Data$CancerType==CancerType3,]
  p = ggplot(Data2, aes(x=rownames(Data2), y=MaxTPM)) + 
    geom_point() + 
    geom_hline(yintercept=median(Data2$MaxTPM), linetype="dashed", color = "blue", size=2) +
    geom_hline(yintercept=2, linetype="dashed", color = "red", size=2) + 
    geom_hline(yintercept=1, linetype="dashed", color = "red", size=2) +
    geom_hline(yintercept=0, linetype="dashed", color = "red", size=2) + 
    xlab('Sample ID') + 
    ylab('Max TPM of not express probability > 0.5')
  
  ggsave(paste(pathout,'/',CancerType3,'.png',sep = ''), width = 10, height = 5)
}

library(ggplot2)
p <- ggplot(Data, aes(x=CancerType, y=MaxTPM)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE)
p + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(paste(pathout,'/','AllCancerTypes_BoxPlot','.png',sep = ''), width = 20, height = 5)



