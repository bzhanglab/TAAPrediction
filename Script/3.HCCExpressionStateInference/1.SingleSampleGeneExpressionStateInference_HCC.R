library("mixtools")
library(plyr)
library(fitdistrplus)
library(reshape)

# HCC
CancerTyptes = c('HCC')
HCC_pathin = 'D:/Project/TumorAntigen/TestData/HCCTumorTissueData/HCCTumorTissueData'
CodingGenes_pathin = 'D:/Project/TumorAntigen/TestData/CodingGeneList/CodingGeneList/CodingGenes_V34.txt'
HouseKeepingGenes_pathin = 'D:/Project/TumorAntigen/TestData/HouseKeepingGeneList/HouseKeepingGeneList/HK_genes.txt'
Pathout = 'D:/Project/TumorAntigen/TestData/Results'

folder_path = paste0(Pathout,'/','HCC_ExpressPosterior')
if (!file.exists(folder_path)){
  dir.create(folder_path)
}


for (c in 1:length(CancerTyptes)){
  CancerTypte = CancerTyptes[c]
  
  hcc_pathin = paste(HCC_pathin,'/','HCC.tumor.ProteinCoding.TPM.txt',sep = '',collapse = '')
  
  HCCData = read.delim(hcc_pathin, header = TRUE, row.names = 1, sep = "\t")
  
  CodingGenesData = read.delim(CodingGenes_pathin, header = TRUE, row.names = 1, sep = "\t")
  
  overlap_HCC_CodingGenes = intersect(rownames(HCCData),rownames(CodingGenesData))
  HCCData_ExpressedGenes = HCCData[overlap_HCC_CodingGenes,]
  HCCData = HCCData_ExpressedGenes
  
  HCCData = HCCData[,-2]
  
  HouseKeepingGenesData = read.delim(HouseKeepingGenes_pathin, header = TRUE, row.names = 1, sep = "\t")
  rownames(HouseKeepingGenesData) = gsub(' ','',rownames(HouseKeepingGenesData))
  overlap_HCC_CodingGenes_HouseKeepingGenes = intersect(HCCData$GeneSymbol2,rownames(HouseKeepingGenesData))
  HCCData_ExpressedGenes_HouseKeepingGenes = HCCData[!is.na(match(HCCData$GeneSymbol2, overlap_HCC_CodingGenes_HouseKeepingGenes)),]
  
  # All the Data
  AllTPM = data.frame(GeneName=rownames(HCCData))
  Allpi0f0 = data.frame(GeneName=rownames(HCCData))
  Allpi1f1 = data.frame(GeneName=rownames(HCCData))
  rownames(AllTPM) = rownames(HCCData)
  rownames(Allpi0f0) = rownames(HCCData)
  rownames(Allpi1f1) = rownames(HCCData)
  
  for (i in 2:ncol(HCCData)){
    #for (i in 2:100){
    data_test = HCCData[,i]
    data_test2 = data_test[data_test>0]
    data_test3 = log(data_test2)
    
    data_test4 = HCCData_ExpressedGenes_HouseKeepingGenes[,i]
    data_test5 = data_test4[data_test4>0]
    data_test6 = log(data_test5)
    
    data_over_0 = data.frame(GeneName=rownames(HCCData)[data_test>0],TPM=data_test3)
    data_0 = data.frame(GeneName=rownames(HCCData)[data_test==0],TPM=0)
    
    # Initialization by using HouseKeeping genes
    fit = fitdist(data_test6,"norm")
    para = fit$estimate
    
    # Guassian Mixture Distribution
    mixmdl <- normalmixEM(data_test3, k = 2,mu = c(-3,para[1]),sigma = c(1,para[2]),lambda = c(0.2,0.8),maxit = 100000)
    #plot(mixmdl, density=TRUE)
    
    # Parameters
    mean0 = mixmdl$mu[1]
    sd0 = mixmdl$sigma[1]
    pi0 = mixmdl$lambda[1]
    
    mean1 = mixmdl$mu[2]
    sd1 = mixmdl$sigma[2]
    pi1 = mixmdl$lambda[2]
    
    # Posterior factors
    sampleName = colnames(HCCData)[i]
    # if x>0
    data_over_0$pi0f0 = pi0*(dnorm(data_test3, mean = mean0, sd = sd0))
    data_over_0$pi1f1 = pi1*(dnorm(data_test3, mean = mean1, sd = sd1))
    
    rownames(data_over_0) = data_over_0$GeneName
    data_over_0 = data_over_0[,-1]
    
    # if x=0
    data_0$pi0f0 = pi0*1
    data_0$pi1f1 = 0
    
    rownames(data_0) = data_0$GeneName
    data_0 = data_0[,-1]
    
    # Combine data
    CombineData = rbind(data_over_0,data_0)
    
    CombineData_TPM = data.frame(sampleName=CombineData$TPM)
    CombineData_pi0f0 = data.frame(sampleName=CombineData$pi0f0)
    CombineData_pi1f1 = data.frame(sampleName=CombineData$pi1f1)
    
    names(CombineData_TPM) = sampleName
    names(CombineData_pi0f0) = sampleName
    names(CombineData_pi1f1) = sampleName
    
    rownames(CombineData_TPM) = rownames(CombineData)
    rownames(CombineData_pi0f0) = rownames(CombineData)
    rownames(CombineData_pi1f1) = rownames(CombineData)
    
    AllTPM$SampleName = CombineData_TPM[rownames(AllTPM),]
    Allpi0f0$SampleName = CombineData_pi0f0[rownames(Allpi0f0),]
    Allpi1f1$SampleName = CombineData_pi1f1[rownames(Allpi1f1),]
    
    AllTPM = rename(AllTPM,c(SampleName=sampleName))
    Allpi0f0 = rename(Allpi0f0,c(SampleName=sampleName))
    Allpi1f1 = rename(Allpi1f1,c(SampleName=sampleName))
    
    print(i)
    
  }
  
  # Posterior
  AllTPM = AllTPM[,-1]
  Allpi0f0 = Allpi0f0[,-1]
  Allpi1f1 = Allpi1f1[,-1]
  
  ExpressPosterior = 1 - Allpi0f0/(Allpi0f0+Allpi1f1)
  nonExpressPosterior = Allpi0f0/(Allpi0f0+Allpi1f1)
  
  write.table(ExpressPosterior, file = paste0(folder_path,'/',CancerTypte,'.ExpressPosterior','.txt',sep=''), sep = "\t",row.names = TRUE,quote = FALSE)
}




