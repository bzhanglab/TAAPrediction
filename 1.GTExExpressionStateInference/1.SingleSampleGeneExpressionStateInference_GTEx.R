library("mixtools")
library(plyr)
library(fitdistrplus)
library(reshape)

# GTEx 54 normal tissues
Tissues = c('Adipose - Subcutaneous','Adipose - Visceral (Omentum)','Adrenal Gland','Artery - Aorta','Artery - Coronary',
            'Artery - Tibial','Bladder','Brain - Amygdala','Brain - Anterior cingulate cortex (BA24)','Brain - Caudate (basal ganglia)',
            'Brain - Cerebellar Hemisphere','Brain - Cerebellum','Brain - Cortex','Brain - Frontal Cortex (BA9)','Brain - Hippocampus',
            'Brain - Hypothalamus','Brain - Nucleus accumbens (basal ganglia)','Brain - Putamen (basal ganglia)',
            'Brain - Spinal cord (cervical c-1)','Brain - Substantia nigra','Breast - Mammary Tissue','Cells - Cultured fibroblasts',
            'Cells - EBV-transformed lymphocytes','Cervix - Ectocervix','Cervix - Endocervix','Colon - Sigmoid','Colon - Transverse',
            'Esophagus - Gastroesophageal Junction','Esophagus - Mucosa','Esophagus - Muscularis','Fallopian Tube',
            'Heart - Atrial Appendage','Heart - Left Ventricle','Kidney - Cortex','Kidney - Medulla',
            'Liver','Lung','Minor Salivary Gland','Muscle - Skeletal','Nerve - Tibial','Ovary','Pancreas',
            'Pituitary','Prostate','Skin - Not Sun Exposed (Suprapubic)','Skin - Sun Exposed (Lower leg)',
            'Small Intestine - Terminal Ileum','Spleen','Stomach','Testis','Thyroid','Uterus','Vagina','Whole Blood'
            )

GTEX_Data_Path = '/Users/xinpeiyi/Library/Mobile Documents/com~apple~CloudDocs/Documents/AssistantProfessor/Project/TumorAntigen/Code/CancerGenesProteins/8.Github/TAAPrediction/Data/GTExTissueData'
CodingGenes_pathin = '/Users/xinpeiyi/Library/Mobile Documents/com~apple~CloudDocs/Documents/AssistantProfessor/Project/TumorAntigen/Code/CancerGenesProteins/8.Github/TAAPrediction/Data/CodingGeneList/CodingGenes.txt'
HouseKeepingGenes_pathin = '/Users/xinpeiyi/Library/Mobile Documents/com~apple~CloudDocs/Documents/AssistantProfessor/Project/TumorAntigen/Code/CancerGenesProteins/8.Github/TAAPrediction/Data/HouseKeepingGeneList/HK_genes.txt'
Pathout = '/Users/xinpeiyi/Library/Mobile Documents/com~apple~CloudDocs/Documents/AssistantProfessor/Project/TumorAntigen/Code/CancerGenesProteins/8.Github/TAAPrediction/1.GTExExpressionStateInference/test'

#for (t in 1:length(Tissues)){
for (t in 2:2){
  
  Tissue = Tissues[t]
  
  print(Tissue)
  
  gtex_pathin = paste0(GTEX_Data_Path,'/',Tissue,'.gene_tpm.gct')
  GtexData = read.delim(gtex_pathin, header = TRUE, row.names = 1, sep = "\t")
  
  CodingGenesData = read.delim(CodingGenes_pathin, header = TRUE, row.names = 1, sep = "\t")
  
  overlap_Gtex_CodingGenes = intersect(rownames(GtexData),rownames(CodingGenesData))
  GtexData_ExpressedGenes = GtexData[overlap_Gtex_CodingGenes,]
  GtexData = GtexData_ExpressedGenes
  
  HouseKeepingGenesData = read.delim(HouseKeepingGenes_pathin, header = TRUE, row.names = 1, sep = "\t")
  rownames(HouseKeepingGenesData) = gsub(' ','',rownames(HouseKeepingGenesData))
  overlap_Gtex_CodingGenes_HouseKeepingGenes = intersect(GtexData$Description,rownames(HouseKeepingGenesData))
  GtexData_ExpressedGenes_HouseKeepingGenes = GtexData[!is.na(match(GtexData$Description, overlap_Gtex_CodingGenes_HouseKeepingGenes)),]
  
  # All the Data
  AllTPM = data.frame(GeneName=rownames(GtexData))
  Allpi0f0 = data.frame(GeneName=rownames(GtexData))
  Allpi1f1 = data.frame(GeneName=rownames(GtexData))
  rownames(AllTPM) = rownames(GtexData)
  rownames(Allpi0f0) = rownames(GtexData)
  rownames(Allpi1f1) = rownames(GtexData)
  
  for (i in 2:ncol(GtexData)){
    data_test = GtexData[,i]
    data_test2 = data_test[data_test>0]
    data_test3 = log(data_test2)
    
    data_test4 = GtexData_ExpressedGenes_HouseKeepingGenes[,i]
    data_test5 = data_test4[data_test4>0]
    data_test6 = log(data_test5)
    
    data_over_0 = data.frame(GeneName=rownames(GtexData)[data_test>0],TPM=data_test3)
    data_0 = data.frame(GeneName=rownames(GtexData)[data_test==0],TPM=0)
    
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
    sampleName = colnames(GtexData)[i]
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
  
  write.table(ExpressPosterior, file = paste0(Pathout,'/',Tissue,'.ExpressPosterior','.txt',sep=''), sep = "\t",row.names = TRUE,quote = FALSE)
  
}












