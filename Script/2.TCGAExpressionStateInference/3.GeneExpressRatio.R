
CodingGenes_pathin = '/Users/xinpeiyi/Library/Mobile Documents/com~apple~CloudDocs/Documents/AssistantProfessor/Project/TumorAntigen/Code/CancerGenesProteins/8.Github/Data/CodingGeneList/CodingGenes.txt'
ThresholdData <- readRDS("/Users/xinpeiyi/Library/Mobile Documents/com~apple~CloudDocs/Documents/AssistantProfessor/Project/TumorAntigen/Code/CancerGenesProteins/8.Github/Data/MaxTPMForCertainProbability_TCGA/MaxTPMForCertainProbability_TCGA.rds")
pathin = '/Users/xinpeiyi/Library/Mobile Documents/com~apple~CloudDocs/Documents/AssistantProfessor/Project/TumorAntigen/Code/CancerGenesProteins/8.Github/Data/TCGATumorTissueData'
Data_TCGA_ExpressMatrix_Path = '/Users/xinpeiyi/Library/Mobile Documents/com~apple~CloudDocs/Documents/AssistantProfessor/Project/TumorAntigen/Code/CancerGenesProteins/8.Github/TAAPrediction/2.TCGAExpressionStateInference/test'
Data_CancerGenes_TCGA_Path = '/Users/xinpeiyi/Library/Mobile Documents/com~apple~CloudDocs/Documents/AssistantProfessor/Project/TumorAntigen/Code/CancerGenesProteins/8.Github/TAAPrediction/2.TCGAExpressionStateInference/test'

CodingGenesData = read.delim(CodingGenes_pathin, header = TRUE, row.names = 1, sep = "\t")

fileNames = list.files(pathin)

GenesCandidates_V = data.frame()
GenesCandidates_F_1 = data.frame()
GenesCandidates_F_2 = data.frame()

GenesCandidates_C_V = data.frame()
GenesCandidates_C_F_1 = data.frame()
GenesCandidates_C_F_2 = data.frame()

totalCancerType = c('GeneName')
for (i in 1:length(fileNames)){
#for (i in 1:10){
  
  if (fileNames[i]=='TGCT.Primary Tumor.bam.gene_tpm.gct'){
    next
  }
  
  name = gsub(".bam.gene_tpm.gct", "", fileNames[i])
  totalCancerType = c(totalCancerType,name)
  
  filePath = paste(pathin,'/',fileNames[i],sep='')
  TCGAData = read.delim(filePath, header = TRUE, row.names = 1, sep = "\t")
  
  overlap_TCGA_CodingGenes = intersect(rownames(TCGAData),rownames(CodingGenesData))
  TCGAData_ExpressedGenes = TCGAData[overlap_TCGA_CodingGenes,]
  TCGAData = TCGAData_ExpressedGenes
  TCGAData2 = TCGAData[,-1]
  
  ThresholdData2 = ThresholdData[ThresholdData$CancerType==name,]
  TCGAData2 = TCGAData[,ThresholdData2$SampleID]
  
  tmp_V = matrix(1, nrow(TCGAData2), ncol(TCGAData2))
  rownames(tmp_V) = rownames(TCGAData2)
  colnames(tmp_V) = colnames(TCGAData2)
  
  tmp_F_1 = matrix(1, nrow(TCGAData2), ncol(TCGAData2))
  rownames(tmp_F_1) = rownames(TCGAData2)
  colnames(tmp_F_1) = colnames(TCGAData2)
  
  tmp_F_2 = matrix(1, nrow(TCGAData2), ncol(TCGAData2))
  rownames(tmp_F_2) = rownames(TCGAData2)
  colnames(tmp_F_2) = colnames(TCGAData2)
  for (j in 1:nrow(ThresholdData2)){
    
    # Variable threshold
    threshold = ThresholdData2$MaxTPM[j]
    tmp_V[TCGAData2[,ThresholdData2$SampleID[j]]<=threshold,ThresholdData2$SampleID[j]] = 0
    
    # Fixed threshold = 1
    threshold = 1
    tmp_F_1[TCGAData2[,ThresholdData2$SampleID[j]]<=threshold,ThresholdData2$SampleID[j]] = 0
    
    # Fixed threshold = 2
    threshold = 2
    tmp_F_2[TCGAData2[,ThresholdData2$SampleID[j]]<=threshold,ThresholdData2$SampleID[j]] = 0
  }
  
  tmp_V_2 = data.frame(apply(tmp_V, 1, sum))
  rownames(tmp_V_2) = rownames(tmp_V)
  colnames(tmp_V_2)[1] = 'Count'
  tmp_V_2$Ratio = tmp_V_2[,1]/ncol(TCGAData2)
  
  tmp_F_1_2 = data.frame(apply(tmp_F_1, 1, sum))
  rownames(tmp_F_1_2) = rownames(tmp_F_1)
  colnames(tmp_F_1_2)[1] = 'Count'
  tmp_F_1_2$Ratio = tmp_F_1_2[,1]/ncol(TCGAData2)
  
  tmp_F_2_2 = data.frame(apply(tmp_F_2, 1, sum))
  rownames(tmp_F_2_2) = rownames(tmp_F_2)
  colnames(tmp_F_2_2)[1] = 'Count'
  tmp_F_2_2$Ratio = tmp_F_2_2[,1]/ncol(TCGAData2)
  
  if (nrow(GenesCandidates_V)==0){
    GenesCandidates_V = data.frame(GeneName=rownames(TCGAData2))
    GenesCandidates_V$CancerType = tmp_V_2$Ratio
    colnames(GenesCandidates_V) = totalCancerType
  }else{
    GenesCandidates_V$CancerType = tmp_V_2$Ratio
    colnames(GenesCandidates_V) = totalCancerType
  }
  
  if (nrow(GenesCandidates_F_1)==0){
    GenesCandidates_F_1 = data.frame(GeneName=rownames(TCGAData2))
    GenesCandidates_F_1$CancerType = tmp_F_1_2$Ratio
    colnames(GenesCandidates_F_1) = totalCancerType
  }else{
    GenesCandidates_F_1$CancerType = tmp_F_1_2$Ratio
    colnames(GenesCandidates_F_1) = totalCancerType
  }
  
  if (nrow(GenesCandidates_F_2)==0){
    GenesCandidates_F_2 = data.frame(GeneName=rownames(TCGAData2))
    GenesCandidates_F_2$CancerType = tmp_F_2_2$Ratio
    colnames(GenesCandidates_F_2) = totalCancerType
  }else{
    GenesCandidates_F_2$CancerType = tmp_F_2_2$Ratio
    colnames(GenesCandidates_F_2) = totalCancerType
  }
  
  if (nrow(GenesCandidates_C_V)==0){
    GenesCandidates_C_V = data.frame(GeneName=rownames(TCGAData2))
    GenesCandidates_C_V$CancerType = tmp_V_2$Count
    colnames(GenesCandidates_C_V) = totalCancerType
  }else{
    GenesCandidates_C_V$CancerType = tmp_V_2$Count
    colnames(GenesCandidates_C_V) = totalCancerType
  }
  
  if (nrow(GenesCandidates_C_F_1)==0){
    GenesCandidates_C_F_1 = data.frame(GeneName=rownames(TCGAData2))
    GenesCandidates_C_F_1$CancerType = tmp_F_1_2$Count
    colnames(GenesCandidates_C_F_1) = totalCancerType
  }else{
    GenesCandidates_C_F_1$CancerType = tmp_F_1_2$Count
    colnames(GenesCandidates_C_F_1) = totalCancerType
  }
  
  if (nrow(GenesCandidates_C_F_2)==0){
    GenesCandidates_C_F_2 = data.frame(GeneName=rownames(TCGAData2))
    GenesCandidates_C_F_2$CancerType = tmp_F_2_2$Count
    colnames(GenesCandidates_C_F_2) = totalCancerType
  }else{
    GenesCandidates_C_F_2$CancerType = tmp_F_2_2$Count
    colnames(GenesCandidates_C_F_2) = totalCancerType
  }
  
  tmp_V_matrix <- data.frame(tmp_V)
  tmp_V_matrix$GeneName <- rownames(tmp_V_matrix)
  tmp_V_matrix$Description <- TCGAData[rownames(tmp_V_matrix),]$Description
  tmp_V_matrix <- tmp_V_matrix[,c(ncol(tmp_V_matrix)-1,ncol(tmp_V_matrix),1:(ncol(tmp_V_matrix)-2))]
  
  write.table(tmp_V_matrix,
              paste0(Data_TCGA_ExpressMatrix_Path,'/',name,'.txt'),
              sep = '\t',
              row.names = FALSE,
              quote = FALSE)
  
  
  print(i)
  
}

pathout = Data_CancerGenes_TCGA_Path

rownames(GenesCandidates_V) = GenesCandidates_V$GeneName
GenesCandidates_V1 = c()
for (i in 1:nrow(GenesCandidates_V)){
  gene = GenesCandidates_V[i,1]
  if (!grepl("_PAR_Y", gene, fixed=TRUE)){
    GenesCandidates_V1 = c(GenesCandidates_V1,gene)
  }
}
GenesCandidates_V2 = GenesCandidates_V[GenesCandidates_V1,]
write.table(GenesCandidates_V2, file = paste(pathout,'/CancerGenesRatio_Pro_0_5.txt',sep=''), sep = "\t",row.names = FALSE,quote = FALSE)

rownames(GenesCandidates_F_1) = GenesCandidates_F_1$GeneName
GenesCandidates_F_11 = c()
for (i in 1:nrow(GenesCandidates_F_1)){
  gene = GenesCandidates_F_1[i,1]
  if (!grepl("_PAR_Y", gene, fixed=TRUE)){
    GenesCandidates_F_11 = c(GenesCandidates_F_11,gene)
  }
}
GenesCandidates_F_12 = GenesCandidates_F_1[GenesCandidates_F_11,]
write.table(GenesCandidates_F_12, file = paste(pathout,'/CancerGenesRatio_TPM_1.txt',sep=''), sep = "\t",row.names = FALSE,quote = FALSE)

rownames(GenesCandidates_F_2) = GenesCandidates_F_2$GeneName
GenesCandidates_F_21 = c()
for (i in 1:nrow(GenesCandidates_F_2)){
  gene = GenesCandidates_F_2[i,1]
  if (!grepl("_PAR_Y", gene, fixed=TRUE)){
    GenesCandidates_F_21 = c(GenesCandidates_F_21,gene)
  }
}
GenesCandidates_F_22 = GenesCandidates_F_2[GenesCandidates_F_21,]
write.table(GenesCandidates_F_22, file = paste(pathout,'/CancerGenesRatio_TPM_2.txt',sep=''), sep = "\t",row.names = FALSE,quote = FALSE)

rownames(GenesCandidates_C_V) = GenesCandidates_C_V$GeneName
GenesCandidates_C_V1 = c()
for (i in 1:nrow(GenesCandidates_C_V)){
  gene = GenesCandidates_C_V[i,1]
  if (!grepl("_PAR_Y", gene, fixed=TRUE)){
    GenesCandidates_C_V1 = c(GenesCandidates_C_V1,gene)
  }
}
GenesCandidates_C_V2 = GenesCandidates_C_V[GenesCandidates_C_V1,]
write.table(GenesCandidates_C_V2, file = paste(pathout,'/CancerGenesCount_Pro_0_5.txt',sep=''), sep = "\t",row.names = FALSE,quote = FALSE)

rownames(GenesCandidates_C_F_1) = GenesCandidates_C_F_1$GeneName
GenesCandidates_C_F_11 = c()
for (i in 1:nrow(GenesCandidates_C_F_1)){
  gene = GenesCandidates_C_F_1[i,1]
  if (!grepl("_PAR_Y", gene, fixed=TRUE)){
    GenesCandidates_C_F_11 = c(GenesCandidates_C_F_11,gene)
  }
}
GenesCandidates_C_F_12 = GenesCandidates_C_F_1[GenesCandidates_C_F_11,]
write.table(GenesCandidates_C_F_12, file = paste(pathout,'/CancerGenesCount_TPM_1.txt',sep=''), sep = "\t",row.names = FALSE,quote = FALSE)

rownames(GenesCandidates_C_F_2) = GenesCandidates_C_F_2$GeneName
GenesCandidates_C_F_21 = c()
for (i in 1:nrow(GenesCandidates_C_F_2)){
  gene = GenesCandidates_C_F_2[i,1]
  if (!grepl("_PAR_Y", gene, fixed=TRUE)){
    GenesCandidates_C_F_21 = c(GenesCandidates_C_F_21,gene)
  }
}
GenesCandidates_C_F_22 = GenesCandidates_C_F_2[GenesCandidates_C_F_21,]
write.table(GenesCandidates_C_F_22, file = paste(pathout,'/CancerGenesCount_TPM_2.txt',sep=''), sep = "\t",row.names = FALSE,quote = FALSE)

