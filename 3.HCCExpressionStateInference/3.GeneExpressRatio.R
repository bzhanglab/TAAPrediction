
CodingGenes_pathin = '/Users/xinpeiyi/Library/Mobile Documents/com~apple~CloudDocs/Documents/AssistantProfessor/Project/TumorAntigen/Code/CancerGenesProteins/8.Github/Data/CodingGeneList/CodingGenes_V34.txt'
ThresholdData <- readRDS("/Users/xinpeiyi/Library/Mobile Documents/com~apple~CloudDocs/Documents/AssistantProfessor/Project/TumorAntigen/Code/CancerGenesProteins/8.Github/Data/MaxTPMForCertainProbability_HCC/MaxTPMForCertainProbability_HCC.rds")
pathin = '/Users/xinpeiyi/Library/Mobile Documents/com~apple~CloudDocs/Documents/AssistantProfessor/Project/TumorAntigen/Code/CancerGenesProteins/8.Github/Data/HCCTumorTissueData'
Data_HCC_ExpressMatrix_Path = '/Users/xinpeiyi/Library/Mobile Documents/com~apple~CloudDocs/Documents/AssistantProfessor/Project/TumorAntigen/Code/CancerGenesProteins/8.Github/Data/Data_HCC_ExpressMatrix'
Data_CancerGenes_HCC_Path = '/Users/xinpeiyi/Library/Mobile Documents/com~apple~CloudDocs/Documents/AssistantProfessor/Project/TumorAntigen/Code/CancerGenesProteins/8.Github/Data/Data_CancerGenes_HCC'

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
  
  if (fileNames[i]=='HCC_proteomics_gene_abundance.tsv'){
    next
  }
  
  name = gsub(".tumor.ProteinCoding.TPM.txt", "", fileNames[i])
  totalCancerType = c(totalCancerType,name)
  
  filePath = paste(pathin,'/',fileNames[i],sep='')
  HCCData = read.delim(filePath, header = TRUE, row.names = 1, sep = "\t")
  
  overlap_HCC_CodingGenes = intersect(rownames(HCCData),rownames(CodingGenesData))
  HCCData_ExpressedGenes = HCCData[overlap_HCC_CodingGenes,]
  HCCData = HCCData_ExpressedGenes
  HCCData2 = HCCData[,-1]
  HCCData2 = HCCData2[,-2]
  
  ThresholdData2 = ThresholdData[ThresholdData$CancerType==name,]
  HCCData2 = HCCData[,ThresholdData2$SampleID]
  
  tmp_V = matrix(1, nrow(HCCData2), ncol(HCCData2))
  rownames(tmp_V) = rownames(HCCData2)
  colnames(tmp_V) = colnames(HCCData2)
  
  tmp_F_1 = matrix(1, nrow(HCCData2), ncol(HCCData2))
  rownames(tmp_F_1) = rownames(HCCData2)
  colnames(tmp_F_1) = colnames(HCCData2)
  
  tmp_F_2 = matrix(1, nrow(HCCData2), ncol(HCCData2))
  rownames(tmp_F_2) = rownames(HCCData2)
  colnames(tmp_F_2) = colnames(HCCData2)
  for (j in 1:nrow(ThresholdData2)){
    
    # Variable threshold
    threshold = ThresholdData2$MaxTPM[j]
    tmp_V[HCCData2[,ThresholdData2$SampleID[j]]<=threshold,ThresholdData2$SampleID[j]] = 0
    
    # Fixed threshold = 1
    threshold = 1
    tmp_F_1[HCCData2[,ThresholdData2$SampleID[j]]<=threshold,ThresholdData2$SampleID[j]] = 0
    
    # Fixed threshold = 2
    threshold = 2
    tmp_F_2[HCCData2[,ThresholdData2$SampleID[j]]<=threshold,ThresholdData2$SampleID[j]] = 0
  }
  
  tmp_V_2 = data.frame(apply(tmp_V, 1, sum))
  rownames(tmp_V_2) = rownames(tmp_V)
  colnames(tmp_V_2)[1] = 'Count'
  tmp_V_2$Ratio = tmp_V_2[,1]/ncol(HCCData2)
  
  tmp_F_1_2 = data.frame(apply(tmp_F_1, 1, sum))
  rownames(tmp_F_1_2) = rownames(tmp_F_1)
  colnames(tmp_F_1_2)[1] = 'Count'
  tmp_F_1_2$Ratio = tmp_F_1_2[,1]/ncol(HCCData2)
  
  tmp_F_2_2 = data.frame(apply(tmp_F_2, 1, sum))
  rownames(tmp_F_2_2) = rownames(tmp_F_2)
  colnames(tmp_F_2_2)[1] = 'Count'
  tmp_F_2_2$Ratio = tmp_F_2_2[,1]/ncol(HCCData2)
  
  if (nrow(GenesCandidates_V)==0){
    GenesCandidates_V = data.frame(GeneName=rownames(HCCData2))
    GenesCandidates_V$CancerType = tmp_V_2$Ratio
    colnames(GenesCandidates_V) = totalCancerType
  }else{
    GenesCandidates_V$CancerType = tmp_V_2$Ratio
    colnames(GenesCandidates_V) = totalCancerType
  }
  
  if (nrow(GenesCandidates_F_1)==0){
    GenesCandidates_F_1 = data.frame(GeneName=rownames(HCCData2))
    GenesCandidates_F_1$CancerType = tmp_F_1_2$Ratio
    colnames(GenesCandidates_F_1) = totalCancerType
  }else{
    GenesCandidates_F_1$CancerType = tmp_F_1_2$Ratio
    colnames(GenesCandidates_F_1) = totalCancerType
  }
  
  if (nrow(GenesCandidates_F_2)==0){
    GenesCandidates_F_2 = data.frame(GeneName=rownames(HCCData2))
    GenesCandidates_F_2$CancerType = tmp_F_2_2$Ratio
    colnames(GenesCandidates_F_2) = totalCancerType
  }else{
    GenesCandidates_F_2$CancerType = tmp_F_2_2$Ratio
    colnames(GenesCandidates_F_2) = totalCancerType
  }
  
  if (nrow(GenesCandidates_C_V)==0){
    GenesCandidates_C_V = data.frame(GeneName=rownames(HCCData2))
    GenesCandidates_C_V$CancerType = tmp_V_2$Count
    colnames(GenesCandidates_C_V) = totalCancerType
  }else{
    GenesCandidates_C_V$CancerType = tmp_V_2$Count
    colnames(GenesCandidates_C_V) = totalCancerType
  }
  
  if (nrow(GenesCandidates_C_F_1)==0){
    GenesCandidates_C_F_1 = data.frame(GeneName=rownames(HCCData2))
    GenesCandidates_C_F_1$CancerType = tmp_F_1_2$Count
    colnames(GenesCandidates_C_F_1) = totalCancerType
  }else{
    GenesCandidates_C_F_1$CancerType = tmp_F_1_2$Count
    colnames(GenesCandidates_C_F_1) = totalCancerType
  }
  
  if (nrow(GenesCandidates_C_F_2)==0){
    GenesCandidates_C_F_2 = data.frame(GeneName=rownames(HCCData2))
    GenesCandidates_C_F_2$CancerType = tmp_F_2_2$Count
    colnames(GenesCandidates_C_F_2) = totalCancerType
  }else{
    GenesCandidates_C_F_2$CancerType = tmp_F_2_2$Count
    colnames(GenesCandidates_C_F_2) = totalCancerType
  }
  
  tmp_V_matrix <- data.frame(tmp_V)
  tmp_V_matrix$GeneName <- rownames(tmp_V_matrix)
  tmp_V_matrix$Description <- HCCData[rownames(tmp_V_matrix),]$GeneSymbol2
  tmp_V_matrix <- tmp_V_matrix[,c(ncol(tmp_V_matrix)-1,ncol(tmp_V_matrix),1:(ncol(tmp_V_matrix)-2))]
  
  write.table(tmp_V_matrix,
              paste0(Data_HCC_ExpressMatrix_Path,'/',name,'.txt'),
              sep = '\t',
              row.names = FALSE,
              quote = FALSE)
  
  
  print(i)
  
}

pathout = Data_CancerGenes_HCC_Path

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


