
CodingGenes_pathin = 'D:/Project/TumorAntigen/TestData/CodingGeneList/CodingGeneList/CodingGenes.txt'
ThresholdData <- readRDS('D:/Project/TumorAntigen/TestData/Results/MaxTPMForCertainProbability_GTEx/MaxTPMForCertainProbability_GTEx.rds')
GTExTissueDataPath <- 'D:/Project/TumorAntigen/TestData/GTExTissueData'

Pathout <- 'D:/Project/TumorAntigen/TestData/Results'

Data_GTEx_ExpressMatrixPath = paste0(Pathout,'/','Data_GTEx_ExpressMatrix')
if (!file.exists(Data_GTEx_ExpressMatrixPath)){
  dir.create(Data_GTEx_ExpressMatrixPath)
}

DormantGenePath = paste0(Pathout,'/','Data_DormantGenes')
if (!file.exists(DormantGenePath)){
  dir.create(DormantGenePath)
}

CodingGenesData = read.delim(CodingGenes_pathin, header = TRUE, row.names = 1, sep = "\t")
fileNames = list.files(GTExTissueDataPath)

GenesCandidates_V = data.frame()
GenesCandidates_F_1 = data.frame()
GenesCandidates_F_2 = data.frame()

totalTissue = c('GeneName')
for (i in 1:length(fileNames)){
#for (i in 1:2){
  
  #if (fileNames[i]=='Testis.gene_tpm.gct'){
  #  next
  #}
  name = gsub(".gene_tpm.gct", "", fileNames[i])
  totalTissue = c(totalTissue,name)
  
  filePath = paste(GTExTissueDataPath,'/',fileNames[i],sep='')
  GtexData = read.delim(filePath, header = TRUE, row.names = 1, sep = "\t")
  
  overlap_Gtex_CodingGenes = intersect(rownames(GtexData),rownames(CodingGenesData))
  GtexData_ExpressedGenes = GtexData[overlap_Gtex_CodingGenes,]
  GtexData = GtexData_ExpressedGenes
  GtexData2 = GtexData[,-1]
  
  ThresholdData2 = ThresholdData[ThresholdData$TissueType==name,]
  GtexData2 = GtexData[,ThresholdData2$SampleID]
  
  tmp_V = matrix(1, nrow(GtexData2), ncol(GtexData2))
  rownames(tmp_V) = rownames(GtexData2)
  colnames(tmp_V) = colnames(GtexData2)
  
  tmp_F_1 = matrix(1, nrow(GtexData2), ncol(GtexData2))
  rownames(tmp_F_1) = rownames(GtexData2)
  colnames(tmp_F_1) = colnames(GtexData2)
  
  tmp_F_2 = matrix(1, nrow(GtexData2), ncol(GtexData2))
  rownames(tmp_F_2) = rownames(GtexData2)
  colnames(tmp_F_2) = colnames(GtexData2)
  for (j in 1:nrow(ThresholdData2)){
    
    # Variable threshold
    threshold = ThresholdData2$MaxTPM[j]
    tmp_V[GtexData2[,ThresholdData2$SampleID[j]]<=threshold,ThresholdData2$SampleID[j]] = 0
    
    # Fixed threshold = 1
    threshold = 1
    tmp_F_1[GtexData2[,ThresholdData2$SampleID[j]]<=threshold,ThresholdData2$SampleID[j]] = 0
    
    # Fixed threshold = 2
    threshold = 2
    tmp_F_2[GtexData2[,ThresholdData2$SampleID[j]]<=threshold,ThresholdData2$SampleID[j]] = 0
  }
  
  tmp_V_2 = data.frame(apply(tmp_V, 1, sum))
  rownames(tmp_V_2) = rownames(tmp_V)
  tmp_V_2$Ratio = tmp_V_2[,1]/ncol(GtexData2)
  
  tmp_F_1_2 = data.frame(apply(tmp_F_1, 1, sum))
  rownames(tmp_F_1_2) = rownames(tmp_F_1)
  tmp_F_1_2$Ratio = tmp_F_1_2[,1]/ncol(GtexData2)
  
  tmp_F_2_2 = data.frame(apply(tmp_F_2, 1, sum))
  rownames(tmp_F_2_2) = rownames(tmp_F_2)
  tmp_F_2_2$Ratio = tmp_F_2_2[,1]/ncol(GtexData2)
  
  if (nrow(GenesCandidates_V)==0){
    GenesCandidates_V = data.frame(GeneName=rownames(GtexData2))
    GenesCandidates_V$Tissue = tmp_V_2$Ratio
    colnames(GenesCandidates_V) = totalTissue
  }else{
    GenesCandidates_V$Tissue = tmp_V_2$Ratio
    colnames(GenesCandidates_V) = totalTissue
  }
  
  if (nrow(GenesCandidates_F_1)==0){
    GenesCandidates_F_1 = data.frame(GeneName=rownames(GtexData2))
    GenesCandidates_F_1$Tissue = tmp_F_1_2$Ratio
    colnames(GenesCandidates_F_1) = totalTissue
  }else{
    GenesCandidates_F_1$Tissue = tmp_F_1_2$Ratio
    colnames(GenesCandidates_F_1) = totalTissue
  }
  
  if (nrow(GenesCandidates_F_2)==0){
    GenesCandidates_F_2 = data.frame(GeneName=rownames(GtexData2))
    GenesCandidates_F_2$Tissue = tmp_F_2_2$Ratio
    colnames(GenesCandidates_F_2) = totalTissue
  }else{
    GenesCandidates_F_2$Tissue = tmp_F_2_2$Ratio
    colnames(GenesCandidates_F_2) = totalTissue
  }
  
  tmp_V_matrix <- data.frame(tmp_V)
  tmp_V_matrix$GeneName <- rownames(tmp_V_matrix)
  tmp_V_matrix$Description <- GtexData[rownames(tmp_V_matrix),]$Description
  tmp_V_matrix <- tmp_V_matrix[,c(ncol(tmp_V_matrix)-1,ncol(tmp_V_matrix),1:(ncol(tmp_V_matrix)-2))]
  
  write.table(tmp_V_matrix,
              paste0(Data_GTEx_ExpressMatrixPath,'/',name,'.txt'),
              sep = '\t',
              row.names = FALSE,
              quote = FALSE)
  
  print(i)
  
}


rownames(GenesCandidates_V) = GenesCandidates_V$GeneName
GenesCandidates_V1 = c()
for (i in 1:nrow(GenesCandidates_V)){
  gene = GenesCandidates_V[i,1]
  if (!grepl("_PAR_Y", gene, fixed=TRUE)){
    GenesCandidates_V1 = c(GenesCandidates_V1,gene)
  }
}
GenesCandidates_V2 = GenesCandidates_V[GenesCandidates_V1,]
write.table(GenesCandidates_V2, file = paste(DormantGenePath,'/DormantGenesRatio_Pro_0_5.txt',sep=''), sep = "\t",row.names = FALSE,quote = FALSE)

rownames(GenesCandidates_F_1) = GenesCandidates_F_1$GeneName
GenesCandidates_F_11 = c()
for (i in 1:nrow(GenesCandidates_F_1)){
  gene = GenesCandidates_F_1[i,1]
  if (!grepl("_PAR_Y", gene, fixed=TRUE)){
    GenesCandidates_F_11 = c(GenesCandidates_F_11,gene)
  }
}
GenesCandidates_F_12 = GenesCandidates_F_1[GenesCandidates_F_11,]
write.table(GenesCandidates_F_12, file = paste(DormantGenePath,'/DormantGenesRatio_TPM_1.txt',sep=''), sep = "\t",row.names = FALSE,quote = FALSE)

rownames(GenesCandidates_F_2) = GenesCandidates_F_2$GeneName
GenesCandidates_F_21 = c()
for (i in 1:nrow(GenesCandidates_F_2)){
  gene = GenesCandidates_F_2[i,1]
  if (!grepl("_PAR_Y", gene, fixed=TRUE)){
    GenesCandidates_F_21 = c(GenesCandidates_F_21,gene)
  }
}
GenesCandidates_F_22 = GenesCandidates_F_2[GenesCandidates_F_21,]
write.table(GenesCandidates_F_22, file = paste(DormantGenePath,'/DormantGenesRatio_TPM_2.txt',sep=''), sep = "\t",row.names = FALSE,quote = FALSE)


