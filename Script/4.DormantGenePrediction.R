# Dormant gene identification
pathin = 'D:/Project/TumorAntigen/TestData/Results/Data_DormantGenes/DormantGenesRatio_Pro_0_5.txt'
DormantGenesRatio = read.delim(pathin, header = TRUE, row.names = 1, sep = "\t")
DormantGenesRatio$Testis <- NULL

Pathout = 'D:/Project/TumorAntigen/TestData/Results/Data_DormantGenes'

names = c()
numrlt = 0
for (i in 1:nrow(DormantGenesRatio)){
  if (!grepl('_PAR_Y', rownames(DormantGenesRatio)[i], fixed = TRUE)){
    numrlt = numrlt + 1
    names[numrlt] = rownames(DormantGenesRatio)[i]
  }
}
DormantGenesRatio = DormantGenesRatio[names,]

names = c()
for (i in 1:nrow(DormantGenesRatio)){
  names[i] = unlist(strsplit(rownames(DormantGenesRatio)[i], ".", fixed = TRUE))[1]
}
rownames(DormantGenesRatio) = names

DormantGenesRatio2 = 1-DormantGenesRatio

DormantGenesRatio3 = data.frame(apply(DormantGenesRatio2, 1, prod))
rownames(DormantGenesRatio3) = rownames(DormantGenesRatio2)

RatioThreshold = 0.9
DormantGenes = rownames(DormantGenesRatio3)[DormantGenesRatio3[,1]>=RatioThreshold]

write.table(DormantGenes, file = paste0(Pathout,'/','DormantGenes_Ratio_Multiply_0_9','.txt',sep=''), sep = "\t",row.names = FALSE,quote = FALSE)

