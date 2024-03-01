##########################################################
CodingGenes_pathin = '/Users/xinpeiyi/Library/Mobile Documents/com~apple~CloudDocs/Documents/AssistantProfessor/Project/TumorAntigen/Code/CancerGenesProteins/2.DormantGeneIdentification/Data/CodingGenes_Gencode/CodingGenes.txt'
CodingGenesData = read.delim(CodingGenes_pathin, header = TRUE, row.names = 1, sep = "\t")
name = c()
for (i in 1:length(rownames(CodingGenesData))){
  name[i] = unlist(strsplit(rownames(CodingGenesData)[i],"[.]"))[[1]]
}
CodingGenesData$gene_name2 = name
##########################################################

pathin = '/Users/xinpeiyi/Library/Mobile Documents/com~apple~CloudDocs/Documents/AssistantProfessor/Project/TumorAntigen/Code/CancerGenesProteins/4.CancerGeneIdentification/Data/Data_CancerGenes_CPTAC_AddHCC/CancerGenesCount_Pro_0_5.txt'
CancerGenesCount = read.delim(pathin, header = TRUE, row.names = 1, sep = "\t")
name = c()
for (i in 1:length(rownames(CancerGenesCount))){
  name[i] = unlist(strsplit(rownames(CancerGenesCount)[i],"[.]"))[[1]]
}
rownames(CancerGenesCount) = name

# Dormant Genes:

pathin = '/Users/xinpeiyi/Library/Mobile Documents/com~apple~CloudDocs/Documents/AssistantProfessor/Project/TumorAntigen/Code/CancerGenesProteins/2.DormantGeneIdentification/Data/Data_DormantGenes/DormantGenesRatio_Pro_0_5.txt'
DormantGenesRatio = read.delim(pathin, header = TRUE, row.names = 1, sep = "\t")
DormantGenesRatio$Testis <- NULL

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

#=======================Store data=======================#
DormantGenes = setdiff(DormantGenes,'ENSG00000213588')
DormantGenes = setdiff(DormantGenes,'ENSG00000124610')
write.table(DormantGenes, file = paste('/Users/xinpeiyi/Library/Mobile Documents/com~apple~CloudDocs/Documents/AssistantProfessor/Project/TumorAntigen/Code/CancerGenesProteins/2.DormantGeneIdentification/Data/Data_DormantGenes/','DormantGenes_Ratio_Multiply_0_9','.txt',sep=''), sep = "\t",row.names = FALSE,quote = FALSE)
write.table(setdiff(rownames(DormantGenesRatio3),DormantGenes), file = paste('/Users/xinpeiyi/Library/Mobile Documents/com~apple~CloudDocs/Documents/AssistantProfessor/Project/TumorAntigen/Code/CancerGenesProteins/2.DormantGeneIdentification/Data/Data_DormantGenes/','RemainingGenes_Ratio_Multiply_0_9','.txt',sep=''), sep = "\t",row.names = FALSE,quote = FALSE)
#========================================================#

pathin = '/Users/xinpeiyi/Library/Mobile Documents/com~apple~CloudDocs/Documents/AssistantProfessor/Project/TumorAntigen/Code/CancerGenesProteins/4.CancerGeneIdentification/Data/Data_CancerGenes_CPTAC_AddHCC/CancerGenesRatio_Pro_0_5.txt'
CancerGenesRatio = read.delim(pathin, header = TRUE, row.names = 1, sep = "\t")

names = c()
numrlt = 0
for (i in 1:nrow(CancerGenesRatio)){
  if (!grepl('_PAR_Y', rownames(CancerGenesRatio)[i], fixed = TRUE)){
    numrlt = numrlt + 1
    names[numrlt] = rownames(CancerGenesRatio)[i]
  }
}
CancerGenesRatio = CancerGenesRatio[names,]

names = c()
for (i in 1:nrow(CancerGenesRatio)){
  names[i] = unlist(strsplit(rownames(CancerGenesRatio)[i], ".", fixed = TRUE))[1]
}
rownames(CancerGenesRatio) = names

RatioThreshold1 = 0.2
Count2 <- function(x){
  y = c()
  for (l in 1:length(x)){
    if (as.double(x[l])>=RatioThreshold1){
      y[l] = 1
    }else{
      y[l] = 0
    }
  }
  return(y)
}

CancerGenesRatio2 = data.frame(apply(CancerGenesRatio, 2, Count2))
rownames(CancerGenesRatio2) = rownames(CancerGenesRatio)
colnames(CancerGenesRatio2) = colnames(CancerGenesRatio)
CancerGenesRatio3 = data.frame(apply(CancerGenesRatio2, 1, sum))
rownames(CancerGenesRatio3) = rownames(CancerGenesRatio2)

CancerGenes = rownames(CancerGenesRatio3)[CancerGenesRatio3[,1]>0]

overlap = intersect(DormantGenes,CancerGenes)
CodingGenesData2 = CodingGenesData[!is.na(match(CodingGenesData$gene_name2, overlap)),]
rownames(CodingGenesData2) = CodingGenesData2$gene_name2

Genes_20 = data.frame(GeneName = overlap, Description = CodingGenesData2[overlap,]$gene_name)
rownames(Genes_20) = Genes_20$GeneName
Genes_20 = cbind(Genes_20,CancerGenesRatio2[rownames(Genes_20),])

Genes_20 = Genes_20[setdiff(rownames(Genes_20),'ENSG00000213588'),]
Genes_20 = Genes_20[setdiff(rownames(Genes_20),'ENSG00000124610'),]
#=======================Store data=======================#
write.table(Genes_20, file = paste('/Users/xinpeiyi/Library/Mobile Documents/com~apple~CloudDocs/Documents/AssistantProfessor/Project/TumorAntigen/Code/CancerGenesProteins/4.CancerGeneIdentification/Data/Data_CancerGenes_CPTAC_AddHCC/CancerGenes/','Ratio_Multiply_0_9_0_2_noHistone','.txt',sep=''), sep = "\t",row.names = TRUE,quote = FALSE)
#======================Draw figure=======================#
CancerGenesCount2 = CancerGenesCount[rownames(Genes_20),]
rownames(CancerGenesCount2) = CodingGenesData2[rownames(CancerGenesCount2),]$gene_name
Antigen <- c(rep(as.character(rownames(CancerGenesCount2)), each = ncol(CancerGenesCount2)))
CancerType  <- c(rep(as.character(names(CancerGenesCount2))[1:ncol(CancerGenesCount2)], times = nrow(CancerGenesCount2)))
Frequency = c()
for (i in 1:nrow(CancerGenesCount2)){
  temp <- as.numeric(CancerGenesCount2[i,1:ncol(CancerGenesCount2)])
  Frequency = append(Frequency,temp)
}
Data <- data.frame(Antigen, CancerType, Frequency)

pdf("/Users/xinpeiyi/Library/Mobile Documents/com~apple~CloudDocs/Documents/AssistantProfessor/Project/TumorAntigen/Code/CancerGenesProteins/4.CancerGeneIdentification/Pics/Ratio_Multiply_0_9_0_2_noHistone.pdf",width=7,height=13,paper='special') 
library(ggplot2)
ggplot(Data, aes(x = reorder(Antigen,Frequency), y = Frequency, fill = CancerType)) +
  geom_bar(stat = "identity") + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text.x = element_text(angle=0, hjust=0.5,vjust = 0.5),axis.title.y=element_blank(),legend.key = element_rect(fill = "white", colour = "black"),legend.title = element_text(face = "bold")) + labs(x = "Cancer genes") + labs(y = "#Sample")+
  coord_flip()
dev.off() 
#========================================================#


RatioThreshold2 = 0.1
Count3 <- function(x){
  y = c()
  for (l in 1:length(x)){
    if (as.double(x[l])>=RatioThreshold2){
      y[l] = 1
    }else{
      y[l] = 0
    }
  }
  return(y)
}
CancerGenesRatio2 = data.frame(apply(CancerGenesRatio, 2, Count3))
rownames(CancerGenesRatio2) = rownames(CancerGenesRatio)
colnames(CancerGenesRatio2) = colnames(CancerGenesRatio)
CancerGenesRatio3 = data.frame(apply(CancerGenesRatio2, 1, sum))
rownames(CancerGenesRatio3) = rownames(CancerGenesRatio2)

CancerGenes = rownames(CancerGenesRatio3)[CancerGenesRatio3[,1]>0]

overlap = intersect(DormantGenes,CancerGenes)
CodingGenesData2 = CodingGenesData[!is.na(match(CodingGenesData$gene_name2, overlap)),]
rownames(CodingGenesData2) = CodingGenesData2$gene_name2

Genes_10 = data.frame(GeneName = overlap, Description = CodingGenesData2[overlap,]$gene_name)
rownames(Genes_10) = Genes_10$GeneName
Genes_10 = cbind(Genes_10,CancerGenesRatio2[rownames(Genes_10),])

Genes_10 = Genes_10[setdiff(rownames(Genes_10),'ENSG00000213588'),]
Genes_10 = Genes_10[setdiff(rownames(Genes_10),'ENSG00000124610'),]
#=======================Store data=======================#
write.table(Genes_10, file = paste('/Users/xinpeiyi/Library/Mobile Documents/com~apple~CloudDocs/Documents/AssistantProfessor/Project/TumorAntigen/Code/CancerGenesProteins/4.CancerGeneIdentification/Data/Data_CancerGenes_CPTAC_AddHCC/CancerGenes/','Ratio_Multiply_0_9_0_1_noHistone','.txt',sep=''), sep = "\t",row.names = TRUE,quote = FALSE)
#======================Draw figure=======================#
CancerGenesCount2 = CancerGenesCount[rownames(Genes_10),]
rownames(CancerGenesCount2) = CodingGenesData2[rownames(CancerGenesCount2),]$gene_name
Antigen <- c(rep(as.character(rownames(CancerGenesCount2)), each = ncol(CancerGenesCount2)))
CancerType  <- c(rep(as.character(names(CancerGenesCount2))[1:ncol(CancerGenesCount2)], times = nrow(CancerGenesCount2)))
Frequency = c()
for (i in 1:nrow(CancerGenesCount2)){
  temp <- as.numeric(CancerGenesCount2[i,1:ncol(CancerGenesCount2)])
  Frequency = append(Frequency,temp)
}
Data <- data.frame(Antigen, CancerType, Frequency)

pdf("/Users/xinpeiyi/Library/Mobile Documents/com~apple~CloudDocs/Documents/AssistantProfessor/Project/TumorAntigen/Code/CancerGenesProteins/4.CancerGeneIdentification/Pics/Ratio_Multiply_0_9_0_1_noHistone.pdf",width=7,height=13,paper='special') 
library(ggplot2)
ggplot(Data, aes(x = reorder(Antigen,Frequency), y = Frequency, fill = CancerType)) +
  geom_bar(stat = "identity") + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text.x = element_text(angle=0, hjust=0.5,vjust = 0.5),axis.title.y=element_blank(),legend.key = element_rect(fill = "white", colour = "black"),legend.title = element_text(face = "bold")) + labs(x = "Cancer genes") + labs(y = "#Sample")+
  coord_flip()
dev.off() 
#========================================================#


RatioThreshold3 = 0.05
Count4 <- function(x){
  y = c()
  for (l in 1:length(x)){
    if (as.double(x[l])>=RatioThreshold3){
      y[l] = 1
    }else{
      y[l] = 0
    }
  }
  return(y)
}
CancerGenesRatio2 = data.frame(apply(CancerGenesRatio, 2, Count4))
rownames(CancerGenesRatio2) = rownames(CancerGenesRatio)
colnames(CancerGenesRatio2) = colnames(CancerGenesRatio)
CancerGenesRatio3 = data.frame(apply(CancerGenesRatio2, 1, sum))
rownames(CancerGenesRatio3) = rownames(CancerGenesRatio2)

CancerGenes = rownames(CancerGenesRatio3)[CancerGenesRatio3[,1]>0]

overlap = intersect(DormantGenes,CancerGenes)
CodingGenesData2 = CodingGenesData[!is.na(match(CodingGenesData$gene_name2, overlap)),]
rownames(CodingGenesData2) = CodingGenesData2$gene_name2

Genes_5 = data.frame(GeneName = overlap, Description = CodingGenesData2[overlap,]$gene_name)
rownames(Genes_5) = Genes_5$GeneName
Genes_5 = cbind(Genes_5,CancerGenesRatio2[rownames(Genes_5),])

Genes_5 = Genes_5[setdiff(rownames(Genes_5),'ENSG00000213588'),]
Genes_5 = Genes_5[setdiff(rownames(Genes_5),'ENSG00000124610'),]
#=======================Store data=======================#
write.table(Genes_5, file = paste('/Users/xinpeiyi/Library/Mobile Documents/com~apple~CloudDocs/Documents/AssistantProfessor/Project/TumorAntigen/Code/CancerGenesProteins/4.CancerGeneIdentification/Data/Data_CancerGenes_CPTAC_AddHCC/CancerGenes/','Ratio_Multiply_0_9_0_0_5_noHistone','.txt',sep=''), sep = "\t",row.names = TRUE,quote = FALSE)
#======================Draw figure=======================#
CancerGenesCount2 = CancerGenesCount[rownames(Genes_5),]
rownames(CancerGenesCount2) = CodingGenesData2[rownames(CancerGenesCount2),]$gene_name
Antigen <- c(rep(as.character(rownames(CancerGenesCount2)), each = ncol(CancerGenesCount2)))
CancerType  <- c(rep(as.character(names(CancerGenesCount2))[1:ncol(CancerGenesCount2)], times = nrow(CancerGenesCount2)))
Frequency = c()
for (i in 1:nrow(CancerGenesCount2)){
  temp <- as.numeric(CancerGenesCount2[i,1:ncol(CancerGenesCount2)])
  Frequency = append(Frequency,temp)
}
Data <- data.frame(Antigen, CancerType, Frequency)

pdf("/Users/xinpeiyi/Library/Mobile Documents/com~apple~CloudDocs/Documents/AssistantProfessor/Project/TumorAntigen/Code/CancerGenesProteins/4.CancerGeneIdentification/Pics/Ratio_Multiply_0_9_0_0_5_noHistone.pdf",width=7,height=13,paper='special') 
library(ggplot2)
ggplot(Data, aes(x = reorder(Antigen,Frequency), y = Frequency, fill = CancerType)) +
  geom_bar(stat = "identity") + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text.x = element_text(angle=0, hjust=0.5,vjust = 0.5),axis.title.y=element_blank(),legend.key = element_rect(fill = "white", colour = "black"),legend.title = element_text(face = "bold")) + labs(x = "Cancer genes") + labs(y = "#Sample")+
  coord_flip()
dev.off() 






















