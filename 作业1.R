setwd("xena")
library(tidyverse)
fpkm1 = read.table(file = 'TCGA-LIHC.htseq_fpkm.tsv', sep = '\t', header = TRUE) 
rownames(fpkm1) <- fpkm1[,1] 
fpkm1 = fpkm1[,-1]
table(substr(colnames(fpkm1),14,16))
fpkm1 <- fpkm1[,substr(colnames(fpkm1),14,16)%in% c("01A","11A")]
table(substr(colnames(fpkm1),14,16))
rownames(fpkm1) <- substr(rownames(fpkm1),1,15)#更改行名为前十五位
fpkm <- fpkm1



Ginfo_0 <- read.table("gene_length_Table.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
Ginfo <- Ginfo_0[which(Ginfo_0$genetype == "protein_coding"),] 
comgene <- intersect(rownames(fpkm),rownames(Ginfo))
fpkm <- fpkm[comgene,]
class(fpkm)#判断数据类型
class(comgene)
Ginfo <- Ginfo[comgene,]
fpkm$Gene <- as.character(Ginfo$genename)   #新增Gene Symbol
fpkm <- fpkm[!duplicated(fpkm$Gene),]   #去重复
rownames(fpkm) <- fpkm$Gene   #将行名变为Gene Symbol
ncol(Ginfo)#计数 列
#nrow 计数 行
fpkm <- fpkm[,-ncol(fpkm)]   #去除最后一列
write.table(fpkm, file = "LIHC_fpkm_mRNA_all.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#保存癌症患者的fpkm
tumor <- colnames(fpkm)[substr(colnames(fpkm),14,16) == "01A"]
fpkm_01A <- fpkm[,tumor]
write.table(fpkm_01A, file = "LIHC_fpkm_mRNA_01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#保存正常样本的fpkm文件
normal <- colnames(fpkm)[substr(colnames(fpkm),14,16) == "11A"]
fpkm_11A <- fpkm[,normal]
write.table(fpkm_11A, file = "LIHC_fpkm_mRNA_11A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#整理完毕#

library(DESeq2)

fpkm_01A <- read.table("LIHC_fpkm_mRNA_01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
fpkm_11A <- read.table("LIHC_fpkm_mRNA_11A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#res需要打开文件夹读取
res_deseq2 <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 3, padj < 0.05)#根据自己需要
#DEG:differentially expressed genes
#读取差异基因文件
#直接在文件夹双击