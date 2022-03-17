library(DESeq2)
# 函数说明：edger差异分析、输入：样本矩阵，矩阵中包含样本名字以及gene_id,以及分组的文件Group
Diff_DESeq2<-function(data,Group){
  #将Treat样本与Control样本分开
  Treat<-paste0(Group$Sample[which(Group$Group == "Treat")],"_count")
  Control<-paste0(Group$Sample[which(Group$Group == "Control")],"_count")
  Treat_Control<-anno_count_rpkm[,c("gene_id",Treat,Control)] %>%
    as.data.frame() %>% 
    column_to_rownames("gene_id") 
  group_info = factor(c(rep("Treat",length(Treat)), rep("Control",length(Control)) ),levels = c("Control","Treat"))
  
  colData <- data.frame(row.names = colnames(Treat_Control),
                        condition = group_info )
  
  dds <- DESeqDataSetFromMatrix(countData = Treat_Control, 
                                colData = colData, 
                                design = ~ condition) 
  dds <- DESeq(dds) 
  res <- results(dds)
  res <- as.data.frame(res) 
  res <- cbind(rownames(res), res) 
  DEGs_res <- res[,c(1,3,6,7)]
  colnames(DEGs_res) <- c("gene_id","logFC","PValue", "FDR")
  return(DEGs_res)
}



