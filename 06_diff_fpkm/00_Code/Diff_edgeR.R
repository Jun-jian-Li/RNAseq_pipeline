library(edgeR)
# 函数说明：edger差异分析、输入：样本矩阵，矩阵中包含样本名字以及gene_id,以及分组的文件Group
Diff_edgeR<-function(data,Group){
  #将Treat样本与Control样本分开
  Treat<-paste0(Group$Sample[which(Group$Group == "Treat")],"_count")
  Control<-paste0(Group$Sample[which(Group$Group == "Control")],"_count")
  Treat_Control<-anno_count_rpkm[,c("gene_id",Treat,Control)] %>%
    as.data.frame() %>% 
    column_to_rownames("gene_id") 
  group_info = factor(c(rep("Treat",length(Treat)), rep("Control",length(Control)) ),levels = c("Control","Treat"))
  # make obj 
  dge.list.obj <- DGEList(counts = Treat_Control, group = group_info)
  #对数据进行过滤，去除低质量的基因
  # keep <- rowSums(cpm(dge.list.obj) > 1 ) >= 2
  # length(which(keep)) 
  # dge.list.obj <- dge.list.obj[keep,,keep.lib.sizes = T] #过滤结束后重新计算样本库size，一般推荐是这样做的。尽管做不做对后续影响都不大.
  # Normalization Normalization method: "TMM","TMMwsp","RLE","upperquartile","none"
  dge.list.obj <- calcNormFactors(dge.list.obj,method = "TMM")
  # make design matrix
  design.mat <- model.matrix(~group_info)
  # estimate dispersion 离散度的评估
  dge.list.obj <- estimateDisp(dge.list.obj,design.mat)
  # test with likelihood ratio test 最大似然估计
  dge.list.res <- exactTest(dge.list.obj)
  DEGs_res <- as.data.frame(topTags(dge.list.res,n=nrow(Treat_Control),sort.by = "logFC"))
  DEGs_res<-data.frame(gene_id=rownames(DEGs_res),DEGs_res)
  
  DEGs_res<-DEGs_res[,c(1,2,4,5)]
  return(DEGs_res)
}

#提取差异表达函数

