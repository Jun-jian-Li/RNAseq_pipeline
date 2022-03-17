#函数说明--------------------------------------------------
# author ：Junjian Li
# time : 20211105 
# description : 该函数对差异表达分析的样本组进行pca分析 
# 参数说明 ：expr ：表达谱
#            ntop ：取前多少个基因进行分析
#            show_name : 是否显示name
#            group : 分组信息
#----------------------------------------------------------


PCA_new <- function(expr, ntop = 500, group, show_name = F){
  library(ggplot2)
  library(ggrepel)
  object <- expr
  rv <- genefilter::rowVars(object)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(object[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  d <- data.frame(PC1 = pca$x[, 1], 
    PC2 = pca$x[, 2], 
    group = group, 
    name = colnames(object))
  attr(d, "percentVar") <- percentVar[1:2]
  if (show_name) {
    ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) + 
      geom_point(size = 2) +
      xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
      ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
      geom_text_repel(aes(label = name),
        size = 3,
        segment.color = "black",
        show.legend = FALSE )
  } else {
    ggplot(data = d, aes_string(x = "PC1", y = "PC2",color = "group")) + 
      geom_point(size = 2) +
      xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
      ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance"))
  }
}
