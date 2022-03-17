library(Hmisc)
library(tidyverse)

#测序类型：single --> rpkm
#测序类型：pair --> fpkm


heatmap_funtion<-function(seq,heatmap_height,heatmap_cellwidth,heatmap_angle_col){
  
  #-参数-------------------
  
  seq = seq 
  heatmap_height = heatmap_height
  heatmap_cellwidth = heatmap_cellwidth
  heatmap_angle_col = heatmap_angle_col
  
  #-----------------------
  file_path_down<- paste0(combat_dir,comp_filename,"/",dir,"/","DEGs_res_diff_down",p_fdr,".txt")
  file_path_up<- paste0(combat_dir,comp_filename,"/",dir,"/","DEGs_res_diff_up",p_fdr,".txt")
  
  Gene_down<-read.table(file_path_down,header = T)
  Gene_up<-read.table(file_path_up,header = T)
  
  Gene_down<-Gene_down[,c(7,grep(paste0(".+",seq),colnames(Gene_down)))]
  Gene_down<-limma::avereps(Gene_down[,-1],ID = Gene_down$gene_name)
  
  Gene_up<-Gene_up[,c(7,grep(paste0(".+",seq),colnames(Gene_up)))]   
  Gene_up<-limma::avereps(Gene_up[,-1],ID = Gene_up$gene_name)
  
  exp<-rbind(Gene_down,Gene_up)
  
  pdf(paste0(combat_dir,comp_filename,"/",dir,"/","heatmap.pdf"),height = 8,width = 8)
  library(pheatmap)
  pheatmap(exp,
           scale = "row",
           color = colorRampPalette(c("#4DBBD5","white", "#E64B35"))(50),
           show_rownames=F,
           border=FALSE,
           cluster_rows = T,
           cluster_col = FALSE,
           cellheight = heatmap_height,
           cellwidth = heatmap_cellwidth,
           angle_col = heatmap_angle_col
  )
  dev.off()
}


