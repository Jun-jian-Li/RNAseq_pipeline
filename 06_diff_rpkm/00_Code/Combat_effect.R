#-----------------------------------------------------------
# author : Junjian Li
# time : 20211105
# description : 对RNAseq数据进行批次矫正的分析
#               It uses a negative binomial regression to model the count matrix, and estimate parameters representing the batch effects. 
#               Then it provides adjusted data by mapping the original data to an expected distribution if there were no batch effects.
# 参数说明：1：count矩阵 
#           2：批次的分组信息 batch
#          3: 生物学的分组信息
#-----------------------------------------------------------

# devtools::install_github("zhangyuqing/sva-devel")


Combat_count<-function(Group,anno_count_fpkm){
  
  #读取比较组的分组
  group<-ifelse(Group$Group == "Treat",1,0)
  
  #读取批次的分组
  batch<-Group$Batch
  
  #读取count的文件数据
  raw_count<-anno_count_fpkm %>% 
              dplyr::select("gene_id" , paste0(Group$Sample,"_count")) %>% 
              column_to_rownames("gene_id") %>% 
              as.matrix()
  
  
  #矫正我们批次的数据
  adjusted_count <- sva::ComBat_seq(raw_count, batch = batch, group = group)
  return(adjusted_count)

}










