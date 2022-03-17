#-----------------------------------------------------------
# author : Junjian Li
# time : 20211105
# description : 对RNAseq数据进行Count值转变为FPKM
#               由于对count表达谱矫正之后,其Count值是矫正后的，FPKM值应该进行相应的矫正 
# 参数说明：1：count矩阵  anno_count_fpkm
#           2：计算每个基因的外显子总长 gene_length
#           3: 生物学的分组信息 Group
# 实验注意：现在用的R = 4.1.1 这个包中某些函数和这个R语言版本不匹配，所以需要对R语言降级
#-----------------------------------------------------------


# 标准化内容：基因的长度（基因中的外显子长度总和）和所有基因上计数后的reads总数(测序深度的)
# BiocManager::install("GenomicFeatures",force = TRUE)
#计算外显子的长度
pacman::p_load(GenomicFeatures)  


  
  if(file.exists("gene_length_mouse.Rdata")){
    load("gene_length_mouse.Rdata")
  }else{
    #载入参考基因组注释文件
    txdb <- makeTxDbFromGFF("./01_Raw_Data/ref/Mus_musculus.GRCm38.99.chr.gtf",format="gtf")
    
    #获取每个gene上的所有外显子的起始位点和终止位点
    exons.list.per.gene <- exonsBy(txdb, by = "gene")
    
    #reduce去除掉重叠冗余的部分，最后计算长度
    exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))})
    
    gene_length <- do.call(rbind,lapply(exonic.gene.sizes, data.frame))
    
    gene_length<-data.frame(gene=rownames(gene_length),length=gene_length)
    colnames(gene_length)<-c("gene_id","length")
    
    save(gene_length,file = "gene_length_mouse.Rdata")
  }
  
  
Count_to_Fpkm<-function( anno_count_fpkm,gene_length,Group  ){
  
  # count 的列选择
  
  #选择交集的gene_id : anno_count_fpkm
  intersect_gene_id<-intersect(gene_length$gene_id,anno_count_fpkm$gene_id)
  
  #外显子长度矩阵
  count_effLen<-anno_count_fpkm[match( intersect_gene_id,  anno_count_fpkm$gene_id  ) ,c("gene_id" , paste0(Group$Sample,"_count"))] %>% 
    as.data.frame() %>% 
    inner_join(gene_length,"gene_id" ) %>% column_to_rownames("gene_id")
  
  
  effLen<-count_effLen$length
  count<-count_effLen %>% as.data.frame() %>% dplyr::select(-length)
  
  
  #构造count文件
  countToFpkm <- function(count, effLen){
    N <- sum(count)
    exp( log(count) + log(1e9) - log(effLen) - log(N) )
  }
  
  # FPKM2TPM <- function(fpkm){
  #   exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
  # }
  
  
  ##count转换为FPKM值
  fpkms <- apply(count, 2, countToFpkm, effLen = effLen) %>% as.data.frame() %>% rownames_to_column("gene_id")
  
  # tpms <- apply(fpkms,2,FPKM2TPM)
  
  return(fpkms)
}



