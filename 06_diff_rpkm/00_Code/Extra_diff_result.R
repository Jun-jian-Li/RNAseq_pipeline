#提取差异表达的结果
#---------------------------------
#参数说明：df:分析后的全部结果
#comp_filename：比较组
#logFC_down_threshold：差异上调的阈值
#threshold：FDR 还是 P
#p_fdr ：差异显著性水平数值
#dir：基因文件夹
#---------------------------------
pacman::p_load(tidyverse)
Extra_diff_result<-function(df, 
                            comp_filename,
                            logFC_up_threshold ,
                            logFC_down_threshold,
                            threshold,
                            p_fdr,
                            dir,
                            combat_dir){
  #提取所有的结果
  
  # df=anno_count_rpkm
  # comp_filename="E12.5_cfp1_wt"
  # logFC_up_threshold = 0.2
  # logFC_down_threshold = -0.2
  # threshold == 
  # FDR_threshold  = 0.001
  # Pval_threshold = 0.05
  # dir ="gene"
  # p_fdr = 0.05
  threshold = threshold
  if( threshold == "FDR"){
    col = "FDR"
  }else{
    col = "PValue"
  }
  #建立比较组文件夹
  if(!exists(file.path(combat_dir,comp_filename))){dir.create(file.path(combat_dir,comp_filename))}
  #建立差异基因文件夹
  if(!exists(file.path(combat_dir,comp_filename,dir))){dir.create(file.path(combat_dir,comp_filename,dir))}
  #建立结果文件夹
  if(!exists(file.path(combat_dir,comp_filename,dir,"All_Gene"))){dir.create(file.path(combat_dir,comp_filename,dir,"All_Gene"))}
  
  #合并差异结果以及表达注释文件夹
  anno_count_rpkm_diff<-df %>% inner_join(DEGs_res,"gene_id")
  
  save(anno_count_rpkm_diff,file=file.path(combat_dir,comp_filename,dir,"All_Gene","anno_count_rpkm_diff.Rdata")) 
  fwrite(anno_count_rpkm_diff,file=file.path(combat_dir,comp_filename,dir,"All_Gene","anno_count_rpkm_diff.txt"),row.names = F,quote = F,sep = "\t")
  
  #提取该比较组下的样本的所有差异结果
  Treat_count<-paste0(Group$Sample[which(Group$Group == "Treat")],"_count")
  Control_count<-paste0(Group$Sample[which(Group$Group == "Control")],"_count")
  
  Treat_rpkm<-paste0(Group$Sample[which(Group$Group == "Treat")],"_rpkm")
  Control_rpkm<-paste0(Group$Sample[which(Group$Group == "Control")],"_rpkm")
  
  
  col_in<-c("seqnames","start",	"end","width","strand","gene_id","gene_name",
    "gene_biotype","logFC","PValue","FDR",Treat_count,Control_count,Treat_rpkm,Control_rpkm)
  
  All_result<-anno_count_rpkm_diff %>%as.data.frame() %>%  dplyr::select(all_of(col_in)) %>% as.data.frame() 
  fwrite(All_result,file=file.path(combat_dir,comp_filename,dir,"All_Gene","All_result.txt"),row.names = F,quote = F,sep = "\t")
  
  #提取差异基因的结果
 
  data<-All_result
  etSig <- data[which(data[,col] < p_fdr & abs(data$logFC) > logFC_up_threshold),] %>% 
         dplyr::mutate(Regulation=if_else(logFC > abs(logFC_up_threshold),"Up","Down"))
   
  if(nrow(etSig) == 0){print("通过阈值筛选后没有差异表达的基因，请重新选取阈值")}else{
    write.table(etSig,file.path(combat_dir,comp_filename,dir,paste0("DEGs_res_diff",p_fdr,".txt")),row.names = F,quote = F,sep = "\t")
  }
  
  #up 基因筛选
  etSig_up<-etSig[which(etSig$Regulation == "Up"),]
  if(nrow(etSig_up) == 0){print("通过阈值筛选后没有Up的差异表达的基因，请重新选取阈值")}else{
    write.table(etSig_up,file.path(combat_dir,comp_filename,dir,paste0("DEGs_res_diff_up",p_fdr,".txt")),row.names = F,quote = F,sep = "\t")
  }
  print(paste0("差异上调基因为：",nrow(etSig_up)))
  #down 基因筛选
  etSig_down<-etSig[which(etSig$Regulation == "Down"),]
  if(nrow(etSig_down) == 0){print("通过阈值筛选后没有Down的差异表达的基因，请重新选取阈值")}else{
    write.table(etSig_down,file.path(combat_dir,comp_filename,dir,paste0("DEGs_res_diff_down",p_fdr,".txt")),row.names = F,quote = F,sep = "\t")
  }
  print(paste0("差异下调基因为：",nrow(etSig_down)))
  
}

