#提取差异表达的结果
#---------------------------------
#参数说明：df:分析后的全部结果
#comp_filename：比较组
#logFC_down_threshold：差异上调的阈值
#threshold：FDR 还是 P
#p_fdr ：差异显著性水平数值
#dir：转录本文件夹
#---------------------------------
pacman::p_load(tidyverse)
trans_Extra_diff_result<-function(df, 
                            comp_filename,
                            logFC_up_threshold ,
                            logFC_down_threshold,
                            threshold,
                            p_fdr,
                            dir,
                            combat_dir){
  #提取所有的结果
  
  # df=anno_count_fpkm
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
  #建立差异转录本文件夹
  if(!exists(file.path(combat_dir,comp_filename,dir))){dir.create(file.path(combat_dir,comp_filename,dir))}
  #建立结果文件夹
  if(!exists(file.path(combat_dir,comp_filename,dir,"All_Trans"))){dir.create(file.path(combat_dir,comp_filename,dir,"All_Trans"))}
  
  #合并差异结果以及表达注释文件夹
  trans_anno_count_fpkm<-df %>% inner_join(trans_DEGs_res,"transcript_id")
  
  save(trans_anno_count_fpkm,file=file.path(combat_dir,comp_filename,dir,"All_Trans","trans_anno_count_fpkm.Rdata")) 
  fwrite(trans_anno_count_fpkm,file=file.path(combat_dir,comp_filename,dir,"All_Trans","trans_anno_count_fpkm.txt"),row.names = F,quote = F,sep = "\t")
  
  #提取该比较组下的样本的所有差异结果
  Treat_count<-paste0(Group$Sample[which(Group$Group == "Treat")],"_count")
  Control_count<-paste0(Group$Sample[which(Group$Group == "Control")],"_count")
  
  Treat_fpkm<-paste0(Group$Sample[which(Group$Group == "Treat")],"_fpkm")
  Control_fpkm<-paste0(Group$Sample[which(Group$Group == "Control")],"_fpkm")
  
  
  col_in<-c("seqnames","start",	"end","width","strand","transcript_id","gene_name",
            "gene_biotype","logFC","PValue","FDR",Treat_count,Control_count,Treat_fpkm,Control_fpkm)
  
  All_result<-trans_anno_count_fpkm %>%as.data.frame() %>%  dplyr::select(all_of(col_in)) %>% as.data.frame() 
  fwrite(All_result,file=file.path(combat_dir,comp_filename,dir,"All_Trans","All_result.txt"),row.names = F,quote = F,sep = "\t")
  
  #提取差异转录本的结果
  
  data<-All_result
  etSig <- data[which(data[,col] < p_fdr & abs(data$logFC) > logFC_up_threshold),] %>% 
    dplyr::mutate(Regulation=if_else(logFC > abs(logFC_up_threshold),"Up","Down"))
  
  if(nrow(etSig) == 0){print("通过阈值筛选后没有差异表达的转录本，请重新选取阈值")}else{
    write.table(etSig,file.path(combat_dir,comp_filename,dir,paste0("trans_DEGs_res_diff",p_fdr,".txt")),row.names = F,quote = F,sep = "\t")
  }
  
  #up 转录本筛选
  etSig_up<-etSig[which(etSig$Regulation == "Up"),]
  if(nrow(etSig_up) == 0){print("通过阈值筛选后没有Up的差异表达的转录本，请重新选取阈值")}else{
    write.table(etSig_up,file.path(combat_dir,comp_filename,dir,paste0("trans_DEGs_res_diff_up",p_fdr,".txt")),row.names = F,quote = F,sep = "\t")
  }
  print(paste0("差异上调转录本为：",nrow(etSig_up)))
  #down 转录本筛选
  etSig_down<-etSig[which(etSig$Regulation == "Down"),]
  if(nrow(etSig_down) == 0){print("通过阈值筛选后没有Down的差异表达的转录本，请重新选取阈值")}else{
    write.table(etSig_down,file.path(combat_dir,comp_filename,dir,paste0("trans_DEGs_res_diff_down",p_fdr,".txt")),row.names = F,quote = F,sep = "\t")
  }
  print(paste0("差异下调转录本为：",nrow(etSig_down)))
  
}

