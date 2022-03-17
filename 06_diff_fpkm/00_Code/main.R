pacman::p_load(rtracklayer,edgeR,dplyr,tidyverse,data.table,genefilter,sva,readr)

#01_合并注释文件表达文件----------------------------------------
mouse_GRC38_99<-rtracklayer::import('./01_Raw_Data/ref/Mus_musculus.GRCm38.99.chr.gtf.gz')
#注释文件
annotation=mouse_GRC38_99 %>% 
                  as.data.frame() %>% 
                  dplyr::filter(type == "gene") %>%                #选择gene进行分析，里面可能有lnc mir snRNA等基因
                  dplyr::select(seqnames,start,end,width,strand,gene_id,gene_name,gene_biotype) %>% 
                  dplyr::filter(gene_biotype == "protein_coding")  #选择编码蛋白质的基因进行分析

#count文件
count_raw_df<-readr::read_csv("./01_Raw_Data/ref_count/gene_ref_count.csv")%>%
        tidyr::separate(col = gene_id, into = c("gene_id", "gene_name"), sep = "\\|", remove = TRUE) %>% 
        rename_at(3:ncol(.),list(~ str_glue("{.}_count"))) %>%  ## 字符串格式化
        select(1,3:ncol(.))


#fpkm文件
fpkm_raw_df<-read_csv("01_Raw_Data/ref_fpkm/gene_ref_fpkm.csv") %>% 
             tidyr::separate(col = gene_id, into = c("gene_id", "gene_name"), sep = "\\|", remove = TRUE) %>% 
             rename_at(3:ncol(.),list(~ str_glue("{.}_fpkm"))) %>% 
             select(1,3:ncol(.))

#注释文件，count文件，fpkm文件合并
anno_count_fpkm<- annotation %>% 
  inner_join(count_raw_df, by="gene_id") %>% 
  inner_join(fpkm_raw_df, by="gene_id") %>% 
  dplyr::distinct(gene_name, .keep_all = TRUE)

save(anno_count_fpkm,file="01_Raw_Data/anno_count_fpkm.RData")
load("01_Raw_Data/anno_count_fpkm.Rdata")


#01.1载入样本分组文件-------------------------------------------------------------------
pacman::p_load(rtracklayer,edgeR,dplyr,tidyverse,data.table,genefilter,sva,readr)

Group_file<-"sample_info/Group_IRI_Sham.txt"
group_name<-str_split(gsub("sample_info/","",Group_file),"\\.")[[1]][1]

#参数设置

p_fdr = 0.05
logFC_up_threshold = 1
logFC_down_threshold = -1                                                                                                                                                                                                          
dir ="gene_protein" #这里有两个选项 gene/gene_protein
threshold = "FDR" #这里有两个选项 FDR/PValue

#读取分组文件
Group<-fread(Group_file)

#02样本的质量控制qc_PCA分析: ---------------------------------------------------
source("00_Code/PCA_sample_qc.R")

pdf( paste0("02_Diff_Analysis/",group_name,"_PCA.pdf" ) )
PCA_new(anno_count_fpkm[,paste0(Group$Sample,"_fpkm")],
  ntop = 500,
  group = Group$Group,
  show_name = T)
dev.off()


pdf( paste0("02_Diff_Analysis/",group_name,"_boxplot.pdf" ) )
boxplot(anno_count_fpkm[,paste0(Group$Sample,"_fpkm")],
  las = 2,
  outline = FALSE,
  col = as.factor(Group$Group))
dev.off()

#03样本的count值进行批次矫正----------------------------------------------------
#提取count的表达文件
# 
# source("00_Code/Combat_effect.R",encoding = "utf-8")
# adjusted_count<-Combat_count(Group,anno_count_fpkm)
# 
# #将原始数据中的count值进行替换
# anno_count_fpkm[,match(paste0(Group$Sample,"_count"),colnames(anno_count_fpkm))]<-adjusted_count
# 
# 
# #04将矫正后的count值进行fpkm值的转换------------------------------------------
# 
# #载入转化函数
# source("00_Code/Count_to_fpkm.R",encoding = "utf-8")
# 
# fpkms<-Count_to_fpkm(anno_count_fpkm,gene_length,Group)
# 
# #将fpkms中的gene_id按照我们的anno_count_fpkm 进行排序
# fpkms<- fpkms[match( anno_count_fpkm$gene_id , fpkms$gene_id ),]
# rownames(fpkms)<-fpkms$gene_id
# fpkms<-fpkms[,-1]
# 
# #将原始数据中的fpkm值进行替换
# anno_count_fpkm[,match(paste0(Group$Sample,"_fpkm"),colnames(anno_count_fpkm))]<-fpkms



#05去除批次之后的质量控制qc_PCA分析:--------------------------------------------

# source("00_Code/PCA_sample_qc.R")
# 
# pdf( "03_Combat_Diff_Analysis/Combat_after.pdf" )
# PCA_new(anno_count_fpkm[,paste0(Group$Sample,"_fpkm")],
#   ntop = 500,
#   group = Group$Group,
#   show_name = T)
# dev.off()

# boxplot(anno_count_fpkm[,paste0(Group$Sample,"_fpkm")],
#   las = 2,
#   outline = FALSE,
#   col = as.factor(Group$Group))

#03_差异分析--------------------------------------------------------------------
source("00_Code/Diff_edgeR.R", encoding = "utf-8")
source("00_Code/Diff_DESeq2.R", encoding = "utf-8")
# edgeR 
#分组的文件表

DEGs_res<-Diff_edgeR(anno_count_fpkm,Group)
# DEGs_res<-Diff_DESeq2(anno_count_fpkm,Group)

#比较组文件名字
source("00_Code/Extra_diff_result.R", encoding = "utf-8")

#03_提取差异结果----------------------------------------------------------------
#提取差异的结果，更改不同的阈值
df=anno_count_fpkm
comp_filename= group_name
combat_dir = "02_Diff_Analysis/"   #选择输出到combat（批次矫正后的文件夹中）
#该文件一旦报错就代表文件中阈值过高，筛选不出来上下调基因
Extra_diff_result(df=df,
                  comp_filename=comp_filename,
                  logFC_up_threshold=logFC_up_threshold,
                  logFC_down_threshold=logFC_down_threshold,
                  p_fdr=p_fdr,
                  dir=dir,
                  threshold=threshold,
                  combat_dir= combat_dir)

# [1] "差异上调基因为：1786"
# [1] "差异下调基因为：1841"

#04_绘制差异分析火山图----------------------------------------------------------
source("00_Code/Volcano.R", encoding = "utf-8")

xaxis = 10
yaxis = 10

volcanno_funtion(xaxis = xaxis,yaxis = yaxis)


#05_绘制上下调基因热图----------------------------------------------------------
source("00_Code/Heatmap.R", encoding = "utf-8")

seq = "fpkm"
heatmap_height = 0.08
heatmap_cellwidth = 20

heatmap_funtion(seq = seq,
                heatmap_height = heatmap_height,
                heatmap_cellwidth = heatmap_cellwidth,
                heatmap_angle_col= "45")


#06_GO_KEGG富集分析---------------------------------------------------------
source("00_Code/Enrich_GO_KEGG.R", encoding = "utf-8")

dir = dir
species = "mouse" #"human"

GO_KEGG_enrichment_funtion(dir,species)

#07_GO富集分析图片----------------------------------------------------------
source("00_Code/GO_BP_CC_MF_barplot.R", encoding = "utf-8")

file="_up" # /"_down" /""
height = 8
width = 10
GO_plot(file = file,height = height,width = width ,dir=dir,p_fdr=p_fdr,combat_dir=combat_dir,comp_filename=comp_filename)


#08_KEGG富集分析图片----------------------------------------------------------
source("00_Code/KEGG_Bubble_plot.R", encoding = "utf-8")

file="_up" # /"_down" /""
height = 8
width = 10

KEGG_plot(file = file,height = height,width = width,dir=dir,p_fdr=p_fdr,combat_dir=combat_dir,comp_filename=comp_filename)
  

#09_对转录本进行分析
pacman::p_load(rtracklayer,edgeR,dplyr,tidyverse,data.table,genefilter,sva,readr)

#01_合并注释文件表达文件----------------------------------------
mouse_GRC38_99<-rtracklayer::import('./01_Raw_Data/ref/Mus_musculus.GRCm38.99.chr.gtf.gz')
#注释文件
annotation_trans=mouse_GRC38_99 %>% 
  as.data.frame() %>% 
  dplyr::filter(type == "transcript") %>%                #选择转录本进行分析
  dplyr::select(seqnames,start,end,width,strand,gene_id,gene_name,gene_biotype,c)

#count文件
trans_count_raw_df<-readr::read_csv("./01_Raw_Data/ref_count/transcript_ref_count.csv")%>%
  rename_at(2:ncol(.),list(~ str_glue("{.}_count"))) %>%  ## 字符串格式化
  select(1,2:ncol(.))

#fpkm文件
trans_fpkm_raw_df<-read_csv("01_Raw_Data/ref_fpkm/transcript_ref_fpkm.csv") %>% 
  rename_at(2:ncol(.),list(~ str_glue("{.}_fpkm"))) %>% 
  select(1,2:ncol(.))

#注释文件，count文件，fpkm文件合并
trans_anno_count_fpkm<- annotation_trans %>% 
  inner_join(trans_count_raw_df, by="transcript_id") %>% 
  inner_join(trans_fpkm_raw_df, by="transcript_id") %>% 
  dplyr::distinct(transcript_id, .keep_all = TRUE)

save(trans_anno_count_fpkm,file="01_Raw_Data/trans_anno_count_fpkm.RData")
load("01_Raw_Data/trans_anno_count_fpkm.Rdata")

source("00_Code/Trans_diff_edgeR.R", encoding = "utf-8")
# edgeR 
#分组的文件表

trans_DEGs_res<-Diff_edgeR(trans_anno_count_fpkm,Group)


#比较组文件名字
source("00_Code/Trans_extra_diff_result.R", encoding = "utf-8")

#03_提取差异结果----------------------------------------------------------------
#提取差异的结果，更改不同的阈值
df=trans_anno_count_fpkm
comp_filename= group_name
combat_dir = "02_Diff_Analysis/"   #选择输出到combat（批次矫正后的文件夹中）

p_fdr = 0.05
logFC_up_threshold = 1
logFC_down_threshold = -1                                                                                                                                                                                                          
dir ="trans" #这里有两个选项 gene/gene_protein
threshold = "FDR" #这里有两个选项 FDR/PValue


#该文件一旦报错就代表文件中阈值过高，筛选不出来上下调基因
trans_Extra_diff_result(df=df,
                  comp_filename=comp_filename,
                  logFC_up_threshold=logFC_up_threshold,
                  logFC_down_threshold=logFC_down_threshold,
                  p_fdr=p_fdr,
                  dir=dir,
                  threshold=threshold,
                  combat_dir= combat_dir)




