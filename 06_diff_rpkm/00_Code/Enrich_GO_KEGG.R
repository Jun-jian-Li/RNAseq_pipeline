#ClusterProfiler富集分析 GO—KEGG#####输入的是GENE—Symbol########################
# BiocManager::install("KEGG.db")
pacman::p_load(clusterProfiler,DO.db,org.Mm.eg.db,ggplot2,stringr,org.Hs.eg.db)

#enrichKEGG使用在线数据速度实在是太慢了，所以可以先使用createKEGGdb生成本地KEGG.db包。
# remotes::install_github("YuLab-SMU/createKEGGdb")
# library(createKEGGdb)
# create_kegg_db("mmu")
# install.packages("KEGG.db_3.2.4.tar.gz",repos=NULL,type="source")
#当使用本地的KEGG.db数据库时，use_internal_data设置为T,使用在线数据库时设置为F.

dir_names = "gene_protein"
species = "mouse" #/ human
GO_KEGG_enrichment_funtion<-function(dir,species){
     
     if(species == "human"){
       species_package = "org.Hs.eg.db"
       organism = "Hs"
     }else{
       species_package = "org.Mm.eg.db"
       organism = "mmu"
     }

    dir_names = dir
    #循环读取文件夹
    for(dir in dir_names){
    
    #获取比较组文件######################################################
    #gene文件富集
    comp_filename<-list.files("02_Diff_Analysis")[ grep( "Group",unlist(list.files("02_Diff_Analysis")) )]
    
    for(i in comp_filename)
    {
      # i = "01_E125_Cfp1_WT"
      #获取富集文件
      files<-list.files(paste0("02_Diff_Analysis/",i,"/",dir,"/"))
      files<-files[grep("*.txt",files)]
    
      #输入文件名字，输出文件在目录下的GO KEGG文件下
      #创建GO文件夹
      if(!dir.exists(paste0("02_Diff_Analysis/",i,"/",dir,"/","GO"))){
        dir.create(paste0("02_Diff_Analysis/",i,"/",dir,"/","GO")) 
      }
      #创建KEGG文件夹
      if(!dir.exists(paste0("02_Diff_Analysis/",i,"/",dir,"/","KEGG"))){
        dir.create(paste0("02_Diff_Analysis/",i,"/",dir,"/","KEGG")) 
      }
     
      enrich_GO_KEGG<-function(file_path){
        
        file = str_split(file_path,"/")[[1]][4]
        diff<-read.table(file_path,header = T)
        Gene<-diff$gene_name
        #转为entrenz_id enrich#####################################################
        #bitr()函数可以转换gene ID -> gene symbo
        #ENSEMBL / SYMBOL /ENTREZID
        #只选已有的id富集分析
        anyid <- bitr( Gene,
                       fromType = "SYMBOL",   #Symbol转entrenz id 
                       toType = c("ENTREZID"),
                       OrgDb = species_package,drop = TRUE)
        
        #BP Biological process########################################################
        go<- enrichGO(anyid$ENTREZID,  #输入symbol
                      OrgDb = species_package, 
                      ont='all',
                      pAdjustMethod = 'BH',
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.2)
        go=DOSE::setReadable(go, OrgDb= species_package,keyType='ENTREZID')
        go<-as.data.frame(go)
        #写出富集结果文件
        if(nrow(go) != 0){
          go<-go[order(go$ONTOLOGY,go$p.adjust),]
          write.table(go,paste0("02_Diff_Analysis/",i,"/",dir,"/GO/",file),quote=F,sep = "\t", row.names = F)
        }
        #排序
        #KEGG#########################################################################
        kegg <- enrichKEGG(anyid$ENTREZID, 
                           organism = organism, 
                           keyType = 'kegg', 
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 0.2, 
                           pAdjustMethod = 'BH'
                           # use_internal_data = T
        )
        #将gene id -> gene name
        kegg<-setReadable(kegg, OrgDb = species_package, keyType="ENTREZID")
        KEGG<-as.data.frame(kegg)
        if(nrow(KEGG) != 0){
          KEGG<-KEGG[order(KEGG$p.adjust),]
          write.table(KEGG,paste0("02_Diff_Analysis/",i,"/",dir,"/KEGG/",file),quote=F,sep = "\t", row.names = F)
        }
      }
      sapply( paste0("02_Diff_Analysis/",i,"/",dir,"/",files) , enrich_GO_KEGG)
    }
  }
}
