#############################################################
#
#  绘制GO富集结果柱状图 BP—CC-MF各自提取 top10
#
#############################################################
library(dplyr)

GO_plot<-function(file = file,height = height,width = width,dir=dir,p_fdr=p_fdr,combat_dir=combat_dir,comp_filename=comp_filename){
   
   
   # file="_up"
   # height = 8
   # width = 15
   
   
   file = paste0("DEGs_res_diff",file)
   GO_diff_gene<-read.table(paste0(combat_dir,comp_filename,"/",dir,"/GO/",file,p_fdr,".txt"),header = T,sep = "\t")
   
   #三种类型各提取前十个作为作图数据
   GO_BP<-GO_diff_gene[which(GO_diff_gene$ONTOLOGY == "BP"),]
   GO_BP_top10<-arrange(GO_BP,desc(Count))[1:10,]
   
   GO_CC<-GO_diff_gene[which(GO_diff_gene$ONTOLOGY == "CC"),]
   GO_CC_top10<-arrange(GO_CC,desc(Count))[1:10,]
   
   GO_MF<-GO_diff_gene[which(GO_diff_gene$ONTOLOGY == "MF"),]
   GO_MF_top10<-arrange(GO_MF,desc(Count))[1:10,]
   
   #合并富集矩阵
   GO_top10<-rbind(GO_BP_top10,GO_CC_top10,GO_MF_top10)
   
   #-------------------------------------------------------------------------------
   colnames(GO_top10)[1]<-"Ontology"
   GO_top10<-na.omit(GO_top10[30:1,])
   
   
   GO_term_order=factor(as.integer(rownames(GO_top10)),
                        labels=GO_top10$Description
                        ,levels=c(30:1))
   COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")
   
   ggplot(data=GO_top10, 
          aes(x=GO_term_order,
              y=Count, 
              fill=Ontology)) +
      geom_bar(stat="identity", width=0.85)  + 
      scale_fill_manual(values = COLS) +
      # theme_bw() +
      coord_flip() + 
      xlab("GO Terms") +
      ylab("Count") + 
      labs(title = "The Most Enriched GO Terms")+ 
      # coord_flip() +
      theme(
         axis.text.x=element_text(face = "bold"
                                  # color="gray50",
                                  # angle = 70,
                                  # vjust = 1, 
                                  # hjust = 1 
         ),
         axis.text.y=element_text(face = "bold"),   #去除背景
         legend.key = element_blank(),
         panel.background = element_blank(),
         panel.border = element_rect(color="black",
                                     fill = "transparent")
      ) 
   
   ggsave(paste0(combat_dir,comp_filename,"/",dir,"/GO/",file,p_fdr,".pdf"),height = height ,width = width)
}

   







