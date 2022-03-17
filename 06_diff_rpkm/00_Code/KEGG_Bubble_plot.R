#########################################
###        KEGG Pathway Plot         ####
#########################################

library(readxl)
library(ggplot2)

KEGG_plot<-function(file = file,height = height,width = width,dir=dir,p_fdr=p_fdr,combat_dir=combat_dir,comp_filename=comp_filename){
   
   file = paste0("DEGs_res_diff",file)
   KEGGPathway<-read.table(paste0(combat_dir,comp_filename,"/",dir,"/KEGG/",file,p_fdr,".txt"),header = T,sep = "\t")
   
   #绘制KEGG前20个通路
   KEGGPathway<-KEGGPathway[1:20,]
   #读入数据
   colnames(KEGGPathway)[c(2,5,6,7)]<-c("KEGG Pathway","Pvalue","Padj","Qvalue")
   
   #将Pathway列转化为因子型
   KEGGPathway$`KEGG Pathway`<-factor(KEGGPathway$`KEGG Pathway`,levels = rev(KEGGPathway$`KEGG Pathway`))
   
   #作图
   #图片背景设定
   mytheme <-
   
   
  ggplot(KEGGPathway,aes(`Pvalue`,`KEGG Pathway`)) +
      geom_point(aes(fill=`Pvalue`,size=Count),alpha=0.9,pch=21,colour="gray25") +  #fill对应点的填充色，colour对应点的边框色
      scale_fill_gradient(low='#F97F51', high='#1B9CFC') + #设定颜色的变化范围
      scale_size_area(max_size = 10, breaks=c(10,20,30,40,50)) + #设定点的大小比例和图例上显示的间隔
      labs( title="The Enriched KEGG Pathways",
            y='KEGG Pathway',x='Pvalue',fill='Pvalue',size='Gene Count' )+ 
      theme(axis.title=element_text(face="bold", size=11,colour = 'gray25'), #坐标轴标题
            axis.text=element_text(face="bold", size=11,colour = 'gray25'), #坐标轴标签
            axis.line = element_line(size=0.5, colour = 'black'),           #轴线
            panel.background = element_rect(color='black'), #绘图区边框
            legend.key = element_blank()                    #关闭图例边框
            ) 
   

   ggsave(paste0(combat_dir,comp_filename,"/",dir,"/KEGG/",file,p_fdr,".pdf"),height = height ,width = width)
}



