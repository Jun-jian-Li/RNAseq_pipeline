#####################################################################
#绘制火山图 
#####################################################################


#纵轴的参数：log10(FDR)
# yaxis = 10

#横轴的参数：logFC
# xaxis = 5

volcanno_funtion<-function(xaxis,yaxis){

Gene_result_exp<-read.table(paste0(combat_dir,comp_filename,"/",dir,"/All_Gene/All_result.txt"),header = T,sep = "\t")
        
#这里用的是pvalue筛选的值
logFC <-Gene_result_exp$logFC
pval <-Gene_result_exp$FDR
        
data <- data.frame(logFC=logFC,pval=pval)
data$sig[(data$pval > 0.05|data$pval=="NA")|(data$logFC < 1)& data$logFC > -1] <- "Not"
        
data$sig[data$pval <= 0.05 & data$logFC >= 1] <- "Up"
data$sig[data$pval <= 0.05 & data$logFC <= -1] <- "Down"
        
length(which(data$sig=="Up"))
length(which(data$sig=="Down"))
        
        # 选最大值作为xlim的上下边界
        #x_lim <- max(logFC,-logFC)
        # 绘制火山图
        library(ggplot2)
        library(RColorBrewer)
        
        # pdf(file = paste0(combat_dir,comp_filename,"/",dir,"/","volcano.pdf"),width=8,height=8)
        theme_set(theme_bw())
        
        ggplot(data,aes(logFC,-1*log10(pval),color = sig))+
               geom_point(alpha=0.5, size=3, 
                          aes(color=sig))+
               xlim(-xaxis,xaxis) +
               ylim(0,yaxis)+
               labs(x="logFC",y="-log10(FDR)")+
               scale_color_manual(values =c("#4DBBD5","grey","#E64B35"))+
               geom_hline(yintercept=-log10(0.05),linetype=1)+
               geom_vline(xintercept=-1,linetype=1,colour="#4DBBD5")+
               geom_vline(xintercept=1,linetype=1,colour="#E64B35")+
               geom_vline(xintercept=0,colour="#990000", linetype="dashed")+
               theme(panel.grid =element_blank())+expand_limits(y=c(0,3))+
               theme(axis.text=element_text(size=15),axis.title=element_text(size=15))+
               labs(colour = "Regulate")+
               theme(legend.position="top")
                     #legend.position="none")
        # dev.off()
        ggsave( paste0(combat_dir,comp_filename,"/",dir,"/","volcano.pdf"),width=8,height=8)
}
