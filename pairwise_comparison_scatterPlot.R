#make scatter plot for all sRNA-mRNA interactions 
library(dplyr)
library(ggplot2)
library(readr)
setwd("~/Desktop/phd/programming/Integrated_analysis/WGCNA/")

lnames = load(file = "WGCNA.RData")
#read metatable_sRNA
#par(mar=c(1,1,1,1))
system("rm -r gene_pairwise_correlation")
system("mkdir gene_pairwise_correlation")
meta=as.data.frame(read_delim("~/Desktop/phd/programming/Integrated_analysis/WGCNA/co_expression_metatable_sRNA.tsv",delim="\t"))
scatter=function(gene1,gene2,df,metatable){
  # if (!file.exists(paste0("gene_pairwise_correlation/",gene1))){
  # system(paste0("mkdir gene_pairwise_correlation/",gene1))
  # }
  
  #  par(mar=c(4,5,4,5))
  print(gene1)
  print(gene2)
  filtered_row=metatable %>% filter(fromNode==gene1,toNode==gene2)
  #filtered_row=
  bicor_coe=filtered_row$correlation
  padj=filtered_row$correlation_padj
  
  #gene1_data=df %>% dplyr::select(gene1)
  #gene2_data=df %>% dplyr::select(gene2)
  
  #  x=sapply(gene1_data,as.numeric)
  # y=sapply(gene2_data,as.numeric)
  
  gene1=gsub(" ","",gene1)
  gene1=gsub("%","",gene1)
  gene1=gsub("\\(","",gene1)
  gene1=gsub(")","",gene1)
  gene1=gsub("'","",gene1)
  gene1=gsub("-","",gene1)
  
  gene2=gsub(" ","",gene2)
  gene2=gsub("%","",gene2)
  gene2=gsub("\\(","",gene2)
  gene2=gsub(")","",gene2)
  gene2=gsub("'","",gene2)
  gene2=gsub("-","",gene2)
  
  colnames(df)=gsub(" ","",colnames(df))
  colnames(df)=gsub("%","",colnames(df))
  colnames(df)=gsub("\\(","",colnames(df))
  colnames(df)=gsub(")","",colnames(df))
  colnames(df)=gsub("'","",colnames(df))
  colnames(df)=gsub("-","",colnames(df))
  
  rownames(df)=gsub(" ","",rownames(df))
  rownames(df)=gsub("%","",rownames(df))
  rownames(df)=gsub("\\(","",rownames(df))
  rownames(df)=gsub(")","",rownames(df))
  rownames(df)=gsub("'","",rownames(df))
  rownames(df)=gsub("-","",rownames(df))
  
  # png(paste0("gene_pairwise_correlation/",gene1,"/",gene2,"_VS_",gene1,".png"))
  title=paste0("bicor = ",signif(bicor_coe,3),"\n padj = ",signif(padj,3))
  
  #  plot(x, y, main=title,
  #    xlab=gene1, ylab=gene2, pch=19,cex=0.5)
  #  abline(lm(y~x), col="red")
  #  dev.off()
  ggplot(df, aes_string(x=gene1, y=gene2)) + 
    geom_point(shape=18, color="blue") + #+
    geom_smooth(method=lm,  linetype="dashed",
                color="darkred", fill="blue")+ ggtitle(title)
  ggsave(paste0("gene_pairwise_correlation/",gene2,"_VS_",gene1,".png"))  
}

apply(meta,1,function(x) scatter(x[1], x[2], as.data.frame(datExpr), meta))