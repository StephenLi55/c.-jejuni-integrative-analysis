### integrate all data. Highlight all sRNA-target interactions
library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(RCy3)
setwd("~/Desktop/phd/programming/Integrated_analysis/WGCNA/")
#start with co-expression metatable with IntaRNA results
co_inta_meta=read_delim("co_expression_metatable_sRNA_Inta.tsv",delim="\t")
#combine with co-expression table from DESeq
DE_meta=read.csv("Cytoscape_DESeq/all_sRNA_mRNA_DEG.csv")
selected_cond_file=read.csv("../selected_pairwise_comparisons.csv")
selected_cond=as.character(selected_cond_file[,1])

#new_df=as.data.frame(matrix(NA,ncol=length(selected_cond),nrow=nrow(co_inta_meta)))
#colnames(new_df)=selected_cond
#for loop to go through co_inta_meta
#unique(DE_meta$sample)
co_pattern=sign(co_inta_meta$correlation)
DEscore=rep(NA,nrow(co_inta_meta))
DE_up_up=rep(NA,nrow(co_inta_meta))
DE_up_down=rep(NA,nrow(co_inta_meta))
DE_down_up=rep(NA,nrow(co_inta_meta))
DE_down_down=rep(NA,nrow(co_inta_meta))
DE_co_match=rep(NA,nrow(co_inta_meta))
DE_co_mismatch=rep(NA,nrow(co_inta_meta))
for (i in 1:nrow(co_inta_meta)){
  gene1=co_inta_meta$fromNode[i]
  gene2=co_inta_meta$toNode[i]
  #add new column for DEup. Paste by ","
  DE_meta_filtered=DE_meta %>% filter(fromNode == gene1 & toNode == gene2) %>%
    filter(sample %in% selected_cond)
  DE_up_up[i]= as.character(DE_meta_filtered %>% filter(from_node_padj <= 0.05 & from_node_log2FoldChange>0 &
                                             to_node_padj <= 0.05 & to_node_log2FoldChange>0) %>%
     summarise(paste(sample,collapse=",")))
  DE_up_down[i]= as.character(DE_meta_filtered %>% filter(from_node_padj <= 0.05 & from_node_log2FoldChange>0 &
                                                          to_node_padj <= 0.05 & to_node_log2FoldChange<0) %>%
                              summarise(paste(sample,collapse=",")))
  DE_down_up[i]= as.character(DE_meta_filtered %>% filter(from_node_padj <= 0.05 & from_node_log2FoldChange<0 &
                                                            to_node_padj <= 0.05 & to_node_log2FoldChange>0) %>%
                                summarise(paste(sample,collapse=",")))
  DE_down_down[i]= as.character(DE_meta_filtered %>% filter(from_node_padj <= 0.05 & from_node_log2FoldChange<0 &
                                                            to_node_padj <= 0.05 & to_node_log2FoldChange<0) %>%
                                summarise(paste(sample,collapse=",")))
  if (co_pattern[i]==1){
    DE_co_match[i]=sum(str_count(c(DE_up_up[i],DE_down_down[i]),"vs"))
    DE_co_mismatch[i]=sum(str_count(c(DE_up_down[i],DE_down_up[i]),"vs"))
  }else{
    DE_co_match[i]=sum(str_count(c(DE_up_down[i],DE_down_up[i]),"vs"))
    DE_co_mismatch[i]=sum(str_count(c(DE_up_up[i],DE_down_down[i]),"vs"))
    #list all conditions of differential upregulation
    #list all conditions of differential downregulation
    
  }
  DEscore[i]=DE_co_match[i] - DE_co_mismatch[i]
}
co_inta_meta_DE=co_inta_meta %>% mutate(DE_up_up,DE_up_down,DE_down_up,DE_down_down,DE_co_match,DE_co_mismatch,DEscore)
#upregulate - downregulate. Need to be positive
#

#find match from crosslinking data
crosslinking_meta=read_delim("~/Desktop/phd/programming/RNA_crosslinking_analysis/Cj_crosslinking/Combined_analysis/metatable_sRNA_mRNA_filtered.tsv",delim="\t")
crosslinking_meta_filtered=crosslinking_meta %>% filter(is.na(OE_VS_AMT_ratio) & is.na(OE_VS_UV_ratio))

match_index=c(which(!is.na(match(paste0(co_inta_meta_DE$fromNode,co_inta_meta_DE$toNode),paste0(crosslinking_meta_filtered$partner1_vec,crosslinking_meta_filtered$partner2_vec)))),
which(!is.na(match(paste0(co_inta_meta_DE$fromNode,co_inta_meta_DE$toNode),paste0(crosslinking_meta_filtered$partner2_vec,crosslinking_meta_filtered$partner1_vec)))))

crosslink_score=rep(NA,nrow(co_inta_meta))
crosslink_score[match_index]=1
co_inta_meta_DE_crosslinking=co_inta_meta_DE %>% mutate(crosslink_score)
#find match in validated sRNA

#also include TSS_UTR_length, module_KEGG_enrichment and validated sRNA name

write.table(co_inta_meta_DE_crosslinking,"final_integration/co_expression_metatable_sRNA_Inta_DE_crosslink.tsv",sep="\t",quote=F,row.names = F )

co_inta_meta_DE_crosslinking=read_delim("final_integration/co_expression_metatable_sRNA_Inta_DE_crosslink.tsv",delim = "\t")
#also include TSS_UTR_length, module_KEGG_enrichment and validated sRNA name
UTR_file=read.csv("../gene_annotation_UTR.csv")[,-1]
fromNode_UTR_length= UTR_file$estimated_5end_length[match(co_inta_meta_DE_crosslinking$fromNode,UTR_file$gene_name)]
toNode_UTR_length=UTR_file$estimated_5end_length[match(co_inta_meta_DE_crosslinking$toNode,UTR_file$gene_name)]

enrichment_file=read.csv("Enrichment_analysis/metatable.csv")
all_mods=co_inta_meta_DE_crosslinking$module

all_mods_unique_KEGG_term=enrichment_file %>% filter( catogory=="KEGG" & pvalue_fdr<=0.05) %>% group_by(sample) %>% summarise(keggs=paste(term_description,collapse=","))

# enrichment_file %>% filter(sample %in% all_mods)
module_enriched_kegg=all_mods_unique_KEGG_term$keggs[match(co_inta_meta_DE_crosslinking$module,all_mods_unique_KEGG_term$sample)]

co_inta_meta_DE_crosslinking_UTR_KEGG = co_inta_meta_DE_crosslinking %>% mutate(fromNode_UTR_length,toNode_UTR_length,module_enriched_kegg)


write.table(co_inta_meta_DE_crosslinking_UTR_KEGG,"final_integration/co_expression_metatable_sRNA_Inta_DE_crosslink_UTR_KEGG.tsv",sep="\t",quote=F,row.names = F )

#create a filtered table
co_inta_meta_DE_crosslinking_highlight = co_inta_meta_DE_crosslinking_UTR_KEGG %>% filter(correlation_padj <=0.05 & DEscore >=2 &
                                                                                   (p.value<=0.05 | crosslink_score==1)) %>%
  filter(!grepl("sRNA_Li_",fromNode) &!grepl("sRNA_Li_",toNode)) 
target_list=apply(co_inta_meta_DE_crosslinking_highlight[,1:2],1,
                  function(x) x[which(!grepl("ANNOgesic",x))])
databases=read_delim("Databases/combined_databases.tsv",delim = "\t")
target_list_STRING=databases$String_annotation[match(target_list,databases$gene_name)]
target_list_KEGG_PATHWAY=databases$DAVID.KEGG_PATHWAY[match(target_list,databases$gene_name)]
target_list_KEGG_Br_C=databases$Br.Level_C_concat[match(target_list,databases$gene_name)]
target_list_KEGG_Br_B=databases$Br.Level_B_concat[match(target_list,databases$gene_name)]
target_list_KEGG_Br_A=databases$Br.Level_A_concat[match(target_list,databases$gene_name)]
target_list_COG=databases$DAVID.COG_ONTOLOGY[match(target_list,databases$gene_name)]
target_list_GO_BP=databases$DAVID.GOTERM_BP_DIRECT[match(target_list,databases$gene_name)]
target_list_GO_CC=databases$DAVID.GOTERM_CC_DIRECT[match(target_list,databases$gene_name)]
target_list_GO_MF=databases$DAVID.GOTERM_MF_DIRECT[match(target_list,databases$gene_name)]
co_inta_meta_DE_crosslinking_highlight = co_inta_meta_DE_crosslinking_highlight %>% 
  mutate(target_list_STRING,target_list_KEGG_PATHWAY,target_list_KEGG_Br_C,
         target_list_KEGG_Br_B,target_list_KEGG_Br_A,
         target_list_COG,target_list_GO_BP,target_list_GO_CC,
         target_list_GO_MF)
                                                                                           
write.table(co_inta_meta_DE_crosslinking_highlight ,"final_integration/co_expression_metatable_sRNA_Inta_DE_crosslink_UTR_KEGG_highlight.tsv",sep="\t",quote=F,row.names = F )
#make a table with only validated sRNA
validated_file=read.csv("~/Desktop/phd/programming/sRNA_prediction/validated_vs_final_prediction.csv")
validated_prediction=as.character(validated_file$predicted_sRNA)
co_inta_meta_DE_crosslinking_highlight_validated=co_inta_meta_DE_crosslinking_highlight %>% filter(fromNode %in% validated_prediction | toNode %in% validated_prediction)
highighted_predicted_sRNA=apply(co_inta_meta_DE_crosslinking_highlight_validated[,1:2],1,
                                function(x) x[which(x %in% validated_prediction )])

validated_prediction_matched=validated_file$validatedsRNA[match(highighted_predicted_sRNA,validated_prediction)]
co_inta_meta_DE_crosslinking_highlight_validated_matched=co_inta_meta_DE_crosslinking_highlight_validated%>%
  mutate(validated_prediction_matched)
write.table(co_inta_meta_DE_crosslinking_highlight_validated_matched,"final_integration/co_expression_metatable_sRNA_Inta_DE_crosslink_UTR_KEGG_highlight_validated.tsv",sep="\t",quote=F,row.names = F )
save(co_inta_meta_DE_crosslinking_highlight_validated_matched,file="final_integration/co_expression_metatable_sRNA_Inta_DE_crosslink_UTR_KEGG_highlight_validated.RData" )

# 
# #make cytoscape image
# nodes <- data.frame(#source=meta_DEG$fromNode,
#   #target=meta_DEG$toNode,
#   id=unique(c(co_inta_meta_DE_crosslinking_highlight$fromNode,co_inta_meta_DE_crosslinking_highlight$toNode))
#  )
# edges <- data.frame(source=co_inta_meta_DE_crosslinking_highlight$fromNode,
#                     target=co_inta_meta_DE_crosslinking_highlight$toNode,
#                   #  weight=meta_DEG$weight,
#                     bicor=co_inta_meta_DE_crosslinking_highlight$correlation,
#                     bicor_sign=sign(co_inta_meta_DE_crosslinking_highlight$correlation),
#                     bicor_padj=co_inta_meta_DE_crosslinking_highlight$correlation_padj) # numeric
# 
# createNetworkFromDataFrames(nodes=nodes,edges=edges,title = "network",collection = "final_integration/sRNA-target_highlight")#,edges)
# all_nodes=unique(c(co_inta_meta_DE_crosslinking_highlight$fromNode,co_inta_meta_DE_crosslinking_highlight$toNode))
# sRNA_index=grep("Li",(all_nodes))
# validated_sRNA_index=which(all_nodes %in% validated_prediction)
# colour=rep("#5577FF",length(all_nodes))
# colour[sRNA_index]="#00FF00"
# colour[validated_sRNA_index]="#FF7755"
# setNodeColorMapping(table.column='id',table.column.values=all_nodes,color=colour,mapping.type="d") #table.column.values=sRNA_names,colour=colour) #c('#5577FF','#FFFFFF','#FF7755'))
# setEdgeLineStyleMapping('bicor_sign',c(1,-1),c('SOLID','LONG_DASH'))
# lockNodeDimensions(TRUE)
# #setNodeShapeDefault('ELLIPSE')
# #setNodeShapeMapping('ABC',c('no','yes'),c('ELLIPSE','Diamond'))
# setNodeShapeDefault('ELLIPSE')
# setNodeSizeDefault(70)
# saveSession(filename ="final_integration/sRNA-target_highlight" )
# exportImage(filename ="final_integration/sRNA-target_highlight")
# deleteAllNetworks() 
# Sys.sleep(3)
for (d in c("",unique(co_inta_meta_DE_crosslinking_highlight$module))){
  #df=co_inta_meta_DE_crosslinking_highlight%>%filter(module ==d)
  if (d==""){
    df=co_inta_meta_DE_crosslinking_highlight
  }else{
    df=co_inta_meta_DE_crosslinking_highlight%>%filter(module ==d)
  }
  nodes <- data.frame(#source=meta_DEG$fromNode,
    #target=meta_DEG$toNode,
    id=unique(c(df$fromNode,df$toNode))
  )
  edges <- data.frame(source=df$fromNode,
                      target=df$toNode,
                      #  weight=meta_DEG$weight,
                      bicor=df$correlation,
                      bicor_sign=sign(df$correlation),
                      bicor_padj=df$correlation_padj) # numeric
  
  createNetworkFromDataFrames(nodes=nodes,edges=edges,title = "network",collection = paste0("final_integration/sRNA-target_highlight_",d))#,edges)
  all_nodes=unique(c(co_inta_meta_DE_crosslinking_highlight$fromNode,co_inta_meta_DE_crosslinking_highlight$toNode))
  sRNA_index=grep("Li",(all_nodes))
  colour=rep("#5577FF",length(all_nodes))
  colour[sRNA_index]="#00FF00"
  colour[validated_sRNA_index]="#FF7755"
  
  setNodeColorMapping(table.column='id',table.column.values=all_nodes,color=colour,mapping.type="d") #table.column.values=sRNA_names,colour=colour) #c('#5577FF','#FFFFFF','#FF7755'))
  setEdgeLineStyleMapping('bicor_sign',c(1,-1),c('SOLID','LONG_DASH'))
  lockNodeDimensions(TRUE)
  #setNodeShapeDefault('ELLIPSE')
  #setNodeShapeMapping('ABC',c('no','yes'),c('ELLIPSE','Diamond'))
  setNodeShapeDefault('ELLIPSE')
  setNodeSizeDefault(70)
  saveSession(filename =paste0("final_integration/Cytoscape/sRNA-target_highlight_",d) )
  exportImage(filename =paste0("final_integration/Cytoscape/sRNA-target_highlight_",d))
  deleteAllNetworks() 
  
}

#make correlation plot for all highlighted sRNA-target interactions
lnames = load(file = "WGCNA.RData")
scatter=function(gene1,gene2,bicor_coe,padj,df,metatable){
  # if (!file.exists(paste0("gene_pairwise_correlation/",gene1))){
  # system(paste0("mkdir gene_pairwise_correlation/",gene1))
  # }
  
  #  par(mar=c(4,5,4,5))
  print(gene1)
  print(gene2)
 # filtered_row=metatable %>% filter(fromNode==gene1,toNode==gene2)
  #filtered_row=
 # bicor_coe=filtered_row$correlation
 # padj=filtered_row$correlation_padj
  
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
  title=paste0("bicor = ",signif(bicor_coe,3),"\npadj = ",signif(padj,3))
  
  #  plot(x, y, main=title,
  #    xlab=gene1, ylab=gene2, pch=19,cex=0.5)
  #  abline(lm(y~x), col="red")
  #  dev.off()
  ggplot(df, aes_string(x=gene1, y=gene2)) + 
    geom_point(shape=18, color="blue") + #+
    geom_smooth(method=lm,  linetype="dashed",
                color="darkred", fill="blue")+ ggtitle(title)
  ggsave(paste0("final_integration/gene_pairwise_correlation/",gene2,"_VS_",gene1,".png"))  
}

apply(co_inta_meta_DE_crosslinking_highlight_validated_matched,1,function(x) scatter(x[1], x[2],as.numeric(x[7]),as.numeric(x[9]), as.data.frame(datExpr)))

#include database annotation for target
 load("final_integration/co_expression_metatable_sRNA_Inta_DE_crosslink_UTR_KEGG_highlight_validated.RData")
 df=co_inta_meta_DE_crosslinking_highlight_validated_matched
 target_list=apply(df[,1:2],1,
                                 function(x) x[which(!grepl("ANNOgesic",x))])
 databases=read_delim("Databases/combined_databases.tsv",delim = "\t")
 target_list_STRING=databases$String_annotation[match(target_list,databases$gene_name)]
 target_list_KEGG_PATHWAY=databases$DAVID.KEGG_PATHWAY[match(target_list,databases$gene_name)]
 target_list_KEGG_Br_C=databases$Br.Level_C_concat[match(target_list,databases$gene_name)]
 target_list_KEGG_Br_B=databases$Br.Level_B_concat[match(target_list,databases$gene_name)]
 target_list_KEGG_Br_A=databases$Br.Level_A_concat[match(target_list,databases$gene_name)]
 target_list_COG=databases$DAVID.COG_ONTOLOGY[match(target_list,databases$gene_name)]
 target_list_GO_BP=databases$DAVID.GOTERM_BP_DIRECT[match(target_list,databases$gene_name)]
 target_list_GO_CC=databases$DAVID.GOTERM_CC_DIRECT[match(target_list,databases$gene_name)]
 target_list_GO_MF=databases$DAVID.GOTERM_MF_DIRECT[match(target_list,databases$gene_name)]
 
 # output_line=paste0(sRNA," (",validated_sRNA,") shows positive co-expression with ",mRNA,". RpoD is a sigma 70 factor. They are in module black, which is enriched which KEGG PATHWAY cje03010 (Ribosome). Both were upregulated mostly in late-staionary phase. That's in contrary to observation in E. coli where RpoS is inhibited in stationary phase")