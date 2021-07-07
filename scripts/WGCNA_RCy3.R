#run RCy3 to automatically generate visualisation for all mRNA-sRNA interactions
#run RCy3
setwd("~/Desktop/phd/programming/Integrated_analysis/WGCNA/")
library(RCy3)
library(readr)
library(dplyr)
lnames = load(file = "ABC_transporters_info.RData");
#The variable lnames contains the names of loaded variables.
ABC_gene

#read metatable in DESeq comparison of interest

#filter out those that are not DEG
#filter out those with insignificant padj
#check whether the pairs are overlaping. Filter out those that are overlapping
Check_overlapping=function(row){
  # print(length(row))
  output="not overlapping"
  partner1_start=as.numeric(row[3])
  partner1_end=as.numeric(row[4])
  partner1_strand=row[5]
  partner2_start=as.numeric(row[8])
  partner2_end=as.numeric(row[9])
  partner2_strand=row[10]
  if (min(partner1_end,partner2_end)+500 >= max(partner1_start,partner2_start) & partner1_strand==partner2_strand){
    output="overlapping"
  }
  return(output)
}
meta=read_delim("co_expression_metatable_sRNA.tsv",delim="\t")
meta_filtered=meta %>% filter(correlation_padj<=0.05)
cytoscapePing()
nodes <- data.frame(#source=meta_DEG$fromNode,
  #target=meta_DEG$toNode,
  id=unique(c(meta_filtered$fromNode,meta_filtered$toNode)))
edges <- data.frame(source=meta_filtered$fromNode,
                    target=meta_filtered$toNode,
                    weight=meta_filtered$weight) # numeric

createNetworkFromDataFrames(nodes=nodes,edges=edges,title = "network",collection = "all_sRNA_mRNA")#,edges)
all_nodes=unique(c(meta_filtered$fromNode,meta_filtered$toNode))
sRNA_index=grep("Li",(all_nodes))
colour=rep("#5577FF",length(all_nodes))
colour[sRNA_index]="#FF7755"
setNodeColorMapping(table.column='id',table.column.values=all_nodes,color=colour,mapping.type="d") #table.column.values=sRNA_names,colour=colour) #c('#5577FF','#FFFFFF','#FF7755'))
# Sys.sleep(3)
system("rm all_sRNA_mRNA.cys")
system("rm all_sRNA_mRNA.png")
saveSession(filename ="all_sRNA_mRNA.cys" )
exportImage(filename ="all_sRNA_mRNA.png","PNG")
deleteAllNetworks() 

# Import into Cytoscape
#colour nodes using 
folder="Cytoscape_DESeq"
system("rm Cytoscape_DESeq/*/*.png")
system("rm Cytoscape_DESeq/*/*.cys")
directory=list.dirs(folder)
deleteAllNetworks() #remove all pre-existing network
for(d in directory[1:length(directory)]){
  if (d == "Cytoscape_DESeq"){
    next()
  }
  if (grepl("Cytoscape_DESeq/LFC=2/",d)){
    next()
  }
  
  input=paste0(d,"/co_expression_metatable_sRNA.tsv.gz")
  if (!file.exists(input)){
    next() #skip those from subdirectory
  }
 # input="Cytoscape_DESeq/iron_lim_M_vs_iron_rep_M/co_expression_metatable_sRNA.tsv.gz"
  system(paste0("gunzip ",input))
  input_gunzip=gsub(".gz","",input)
  meta=read_delim(input_gunzip,delim="\t")
  #input_gunzip_DEG=
  system(paste0("gzip ",input_gunzip))
  #folder="results_pairwise_combination"
  #files=list.files(folder)
 # meta_DEG=meta %>% filter(correlation_padj<=0.05 ,abs(from_node_log2FoldChange)>2 & abs(to_node_log2FoldChange)>2 & from_node_padj<=0.05 & to_node_padj<=0.05)
  meta_DEG=meta %>% filter(correlation_padj<=0.05 ,from_node_padj<=0.05 & to_node_padj<=0.05)
  if (dim(meta_DEG)[1]==0){
    next()
  }
  cytoscapePing()
  unique_genes=unique(c(meta_DEG$fromNode,meta_DEG$toNode))
  unique_index=which(!duplicated(c(meta_DEG$fromNode,meta_DEG$toNode)))
  unique_LFC=c(meta_DEG$from_node_log2FoldChange,meta_DEG$to_node_log2FoldChange)[unique_index]
  unique_genes_ABC=rep("no",length(unique_genes))
  unique_genes_ABC[which(unique_genes %in% ABC_gene_name)]="yes"
  nodes <- data.frame(#source=meta_DEG$fromNode,
    #target=meta_DEG$toNode,
    id=unique_genes,
    LFC=unique_LFC,
    LGC_sign=sign(unique_LFC),
    ABC=unique_genes_ABC)
  edges <- data.frame(source=meta_DEG$fromNode,
                      target=meta_DEG$toNode,
                      weight=meta_DEG$weight,
                      bicor=meta_DEG$correlation,
                      bicor_sign=sign(meta_DEG$correlation),
                      bicor_padj=meta_DEG$correlation_padj) # numeric
  
  createNetworkFromDataFrames(nodes=nodes,edges=edges,title = "network",collection = gsub("Cytoscape_DESeq/","",d))#,edges)
  all_nodes=unique(c(meta_DEG$fromNode,meta_DEG$toNode))
  sRNA_index=grep("Li",(all_nodes))
  colour=rep("#5577FF",length(all_nodes))
  colour[sRNA_index]="#FF7755"
  setNodeColorMapping(table.column='id',table.column.values=all_nodes,color=colour,mapping.type="d") #table.column.values=sRNA_names,colour=colour) #c('#5577FF','#FFFFFF','#FF7755'))
 # Sys.sleep(3)
  setEdgeLineStyleMapping('bicor_sign',c(1,-1),c('SOLID','LONG_DASH'))
  setNodeBorderColorMapping('LGC_sign', c(1,-1), c('#FF0000','#0000A0'),mapping.type = "d")
  setNodeBorderWidthDefault(10)
  lockNodeDimensions(TRUE)
  #setNodeShapeDefault('ELLIPSE')
  setNodeShapeMapping('ABC',c('no','yes'),c('ELLIPSE','Diamond'))
  setNodeSizeMapping('LFC', c(min(abs(unique_LFC)),max(abs(unique_LFC))), c(50,70))
  setNodeFontSizeDefault(4)
  saveSession(filename =paste0(d,"/","sRNA-mRNA.cys") )
  exportImage(filename =paste0(d,"/","sRNA-mRNA.png"),"PNG")
  deleteAllNetworks() 
  
  
}

all_meta=matrix(,nrow=0, ncol=18)
for(d in directory[!directory %in% "Cytoscape_DESeq"]){
  if (d == "Cytoscape_DESeq"){
    next()
  }
  input=paste0(d,"/co_expression_metatable_sRNA.tsv.gz")
  if (!file.exists(input)){
    next() #skip those from subdirectory
  }
  # input="Cytoscape_DESeq/iron_lim_M_vs_iron_rep_M/co_expression_metatable_sRNA.tsv.gz"
  system(paste0("gunzip ",input))
  input_gunzip=gsub(".gz","",input)
  meta=read_delim(input_gunzip,delim="\t")
  #input_gunzip_DEG=
  system(paste0("gzip ",input_gunzip))
  #folder="results_pairwise_combination"
  #files=list.files(folder)
  #meta_DEG=meta %>% filter(correlation_padj<=0.05 ,abs(from_node_log2FoldChange)>2 & abs(to_node_log2FoldChange)>2 & from_node_padj<=0.05 & to_node_padj<=0.05)
  meta_DEG=meta %>% filter(correlation_padj<=0.05 ,abs(from_node_log2FoldChange)>0 & abs(to_node_log2FoldChange)>0 & from_node_padj<=0.05 & to_node_padj<=0.05)
  if (dim(meta_DEG)[1]==0){
    next()
  }
  meta_DEG_cbind=cbind(meta_DEG,gsub("Cytoscape_DESeq/","",d))
  all_meta=rbind(all_meta,meta_DEG_cbind)
  
}
colnames(all_meta)[ncol(all_meta)]="sample"
write.csv(all_meta,"Cytoscape_DESeq/all_sRNA_mRNA_DEG.csv")


#Label edge based on correlation sign. And their width by bicor


#=================================
#create Cytoscape image for sRNA of interest
setwd("~/Desktop/phd/programming/Integrated_analysis/WGCNA/")
library(RCy3)
library(readr)
library(dplyr)

sRNA_of_interest="Li_ANNOgesic_110"
folder="Cytoscape_DESeq"
cytoscapePing()
#system("rm Cytoscape_DESeq/*/*.png")
#system("rm Cytoscape_DESeq/*/*.cys")
directory=list.dirs(folder)
deleteAllNetworks() #remove all pre-existing network
for(d in directory[1:length(directory)]){
  if (d == "Cytoscape_DESeq"){
    next()
  }
  if (grepl("Cytoscape_DESeq/LFC=2/",d)){
    next()
  }
  
  input=paste0(d,"/co_expression_metatable_sRNA.tsv.gz")
  if (!file.exists(input)){
    next() #skip those from subdirectory
  }
  # input="Cytoscape_DESeq/iron_lim_M_vs_iron_rep_M/co_expression_metatable_sRNA.tsv.gz"
  system(paste0("gunzip ",input))
  input_gunzip=gsub(".gz","",input)
  meta=read_delim(input_gunzip,delim="\t")
  #input_gunzip_DEG=
  system(paste0("gzip ",input_gunzip))
  #folder="results_pairwise_combination"
  #files=list.files(folder)
  # meta_DEG=meta %>% filter(correlation_padj<=0.05 ,abs(from_node_log2FoldChange)>2 & abs(to_node_log2FoldChange)>2 & from_node_padj<=0.05 & to_node_padj<=0.05)
  meta_DEG=meta %>% filter(correlation_padj<=0.05 ,from_node_padj<=0.05 & to_node_padj<=0.05) %>%
    filter(fromNode == sRNA_of_interest | toNode==sRNA_of_interest)
  if (dim(meta_DEG)[1]==0){
    next()
  }
  cytoscapePing()
  unique_genes=unique(c(meta_DEG$fromNode,meta_DEG$toNode))
  unique_index=which(!duplicated(c(meta_DEG$fromNode,meta_DEG$toNode)))
  unique_LFC=c(meta_DEG$from_node_log2FoldChange,meta_DEG$to_node_log2FoldChange)[unique_index]
  unique_genes_ABC=rep("no",length(unique_genes))
  unique_genes_ABC[which(unique_genes %in% ABC_gene_name)]="yes"
  nodes <- data.frame(#source=meta_DEG$fromNode,
    #target=meta_DEG$toNode,
    id=unique_genes,
    LFC=unique_LFC,
    LGC_sign=sign(unique_LFC),
    ABC=unique_genes_ABC)
  edges <- data.frame(source=meta_DEG$fromNode,
                      target=meta_DEG$toNode,
                      weight=meta_DEG$weight,
                      bicor=meta_DEG$correlation,
                      bicor_sign=sign(meta_DEG$correlation),
                      bicor_padj=meta_DEG$correlation_padj) # numeric
  
  createNetworkFromDataFrames(nodes=nodes,edges=edges,title = "network",collection = gsub("Cytoscape_DESeq/","",d))#,edges)
  all_nodes=unique(c(meta_DEG$fromNode,meta_DEG$toNode))
  sRNA_index=grep("Li",(all_nodes))
  colour=rep("#5577FF",length(all_nodes))
  colour[sRNA_index]="#FF7755"
  setNodeColorMapping(table.column='id',table.column.values=all_nodes,color=colour,mapping.type="d") #table.column.values=sRNA_names,colour=colour) #c('#5577FF','#FFFFFF','#FF7755'))
  # Sys.sleep(3)
  setEdgeLineStyleMapping('bicor_sign',c(1,-1),c('SOLID','LONG_DASH'))
  setNodeBorderColorMapping('LGC_sign', c(1,-1), c('#FF0000','#0000A0'),mapping.type = "d")
  setNodeBorderWidthDefault(10)
  lockNodeDimensions(TRUE)
  #setNodeShapeDefault('ELLIPSE')
  setNodeShapeMapping('ABC',c('no','yes'),c('ELLIPSE','Diamond'))
  setNodeSizeMapping('LFC', c(min(abs(unique_LFC)),max(abs(unique_LFC))), c(50,70))
  setNodeFontSizeDefault(4)
  saveSession(filename =paste0(d,"/","sRNA-mRNA_",sRNA_of_interest,".cys") )
  exportImage(filename =paste0(d,"/","sRNA-mRNA_",sRNA_of_interest,".png"),"PNG")
  deleteAllNetworks() 
  
  
}

