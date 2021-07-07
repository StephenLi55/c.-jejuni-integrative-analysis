### for each cytoscape input file, add columns about ABC transporter genes
setwd("~/Desktop/phd/programming/Integrated_analysis/WGCNA/")
library(readr)
library(dplyr)
#read the cytoscape input file

lnames = load(file = "WGCNA.RData");
#The variable lnames contains the names of loaded variables.
lnames
#all_mods
membrane_protein=read.csv("~/Desktop/phd/programming/Integrated_analysis/Jenna_pathway/ABC_transporters_all_databases.csv")
#ABC=membrane_protein %>% filter(system == "ABC transporter")
#ABC_system=as.character(ABC$pathway)
#ABC_system=as.character(ABC$pathway)
ABC=membrane_protein %>% distinct(Gene_id, .keep_all = TRUE)
ABC_system=as.character(ABC$KEGG_brite_Substrate)
ABC_system[which(is.na(ABC_system))]="unknwon"
#ABC_gene=as.character(ABC$Gene)
ABC_gene=as.character(ABC$Gene_id)
gene_annotation=read.csv("../gene_annotation.csv")
ABC_gene_name=gene_annotation$gene_name[match(ABC_gene,gene_annotation$gene_id)]

for (mods in all_mods){
  file_name=paste0("CytoscapeInput-edges-",mods,"_adjacency.txt")
  cyto_input=as.data.frame(read_delim(file_name,delim = "\t"))
  from_node=as.character(cyto_input$fromNode)
  to_node=as.character(cyto_input$toNode)
  #read the membrane protein file
  
  #match membrane protein name with cytoscape gene name
  from_node_output=rep(NA,length(from_node))
  to_node_output=rep(NA,length(to_node))
  #create a NA vector with length of the row number of the cytoscape file
  from_node_output[which(!is.na(match(from_node,ABC_gene_name)))]=as.character(ABC_system[match(from_node[which(!is.na(match(from_node,ABC_gene_name)))],ABC_gene_name)])
  from_node_output[which(is.na(match(from_node,ABC_gene_name)))]="null"
  
  to_node_output[which(!is.na(match(to_node,ABC_gene_name)))]=as.character(ABC_system[match(to_node[which(!is.na(match(to_node,ABC_gene_name)))],ABC_gene_name)])
  to_node_output[which(is.na(match(to_node,ABC_gene_name)))]="null"
  #add new vector 
  cyto_input_ABC=cyto_input %>% mutate(from_node_ABC=from_node_output, to_node_ABC=to_node_output)
  #View(cyto_input_ABC %>% filter(from_node_ABC != "null"))
 # View(cyto_input_ABC %>% filter(to_node_ABC != "null"))
  #write new file
  output_file=gsub(".txt","_ABC.txt",file_name)
  write.table(cyto_input_ABC,output_file,sep="\t",row.names = F , quote = F)
}
save(ABC_system,ABC_gene,ABC_gene_name,all_mods,file="ABC_transporters_info.RData")
#===============================================
setwd("~/Desktop/phd/programming/Integrated_analysis/WGCNA/")
library(readr)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(reshape2)
###identify condition and module of interest
#check which module does each ABC transporters belong to 
lnames = load(file = "ABC_transporters_info.RData");
#The variable lnames contains the names of loaded variables.
ABC_gene

#check the edge it forms with sRNA, and their respectively weight
#go through all nodes file
for (mods in all_mods){
  file_name=paste0("CytoscapeInput-nodes-",mods,".txt")
  cyto_input=as.data.frame(read_delim(file_name,delim = "\t"))
#  cyto_input=cbind(cyto_input,mods)
  if (mods == all_mods[1]){
    merged_cyto_input=cyto_input
  }
  else{
    merged_cyto_input=rbind(merged_cyto_input,cyto_input)
  }
}
#put them all in one data frame. Add one colour of their repectively modules
#filter out rows that have those ABC genes
colnames(merged_cyto_input)[3]="mods"
merged_cyto_input_filtered=merged_cyto_input %>% filter(nodeName %in% ABC_gene_name) %>% dplyr::select(nodeName,mods)
dat=as.data.frame(melt(table(merged_cyto_input_filtered)))
#dat$value=factor(dat$value)]
#cols=c("Absent"="red","Present"="orange")
ggplot(dat , aes(x = mods , y = nodeName ,fill=factor(value)))+
         geom_tile()+labs(x="modules", y="ABC transporter genes")+
  scale_fill_discrete(name="",labels=c("Absent","Present"))+ theme(axis.text.x = element_text(angle =45, vjust = 0.8, hjust=1))
ggsave("ABC_genes_modules_heatmaps.png",height = 10,width=7)
  
#visualise in heatmap and barplot. 
#also check DESeq2 output for each ABC transporter for each comparison


#===============================================
#include DESeq analysis
setwd("~/Desktop/phd/programming/Integrated_analysis/WGCNA/")

gene_name_file=read.csv("~/Desktop/phd/programming/Integrated_analysis/gene_annotation.csv")
system("rm -r Cytoscape_DESeq/")
system("mkdir Cytoscape_DESeq")
all_pairwise_input=list.files("~/Desktop/phd/programming/Integrated_analysis/Stephen_DESeq/results_pairwise_combination/")
for (input in all_pairwise_input){
  cond=gsub("control.csv","",input)
  
  system(paste0("mkdir Cytoscape_DESeq/",cond,"/"))
  DESeq_file=read.csv(paste0("~/Desktop/phd/programming/Integrated_analysis/Stephen_DESeq/results_pairwise_combination/",input))
  
  #DESeq_file[,1]
    #DESeq_extracted=DESeq_file[match(gene_name_file[,1],DESeq_file[,1]),]
  
  #DESeq_extracted[,1]=gene_name_file[,2]
  #View(DESeq_extracted)
  #which(is.na(DESeq_extracted[,2])) #check any mismatch
  #pick one pairwise comparison (L37_vs_S37control)
  #include log2FoldChange and padj for each node
  #save new file
  for (mods in all_mods){
    file_name=paste0("CytoscapeInput-edges-",mods,"_adjacency_ABC.txt")
    cyto_input=as.data.frame(read_delim(file_name,delim = "\t"))
    from_node=as.character(cyto_input$fromNode)
    to_node=as.character(cyto_input$toNode)
    #read the membrane protein file
    from_node_output=DESeq_file[match(from_node,DESeq_file[,1]),c(3,7)]
    colnames(from_node_output)=c("from_node_log2FoldChange","from_node_padj")
    to_node_output=DESeq_file[match(to_node,DESeq_file[,1]),c(3,7)]
    colnames(to_node_output)=c("to_node_log2FoldChange","to_node_padj")
    cyto_input_DE=cbind(cyto_input,from_node_output,to_node_output)
    #View(cyto_input_ABC %>% filter(from_node_ABC != "null"))
    # View(cyto_input_ABC %>% filter(to_node_ABC != "null"))
    #write new file
    
    output_path=paste0("Cytoscape_DESeq/",cond,"/")
    output_file=paste0(output_path,gsub(".txt",paste0("_",cond,".txt"),file_name))
    write.table(cyto_input_DE,output_file,sep="\t",row.names = F , quote = F)
    system(paste0("gzip ",output_file))
    system(paste0("rm ",output_file))
  }
}

#input=paste0("~/Desktop/phd/programming/Integrated_analysis/Stephen_DESeq/results_pairwise_combination/",cond,"control.csv")


#===============================================
#include GO term and KEGG info

