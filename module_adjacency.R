#go through all edge files. Find out all sRNA-mediated interactions
all_input=list.files(pattern = "_adjacency_ABC.txt")
output=matrix(, nrow = 0, ncol = 13)
sRNA_output=matrix(, nrow = 0, ncol = 13)
for (file in all_input){
  print(file)
  module=gsub("_adjacency_ABC.txt","",gsub("CytoscapeInput-edges-","",file))
  input=paste0("./",file)
  df=as.data.frame(read_delim(input,delim = "\t"))
  extracted_row= df %>% filter(grepl("ANNOgesic",fromNode) | grepl("ANNOgesic",toNode) ) #pick out rows with the first column consisting "Li_ANNOgesic"
  output=rbind(output,cbind(df,module))
  sRNA_output=rbind(sRNA_output,cbind(extracted_row,module))
}
output_name= "co_expression_metatable.tsv"
write.table(output,output_name,sep="\t",col.names = TRUE,row.names = FALSE)

sRNA_output_name= "co_expression_metatable_sRNA.tsv"
write.table(sRNA_output,sRNA_output_name,sep="\t",col.names = TRUE,row.names = FALSE)

all_pairwise_input=list.files("~/Desktop/phd/programming/Integrated_analysis/Stephen_DESeq/results_pairwise_combination/")
for (input in all_pairwise_input){
  cond=gsub("control.csv","",input)
  # system(paste0("mkdir Cytoscape_DESeq/",cond,"/"))
  DESeq_file=read.csv(paste0("~/Desktop/phd/programming/Integrated_analysis/Stephen_DESeq/results_pairwise_combination/",input))
  cyto_input=output
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
  output_file=paste0(output_path,output_name)
  write.table(cyto_input_DE,output_file,sep="\t",row.names = F , quote = F)
  system(paste0("rm ",output_file,".gz"))
  system(paste0("gzip ",output_file))
  #system(paste0("rm ",output_file))
  sRNA_cyto_input_DE=cyto_input_DE %>% filter(grepl("ANNOgesic",fromNode) | grepl("ANNOgesic",toNode))
  output_file=paste0(output_path,sRNA_output_name)
  write.table(sRNA_cyto_input_DE,output_file,sep="\t",row.names = F , quote = F)
  system(paste0("rm ",output_file,".gz"))
  system(paste0("gzip ",output_file))
  #  system(paste0("rm ",output_file))
  #DESeq_file[,1]
  #DESeq_extracted=DESeq_file[match(gene_name_file[,1],DESeq_file[,1]),]
  
  #DESeq_extracted[,1]=gene_name_file[,2]
  #View(DESeq_extracted)
  #which(is.na(DESeq_extracted[,2])) #check any mismatch
  #pick one pairwise comparison (L37_vs_S37control)
  #include log2FoldChange and padj for each node
  #save new file
  #  for (mods in all_mods){
  #  file_name=paste0("CytoscapeInput-edges-",mods,"_adjacency_ABC.txt")
  
  # }
}

#input=paste0("~/Desktop/phd/programming/Integrated_analysis/Stephen_DESeq/results_pairwise_combination/",cond,"control.csv")