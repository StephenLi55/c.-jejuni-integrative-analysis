###link co-expressed genes to DAVID , Kegg brite and String
library(dplyr)
library(readr)
library(readxl)
setwd("~/Desktop/phd/programming/Integrated_analysis/WGCNA/Databases/")
#read DAVID, Kegg brite file and string file
DAVID=as.data.frame(read_delim("Cj_DAVID_list.txt",delim = "\t"))
Br=as.data.frame(read_excel("KEGG_brite.xlsx"))
String_cluster_protein=as.data.frame(read_delim("192222.clusters.proteins.v11.0.txt",delim = "\t"))
String_cluster=as.data.frame(read_delim("192222.clusters.info.v11.0.txt",delim = "\t"))
String_protein_info=as.data.frame(read_delim("192222.protein.info.v11.0.txt",delim = "\t"))


DAVID_id=unique(DAVID[,1])

String_protein_info_id= gsub("192222.","",String_protein_info[,1])
#use the DAVID file as a start
#DAVID_missing=String_protein_info_id[which(!String_protein_info_id %in% DAVID_id)]
#ensure all gene names match
combined_df_String=data.frame(gene_id=String_protein_info_id, gene_name=String_protein_info[,2],
                       String_annotation=String_protein_info[,4])
combined_df_String_DAVID = combined_df_String %>% mutate(DAVID=DAVID[match(String_protein_info_id,DAVID_id),4:ncol(DAVID)])

Br_concat=as.data.frame(Br %>% group_by(Locus_tag) %>% mutate(Level_C_concat=paste0(`Level-C`,collapse = ","),
                                                              Path_concat=paste0(Path,collapse = ","),
                                                              Level_B_concat=paste0(`Level-B`,collapse = ","),
                                                              Level_A_concat=paste0(`Level-A`,collapse = ",")))
Br_concat_unique= Br_concat %>% distinct(Locus_tag, .keep_all = TRUE)
combined_df_String_DAVID_Br = combined_df_String_DAVID %>% mutate(Br=Br_concat_unique[match(String_protein_info_id,Br_concat_unique[,1]),(ncol(Br_concat_unique)-3):ncol(Br_concat_unique)])


#now include string cluster info 

String_cluster_protein_info= String_cluster_protein %>% mutate(gene_id= gsub("192222.","",String_cluster_protein_info$protein_id) ,
                                                               best_described_by=String_cluster$best_described_by[match(String_cluster_protein$cluster_id,String_cluster$cluster_id)])
String_cluster_concat=as.data.frame(String_cluster_protein_info %>% group_by(gene_id) %>% mutate(cluster_id_concat= paste0(cluster_id,collapse = ";"),
                                                                                                 best_described_by_concat=paste0(best_described_by,collapse = ";")))
String_cluster_concat_unique= String_cluster_concat %>% distinct(gene_id, .keep_all = TRUE)
combined_df_String_DAVID_Br_cluster = as.data.frame(combined_df_String_DAVID_Br %>% mutate(String_cluster=String_cluster_concat_unique[match(String_protein_info_id,String_cluster_concat_unique[,4]),(ncol(String_cluster_concat_unique)-1):ncol(String_cluster_concat_unique)]))

#now identify genes modules

#read WGCNA metatable
meta=as.data.frame(read_delim("../co_expression_metatable_sRNA.tsv",delim = "\t"))
temp_df=data.frame(genes=c(meta$fromNode,meta$toNode), modules=c(meta$module,meta$module))
temp_df_unique= temp_df %>% distinct(genes, .keep_all = T)
write.csv(temp_df_unique,"../all_nodes_modules.csv")

gene_annotation=read.csv("../../gene_annotation.csv")
temp_df_unique$gene_id=gene_annotation[match(temp_df_unique$genes,gene_annotation[,2]),1]

#cbind(as.character(temp_df_unique$genes),as.character(test))
combined_df_String_DAVID_Br_cluster$WGCNA_module=temp_df_unique$modules[match(combined_df_String_DAVID_Br_cluster$gene_id,temp_df_unique$gene_id)]
#temp_df_unique$genes[which(!temp_df_unique$genes %in% combined_df_String_DAVID_Br_cluster$gene_name)]

#all_genes=c(meta$fromNode,meta$toNode)
#all_modules=c(meta$module,meta$module)
#find out which modules each mRNA belongs to 




write.table(combined_df_String_DAVID_Br_cluster,"combined_databases.tsv",sep="\t")
#String_cluster_protein_info_id=
combined_databases=as.data.frame(read_delim("combined_databases.tsv",delim = "\t"))


