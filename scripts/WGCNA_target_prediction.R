#integrate target prediction results into my crosslinking model
library(dplyr)
library(plyr)
library(seqinr)
library(readr)
library(GenomicRanges)
library(Biostrings)
# path="~/Desktop/phd/programming/sRNA_prediction/annogesic/output/sRNA_targets/IntaRNA_output/"
# setwd(path)
# sense=read.csv("IntaRNA_output_sense.csv",sep = ";")
# antisense=read.csv("IntaRNA_output_antisense.csv",sep = ";")
# Inta=rbind(sense,antisense)
# Inta=Inta[-which(gsub("\\|.*","",Inta[,4]) %in% sRNA_to_remove),]
# rm(sense)
# rm(antisense)
# sRNA_to_remove=as.character(read.csv("sRNA_to_remove.csv")[,1])
# #mRNA_coordinate=gsub("0_Cj0001_cds-CAL34182.1_","",Inta[,1])
# mRNA_name=as.character(sapply(as.character(Inta[,1]),function(x) substr(x,3,nchar(x))))
# #mRNA_name=gsub("_cds.*","",mRNA_name2)
# mRNA_name_split=strsplit(mRNA_name,split = "_")
# mRNA_name_split_coordinates=sapply(mRNA_name_split,function(x) x[3])
# mRNA_name_split_start=gsub("-.*","",mRNA_name_split_coordinates)
# mRNA_name_split_end=gsub(".*-","",mRNA_name_split_coordinates)
# mRNA_name_split_strand=sapply(mRNA_name_split,function(x) x[4])
# inta_mRNA_df=data.frame(chr="NC_002163.1",start=mRNA_name_split_start,end=mRNA_name_split_end,strand=mRNA_name_split_strand)
# inta_mRNA_df_Grange=makeGRangesFromDataFrame(inta_mRNA_df)
 annotation=as.data.frame(read_delim("~/Desktop/phd/programming/RNA_crosslinking_analysis/Cj_crosslinking/genomeDir/11168_genome_sRNAs_exon_updated_shortened_visualised.gtf", delim="\t",col_names = F))
 annotation_remove_start=annotation %>% filter(grepl("antisense",X9)) %>% select(X4)
 gene_id=gsub("gene_id \"","",gsub("\";.*","",annotation[,9]))
# colnames(annotation)[c(1,4,5,7)]=c("chr","start","end","strand")
# annotation_Grange=makeGRangesFromDataFrame(annotation)
# 
# Inta_mRNA_overlaps=findOverlaps(inta_mRNA_df_Grange,annotation_Grange,type = 'equal')
# Inta_mRNA=gene_id[subjectHits(Inta_mRNA_overlaps)]
# 
# sRNA_name_split=strsplit(as.character(Inta[,4]),split = "\\|")
# sRNA_name_split_start=sapply(sRNA_name_split,function(x) x[3])
# 
# sRNA_name_split_end=sapply(sRNA_name_split,function(x) x[4])
# sRNA_name_split_strand=sapply(sRNA_name_split,function(x) x[5])
# Inta_sRNA_df=data.frame(chr="NC_002163.1",start=sRNA_name_split_start,end=sRNA_name_split_end,
#                         strand=sRNA_name_split_strand)
# Inta_sRNA_df_Grange=makeGRangesFromDataFrame(Inta_sRNA_df)
# Inta_sRNA_overlaps=findOverlaps(Inta_sRNA_df_Grange,annotation_Grange,type = 'equal')
# Inta_sRNA=gene_id[subjectHits(Inta_sRNA_overlaps)]
# 
# Inta_mRNA_sRNA= Inta %>% mutate(mRNA=Inta_mRNA , sRNA=Inta_sRNA)
# #and for the interaction energy, use IntaRNA
# #pick out all predicted sRNA. Pick out all elements with "ANNOgesic"
# #sRNA_list=unique(c(as.character(df$gene1[grep("ANNOgesic",df$gene1)]),as.character(df$gene2[grep("ANNOgesic",df$gene2)])))


#====================================================
#df=read.csv("~/Desktop/phd/programming/sRNA_prediction/sRNA_final_prediction.csv")
#sRNA_list=unique(c(as.character(df$gene1[grep("ANNOgesic",df$gene1)]),as.character(df$gene2[grep("ANNOgesic",df$gene2)])))
#sRNA_list=as.character(df$name)
sRNA_list=gene_id[grep("ANNOgesic",gene_id)]
#defube a function. Use biostring to extract the DNA sequence
sRNA_blast=data.frame()
#run blastSequences. Store output in a list
#finalised_sRNA=read.csv("~/Desktop/phd/programming/sRNA_prediction/sRNA_final_prediction.csv")
genome=readDNAStringSet("~/Desktop/phd/programming/RNA_crosslinking_analysis/Cj_crosslinking/genomeDir/NCTC11168_genome.fasta.txt")
meta=as.data.frame(read_delim("~/Desktop/phd/programming/Integrated_analysis/WGCNA/co_expression_metatable_sRNA.tsv",delim="\t"))
sRNA_mRNA=meta %>% filter((grepl("Li_",fromNode) & !grepl("Li_",toNode)) | 
                                (!grepl("Li_",fromNode) & grepl("Li_",toNode)) )
sRNA=apply(sRNA_mRNA,1,function(x) if (grepl("Li_",x[1])){
  return (x[1])
}else{
  return(x[2])
})
mRNA=apply(sRNA_mRNA,1,function(x) if (grepl("Li_",x[1])){
  return (x[2])
}else{
  return(x[1])
})
sRNA_mRNA_df= data.frame(sRNA=as.character(sRNA), mRNA=as.character(mRNA)) %>% filter(!grepl("rna",mRNA))
write.table(sRNA_mRNA_df,"sRNA_mRNA_df.tsv",quote = F,row.names = F,sep="\t")
sRNA_seq=function(sRNA){
  print(sRNA)
  index=which(gene_id==as.character(sRNA))
  print(index)
  sequence=substring(genome,annotation[index,4],annotation[index,5])
  if (annotation[index,7]=="-"){
    sequence=reverseComplement(DNAString(sequence))
  }
  return(RNAString(DNAString(sequence)))
}
TSS_file=read.csv("~/Desktop/phd/programming/Jenna/categorised_TSS.csv")
TSS_file_interest=TSS_file %>% filter(TSS.type=="Primary" | TSS.type=="Secondary")
TSS_file_interest$Locus.tag_modified=gsub("_([0-9])","",TSS_file_interest$Locus.tag)

TSS_aggregated=as.data.frame(TSS_file_interest) %>% group_by(Locus.tag_modified) %>% 
  dplyr::summarise(TSS.pos= paste(TSS.position,collapse = ","))
  


all_genes=read.csv("~/Desktop/phd/programming/Integrated_analysis/gene_annotation.csv")
all_genes$TSS=TSS_aggregated$TSS.pos[match(all_genes$gene_id , TSS_aggregated$Locus.tag_modified)]
annotation_all_genes=annotation[match(all_genes$gene_name,gene_id),]
#annotation
df=all_genes %>% mutate(start=annotation_all_genes$X4,
                                         end=annotation_all_genes$X5,
                                         strand=annotation_all_genes$X7)
df$estimated_5end=NA
df$estimated_5end_length=NA
for ( i in 1:nrow(df)){
  if (df$strand[i] == "+" & !is.na(df$TSS[i])){
    df$estimated_5end[i]=min(as.numeric(unlist(strsplit(df$TSS[i],","))))
    df$estimated_5end_length[i]=df$start[i]-df$estimated_5end[i]
  }else if(df$strand[i] == "-" & !is.na(df$TSS[i])){
    df$estimated_5end[i]=max(as.numeric(unlist(strsplit(df$TSS[i],","))))
    df$estimated_5end_length[i]=df$estimated_5end[i]-df$end[i]
  }
  else if(df$strand[i] == "+" & is.na(df$TSS[i])){
    df$estimated_5end[i]=df$start[i]
    df$estimated_5end_length[i]=0
  }else{
    df$estimated_5end[i]=df$end[i]
    df$estimated_5end_length[i]=0
  }
}
#all_genes$TSS.type=TSS_file_interest$TSS.type[match(all_genes$gene_id , TSS_file_interest$Locus.tag_modified)]
#all_genes$TSS.identified=TSS_file_interest$TSS_identifier[match(all_genes$gene_id , TSS_file_interest$Locus.tag_modified)]
write.csv(df,"~/Desktop/phd/programming/Integrated_analysis/gene_annotation_UTR.csv")
#improve this function. Use Jenna's file. Find associated primary or secondary TSS using Jenna's file
#don't change 3' UTR for now
mRNA_seq=function(mRNA,df){
 # index=which(gene_id==mRNA)
  print(mRNA)
  index=which(df$gene_name == as.character(mRNA))
  #UTR_length=200 #default length
  #if mRNA is in Jenna's file, find its primary and secondary TSS
    #find the maximum abs distance between mRNA 5' end and TSS
    #cvhange the UTR_length variable accordingly
  
  sequence=substring(genome,df$estimated_5end[index],df$end[index])
  if (df$strand[index]=="-"){
    sequence=substring(genome,df$start[index],df$estimated_5end[index])
    sequence=reverseComplement(DNAString(sequence))
  }
  return(RNAString(DNAString(sequence)))
         
}
Run_Inta=function(sRNA_id,mRNA_id){
  sRNA=sRNA_seq(sRNA_id)
  mRNA=mRNA_seq(mRNA_id,df)
  output=paste0(sRNA_id,"_",gsub(" ","",gsub("\\(","",gsub("\\)","",gsub("\'","",mRNA_id)))),".csv")
  #system(paste0("~/miniconda3/bin/IntaRNA -q ~/Desktop/phd/programming/RNA_crosslinking_analysis/Cj_crosslinking/IntaRNA_output/sRNA.fa -t ~/Desktop/genomeDir/all_gene.fasta --personality=IntaRNAsTar --threads 4 -n 10 --outMode C --out ",output ))
  system(paste0("~/anaconda3/bin/IntaRNA -q ",as.character(sRNA)," -t ",as.character(mRNA)," --threads 4 -n 1 --outMode C --out ",output,"  --personality=IntaRNAsTar" ))
}
Run_Inta_fig=function(sRNA_id,mRNA_id){
  sRNA=sRNA_seq(sRNA_id)
  print(sRNA)
  mRNA=mRNA_seq(mRNA_id,df)
  print(mRNA)
  output=paste0(sRNA_id,"_",gsub(" ","",gsub("\\(","",gsub("\\)","",gsub("\'","",mRNA_id)))),"_fig")
  #system(paste0("~/miniconda3/bin/IntaRNA -q ~/Desktop/phd/programming/RNA_crosslinking_analysis/Cj_crosslinking/IntaRNA_output/sRNA.fa -t ~/Desktop/genomeDir/all_gene.fasta --personality=IntaRNAsTar --threads 4 -n 10 --outMode C --out ",output ))
  system(paste0("~/anaconda3/bin/IntaRNA -q ",as.character(sRNA)," -t ",as.character(mRNA)," --threads 4 -n 1 --outMode N --personality=IntaRNAsTar" ))
}
#apply(sRNA_mRNA_df,1,function(x) Run_Inta(x[1],x[2]))
for (i in 1:nrow(sRNA_mRNA_df)){

  Run_Inta(sRNA_mRNA_df[i,1],sRNA_mRNA_df[i,2])
}
# for (i in 1:nrow(sRNA_mRNA_df)){
#  if (sRNA_mRNA_df[i,1]=="Li_ANNOgesic_110"){
#    Run_Inta_fig(sRNA_mRNA_df[i,1],sRNA_mRNA_df[i,2])
#  }else{
#    next()
#  }
# }
setwd("~/Desktop/phd/programming/Integrated_analysis/WGCNA/IntaRNA/output/")
#dataset=ldply(list.files(pattern = ".csv"),
             #  read.csv,header=T)
#repetition=nrow(dataset)/length(list.files(pattern = ".csv"))
meta=matrix(,nrow=0,ncol=8)
no_interaction_meta=matrix(,nrow=0,ncol=1)
for (file in list.files(pattern = ".csv")){
  input=read.csv(file,sep = ";",header=T)
  if (nrow(input)==0){
    print(file)
    no_interaction_meta=rbind(no_interaction_meta,file)
    next()
  }
  input_cbind=cbind(input,file)
  meta=rbind(meta,input_cbind)
}
#write.table(meta,"metatable.tsv",sep = "\t",col.names = NA)
write.table(no_interaction_meta,"no_interactions.tsv",sep = "\t",row.names = F)
#sRNA_mRNA_df_seq=sRNA_mRNA_df %>% mutate(sRNA_seq=sapply(sRNA_mRNA_df[,1],sRNA_seq),mRNA_seq=sapply(sRNA_mRNA_df[,2],mRNA_seq))
#write.fasta(as.character(unlist(sRNA_mRNA_df_seq[,3])),as.character(sRNA_mRNA_df_seq[,1]),"~/Desktop/phd/programming/Integrated_analysis/WGCNA/IntaRNA/sRNA.fa")
#write.csv(sRNA_mRNA_df_seq,"~/Desktop/phd/programming/Integrated_analysis/WGCNA/IntaRNA/sRNA_mRNA_sequences.csv")
df=read.csv("~/Desktop/phd/programming/Integrated_analysis/gene_annotation_UTR.csv")
#for 5' UTR, use the bp between UTR and the annotated 5'end
binding_catogory=function(target,target_start,target_end,df){
 # annotation=df
  target_id=gsub("Li_ANNOgesic_[0-999]+_","",gsub(".csv","",target))
  print(target_id)
  print(target_start)
  print(target_end)
 # gene_id_gsubbed=gsub(" ","",gsub("\\(","",gsub("\\)","",gsub("\'","",gene_id))))
  gene_id_gsubbed=gsub(" ","",gsub("\\(","",gsub("\\)","",gsub("\'","",df$gene_name))))
  index=which(gene_id_gsubbed==target_id)
  print(index)
  target_strand=as.character(df$strand[index])
  print(target_strand)
  target_df=data.frame(chr="NC002163.1",start=target_start,end=target_end,strand=target_strand)
  target_Grange=makeGRangesFromDataFrame(target_df)
  
  mRNA_start=as.numeric(df$start[index])
  mRNA_end=as.numeric(df$end[index])
  mRNA_length=mRNA_end-mRNA_start
  UTR_length=as.numeric(df$estimated_5end_length[index])
  
  if (target_strand=="+"){
    five_UTR_Grange=makeGRangesFromDataFrame(data.frame(chr="NC002163.1",start=0,end=UTR_length,strand=target_strand))
    three_UTR_Grange=makeGRangesFromDataFrame(data.frame(chr="NC002163.1",start=mRNA_length,end=mRNA_length + 150 ,strand=target_strand))
  }else{
    five_UTR_Grange=makeGRangesFromDataFrame(data.frame(chr="NC002163.1",start=mRNA_length,end=mRNA_length + UTR_length ,strand=target_strand))
    three_UTR_Grange=makeGRangesFromDataFrame(data.frame(chr="NC002163.1",start=0,end=150,strand=target_strand))
  }
  five_UTR_Grange_overlap=findOverlaps(five_UTR_Grange,target_Grange,type = "any")
  three_UTR_Grange_overlap=findOverlaps(three_UTR_Grange,target_Grange,type = "any")
  if (length(five_UTR_Grange_overlap)!=0){
    output="5-UTR"
  }else if(length(three_UTR_Grange_overlap)!=0){
    output="3-UTR"
  }else{
    output="CDs"
  }
  return(output)
}
binding_cat=apply(meta,1,function(x) binding_catogory(x[8],as.numeric(x[3])
                                          ,as.numeric(x[4]),df))
meta_binding = meta %>% mutate(target_binding_site=binding_cat)
write.table(meta_binding,"metatable.tsv",sep = "\t",col.names = NA)

#write the table of interest as latex input
system("mv *.tsv co_expressed/")
system("mv Li_ANNOgesic_1* co_expressed/")
system("mv Li_ANNOgesic_2* co_expressed/")
system("mv Li_ANNOgesic_3* co_expressed/")
system("mv Li_ANNOgesic_4* co_expressed/")
system("mv Li_ANNOgesic_5* co_expressed/")
system("mv Li_ANNOgesic_6* co_expressed/")
system("mv *.csv co_expressed/")
#====================================
library(taRifx)
library(xtable)
library(readr)
library(dplyr)
gene_of_interest="Li_ANNOgesic_110"
write.table(meta_binding %>% filter(grepl(gene_of_interest,file)),paste0("metatable",gene_of_interest,".tsv"),sep = "\t",col.names = NA)

#gene_list_input=read.csv(paste0("~/Desktop/phd/programming/Integrated_analysis/Stephen_DESeq/Enrichment_analysis/mapped_",direction,"_",pair,"control.csv"))
gene_list_input_extract=as.data.frame(meta_binding %>% filter(grepl(gene_of_interest,file)) %>% select(file,target_binding_site,E))
gene_list_input_extract[,1]=gsub(paste0(gene_of_interest,"_"),"",gsub(".csv","",gene_list_input_extract[,1]))
colnames(gene_list_input_extract)[c(1,2,3)]=c("target","target_binding_site","minimal_binding_energy (kcal/mol)")
gene_list_input_extract.big=xtable(gene_list_input_extract, caption = paste0("IntaRNA targets of ",gene_of_interest)
                                #   ,label=paste0('IntaRNA ',gene_of_interest),digits=-2)
                                ,label=paste0('IntaRNA ',gene_of_interest))
#digits(gene_list_input_extract.big) <- c(0,0,0,8,8,8)
#align(gene_list_input_extract.big) <-"|p{.10\\textwidth}|p{.20\\textwidth}|p{.20\\textwidth}|p{.20\\textwidth}|p{.10\\textwidth}|p{.10\\textwidth}|"
align(gene_list_input_extract.big) <-"|l|l|l|l|"

output=print(gene_list_input_extract.big,math.style.exponents = TRUE, hline.after=c(-1, 0), tabular.environment = "longtable",include.rownames = F)
write(output,file=paste0("~/Desktop/phd/programming/Integrated_analysis/WGCNA/IntaRNA/output/",gene_of_interest,".tex"))
system("mv *.tex tex/")

#=======================================
setwd("~/Desktop/phd/programming/Integrated_analysis/WGCNA/IntaRNA/output/random/")

sRNA_mRNA_df=read_delim("../sRNA_mRNA_df.tsv",delim = "\t")
# all_sRNA=sRNA_mRNA_df$sRNA
# all_mRNA=sRNA_mRNA_df$mRNA
# sRNA_mRNA_pair=paste0(all_sRNA,"_",all_mRNA)
# i=0
# sRNA_rand_list=c()
# mRNA_rand_list=c()
# while (i < 10000){
#   i=i+1
#   print(i)
#   sRNA_rand=sample(unique(all_sRNA),1)
#   mRNA_rand=sample(unique(all_mRNA),1)
#   rand_pair=paste0(sRNA_rand,"_",mRNA_rand)
#   if (rand_pair %in% sRNA_mRNA_pair){
#     i=i-1
#     next()
#   }
#   else{
#     sRNA_rand_list=c(sRNA_rand_list,sRNA_rand)
#     mRNA_rand_list=c(mRNA_rand_list,mRNA_rand)
#   }
# }
#rand_sRNA_mRNA_df=data.frame(sRNA=sRNA_rand_list,mRNA=mRNA_rand_list)
#write.table(rand_sRNA_mRNA_df,"rand_sRNA_mRNA_df.tsv",sep="\t",row.names = F, quote = F)
rand_sRNA_mRNA_df=read_delim("rand_sRNA_mRNA_df.tsv",delim = "\t")
annotation=as.data.frame(read_delim("~/Desktop/phd/programming/RNA_crosslinking_analysis/Cj_crosslinking/genomeDir/11168_genome_sRNAs_exon_updated_shortened_visualised.gtf", delim="\t",col_names = F))
annotation_remove_start=annotation %>% filter(grepl("antisense",X9)) %>% select(X4)
gene_id=gsub("gene_id \"","",gsub("\";.*","",annotation[,9]))
# Run_Inta_rand=function(sRNA_id,mRNA_id){
#   sRNA=sRNA_seq(sRNA_id)
#   mRNA=mRNA_seq(mRNA_id,df)
#   output=paste0("~/Desktop/phd/programming/Integrated_analysis/WGCNA/IntaRNA/output/random/",sRNA_id,"_",gsub(" ","",gsub("\\(","",gsub("\\)","",gsub("\'","",mRNA_id)))),".csv")
#   #system(paste0("~/miniconda3/bin/IntaRNA -q ~/Desktop/phd/programming/RNA_crosslinking_analysis/Cj_crosslinking/IntaRNA_output/sRNA.fa -t ~/Desktop/genomeDir/all_gene.fasta --personality=IntaRNAsTar --threads 4 -n 10 --outMode C --out ",output ))
#   system(paste0("~/anaconda3/bin/IntaRNA -q ",as.character(sRNA)," -t ",as.character(mRNA)," --threads 4 -n 1 --outMode C --out ",output,"  --personality=IntaRNAsTar" ))
# }
for (i in 1:nrow(rand_sRNA_mRNA_df)){
  #  print(paste0("sRNA:",sRNA_mRNA_df[i,1]))
  #  print(paste0("mRNA:",sRNA_mRNA_df[i,2]))
  Run_Inta(rand_sRNA_mRNA_df[i,1],rand_sRNA_mRNA_df[i,2])
}
#run IntaRNA
#setwd("")
meta_rand=matrix(,nrow=0,ncol=8)
no_interaction_meta_rand=matrix(,nrow=0,ncol=1)
for (file in list.files(pattern = ".csv")){
  input=read.csv(file,sep = ";",header=T)
  if (nrow(input)==0){
    print(file)
    no_interaction_meta_rand=rbind(no_interaction_meta_rand,file)
    next()
  }
  input_cbind=cbind(input,file)
  meta_rand=rbind(meta_rand,input_cbind)
}
write.table(no_interaction_meta_rand,"no_interactions_rand.tsv",sep = "\t",row.names = F)
binding_cat_rand=apply(meta_rand,1,function(x) binding_catogory(x[8],as.numeric(x[3])
                                                      ,as.numeric(x[4]),df))
meta_binding_rand = meta_rand %>% mutate(target_binding_site=binding_cat_rand)
write.table(meta_binding_rand,"metatable_rand.tsv",sep = "\t",col.names = NA)
system("mv *.tsv random/")
system("mv Li_ANNOgesic_1* random/")
system("mv Li_ANNOgesic_2* random/")
system("mv Li_ANNOgesic_3* random/")
system("mv Li_ANNOgesic_4* random/")
system("mv Li_ANNOgesic_5* random/")
system("mv Li_ANNOgesic_6* random/")
system("mv *.csv random/")
#=======================================
###make energy distribution curve
library(readr)
library(dplyr)
library(ggplot2)

setwd("~/Desktop/phd/programming/Integrated_analysis/WGCNA/IntaRNA/output/co_expressed/")

#read metatable
meta=as.data.frame(read_delim("metatable.tsv",delim = "\t"))
#read the energy columns. Plot curve for CDs, 5'UTR, 3'UTR
all_E=data.frame(Binding_energy=meta$E,target_binding_site=meta$target_binding_site,type="co-expressed")
 density_plot_co=ggplot(all_E,aes(x=-(Binding_energy),colour=target_binding_site))+
   geom_density()+labs(x="-(binding energy) kcal/mol",y="Density of sRNA-mRNA pairs",colour="targeted regions")+
   scale_colour_manual(values = c("darkgreen","red","blue"))+
   theme(axis.title.y = element_text(size = 7))
#   theme(axis.title.y = element_text(size = 7),legend.title = element_blank())
 density_plot_co
#bar plot for numbers between these 3 features
#compare numbers between those with interaction and without interactions

no_interaction_meta=read_delim("no_interactions.tsv",delim = "\t")

bar_df=as.data.frame(table(meta$target_binding_site))
bar_df_added=rbind(bar_df,data.frame(Var1="no_intereactions",Freq=nrow(no_interaction_meta)))
bar_df_added = bar_df_added%>% mutate(type="co-expressed",proportion=Freq/sum(Freq))
# bar_plot_rand=ggp
 bar_plot_co=ggplot(bar_df_added,aes(x=Var1,y=proportion,fill=Var1))+theme_minimal()+
   geom_bar(stat="identity")+labs(x="binding position",y="Proportion of sRNA-mRNA pairs")+
   geom_text(aes(label=signif(proportion,3)),vjust=-0.2,size=2)+theme(legend.position="none",axis.title.y = element_text(size = 7))
 bar_plot_co  
#Randomly generate 10000 combination between sRNA and mRNA
#for loop
#randon samples for sRNA and mRNA
#if any row in metatable match the combination, i=i-1
#next()
setwd("~/Desktop/phd/programming/Integrated_analysis/WGCNA/IntaRNA/output/random/")

#read metatable
meta_rand=as.data.frame(read_delim("metatable_rand.tsv",delim = "\t"))
#read the energy columns. Plot curve for CDs, 5'UTR, 3'UTR
all_E_rand=data.frame(Binding_energy=meta_rand$E,target_binding_site=meta_rand$target_binding_site,type="not co-expressed")
 density_plot_rand=ggplot(all_E_rand,aes(x=-(Binding_energy),colour=target_binding_site))+
   geom_density(linetype="dashed")+labs(x="-(binding energy) kcal/mol",y="Density of sRNA-mRNA pairs",colour="targeted regions")+
   scale_colour_manual(values = c("darkgreen","red","blue"))+
   theme(axis.title.y = element_text(size = 7))
 density_plot_rand
#bar plot for numbers between these 3 features
#compare numbers between those with interaction and without interactions

no_interaction_meta=read_delim("no_interactions_rand.tsv",delim = "\t")

bar_df_rand=as.data.frame(table(meta_rand$target_binding_site))
bar_df_added_rand=rbind(bar_df_rand,data.frame(Var1="no_intereactions",Freq=nrow(no_interaction_meta)))
#bar_df_added_rand$type="not co-expressed"
bar_df_added_rand = bar_df_added_rand%>% mutate(type="not co-expressed",proportion=Freq/sum(Freq))
bar_plot_rand=ggplot(bar_df_added_rand,aes(x=Var1,y=proportion,fill=Var1))+theme_minimal()+
  geom_bar(stat="identity")+labs(x="binding position",y="Proportion of sRNA-mRNA pairs")+
  geom_text(aes(label=signif(proportion,3)),vjust=-0.2,size=2)+theme(legend.position="none",axis.title.y = element_text(size = 7))
bar_plot_rand  

#restructure df
#add column, state whether it is co-expressed or random
#random = dotted line
#
combined_all_E_df=rbind(all_E,all_E_rand)
combined_bar_df=rbind(bar_df_added,bar_df_added_rand)


 density_plot_5UTR=ggplot(combined_all_E_df %>% filter(target_binding_site=="5-UTR"),aes(x=-(Binding_energy),colour=target_binding_site,linetype=type))+
   geom_density()+labs(x="-(binding energy) kcal/mol",y="Density of sRNA-mRNA pairs",colour="targeted regions",linetype="sRNA-mRNA co-expression")+
   scale_color_manual(values=c("red"))+
   theme(axis.title.y = element_text(size = 7))+
   guides(colour=guide_legend(order=1),linetype=guide_legend(order=2))
 density_plot_5UTR
 
 density_plot_3UTR=ggplot(combined_all_E_df %>% filter(target_binding_site=="3-UTR"),aes(x=-(Binding_energy),colour=target_binding_site,linetype=type))+
   geom_density()+labs(x="-(binding energy) kcal/mol",y="Density of sRNA-mRNA pairs",colour="targeted regions",linetype="sRNA-mRNA co-expression")+
   scale_color_manual(values=c("darkgreen"))+
   theme(axis.title.y = element_text(size = 7))+
   guides(colour=guide_legend(order=1),linetype=guide_legend(order=2))
 
 density_plot_3UTR
 
 density_plot_CDs=ggplot(combined_all_E_df %>% filter(target_binding_site=="CDs"),aes(x=-(Binding_energy),colour=target_binding_site,linetype=type))+
   geom_density()+labs(x="-(binding energy) kcal/mol",y="Density of sRNA-mRNA pairs",colour="targeted regions",linetype="sRNA-mRNA co-expression")+
   scale_color_manual(values=c("blue"))+
   theme(axis.title.y = element_text(size = 7))+
   guides(colour=guide_legend(order=1),linetype=guide_legend(order=2))
 
 density_plot_CDs
 
 bar_plot_combined=ggplot(combined_bar_df,aes(x=Var1,y=proportion,fill=type))+theme_minimal()+
   geom_bar(stat="identity",position = position_dodge())+labs(x="binding position",y="Proportion of sRNA-mRNA pairs")+
   geom_text(aes(label=signif(proportion,3)),vjust=-0.2,position=position_dodge(0.9),size=2)+
   theme(axis.title.y = element_text(size = 7))
 
 bar_plot_combined 
 
 setwd("~/Desktop/phd/programming/Integrated_analysis/WGCNA/IntaRNA/output/figures/")
 density_plot_co
 ggsave("Energy_distribution_density_co-expressed.png")
 density_plot_rand
 ggsave("Energy_distribution_density_random.png")
 density_plot_5UTR
 ggsave("Energy_distribution_density_5UTR.png")
 density_plot_3UTR
 ggsave("Energy_distribution_density_3UTR.png")
 density_plot_CDs
 ggsave("Energy_distribution_density_CDs.png")
 bar_plot_co
 ggsave("Energy_distribution_bar_co-expressed.png")
 bar_plot_rand
 ggsave("Energy_distribution_bar_random.png")
 bar_plot_combined
 ggsave("Energy_distribution_bar_combined.png")

 #=====================================
 library(readr)
 library(dplyr)
 library(ggplot2)
 ###now include correlation coefficient . See if there is a difference in positive and negative coefficient
 
 ###or any difference between modules
 
 
 #read co-expression table
 co_meta=as.data.frame(read_delim("~/Desktop/phd/programming/Integrated_analysis/WGCNA/co_expression_metatable_sRNA.tsv",delim = "\t"))
 co_meta_sRNA=apply(co_meta,1,function(x) if (grepl("Li_",x[1])){
   return (x[1])
 }else{
   return(x[2])
 })
 co_meta_mRNA=apply(co_meta,1,function(x) if (grepl("Li_",x[1])){
   return (x[2])
 }else{
   return(x[1])
 })
 co_meta_mRNA_gsubbed=gsub(" ","",gsub("\\(","",gsub("\\)","",gsub("\'","",co_meta_mRNA))))
 co_meta_sRNA_mRNA=paste0(co_meta_sRNA,"_",co_meta_mRNA_gsubbed)
 #read IntaRNA metatable
 Inta_meta=as.data.frame(read_delim("~/Desktop/phd/programming/Integrated_analysis/WGCNA/IntaRNA/output/co_expressed/metatable.tsv",delim = "\t"))
 Inta_meta_sRNA_mRNA=gsub(".csv","",Inta_meta$file)
 #match columns
 index=match(Inta_meta_sRNA_mRNA,co_meta_sRNA_mRNA)
 co_meta_matched=co_meta[index,]
 #fit in correlation coefficient
 #fit in module
 Inta_meta_added=Inta_meta %>% mutate(bicor=co_meta_matched$correlation,bicor_sign=sign(co_meta_matched$correlation),bicor_padj=co_meta_matched$correlation_padj,
                      module=co_meta_matched$module,weight=co_meta_matched$weight)
 #make bar plot and density plot.
 setwd("~/Desktop/phd/programming/Integrated_analysis/WGCNA/IntaRNA/output/figures/")
 
 Inta_meta_added_density=ggplot( Inta_meta_added,aes(x=-(E),colour=as.factor(bicor_sign)))+
   geom_density()+labs(x="-(binding energy) kcal/mol",y="Density of sRNA-mRNA pairs",colour="Correlation coefficient sign")+
   scale_colour_manual(values = c("darkgreen","red"))+
   theme(axis.title.y = element_text(size = 7))
 Inta_meta_added_density
 ggsave("Energy_distribution_density_co-expressed_bicor.png")
  
 Inta_meta_added_df=as.data.frame(Inta_meta_added %>% group_by(bicor_sign) %>%
   count(target_binding_site) )#%>% 
 #  mutate(proportion=n/sum(n)))
 
 no_interaction=as.data.frame(read_delim("~/Desktop/phd/programming/Integrated_analysis/WGCNA/IntaRNA/output/co_expressed/no_interactions.tsv",delim = "\t"))[,1]
 no_interaction_sRNA_mRNA=gsub(".csv","",no_interaction)
 no_interaction_index=match(no_interaction_sRNA_mRNA,co_meta_sRNA_mRNA)
 no_interaction_co_meta_matched=co_meta[no_interaction_index,]
 #fit in correlation coefficient
 #fit in module
 no_interaction_added=as.data.frame(no_interaction) %>% mutate(bicor=no_interaction_co_meta_matched$correlation,bicor_sign=sign(no_interaction_co_meta_matched$correlation),bicor_padj=no_interaction_co_meta_matched$correlation_padj,
                                      module=no_interaction_co_meta_matched$module,weight=no_interaction_co_meta_matched$weight)
 no_interaction_added_counted=as.data.frame(no_interaction_added  %>% group_by(bicor_sign)%>% 
                                              count() %>% mutate(target_binding_site="no interactions"))
 Inta_meta_combined_df=rbind(Inta_meta_added_df,no_interaction_added_counted) %>% group_by(bicor_sign) %>% mutate(proportion=n/sum(n))
 
 Inta_meta_bar_plot=ggplot(Inta_meta_combined_df,aes(x=target_binding_site,y=proportion,fill=as.factor(bicor_sign)))+
   geom_bar(stat="identity",position = position_dodge())+labs(x="binding position",y="Proportion of sRNA-mRNA pairs",fill="Correlation coefficient sign")+
   geom_text(aes(label=signif(proportion,3)),vjust=-0.2,position=position_dodge(0.9),size=2)+
   theme(axis.title.y = element_text(size = 7))
 
 Inta_meta_bar_plot
 ggsave("Energy_distribution_bar_co-expressed_bicor.png")
 
 #cor(Inta_meta_added %>% filter(bicor_sign==-1) %>% select(bicor),Inta_meta_added %>% filter(bicor_sign==-1) %>% select(E))
# cor(Inta_meta_added %>% filter(bicor_sign==1) %>% select(bicor),Inta_meta_added %>% filter(bicor_sign==1) %>% select(E))
 #colour with different coefficient and module
 #make pearson correlation plot between bicor correlation coefficient and binding energy
 negative_df=as.data.frame(Inta_meta_added %>% filter(bicor_sign==-1))
# title=paste0("bicor = ",signif(cor(Inta_meta_added %>% filter(bicor_sign==-1) %>% select(bicor),Inta_meta_added %>% filter(bicor_sign==-1) %>% select(E))
                       #         ,3),"\n padj = ",signif(padj,3))
 
 ggplot(negative_df, aes_string(x= "bicor", y="E")) + 
   labs(y="binding energy (kcal/mol)", x="co-expression correlation coeffiecient")+
   geom_point(shape=18, color="blue") + #+
   geom_smooth(method=lm,  linetype="dashed",
               color="darkred", fill="blue")+ ggtitle("IntaRNA binding energy VS bicor coefficient \n (sRNA-mRNA pairs with negative correlation)")
# ggsave("correlation_positive.png")
 ggsave("correlation_negative.png")
 positive_df=as.data.frame(Inta_meta_added %>% filter(bicor_sign==1))
 ggplot(positive_df, aes_string(x= "bicor", y="E")) + 
   labs(y="binding energy (kcal/mol)", x="co-expression correlation coeffiecient")+
   geom_point(shape=18, color="blue") + #+
   geom_smooth(method=lm,  linetype="dashed",
               color="darkred", fill="blue")+ ggtitle("IntaRNA binding energy VS bicor coefficient \n (sRNA-mRNA pairs with positive correlation)")
 ggsave("correlation_positive.png")
 
 #==================================================
 
 path="~/Desktop/phd/programming/Integrated_analysis/WGCNA/IntaRNA/output/"
 setwd(path)
 library(dplyr)
 library(readr)
 library(annotate)
 library(Biostrings)
 library(rBLAST)
 library(seqinr)
 
 #design a function to pick out all coding genes .
 fasta_vector=read.fasta("~/Desktop/phd/programming/RNA_crosslinking_analysis/Cj_crosslinking/genomeDir/NCTC11168_genome.fasta.txt")[[1]]
 fasta_seq=paste(fasta_vector,collapse="")
 
 df=read.csv("~/Desktop/phd/programming/Integrated_analysis/gene_annotation_UTR.csv")
 # mRNA_seq=function(mRNA,df){
 #   # index=which(gene_id==mRNA)
 #   print(mRNA)
 #   index=which(df$gene_name == as.character(mRNA))
 #   #UTR_length=200 #default length
 #   #if mRNA is in Jenna's file, find its primary and secondary TSS
 #   #find the maximum abs distance between mRNA 5' end and TSS
 #   #cvhange the UTR_length variable accordingly
 #   
 #   sequence=substring(genome,df$estimated_5end[index],df$end[index])
 #   if (df$strand[index]=="-"){
 #     sequence=substring(genome,df$start[index],df$estimated_5end[index])
 #     sequence=reverseComplement(DNAString(sequence))
 #   }
 #   return(RNAString(DNAString(sequence)))
 #   
 # }
 sequence_extract=function(start,end,strand,gene_id,df){
   # fasta_seq=read.fasta("~/Desktop/genomeDir/NCTC11168_genome.fasta.txt")[[1]]
   start=as.numeric(start)
   end=as.numeric(end)
   print(gene_id)
   index=which(df$gene_name == as.character(gene_id))
   if (strand=="+" & grepl("Li_",gene_id)==FALSE){
     #to include UTR, include genomic subsequence of eg from 200 up- to 100 downstream of the annotated start codon.
     #record their start, stop and strand
     #do so for sense and antisense strand
     # start=start-200
     start=df$estimated_5end[index]
     # end=end+100
     # end=start+300  #for looking into the 5 UTR
     end=df$end[index]
     if (start <=0){
       start=1
     }
     if (end > nchar(fasta_seq)){
       end=nchar(fasta_seq)
     }
   }
   if (strand=="-" & grepl("Li_",gene_id)==FALSE){
     # start=start-100
     # end=end+200
     end=df$estimated_5end[index]
     # start=end-300 #for looking into the 5UTR
     start=df$start[index]
     
     if (start <=0){
       start=1
     }
     if (end > nchar(fasta_seq)){
       end=nchar(fasta_seq)
     }
   }
   
   output=substring(fasta_seq,start,end)
   if (strand=="-"){
     output=as.character(reverseComplement(DNAString(output)))
   }
   return (RNAString(DNAString(output)))
 }
 
 gff=read.table("~/Desktop/phd/programming/RNA_crosslinking_analysis/Cj_crosslinking/genomeDir/11168_genome_sRNAs_exon_updated_shortened_visualised.gtf",sep="\t")
 gene_name=gsub(";.*","",gsub("gene_id ","",gff$V9))
 gff_cbind=cbind(gff,gene_name) 
 gff_cbind_filtered=gff_cbind %>% filter(!grepl("rna",gene_name))
 gene_seq=as.list(apply(gff_cbind_filtered,1,function(x) sequence_extract(x[4],x[5],x[7],x[10],df)))
 #write.fasta(gene_seq,gene_name,"~/Desktop/genomeDir/all_gene.fasta")
 #gsub(" ","",gsub("\\(","",gsub("\\)","",gsub("\'","",gene_name[!grepl("rna",gene_name)]))))
 write.fasta(gene_seq,gsub(" ","",gsub("\\(","",gsub("\\)","",gsub("\'","",gene_name[!grepl("rna",gene_name)])))),"~/Desktop/phd/programming/Integrated_analysis/WGCNA/IntaRNA/all_gene.fasta")
 #base on that, fine the fasta substring for all rows
 #put those sequences and names into a list
 #write them into a fasta file
 
 
 
 
 
 
 #pick out all predicted sRNA. Pick out all elements with "ANNOgesic"
 #sRNA_list=unique(c(as.character(df$gene1[grep("ANNOgesic",df$gene1)]),as.character(df$gene2[grep("ANNOgesic",df$gene2)])))
 finalised_sRNA=gff %>% filter(grepl("ANNOgesic",gene_name))
 #sRNA_list=as.character(df$name)
 sRNA_list=gene_name[grepl("ANNOgesic",gene_name)]
 #defube a function. Use biostring to extract the DNA sequence
 # sRNA_blast=data.frame()
 #run blastSequences. Store output in a list
 # finalised_sRNA=read.csv("~/Desktop/phd/programming/sRNA_prediction/sRNA_final_prediction.csv")
 genome=readDNAStringSet("~/Desktop/phd/programming/RNA_crosslinking_analysis/Cj_crosslinking/genomeDir/NCTC11168_genome.fasta.txt")
 sRNA_seq=c()
 for (sRNA in sRNA_list){
   print(sRNA)
   #for (sRNA in sRNA_list[14:length(sRNA_list)]){ #the algorithm terminated after 11 outputs
   # finalised_sRNA=read.csv("~/Desktop/phd/programming/sRNA_prediction/sRNA_final_prediction.csv")
   # genome=readDNAStringSet("~/Desktop/genomeDir/NCTC11168_genome.fasta.txt")
   row=which(sRNA_list==sRNA)
   print(row)
   sequence=substring(genome,finalised_sRNA$V4[row],finalised_sRNA$V5[row])
   # as.character(substring(genome,1,10))
   if (finalised_sRNA$V7[row]=="-"){
     sequence=reverseComplement(DNAString(sequence))
   }
   sRNA_seq=c(sRNA_seq,RNAString(DNAString(sequence)))
   
  
   
}
 write.fasta(sRNA_seq, sRNA_list,"~/Desktop/phd/programming/Integrated_analysis/WGCNA/IntaRNA/all_sRNA.fa")
 
    
 output="~/Desktop/phd/programming/Integrated_analysis/WGCNA/IntaRNA/output/genome/all_sRNA_genome.csv"
 system(paste0("~/anaconda3/bin/IntaRNA -q ~/Desktop/phd/programming/Integrated_analysis/WGCNA/IntaRNA/all_sRNA.fa -t ~/Desktop/phd/programming/Integrated_analysis/WGCNA/IntaRNA/all_gene.fasta --personality=IntaRNAsTar --threads 4 -n 1 --outMode C --out ",output,
               " --outCsvCols 'id1,id2,seq1,seq2,subseq1,subseq2,subseqDP,subseqDB,start1,end1,start2,end2,hybridDP,hybridDPfull,hybridDB,hybridDBfull,E'" ))
 
# output=paste0("~/Desktop/phd/programming/RNA_crosslinking_analysis/Cj_crosslinking/IntaRNA_output/",sRNA,"_genome")
# system(paste0("~/miniconda3/bin/IntaRNA -q ~/Desktop/phd/programming/RNA_crosslinking_analysis/Cj_crosslinking/IntaRNA_output/sRNA.fa -t ~/Desktop/genomeDir/all_gene.fasta --personality=IntaRNAsTar --threads 4 -n 10 --outMode C --out ",output," --outCsvCols '*'" ))
 
#system("Rscript --vanilla ../IntaRNA_CSV_p-value.R genome/all_sRNA_genome.csv")
 
 #compute p-value for the output energy using IntaRNA_CSV_p-value.R. That was ran in the CLIMB server (somehow it didnt work in my current laptop)
 
#Inta_fdr=read.csv("genome/all_sRNA_genome.csv",sep=";")
all_output=read.csv("genome/all_sRNA_genome.csv",sep=";")
 for (sRNA in sRNA_list){
   all_output_sRNA=all_output %>% filter(id2==sRNA)
   print(sRNA)
   print(nrow(all_output_sRNA))
   csv_file=paste0("genome/all_sRNA_genome_",sRNA,".csv")
   write.table(all_output_sRNA,csv_file,sep=";",row.names = F)
   system(paste0("Rscript --vanilla ../IntaRNA_CSV_p-value.R ", csv_file))
   
 }
files=list.files(path = "genome/", pattern="ANNOgesic")
meta=matrix(, nrow = 0, ncol = 8)
colnames(meta)=c(colnames(all_output),"fdr")
for (file in files){
  input=paste0("genome/",file)
  meta=rbind(meta,read_delim(input,delim=";"))
}
#meta_df=as.data.frame(meta)
write.table(as.data.frame(meta),"genome/metatable.tsv",sep="\t")
View(as.data.frame(meta) %>% filter(fdr <= 0.05))
#system("Rscript --vanilla ../IntaRNA_CSV_p-value.R genome/all_sRNA_genome.csv")

#Venn diagram for co-expression, crosslinkinng and IntaRNA

#==========================================
#for genome prediction, integrate with co-expression analysis

#run blastSequences. Store output in a list
#finalised_sRNA=read.csv("~/Desktop/phd/programming/sRNA_prediction/sRNA_final_prediction.csv")
co_meta=as.data.frame(read_delim("~/Desktop/phd/programming/Integrated_analysis/WGCNA/co_expression_metatable_sRNA.tsv",delim="\t"))
#sRNA_mRNA=co_meta %>% filter((grepl("ANNOgesic",fromNode) & !grepl("ANNOgesic",toNode)) | 
                 #           (!grepl("ANNOgesic",fromNode) & grepl("ANNOgesic",toNode)) )
co_sRNA_duplex=apply(co_meta,1,function(x) if (grepl("ANNOgesic",x[1])){
  return (paste0(x[1],"_",x[2]))
}else{
  return (paste0(x[2],"_",x[1]))
})
#mRNA=apply(sRNA_mRNA,1,function(x) if (grepl("ANNOgesic",x[1])){
 # return (x[2])
#}else{
 # return(x[1])
#})
#sRNA_duplex_pasted=
#read metatable for IntaRNA
#read metatable for co-expression
Inta_meta=Inta_meta=as.data.frame(read.csv("~/Desktop/phd/programming/Integrated_analysis/WGCNA/IntaRNA/output/genome/metatable.tsv",sep="\t"))
Inta_sRNA_duplex=paste0(Inta_meta$id2,"_",Inta_meta$id1)
#in co-expression metatable, include binding energy, p-values and fdr

co_meta_Inta=cbind(co_meta,Inta_meta[match(co_sRNA_duplex,Inta_sRNA_duplex),3:ncol(Inta_meta)])
write.table(co_meta_Inta,"~/Desktop/phd/programming/Integrated_analysis/WGCNA/co_expression_metatable_sRNA_Inta.tsv",sep="\t",row.names = F,quote = F)
write.table(co_meta_Inta %>% filter(p.value <= 0.05),"~/Desktop/phd/programming/Integrated_analysis/WGCNA/co_expression_metatable_sRNA_Inta_filtered.tsv",sep="\t",row.names = F,quote = F)

#highlight those that are significant
#system(paste0("~/anaconda3/bin/IntaRNA -q ",as.character(sRNA)," -t ",as.character(mRNA)," --threads 4 -n 1 --outMode C --out ",output,"  --personality=IntaRNAsTar " ))
#=======
