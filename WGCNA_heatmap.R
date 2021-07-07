library(pheatmap)
library(dplyr)
library(DESeq2)
library(readr)
library(RColorBrewer)
### create heatmap for all genes
setwd("~/Desktop/phd/programming/Integrated_analysis/WGCNA/")
lnames = load(file = "WGCNA.RData");
#The variable lnames contains the names of loaded variables.
lnames
#load read count data and coldata
samples=read.csv("~/Desktop/phd/programming/Integrated_analysis/samples.csv")  #read sample table (with batch, conditions...)
kable(samples)

raw_counts <- as.data.frame(read_delim('~/Desktop/phd/programming/Integrated_analysis/Count_table.tsv',delim = "\t")) #read raw count table
rownames(raw_counts)=raw_counts[,1]
raw_counts=raw_counts[,-1]
head(raw_counts)
samples=samples[which(samples$sample_id %in% colnames(raw_counts)),] 
raw_counts=raw_counts[,match(as.character(samples$sample_id),colnames(raw_counts))]
gene_ids= rownames(raw_counts)#read all gene ids
### Data Preparation
num_conditions=nlevels(samples$condition)
#pal=colorRampPalette(brewer.pal(num_conditions,"Set1"))(num_conditions)
pal=colorRampPalette(brewer.pal(9,"Set1"))(num_conditions) #maximum number of colour in set 1 is 9
cond_colors=pal[as.integer(samples$condition)]
cond_df=data.frame(cond=samples$condition)
rownames(cond_df)=rownames(cor(raw_counts))
#par(mar=c(1,1,1,1))
samples_filtered=samples %>% filter(Author=="Jenna")
raw_counts=raw_counts[,match(as.character(samples_filtered$sample_id),colnames(raw_counts))]
read_counts=raw_counts
#read the TPM file
######   calculate TPM from read counts author - Jenna
###### need raw count files which can be generated using combining_read_counts_Jenna.R from filtered text files taken off from CLIMB which is basically from coverageBed but filtered

#setwd("~/OneDrive - University of Warwick/OneDrive/PhD/RNAtagSeq (differential gene expression)/calculate_TPM/")
#read_counts = read.csv("counts.all.csv")
#which(is.na(read_counts)==TRUE)
#read_counts$X = gsub('"','',read_counts$X)
#rownames(read_counts)= read_counts$X
#read_counts = read_counts[,-1]
GFF_file = read.delim("~/Desktop/phd/programming/RNA_crosslinking_analysis/Cj_crosslinking/genomeDir/11168_genome_sRNAs_exon_updated_shortened_visualised.gtf", header = F)
#GFF_file = GFF_file[-1,]#be careful of running these one after the other but it gets rid of the two lines at the top
#GFF_file = GFF_file[-1,]
#### STOP #### don't repeat above
##gene_feature = as.character(GFF_file$V3) #i.e. gene, ncRNA, pseudogene, etc
#ID = as.character(GFF_file$V9)


gene_length = data.frame(Gene_name = rownames(read_counts)) #### <- start from here, make a data_frame of the gene lengths
G = gene_length$Gene_name
start = as.numeric(GFF_file$V4) #start of gene feature as a number
end = as.numeric(GFF_file$V5)
Gene_Length = end - start#empty NA vector for combining to dataframe later


gene_length = cbind(gene_length, Gene_Length)
#which(is.na(gene_length)==TRUE)
write.csv(gene_length, file="Gene_Length.csv")


gene_length = read.csv("Gene_Length.csv") ### <- start here if you have already created Gene_length.csv

### same order as read counts
##################### now do the TPM calculation
RPK = data.frame(matrix(NA, nrow=nrow(read_counts), ncol = ncol(read_counts)) )#RPK is dividing the read count by gene length, make an empty dataframe same dimensions as read_counts
names(RPK) = names(read_counts)
rownames(RPK)= rownames(read_counts)
for(i in 1:nrow(read_counts)){
  
  for(j in 1:ncol(read_counts)){
    #if(as.character(read_counts$X)[i] == as.character(gene_length$Gene_name)[i]){#same order as read counts but might as well just in case
    #print("yes")
    new_num = as.numeric(read_counts[i,j])/gene_length$Gene_Length[i] #for every row i do the division
    RPK[i,j] = new_num
    
    #}
  }
  
}

############ scaling factor is the sum of each column of RPK and then divide by 1 million
scaling = c()
for(i in 1:ncol(RPK)){
  sum_col = sum(RPK[,i])/1000000 #sum divided by million
  scaling = c(scaling, sum_col)
}

########### now divide the RPK value by scaling factor
TPM = data.frame(matrix(NA, nrow=nrow(read_counts), ncol = ncol(read_counts))) #another empty dataframe but will be final output
names(TPM) = names(read_counts)
rownames(TPM)= rownames(read_counts)
for(j in 1:ncol(RPK)){#do it by column first since the scaling factor is per column
  for(i in 1:nrow(RPK)){
    tpm_value = RPK[i,j]/scaling[j]
    TPM[i,j] = tpm_value
    
  }
}
##### check if the sum adds up to 1 million for each column to see if it has worked
i = 1
while(i <length(scaling)){
  sum_final = sum(TPM[,i])
  print(sum_final)
  i=i+1
}

which(is.na(TPM) == TRUE)

write.csv(TPM, file="TPM.csv")
system("rm -r ../gene_heatmap/")
system("mkdir ../gene_heatmap")
colnames(cond_df)="condition"
vst_counts=vst(as.matrix(raw_counts))
png(file="~/Desktop/phd/programming/Integrated_analysis/gene_heatmap/all_log2TPM.png", width = 1200, height = 1100)
#pheatmap(log2(TPM+1) %>% filter(grepl("Li",row.names(TPM))))
#input=log2(TPM+1)

#pheatmap(log2(TPM+1),annotation_col =cond_df,show_rownames=F,show_colnames =F,col=brewer.pal(9,"Blues"))
#pheatmap(log2(TPM+1),annotation_col =cond_df,show_rownames=F,show_colnames =T,scale="row",col=rev(brewer.pal(9,"RdBu")))
pheatmap(log2(TPM+1),annotation_col =cond_df,show_rownames=F,show_colnames =T,scale="row",color = colorRampPalette(c("green","black", "red"))(100))

#pheatmap(log2(TPM+1),annotation_col =cond_df,show_rownames=F,show_colnames =T,scale="row",col=rev(brewer.pal(9,"RdYlGn")))
#setHook("grid.newpage", NULL, "replace")
#grid.text("xlabel example", y=-0.07, gp=gpar(fontsize=16))
#grid.text("ylabel example", x=-0.07, rot=90, gp=gpar(fontsize=16))
#dev.off()
dev.off()
png(file="~/Desktop/phd/programming/Integrated_analysis/gene_heatmap/all_vst.png", width = 1200, height = 1100)
pheatmap(vst_counts,annotation_col =cond_df,show_rownames=F,show_colnames =T,scale="row",color = colorRampPalette(c("green","black", "red"))(100))
dev.off()
png(file="~/Desktop/phd/programming/Integrated_analysis/gene_heatmap/sRNA_vst.png", width = 1200, height = 1100)
#pheatmap(as.data.frame(vst_counts)%>% filter(grepl("ANNOgesic",row.names(vst_counts))),annotation_col =cond_df,show_rownames=T,show_colnames =T,col=brewer.pal(9,"Blues"),scale="row")
pheatmap(as.data.frame(vst_counts)%>% filter(grepl("ANNOgesic",row.names(vst_counts))),annotation_col =cond_df,show_rownames=T,show_colnames =T,scale="row",color = colorRampPalette(c("green","black", "red"))(100))

dev.off()


for (mods in all_mods){
  system(paste0("rm -r ../gene_heatmap/",mods,"/"))
  system(paste0("mkdir ../gene_heatmap/",mods,"/"))
  nodes_file=paste0("CytoscapeInput-nodes-",mods,".txt")
  folder=paste0("../gene_heatmap/",mods,"/")
  nodes=as.data.frame(read_delim(nodes_file,delim="\t"))[,1]
  png(file=paste0(folder,"sRNA_vst.png"), width = 1200, height = 1100)
  pheatmap(as.data.frame(vst_counts)%>% filter( row.names(vst_counts) %in% nodes),annotation_col =cond_df,show_rownames=T,show_colnames =T,scale="row",color = colorRampPalette(c("green","black", "red"))(100))
  dev.off()
}
#png(file="~/Desktop/phd/programming/Integrated_analysis/WGCNA/heaptmap_log_author.png", width = 1200, height = 1100)
#heatmap(cor(log_counts),RowSideColors = cond_colors,ColSideColors = cond_colors,main="Sample correlations (log transformed)",
#        scale = "none" , labRow=FALSE, labCol = FALSE,col=brewer.pal(9,"Blues"))
lnames = load(file = "ABC_transporters_info.RData");
#The variable lnames contains the names of loaded variables.
#subset for all ABC transporters
#ABC_gene
#ABC=ABC_gene
ABC_system[which(ABC_system=="Iron(III)")]="ferric"
ABC_df=data.frame(ABC_system=ABC_system)
rownames(ABC_df)=ABC_gene_name
png(file="~/Desktop/phd/programming/Integrated_analysis/gene_heatmap/ABC_vst.png", width = 1200, height = 1100)
pheatmap(as.data.frame(vst_counts) %>% filter(row.names(vst_counts) %in% ABC_gene_name),annotation_col =cond_df, annotation_row = ABC_df,show_rownames=T,show_colnames =T,scale = "row",color = colorRampPalette(c("green","black", "red"))(100))
dev.off()

ABC_df_iron=ABC_df %>% filter(ABC_system =="ferric" | row.names(ABC_df) %in% c("Cj1661","Cj1662","Cj1663"))
iron_system=c(rep("Rhodotorulic acid",2),rep("Enterochelin",4),
              rep("Haem",3),rep("Transferrins",3),"Rhodotorulic acid")
iron_system_df=data.frame(iron_system=iron_system)
rownames(iron_system_df)=rownames(ABC_df_iron)
png(file="~/Desktop/phd/programming/Integrated_analysis/gene_heatmap/ABC_iron_vst.png", width = 1200, height = 1100)
pheatmap(as.data.frame(vst_counts) %>% filter(row.names(vst_counts) %in% row.names(ABC_df_iron)),annotation_col =cond_df, annotation_row = iron_system_df,show_rownames=T,show_colnames =T,,scale = "row",color = colorRampPalette(c("green","black", "red"))(100))
dev.off()
#====================
#list all pairwise comparisons
library(pheatmap)
raw_counts <- as.data.frame(read_delim('~/Desktop/phd/programming/Integrated_analysis/Count_table.tsv',delim = "\t")) #read raw count table
rownames(raw_counts)=raw_counts[,1]
raw_counts=raw_counts[,-1]

system("rm ../Stephen_DESeq/results_pairwise_combination/metatable*")
folder="../Stephen_DESeq/results_pairwise_combination/"
all_files=list.files(folder)
table_binary=as.data.frame(matrix(,ncol=0,nrow=nrow(raw_counts)))
table_LFC=as.data.frame(matrix(,ncol=0,nrow=nrow(raw_counts)))
for(f in all_files){
 # if (f == "Cytoscape_DESeq"){
#  next()
  #}
  input=paste0(folder,"/",f)
  # input="Cytoscape_DESeq/iron_lim_M_vs_iron_rep_M/co_expression_metatable_sRNA.tsv.gz"
 # system(paste0("gunzip ",input))
 # input_gunzip=gsub(".gz","",input)
  re=read.csv(input)
  re$padj[which(is.na(re$padj))]=1 #ignore outliers
  #DE=as.numeric(apply(re,1,function(x) abs(as.numeric(x[3]))>=2 & as.numeric(x[7])<=0.05))
  DE=as.numeric(apply(re,1,function(x) abs(as.numeric(x[3]))>=0 & as.numeric(x[7])<=0.05))
 # which_down=which(apply(re,1,function(x) as.numeric(x[3])<=-2 & as.numeric(x[7])<=0.05))
  which_down=which(apply(re,1,function(x) as.numeric(x[3])<0 & as.numeric(x[7])<=0.05))
  DE[which_down]=-1
  output_binary=DE
 # output_LFC=apply(re,1,function(x) if(abs(as.numeric(x[3]))>=2 & as.numeric(x[7])<=0.05){return(as.numeric(x[3]))}else{return(0)})
  output_LFC=apply(re,1,function(x) if(abs(as.numeric(x[3]))>=0 & as.numeric(x[7])<=0.05){return(as.numeric(x[3]))}else{return(0)})
  
  name=gsub("control.csv","",gsub("../Stephen_DESeq/results_pairwise_combination//","",input))
  table_binary=cbind(table_binary,output_binary)
  colnames(table_binary)[ncol(table_binary)]=name
  
  #name=gsub("control.csv","",gsub("../Stephen_DESeq/results_pairwise_combination//","",input))
  table_LFC=cbind(table_LFC,output_LFC)
  colnames(table_LFC)[ncol(table_LFC)]=name
}
rownames(table_binary)=rownames(raw_counts)
rownames(table_LFC)=rownames(raw_counts)

write.csv(table_binary,"../Stephen_DESeq/results_pairwise_combination/metatable_binary.csv")
write.csv(table_LFC,"../Stephen_DESeq/results_pairwise_combination/metatable_LFC.csv")
png(file="~/Desktop/phd/programming/Integrated_analysis/gene_heatmap/all_DESeq_binary.png", width = 1200, height = 1100)
#row_to_remove=which(apply(table_binary,1,function(x) sd(x)==0))
pheatmap(table_binary,show_rownames=F,show_colnames =T,col=rev(brewer.pal(3,"RdBu")))
#pheatmap(log2(TPM+1),annotation_col =cond_df,show_rownames=F,show_colnames =F,col=brewer.pal(9,"Blues"))
dev.off()
png(file="~/Desktop/phd/programming/Integrated_analysis/gene_heatmap/all_DESeq_LFC.png", width = 1200, height = 1100)
#pheatmap(table_LFC,show_rownames=F,show_colnames =F, breaks=seq(-10, 10, by = 0.2))
range=abs(max(table_LFC))
pheatmap(table_LFC,show_rownames=F,show_colnames =F,color = colorRampPalette(c("blue","white", "red"))(100),
         breaks = seq(-range, range, length.out = 100))
#pheatmap(log2(TPM+1),annotation_col =cond_df,show_rownames=F,show_colnames =F,col=brewer.pal(9,"Blues"))
dev.off()

View(all_cond)
all_cond=gsub("control.csv","",gsub("Cytoscape_DESeq/","",directory))
selected_cond_file=read.csv("../selected_pairwise_comparisons.csv")
selected_cond=as.character(selected_cond_file[,1])

png(file="~/Desktop/phd/programming/Integrated_analysis/gene_heatmap/selected_DESeq_LFC.png", width = 1200, height = 1100)
range=abs(max(table_LFC))
pheatmap(table_LFC %>% select(selected_cond),show_rownames=F,show_colnames =T,color = colorRampPalette(c("blue","white", "red"))(100),
         breaks = seq(-range, range, length.out = 100))
#pheatmap(log2(TPM+1),annotation_col =cond_df,show_rownames=F,show_colnames =F,col=brewer.pal(9,"Blues"))
dev.off()





temp_sRNA=c()
temp_modules=c()
for (mods in all_mods){
 # system(paste0("rm -r ../gene_heatmap/",mods,"/"))
 # system(paste0("mkdir ../gene_heatmap/",mods,"/"))
  nodes_file=paste0("CytoscapeInput-nodes-",mods,".txt")
  #folder=paste0("../gene_heatmap/",mods,"/")
  nodes=as.data.frame(read_delim(nodes_file,delim="\t"))[,1]
  #png(file=paste0(folder,"sRNA_vst.png"), width = 1200, height = 1100)
#  pheatmap(as.data.frame(vst_counts)%>% filter( row.names(vst_counts) %in% nodes),annotation_col =cond_df,show_rownames=T,show_colnames =T,col=brewer.pal(9,"Blues"))
 # dev.off()
  sRNA_nodes=nodes[grep("ANNOgesic",nodes)]
  modules=rep(mods,length(sRNA_nodes))
  temp_sRNA=c(temp_sRNA,sRNA_nodes)
  temp_modules=c(temp_modules,modules)
}
sRNA_module_df=data.frame(module=temp_modules)
rownames(sRNA_module_df)=temp_sRNA
#Var1        <-unique(sRNA_module_df[,1])
#names(Var1) <- unique(sRNA_module_df[,1])
#anno_colors <- list(module = Var1)
png(file="~/Desktop/phd/programming/Integrated_analysis/gene_heatmap/selected_sRNA_DESeq_LFC.png", width = 1200, height = 1100)
range=abs(max(table_LFC))
pheatmap(table_LFC %>% select(selected_cond) %>%
           filter(grepl("ANNOgesic",row.names(table_LFC))),show_rownames=T,show_colnames =T,color = colorRampPalette(c("blue","white", "red"))(100),
         breaks = seq(-range, range, length.out = 100))
# color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), breaks=seq(-7.5, 7.5, by = 1))
#pheatmap(log2(TPM+1),annotation_col =cond_df,show_rownames=F,show_colnames =F,col=brewer.pal(9,"Blues"))
dev.off()

png(file="~/Desktop/phd/programming/Integrated_analysis/gene_heatmap/ABC_DESeq_LFC.png", width = 1200, height = 1100)
pheatmap(table_LFC %>% select(selected_cond) %>%
           filter(row.names(table_LFC) %in% ABC_gene_name),show_rownames=T,show_colnames =T,color = colorRampPalette(c("blue","white", "red"))(100),
         breaks = seq(-range, range, length.out = 100))
# color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), breaks=seq(-7.5, 7.5, by = 1))
#pheatmap(log2(TPM+1),annotation_col =cond_df,show_rownames=F,show_colnames =F,col=brewer.pal(9,"Blues"))
dev.off()

png(file="~/Desktop/phd/programming/Integrated_analysis/gene_heatmap/ABC_iron_DESeq_LFC.png", width = 1200, height = 1100)
pheatmap(table_LFC %>% select(selected_cond) %>%
           filter(row.names(table_LFC) %in% row.names(ABC_df_iron)),show_rownames=T,show_colnames =T,annotation_row = iron_system_df,color = colorRampPalette(c("blue","white", "red"))(100),
         breaks = seq(-range, range, length.out = 100))
#pheatmap(as.data.frame(vst_counts) %>% filter(row.names(vst_counts) %in% row.names(ABC_df_iron)),annotation_col =cond_df, annotation_row = iron_system_df,show_rownames=T,show_colnames =T,col=brewer.pal(9,"Blues"))

# color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), breaks=seq(-7.5, 7.5, by = 1))
#pheatmap(log2(TPM+1),annotation_col =cond_df,show_rownames=F,show_colnames =F,col=brewer.pal(9,"Blues"))
dev.off()
#all_sRNA_meta=as.data.frame(read_delim("co_expression_metatable_sRNA.tsv",delim = "\t"))
#sRNA_module_df=all_sRNA_meta %>% filter(grepl("ANNOgesic",fromNode))
#select pairwise compairsons of interest
#create Jenna's binary table
#for differentially expressed genes, include LGC
#create heatmap by LGC

#========================================
#pick out all DEG for Li_ANNOgesic_110
DEG_meta=read.csv("Cytoscape_DESeq/all_sRNA_mRNA_DEG.csv")
selected_sRNA="Li_ANNOgesic_110"
selected_df=DEG_meta %>% filter(fromNode==selected_sRNA | toNode==selected_sRNA) 
DE_cond=as.character(unique(selected_df$sample))
DE_cond_filter=DE_cond#[!grepl("vs_starv|vs_nacl|37_ES_vs|37_M_vs",DE_cond)]
selected_mod=as.character(unique(selected_df$module))
selected_gene_df=DEG_meta %>% filter(module==selected_mod) %>% select(fromNode, toNode)
selected_gene=unique(c(as.character(selected_gene_df[,1]),as.character(selected_gene_df[,2])))
png_name=paste0("~/Desktop/phd/programming/Integrated_analysis/gene_heatmap/",selected_sRNA,"_DESeq_LFC.png")
png(file=png_name, width = 1200, height = 1100)
pheatmap(table_LFC %>% select(DE_cond_filter) %>%
           filter(row.names(table_LFC) %in% selected_gene),show_rownames=T,show_colnames =T,color = colorRampPalette(c("blue","white", "red"))(100),
         breaks = seq(-range, range, length.out = 100))
# color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), breaks=seq(-7.5, 7.5, by = 1))
#pheatmap(log2(TPM+1),annotation_col =cond_df,show_rownames=F,show_colnames =F,col=brewer.pal(9,"Blues"))
dev.off()
