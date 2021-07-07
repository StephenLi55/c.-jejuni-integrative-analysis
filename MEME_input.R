###extract sequnece to make input for meme suite
library("Biostrings")
library(seqinr)
library(readr)
library(dplyr)
library(readr)
library(GenomicRanges)
#read all_gene_fasta
fa=DNAStringSet(readRNAStringSet("~/Desktop/phd/programming/Integrated_analysis/WGCNA/IntaRNA/all_gene.fasta"))
fa_name=names(fa)
fa_seq=paste(fa)
#for each module, find sequence
setwd("~/Desktop/phd/programming/Integrated_analysis/WGCNA/")
#write into new input fasta. Name by
all_mods= gsub(".txt","",gsub("CytoscapeInput-nodes-","",list.files(pattern = "CytoscapeInput-nodes-")))
for (file in list.files(pattern = "CytoscapeInput-nodes-")){
  input=read_delim(file,delim="\t")
  genes= input$nodeName
  fa_extract=fa[which(fa_name %in% genes)]
  write.fasta(as.list(paste(fa_extract)),gsub(" ","",gsub("\\(","",gsub("\\)","",gsub("\'","",names(fa_extract))))),paste0("~/Desktop/phd/programming/Integrated_analysis/WGCNA/MEME/",file))
}
#do that again for sRNA of interest. Find all its co-expressed targets
#make new file

#============================
###search for promoter motif or just CDS
#sRNA_list=unique(c(as.character(df$gene1[grep("ANNOgesic",df$gene1)]),as.character(df$gene2[grep("ANNOgesic",df$gene2)])))
#sRNA_list=as.character(df$name)

all_genes=read.csv("~/Desktop/phd/programming/Integrated_analysis/gene_annotation.csv")
all_genes$TSS=TSS_aggregated$TSS.pos[match(all_genes$gene_id , TSS_aggregated$Locus.tag_modified)]
annotation_all_genes=annotation[match(all_genes$gene_name,gene_id),]

annotation=as.data.frame(read_delim("~/Desktop/phd/programming/RNA_crosslinking_analysis/Cj_crosslinking/genomeDir/11168_genome_sRNAs_exon_updated_shortened_visualised.gtf", delim="\t",col_names = F))
#annotation_remove_start=annotation %>% filter(grepl("antisense",X9)) %>% select(X4)
gene_id=gsub("gene_id \"","",gsub("\";.*","",annotation[,9]))

sRNA_list=gene_id[grep("ANNOgesic",gene_id)]
anno_sRNA=annotation[gene_id %in% sRNA_list,]
colnames(anno_sRNA)[c(1,4,5,7)]=c("chr","start","end","strand")

TSS_file=read.csv("~/Desktop/phd/programming/Jenna/categorised_TSS.csv")
#TSS_file_interest=TSS_file %>% filter(TSS.type=="Primary" | TSS.type=="Secondary")
TSS_file$Locus.tag_modified=gsub("_int","",gsub("as_","",gsub("_([0-9])","",TSS_file$Locus.tag)))
#TSS_file$Locus.tag_modified=gsub(" ","",gsub("\\(","",gsub("\\)","",gsub("\'","",TSS_file$Locus.tag_modified))))
TSS_file$Locus.tag_modified_name=all_genes$gene_name[match(TSS_file$Locus.tag_modified,all_genes$gene_id)]
TSS_file$Locus.tag_modified_name=gsub(" ","",gsub("\\(","",gsub("\\)","",gsub("\'","",TSS_file$Locus.tag_modified_name))))
#TSS_aggregated=as.data.frame(TSS_file_interest) %>% group_by(Locus.tag_modified) %>% 
  #dplyr::summarise(TSS.pos= paste(TSS.position,collapse = ","))
co_meta=as.data.frame(read_delim("~/Desktop/phd/programming/Integrated_analysis/WGCNA/co_expression_metatable_sRNA.tsv",delim = "\t"))
#TSS_file$Locus.tag_modified_name_module=mat
TSS_file$Locus.tag_modified_module=co_meta$module[match(TSS_file$Locus.tag_modified_name,co_meta$fromNode)]
#now try to match sRNA
anno_sRNA_TSS=anno_sRNA
for (i in 1:nrow(anno_sRNA_TSS)){
  if (anno_sRNA_TSS[i,7]=="+"){
    anno_sRNA_TSS[i,5]=anno_sRNA_TSS[i,4]+5
    anno_sRNA_TSS[i,4]=anno_sRNA_TSS[i,4]-5
  }else{
    anno_sRNA_TSS[i,4]=anno_sRNA_TSS[i,5]-5
    anno_sRNA_TSS[i,5]=anno_sRNA_TSS[i,5]+5
  }
}
anno_sRNA_TSS_grange=makeGRangesFromDataFrame(anno_sRNA_TSS)
cat_TSS_grange=GRanges(seqnames = unique(anno_sRNA$chr),ranges = IRanges(start=TSS_file$TSS.position,end=TSS_file$TSS.position) , strand = as.character(TSS_file$Strand))
#TSS_file$sRNA=match(TSS_file$TSS.position,anno_sRNA$X4)
overlap=findOverlaps(cat_TSS_grange,anno_sRNA_TSS_grange)
TSS_file$sRNA=NA
TSS_file$sRNA[queryHits(overlap)]=sRNA_list[subjectHits(overlap)]
TSS_file$sRNA_module=co_meta$module[match(TSS_file$sRNA,co_meta$fromNode)]
write.csv(TSS_file,"~/Desktop/phd/programming/Integrated_analysis/WGCNA/categorised_TSS_WGCNA.csv")
#annotation
#TSS_file_interest=TSS_file %>% filter(TSS.type=="Primary" | TSS.type=="Secondary" )

#for loop to go through all Cytoscape node

setwd("~/Desktop/phd/programming/Integrated_analysis/WGCNA/")
#write into new input fasta. Name by
all_mods= gsub(".txt","",gsub("CytoscapeInput-nodes-","",list.files(pattern = "CytoscapeInput-nodes-")))
genome=readDNAStringSet("~/Desktop/phd/programming/RNA_crosslinking_analysis/Cj_crosslinking/genomeDir/NCTC11168_genome.fasta.txt")

for (mods in all_mods){
  TSS_file_filter=TSS_file %>% filter((Locus.tag_modified_module==mods & (TSS.type=="Primary" | TSS.type=="Secondary"))| sRNA_module==mods) 
  fa_extract=TSS_file_filter$Promoter.sequence
  
  fa_names=as.character(TSS_file_filter$TSS_identifier)
  sRNA_index=which(!is.na(TSS_file_filter$sRNA))
  fa_names[sRNA_index]=TSS_file_filter$sRNA[sRNA_index]
  write.fasta(as.list(paste(fa_extract)),fa_names,paste0("~/Desktop/phd/programming/Integrated_analysis/WGCNA/MEME/CytoscapeInput-nodes-",mods,"_TSS.fasta"))
  
  #now create fasta files for CDS
  node_file=read_delim(paste0("~/Desktop/phd/programming/Integrated_analysis/WGCNA/CytoscapeInput-nodes-",mods,".txt"),delim = "\t")
  annotation_nodes=annotation[match(node_file$nodeName,gene_id),]
  gene_id_modes=gene_id[match(node_file$nodeName,gene_id)]
  sequence=substring(genome,annotation_nodes[,4],annotation_nodes[,5])
  rev_index=which(annotation_nodes$X7=="-")
  sequence[rev_index]=reverseComplement(DNAStringSet(sequence[rev_index]))
  write.fasta(as.list(paste(sequence)),gene_id_modes,paste0("~/Desktop/phd/programming/Integrated_analysis/WGCNA/MEME/CytoscapeInput-nodes-",mods,"_CDS.fasta"))
  
  
}