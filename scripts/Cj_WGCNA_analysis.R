library(grid)
library('ggplot2')
library('knitr')
library('limma')
library('reshape2')
library('RColorBrewer')
library('WGCNA')
library(readr)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(doParallel)
##=================================
#write a function which remove genes whose counts are consistently low (for example, removing all features that have a count of less than say 1th decile in more than 90% of the samples)
data_filtering=function(raw_counts){
  low_reads_index=apply(raw_counts,2, function(x) which(ntile(x,10)==1))
  low_reads_table=table(as.numeric(low_reads_index))
  filter_index=which(as.numeric(low_reads_table)>0.9*ncol(raw_counts))
  rows_to_remove=as.numeric(names(low_reads_table)[filter_index])
  return(rows_to_remove)
}
source("~/Desktop/phd/programming/Integrated_analysis/WGCNA/csuWGCNA-master/source code/Hadjacency.r")
source("~/Desktop/phd/programming/Integrated_analysis/WGCNA/csuWGCNA-master/source code/hpickSoftThreshold.r")

#=================================
#data input
#=================================

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

png(file="~/Desktop/phd/programming/Integrated_analysis/heaptmap_raw_cond.png", width = 1200, height = 1100)
#heatmap(cor(raw_counts),RowSideColors = cond_colors,ColSideColors = cond_colors,main="Sample correlations (raw)",
#  scale = "none" , labRow=FALSE, labCol = FALSE,col=brewer.pal(9,"Blues"))
#legend(x="topright", legend=as.character(unique(samples$condition_original)), 
#      fill=unique(cond_colors),ncol=4,cex=0.8)
pheatmap(cor(raw_counts),annotation_row =cond_df ,annotation_col =cond_df)
dev.off()
num_Author=nlevels(samples$Author)
#pal=colorRampPalette(brewer.pal(num_conditions,"Set1"))(num_conditions)
pal=colorRampPalette(brewer.pal(num_Author,"Set1"))(num_Author) #maximum number of colour in set 1 is 9
author_colors=pal[as.integer(samples$Author)]
author_df=data.frame(author=samples$Author)
rownames(author_df)=rownames(cor(raw_counts))
#par(mar=c(1,1,1,1))
png(file="~/Desktop/phd/programming/Integrated_analysis/WGCNA/heaptmap_raw_author.png", width = 1200, height = 1100)
#heatmap(cor(raw_counts),RowSideColors = author_colors,ColSideColors = author_colors,main="Sample correlations (raw)",
#       scale = "none" , labRow=FALSE, labCol = FALSE,col=brewer.pal(9,"Blues"))

#legend(x="topright", legend=as.character(unique(samples$Author)), 
#     fill=unique(author_colors),ncol=4,cex=1.2)
pheatmap(cor(raw_counts),annotation_row =author_df ,annotation_col =author_df)

dev.off()

png(file="~/Desktop/phd/programming/Integrated_analysis/WGCNA/heaptmap_log_author.png", width = 1200, height = 1100)
#heatmap(cor(log_counts),RowSideColors = cond_colors,ColSideColors = cond_colors,main="Sample correlations (log transformed)",
#        scale = "none" , labRow=FALSE, labCol = FALSE,col=brewer.pal(9,"Blues"))

#legend(x="topright", legend=as.character(unique(samples$condition_original)), 
#      fill=unique(cond_colors),ncol=4,cex=0.8)
pheatmap(cor(log_counts),annotation_row =author_df ,annotation_col =author_df)
dev.off()
png(file="~/Desktop/phd/programming/Integrated_analysis/WGCNA/heaptmap_log_cond.png", width = 1200, height = 1100)
#heatmap(cor(log_counts),RowSideColors = cond_colors,ColSideColors = cond_colors,main="Sample correlations (log transformed)",
#        scale = "none" , labRow=FALSE, labCol = FALSE,col=brewer.pal(9,"Blues"))

#legend(x="topright", legend=as.character(unique(samples$condition_original)), 
#      fill=unique(cond_colors),ncol=4,cex=0.8)
pheatmap(cor(log_counts),annotation_row =cond_df ,annotation_col =cond_df)
dev.off()


#run hclust on data.
sampleTree = hclust(dist(t(log_counts)), method = "average");
#plot dendrogram, and heatmap
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
png(file="~/Desktop/phd/programming/Integrated_analysis/WGCNA/samples_clustering.png", width = 1200, height = 1100)

plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
dev.off()
#cols_keep <- samples %>% filter(Author == "Jenna" | Author == "Emily") %>% select(sample_id)
#cols_keep <- samples %>% filter(Author == "Jenna") %>% select(sample_id)
#raw_counts=raw_counts[, (colnames(raw_counts) %in% cols_keep[,1])]
#low_count_mask=rowSums(raw_counts) < ncol(raw_counts)
low_count_mask=data_filtering(raw_counts)
sprintf("Removing %d low-count genes.", length(low_count_mask))
#raw_counts_new=raw_counts[-low_count_mask,]
raw_counts_new=raw_counts
#removed_genes=raw_counts[low_count_mask,]
#write.csv(removed_genes,"~/Desktop/phd/programming/Integrated_analysis/WGCNA/removed_genes.csv")
#print(which(low_count_mask==T))
#log transform data
#WGCNA was initially designed for microarray data, not RNA-Seq
#log2 transform data. To make RNA-Seq data look like microarray data

#log_counts=log2(raw_counts_new+1)
log_counts=vst(as.matrix(raw_counts))
x=melt(as.matrix(log_counts))
colnames(x)=c("gene_id","sample","value")
#png(file="~/Desktop/test.png")
ggplot(x,aes(x=value,color=sample))+geom_density()# theme(legend.position = "none") 
#ggsave("~/Desktop/phd/programming/Integrated_analysis/WGCNA/read_distribution_log.png", width = 20, height = 10)
ggsave("~/Desktop/phd/programming/Integrated_analysis/WGCNA/read_distribution_vst.png", width = 20, height = 10)






#identify outliers. Remove them
#cols_keep <- samples %>% filter(Author == "Jenna" | Author == "Emily") %>% select(sample_id)

#filter by mean instead?

#low_count_mask=rowSums(raw_counts[, (colnames(log_counts) %in% cols_keep[,1])]) < ncol(raw_counts[, (colnames(log_counts) %in% cols_keep[,1])])

#sprintf("Removing %d low-count genes (%d remaining).", sum(low_count_mask),sum(!low_count_mask))

#log_counts[, (colnames(log_counts) %in% cols_keep[,1])]

#log_counts=log_counts[, (colnames(log_counts) %in% cols_keep[,1])]

#log_counts=log_counts[which(low_count_mask==F),]
#print(which(low_count_mask==T))
#=================================
#module detection
#=================================
datExpr=t(log_counts)
#calculate soft-threshold power . Scale independence
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
#enableWGCNAThreads()
allowWGCNAThreads()
# Load the data saved in the first part
#make topological overlap matrix
powers = c(c(1:10), seq(from = 12, to=80, by=2))
#powers = c( seq(from = 60, to=100, by=2))
# Call the network topology analysis function
#sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5 , corFnc="bicor", networkType="unsigned")
#sft<-hpickSoftThreshold(datExpr,powerVector = powers, corFnc = bicor, networkType = "hybrid2", verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
png(file="~/Desktop/phd/programming/Integrated_analysis/WGCNA/scale_independence.png", width = 1200, height = 1100)

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.80,col="red")
dev.off()
png(file="~/Desktop/phd/programming/Integrated_analysis/WGCNA/mean_connectivity.png", width = 1200, height = 1100)

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()



softPower = 5;
#adjacency = adjacency(datExpr, power = softPower);


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# Turn adjacency into topological overlap
sim_matrix=bicor(datExpr)
write.csv(sim_matrix,"~/Desktop/phd/programming/Integrated_analysis/WGCNA/sim_matrix.csv")
sim_matrix_p=bicorAndPvalue(datExpr)
write.csv(sim_matrix_p$p,"~/Desktop/phd/programming/Integrated_analysis/WGCNA/sim_matrix_p.csv")
sim_matrix_p.padj = matrix(p.adjust(as.numeric(unlist(sim_matrix_p$p)),"BH"),nrow=nrow(sim_matrix_p$p),ncol=ncol(sim_matrix_p$p))
write.csv(sim_matrix_p.padj,"~/Desktop/phd/programming/Integrated_analysis/WGCNA/sim_matrix_padj.csv")

#sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
##par(cex = 0.6);
#par(mar = c(0,4,2,0))
#png(file="~/Desktop/phd/programming/Integrated_analysis/WGCNA/sim_matrix.png", width = 1200, height = 1100)
##pheatmap::pheatmap(t(sim_matrix),
   #             labRow = NA, labCol = NA,
    #             xlab="Gene",ylab="Gene",
    #             main="Similarity matrix",revC = TRUE,col=brewer.pal(11,"RdBu"), show_rownames=F,show_colnames=F,)
#dev.off()

#adj_matrix<-Hadjacency(datExpr=datExpr,type='hybrid2',power=softPower, corFnc='bicor')
adj_matrix=adjacency.fromSimilarity(sim_matrix,power =softPower)
rm(sim_matrix) #remove varialbes to free up memory
gene_ids=rownames(adj_matrix)
#rownames(adj_matrix)=gene_ids
#colnames(adj_matrix)=gene_ids

write.table(t(adj_matrix),"~/Desktop/phd/programming/Integrated_analysis/WGCNA/adj_mat.tsv",sep="\t")
#sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
#par(cex = 0.6);
#par(mar = c(0,4,2,0))
#png(file="~/Desktop/phd/programming/Integrated_analysis/WGCNA/adj_matrix.png", width = 1200, height = 1100)
#pheatmap::pheatmap(t(adj_matrix),
     #             labRow = NA, labCol = NA,
   #              xlab="Gene",ylab="Gene",
    #              main="Adjacency matrix",revC = TRUE,col=rev(brewer.pal(11,"RdBu")), show_rownames=F,show_colnames=F)
#dev.off()


#TOM = TOMsimilarityFromExpr(datExpr  ,power = softPower,TOMType = "signed",corType="bicor", networkType="unsigned", maxPOutliers = 0.05);
TOM=TOMsimilarity(adj_matrix)
dissTOM = 1-TOM


#=====================================================================================
#
#  Code chunk 5adjacency
#
#=====================================================================================


# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
png(file="~/Desktop/phd/programming/Integrated_analysis/WGCNA/TOM_gene_clustering.png", width = 1200, height = 1100)

plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
dev.off()

#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(12, 9)
png(file="~/Desktop/phd/programming/Integrated_analysis/WGCNA/gene_dendrogram.png", width = 1200, height = 1100)

plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
png(file="~/Desktop/phd/programming/Integrated_analysis/WGCNA/eigengenes_clustering.png", width = 1200, height = 1100)

plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()
#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================



# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;


#=====================================================================================
#
#  Code chunk 10
#
#=====================================================================================
sizeGrWindow(12, 9)
png(file="~/Desktop/phd/programming/Integrated_analysis/WGCNA/gene_dendrogram_merged.png", width = 1200, height = 1100)


#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


#=====================================================================================
#
#  Code chunk 11
#
#=====================================================================================


# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
#save(MEs, moduleLabels, moduleColors, geneTree, file = "FemaleLiver-02-networkConstruction-stepByStep.RData")
#=================================
#Identfy eigengenes. Correlate with conditions
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs= orderMEs(MEs0)
write.csv(MEs,"~/Desktop/phd/programming/Integrated_analysis/WGCNA/moduleEigengenes.csv")
#moduleTraitCor = cor(MEs, datTraits, use = "p");
#moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
png(file="~/Desktop/phd/programming/Integrated_analysis/WGCNA/eigengenes_heatmap.png", width = 1200, height = 1100)
pheatmap::pheatmap(MEs,
                   labRow = NA, labCol = NA,
                   xlab="Gene",ylab="Gene",
                   main="eigengenes heatmap",revC = TRUE,col=brewer.pal(9,"RdBu"))
dev.off()


#plot cluster dendrogram using dynamic tree cut
#=================================
#visualise similarity matrix and adjacency matrix
#=================================

#cordist = function(dat){
 # cor_matrix = cor(t(dat))
 # dist_matrix=as.matrix(dist(dat,diag = TRUE,upper=TRUE))
 # dist_matrix=log1p(dist_matrix)
 # dist_matrix=1-(dist_matrix/max(dist_matrix))
  
#  sign(cor_matrix)*((abs(cor_matrix)+dist_matrix)/2)
#}
#sim_matrix=cordist(log_counts)
#write.csv(sim_matrix,"~/Desktop/phd/programming/Integrated_analysis/sim_matrix.csv")
#take a subsample, just to save time
#heatmap_indices=sample(nrow(sim_matrix),500)
#png(file="~/Desktop/phd/programming/Integrated_analysis/sim_matrix.png", width = 1200, height = 1100)
#pheatmap::pheatmap(t(sim_matrix),
   #                labRow = NA, labCol = NA,
  #                 xlab="Gene",ylab="Gene",
  #                 main="Similarity matrix \n (by Pearson correlation coefficient and Euclidean distance)",revC = TRUE,col=brewer.pal(11,"RdBu"), show_rownames=F,show_colnames=F,)
#dev.off()

#now construct adjacency matrix by using power transformation
#adj_matrix=adjacency.fromSimilarity(sim_matrix,power = 22)
#sim_matrix=bicor(datExpr)
#rm(sim_matrix) #remove varialbes to free up memory
#gc()
#gene_ids=rownames(adj_matrix)
#rownames(adj_matrix)=gene_ids
#colnames(adj_matrix)=gene_ids
#write.table(t(adj_matrix),"~/Desktop/phd/programming/Integrated_analysis/adj_mat.tsv",sep="\t")
#png(file="~/Desktop/phd/programming/Integrated_analysis/adj_matrix.png", width = 1200, height = 1100)
#pheatmap::pheatmap(t(adj_matrix),
 #                  labRow = NA, labCol = NA,
  #                 xlab="Gene",ylab="Gene",
 #                  main="Adjacency matrix",revC = TRUE,col=rev(brewer.pal(11,"RdBu")), show_rownames=F,show_colnames=F)
#dev.off()




#=================================
#Find key modules. enrichment analysis
# Read in the probe annotation
annot = read.csv(file = "~/Desktop/phd/programming/GeneAnnotation.csv");
# Match probes in the data set to the probe IDs in the annotation file
probes = names(datExpr)
probes2annot = match(probes, annot$substanceBXH)
# Get the corresponding Locuis Link IDs
allLLIDs = annot$LocusLinkID[probes2annot];
# $ Choose interesting modules
intModules = c("brown", "red", "salmon")
for (module in intModules)
{
  # Select module probes
  modGenes = (moduleColors==module)
  # Get their entrez ID codes
  modLLIDs = allLLIDs[modGenes];
  # Write them into a file
  fileName = paste("LocusLinkIDs-", module, ".txt", sep="");
  write.table(as.data.frame(modLLIDs), file = fileName,
              row.names = FALSE, col.names = FALSE)
}
# As background in the enrichment analysis, we will use all probes in the analysis.
fileName = paste("LocusLinkIDs-all.txt", sep="");
write.table(as.data.frame(allLLIDs), file = fileName,
            row.names = FALSE, col.names = FALSE)
GOenr = GOenrichmentAnalysis(moduleColors, allLLIDs, organism = "mouse", nBestP = 10);
tab = GOenr$bestPTerms[[4]]$enrichment
names(tab)
write.table(tab, file = "GOEnrichmentTable.csv", sep = ",", quote = TRUE, row.names = FALSE)

#do similar thing for Kegg

#output gene lists for genes from the module of interest. Input into DAVID/Kegg

#=================================
#visualise gene network
#=================================
#make heatmap for TOM
#visualise eigengene network
#make Eigengene adjacency heatmapo


#=================================
#integrate with DE values
#=================================

#=================================
#export to cytoscape
#=================================
# Recalculate topological overlap if needed
#TOM = TOMsimilarityFromExpr(datExpr,TOMType = "signed", networkType = "signed"  ,power = 22);
#TOM = TOMsimilarityFromExpr(datExpr  ,power = softPower,TOMType = "signed",corType="bicor", maxPOutliers = 0.05);
#TOM = TOMsimilarityFromExpr(datExpr  ,power = softPower,TOMType = "unsigned",corType="bicor", networkType="unsigned", maxPOutliers = 0.05);

# Read in the annotation file
#annot = read.csv(file = "GeneAnnotation.csv");
# Select modules
all_mods=unique(moduleColors)
setwd("~/Desktop/phd/programming/Integrated_analysis/WGCNA/")
system("rm VisANTInput-*")
system("rm CytoscapeInput*")
for (modules in all_mods){
 # modules = c("green");
  # Select module probes
  probes = colnames(datExpr)
  inModule = is.finite(match(moduleColors, modules));
  modProbes = probes[inModule];
  modGenes = modProbes
  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule];
  nTop = 30;
  IMConn = softConnectivity(datExpr[, modProbes]);
  top = (rank(-IMConn) <= nTop)
 # vis = exportNetworkToVisANT(modTOM[top, top],
              #                file = paste("VisANTInput-", modules, "-top30.txt", sep=""),
               #               weighted = TRUE,
                  #            threshold = 0,
                 #             probeToGene = data.frame(probes, probes))
  
  dimnames(modTOM) = list(modProbes, modProbes)
  # Export the network into edge and node list files Cytoscape can read
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE,
                               #  threshold = 0.02,
                                 threshold = 0.02,
                                 nodeNames = modProbes,
                                 altNodeNames = modGenes,
                                 nodeAttr = moduleColors[inModule])
  
}
save(datExpr,MEs, moduleLabels, moduleColors, geneTree,TOM ,all_mods,file = "WGCNA.RData")


#==================================

lnames = load(file = "WGCNA.RData");
#The variable lnames contains the names of loaded variables.
lnames
traitData = read.csv("~/Desktop/phd/programming/Integrated_analysis/WGCNA/samples_trait.csv");
dim(traitData)
names(traitData)=gsub("X","",names(traitData))
names(traitData)=gsub("5.","5%",names(traitData))
# remove columns that hold information we do not need.
# allTraits = traitData[, -c(31, 16)];
# allTraits = allTraits[, c(2, 11:36) ];
allTraits = traitData
dim(allTraits)
names(allTraits)
# Form a data frame analogous to expression data that will hold the clinical traits.
Samples = rownames(datExpr);
traitRows = match(Samples, allTraits$sample_id);
datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 1];

collectGarbage();



nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
#png("~/Desktop/phd/programming/Integrated_analysis/WGCNA/Module-trait_relationships.png", width = 550, height = 400)
png("~/Desktop/phd/programming/Integrated_analysis/WGCNA/Module-trait_relationships.png")
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = gsub("ME","",names(MEs)),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"),
               cex.lab.y=0.7)
dev.off()
# Re-cluster samples
png("~/Desktop/phd/programming/Integrated_analysis/WGCNA/sample_tree2.png",height=1000, width=500)

sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()


#now check the gene membership
#======================
#now check the gene membership
system("rm -r module_membership/*")
system("mkdir module_membership")
for (module in  modNames){
  
  system(paste0("mkdir module_membership/",module,"/"))
}
for (trait_selected in colnames(datTraits)){
  trait_selected_df = datTraits %>% dplyr::select(trait_selected)
  names(trait_selected_df ) = trait_selected
  # names (colors) of the modules
  modNames = substring(names(MEs), 3)
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  geneTraitSignificance = as.data.frame(cor(datExpr, trait_selected_df, use = "p"));
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
  names(geneTraitSignificance) = paste("GS.", names(trait_selected_df), sep="");
  names(GSPvalue) = paste("p.GS.", names(trait_selected_df), sep="")
  #module = "greenyellow"
  for (module in  modNames){
    
    column = match(module, modNames);
    moduleGenes = moduleColors==module;
    sizeGrWindow(7, 7);
    par(mfrow = c(1,1));
    trait_selected=gsub("5%","5percents",trait_selected)
    png(paste0("~/Desktop/phd/programming/Integrated_analysis/WGCNA/module_membership/",module,"/module_membership_",trait_selected,".png"),height=300, width=600)
    
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                       abs(geneTraitSignificance[moduleGenes, 1]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = paste0("Gene significance for ",trait_selected),
                       main = paste("Module membership vs. gene significance\n"),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 0.8, col = "black")
    dev.off()
    df=cbind(abs(geneModuleMembership[moduleGenes, column]),abs(geneTraitSignificance[moduleGenes, 1]))
    colnames(df)=c("Module membership","Gene-trait significance")
    rownames(df)=rownames(geneModuleMembership[moduleGenes, ])
    write.csv(df,paste0("~/Desktop/phd/programming/Integrated_analysis/WGCNA/module_membership/",module,"/module_membership_",trait_selected,".csv"))
  }
  
}

#=================
#=================================
#visualise eigenegene adjacency
# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
#TOM = TOMsimilarityFromExpr(datExpr  ,power = softPower,TOMType = "signed",corType="bicor", maxPOutliers = 0.05);
#TOM = TOMsimilarityFromExpr(datExpr  ,power = softPower,TOMType = "unsigned",corType="bicor", networkType="unsigned", maxPOutliers = 0.05);

# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
par(mfrow = c(1,1));
sizeGrWindow(9,9)
png(file="~/Desktop/phd/programming/Integrated_analysis/WGCNA/network_heatmap.png", width = 1200, height = 1100)

TOMplot(plotTOM, geneTree, moduleColors, main = "Topological Overlap Matrix, all genes")
dev.off()

# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# Isolate weight from the clinical traits
system("rm -r eigengene_trait/")
system("mkdir eigengene_trait")
for (Selected_trait in colnames(datTraits)){
  
  #assign("x",Selected_trait)
  Selected_trait_df = datTraits %>% dplyr::select(Selected_trait)
  names(Selected_trait_df) = Selected_trait
  # Add the weight to existing module eigengenes
  MET = orderMEs(cbind(MEs, Selected_trait_df))
  # Plot the relationships among the eigengenes and the trait
  sizeGrWindow(5,7.5);
  par(cex = 0.9)
  plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                        = 90)
  # Plot the dendrogram
  sizeGrWindow(6,6);
  Selected_trait=gsub("5%","5percents",Selected_trait)
  png(paste0("~/Desktop/phd/programming/Integrated_analysis/WGCNA/eigengene_trait/Eigengene_dendrogram_",Selected_trait,".png"))
  par(cex = 1.0)
  plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                        plotHeatmaps = FALSE)
  # Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
  dev.off()
  png(paste0("~/Desktop/phd/programming/Integrated_analysis/WGCNA/eigengene_trait/Eigengene_adjacency_heatmap_",Selected_trait,".png"),width = 700, height=700)
  par(cex = 1.0)
  plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(6,8,2,2),
                        plotDendrograms = FALSE, xLabelsAngle = 90, cex.lab.x = 0.8,cex.lab.y = 0.8)
  dev.off()
}
#Selected_trait="iron_lim_M"

#====================
rm(sim_matrix)
rm(sim_matrix_p)
rm(sim_matrix_p.padj)
rm(adj_matrix)
rm(datExpr)
rm(plotTOM)
rm(dissTOM)
rm(log_counts)
rm(raw_counts)
rm(raw_counts_new)
rm(TOM)
#connect Cytoscape output edge with simialarity, adjacency and topology
#read all 3 files
#define a function that take 2 genes, and return their value
sim=read.csv("sim_matrix.csv")[,-1]
rownames(sim)=gene_ids
colnames(sim)=gene_ids
sim_p_list=bicorAndPvalue(datExpr)
sim_p=sim_p_list$p
rm(sim_p_list)
rownames(sim_p)=gene_ids
colnames(sim_p)=gene_ids
sim_padj=matrix(p.adjust(as.numeric(unlist(sim_p)),"BH"),nrow=nrow(sim_p),ncol=ncol(sim_p))
rownames(sim_padj)=gene_ids
colnames(sim_padj)=gene_ids
adj=as.data.frame(read.table("adj_mat.tsv",sep = "\t"))
rownames(adj)=gene_ids
colnames(adj)=gene_ids
#TOM=TOMsimilarity(sapply(adj,as.numeric))
#rm(adj_matrix)
#top=TOM
#rownames(top)=gene_ids
#colnames(top)=gene_ids
#rm(TOM)
#for loop to go through all cytoscape edge file
#apply the function to all rows
cytoscape_details=function(input_row){
  gene1=input_row[1]
  gene2=input_row[2]
  sim_value=as.numeric(as.data.frame(sim) %>% dplyr::select(gene2) %>% filter(row.names(sim)==gene1))
  sim_p_value=as.numeric(as.data.frame(sim_p) %>% dplyr::select(gene2) %>% filter(row.names(sim_p)==gene1))
  sim_padj_value=as.numeric(as.data.frame(sim_padj) %>% dplyr::select(gene2) %>% filter(row.names(sim_padj)==gene1))
  adj_value=as.numeric(as.data.frame(adj) %>% dplyr::select(gene2) %>% filter(row.names(adj)==gene1))
 # top_value=as.numeric(as.data.frame(top) %>% dplyr::select(gene2) %>% filter(row.names(top)==gene1))
  output=c(input_row,sim_value,sim_p_value,sim_padj_value,adj_value )
  return(output)
}
rm("rm *_adjacency.txt")
files = list.files(pattern="CytoscapeInput-edges")
for (file in files){
  input=as.data.frame(read_delim(file , delim="\t"))
  row_index=match(input[,1],gene_ids)
  col_index=match(input[,2],gene_ids)
  sim_list=c()
  sim_p_list=c()
  sim_padj_list=c()
  adj_list=c()
  for (i in 1:length(row_index)){
    sim_list=c(sim_list,sim[row_index[i],col_index[i]])
    sim_p_list=c(sim_p_list,sim_p[row_index[i],col_index[i]])
    sim_padj_list=c(sim_padj_list,sim_padj[row_index[i],col_index[i]])
    adj_list=c(adj_list,adj[row_index[i],col_index[i]])
  }
  input_added=input %>% mutate(correlation=sim_list,
                               correlation_p=sim_p_list,
                               correlation_padj=sim_padj_list,
                               adjacency=adj_list)
 # colnames(input_added)=c(colnames(input),"correlation","correlation_p","correlation_padj","adjacency")
  name=gsub(".txt","_adjacency.txt",file)
  write.table(input_added,name,sep="\t",quote=F,row.names = F)
  
}
#=================================
source("cytoscape_ABC_transporter.R")
#================================
source("module_adjacency.R")
#=================================
source(file="pairwise_comparison_scatterPlot.R")
#==========================
#combine DESeq output and co-expression output. Check how well they agree
source(file = " WGCNA_RCy3.R")
#==========================
source(file = " WGCNA_heatmap.R")
#==========================
#make another script that compare the gene significance to all traits
#for each sRNA
source(file = "sRNA_gene_significance.R")

#=============================
#match against crosslinking and target prediction results
source(file = "WGCNA_crosslinking.R")
source(file = "WGCNA_target_prediction.R")
#=============================
#annotate my data with DAVID. Visualise the pathway / ontology that has been more enriched
source(file = "WGCNA_databases_annotation.R")
#==================================================
#create input file for MEME analysis
source(file= "MEME_input.R")
#==================================================
#final intergration
source(file= "integrate_sRNA_target.R")
