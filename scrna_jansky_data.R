


#Script I wrote for Jansky dataset that can be run on my laptop (actually turns out the integrated object is small enough to run on my laptop)


library(readxl)
library(scran)
library(scater)
library(bluster)
library(tidyverse)
library(Seurat)
library(openxlsx)
library(pheatmap)
library(RColorBrewer)
library(SeuratData)
library(batchelor)
library(scCustomize)
library(viridis)
library(scales) #shows full palette
library(extrafont)


loadfonts(device = c("all", "pdf", "postscript", "win"), quiet = FALSE) 

windowsFonts("Papyrus" = windowsFont("Papyrus")) #Papyrus is a good font to check with
windowsFonts("Lato" = windowsFont("Lato"))
windowsFonts("Gadugi" = windowsFont("Gadugi"))
windowsFonts("Roboto" = windowsFont("Roboto"))

##----------Load list of marker genes----------##

setwd("/Users/jenna/Downloads/markers")
celltype_markers <- read_excel('cell_markers.xlsx',col_names = TRUE)
mes_adr_markers <- read_excel('mes_adr.xlsx',col_names = TRUE)

celltype_list <- celltype_markers %>%
  gather(key = "celltype", value = "gene")
celltype_list <- na.omit(celltype_list)
celltype_list <- celltype_list$gene
celltype_list <- unique(celltype_list)



##----------Read in the sce object (already annotated)--------------##

setwd("/Users/jenna/Downloads/project")

sce.seurat.named <- readRDS(file = 'jansky_annotated_seurat_2paths')
sce.seurat.named  <- RenameIdents(object = sce.seurat.named , `Fibroblasts` = "Mesenchyme")
sce.seurat.named  <- RenameIdents(object = sce.seurat.named , `Erythrocyte` = "Erythrocytes")


#rearrange the levels so that they fit with the colour scheme I want
levels(sce.seurat.named) <- c("Sympathoblasts", "Chromaffin cells", "Schwann cell precursors", "Adrenal cortex", "Mesenchyme", "Macrophages" , "Endothelium", "Erythrocytes")
table <- sce.seurat.named@meta.data



palette <- c("#F8766D",  "#00B8E7", "#C77CFF", "#E68613", "#CD9600", "#00C19A", "#8494FF", "#FF68A1")             
show_col(palette)





##----------Read in the sce object (not annotated)--------------##

sce <- readRDS(file = 'jansky_integrated')
sce <- as.SingleCellExperiment(sce)

sce$batch <- factor(sce$SampleName)
rownames(colData(sce)) <- make.unique(as.character(rownames(colData(sce))))

#check the labels
sce_label <- as.data.frame(colData(sce))
write.xlsx(sce_label,'plots/coldata2.xlsx',colNames = TRUE)

#I duplicate the ident column so it is not lost later
sce$group <- sce$ident





patientdata <- read_excel('jansky_patientdata.xlsx',col_names = TRUE)

#genes <- as.data.frame(rownames(sce))


##----------Add in patient data--------------##


sce$treatment <- sce$SampleName
sce$risk <- sce$SampleName
sce$INSS <- sce$SampleName
sce$MYCN_status <- sce$SampleName

for(i in 1:28 ){
  sce$treatment[sce$treatment == patientdata[[i,1]]] <- patientdata[[i,2]]
  sce$risk[sce$risk == patientdata[[i,1]]] <- patientdata[[i,3]]
  sce$INSS[sce$INSS == patientdata[[i,1]]] <- patientdata[[i,4]]
  sce$MYCN_status[sce$MYCN_status == patientdata[[i,1]]] <- patientdata[[i,5]]
}  




sce_label <- as.data.frame(colData(sce))
write.xlsx(sce_label,'plots/coldata2.xlsx')


#saveRDS(sce, file = 'jansky_integrated2')


##-----------Dimensionality reduction------##




#Plot TSNE (stochastic)
#Set seed for reproducible analyses 
set.seed(100)
sce <- runTSNE(sce, 
               perplexity = 50, 
               dimred = "PCA",
               n_dimred = 10)


#Plot UMAP (stochastic)
#Set seed for reproducible analyses
set.seed(100)
sce <- runUMAP(sce, 
               dimred = "PCA",
               n_neighbors = 50,
               n_dimred = 10)



#For Seurat object
set.seed(100)
sce <- runUMAP(sce, 
               dimred = "HARMONY",
               n_neighbors = 50,
               n_dimred = 10)

##-----------Clustering------##




#The following plots are used to test how accurate the clustering was
plotSilBeeswarm <- function(silDat){
  silTab <- silDat %>% 
    as.data.frame() %>% 
    mutate(closestCluster = ifelse(width > 0, cluster, other) %>% factor())
  
  plt <- silTab %>% 
    ggplot(aes(x=cluster, y=width, colour=closestCluster)) +
    ggbeeswarm::geom_quasirandom(method="smiley", alpha=0.6) +
    theme_bw()
  
  plt <- scater:::.resolve_plot_colours(plt, silTab$closestCluster, "closestCluster")
  plt
}


plotSilGrid <- function(silDat){
  silDat %>% 
    as.data.frame() %>% 
    mutate(closestCluster = ifelse(width > 0, cluster, other) %>% factor()) %>% 
    count(cluster, closestCluster,  name="olap") %>% 
    group_by(cluster) %>% 
    mutate(total  = sum(olap)) %>% 
    mutate(proportion = olap / total) %>% 
    mutate(proportion = ifelse(cluster == closestCluster, proportion, -proportion)) %>% 
    ggplot(aes(x = cluster, y = closestCluster)) +
    geom_tile(aes(fill = proportion)) +
    geom_text(aes(label = olap), size=5) +
    scale_fill_gradientn(colors = c("#fc8d59", "#ffffbf", "#91cf60"),
                         limits = c(-1, 1)) +
    geom_vline(xintercept=seq(0.5, 30.5, by=1)) +
    geom_hline(yintercept=seq(0.5, 30.5, by=1), colour="lightgrey", linetype=2) +
    guides(fill = "none") +
    theme(
      aspect.ratio = 1,
      panel.background = element_blank())
}



#Cluster with leiden method (stochastic)
set.seed(100)
sce$leiden <- clusterCells(sce, 
                           use.dimred = "PCA", 
                           BLUSPARAM = SNNGraphParam(k = 60, 
                                                     cluster.fun = "leiden"))

set.seed(100)
sce$leiden <- clusterCells(sce, 
                           use.dimred = "HARMONY", 
                           BLUSPARAM = SNNGraphParam(k = 60, 
                                                     cluster.fun = "leiden"))


table <- table(sce$leiden)
write.xlsx(table,'plots/clusters.xlsx',colNames = TRUE)

sil.approx <- approxSilhouette(reducedDim(sce, "PCA"),
                               clusters=sce$leiden)




pdf("plots/leiden.pdf")

table(sce$leiden, sce$batch) %>% 
  pheatmap(cluster_rows = FALSE, cluster_cols = FALSE)

plotSilBeeswarm(sil.approx)

plotSilGrid(sil.approx)

dev.off()







##-----------Conversion to Seurat object------##

sce.seurat <- as.Seurat(sce, counts = "counts", data = "logcounts")

#Set main Seurat label to the cluster numbers
sce.seurat <- SetIdent(sce.seurat, value = "leiden")

##-----------Seurat plots------##

#PCA plot coloured by batch
pdf("plots/PCA_batch.pdf")
DimPlot(sce.seurat, reduction = "PCA", group.by = "batch")
dev.off()

#PCA plot coloured by cluster, not yet annotated
pdf("plots/PCA_cluster.pdf")
DimPlot(sce.seurat, reduction = "PCA", label = TRUE, repel = TRUE)
dev.off()


#tSNE plot coloured by batch
pdf("plots/tSNE_batch.pdf")
DimPlot(sce.seurat, reduction = "TSNE", group.by = "batch")
dev.off()

#tSNE plot coloured by clusters, not yet annotated
pdf("plots/tSNE_cluster.pdf")
DimPlot(sce.seurat, reduction = "TSNE", label = TRUE, repel = TRUE)
dev.off()



#UMAP plot coloured by batch
pdf("plots/UMAP_batch.pdf")
DimPlot(sce.seurat, reduction = "UMAP", group.by = "batch")
dev.off()

#UMAP plot coloured by clusters, not yet annotated
pdf("plots/UMAP_cluster.pdf")
DimPlot(sce.seurat, reduction = "UMAP", label = TRUE, repel = TRUE)
dev.off()

#We can plot by patient but I don't think we have a need to yet
#DimPlot(sce.seurat, reduction = "UMAP", group.by = "label", split.by = "batch")


DimPlot(sce.seurat.named, reduction = "UMAP", label = TRUE, repel = TRUE)



#Plotting cell types all together

n <- 18

pdf("plots/allplots.pdf")

for(i in 1:n ){
  
  celltype <- unlist(celltype_markers[[i]])
  celltype <- na.omit(celltype)
  sce.seurat <- AddModuleScore(sce.seurat,
                               features = list(celltype),
                               name= colnames(celltype_markers[i]))
  
  name <- (paste0(colnames(celltype_markers[i]),"1"))
  print(name)
  
  
  
  #must specify print() or the file will be corrupted
  print(FeaturePlot(sce.seurat,
                    features = name) +
          scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))))
  
  
}

dev.off()


pdf("plots/mes_adr.pdf")

for(i in 1:3 ){
  
  celltype <- unlist(mes_adr_markers[[i]])
  celltype <- na.omit(celltype)
  sce.seurat <- AddModuleScore(sce.seurat,
                               features = list(celltype),
                               name= colnames(mes_adr_markers[i]))
  
  name <- (paste0(colnames(mes_adr_markers[i]),"1"))
  print(name)
  
  
  
  #must specify print() or the file will be corrupted
  print(FeaturePlot(sce.seurat,
                    features = name) +
          scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))))
  
  
}

dev.off()



##-----------Checking the 14 gene signature for ALL cell types-----##

#High expression of these genes associated with poor prognosis:
#  HSD17B12, NXT2, TXNDC5, RACK1, CD147, CCDC125, HAPLN4, HEBP2

#High expression of these genes associated with good prognosis:
#  B3GAT1, FNBP1, IGSF10, IQCE, KCNQ3, TOX2

pal <- viridis(n = 10, option = "D")
#https://samuel-marsh.github.io/scCustomize/reference/FeaturePlot_scCustom.html            #colors_use = viridis_plasma_dark_high



pdf("plots/gene_sig.pdf")


poor_prognosis <- c("HSD17B12", "NXT2", "TXNDC5", "RACK1", "CD147", "CCDC125", "HAPLN4", "HEBP2")
sce.seurat <- AddModuleScore( sce.seurat, features = list(poor_prognosis), name="poor_prognosis")
FeaturePlot_scCustom(sce.seurat, features = "poor_prognosis1", pt.size=0.5,na_color = "lightgray",colors_use = pal) + ggtitle('Poor prognosis genes')



good_prognosis <- c("B3GAT1", "FNBP1", "IGSF10", "IQCE", "KCNQ3", "TOX2")
sce.seurat <- AddModuleScore( sce.seurat.named,features = list(good_prognosis),name="good_prognosis")
FeaturePlot_scCustom(sce.seurat, features = "good_prognosis1", pt.size=0.5,na_color = "lightgray",colors_use = pal) + ggtitle('Good prognosis genes')


dev.off()


#RACK1, CD147 not present
all_genes <- c("NTRK1","TP53","PTPN6","HSD17B12", "NXT2", "TXNDC5", "CCDC125", "HAPLN4", "HEBP2","B3GAT1", "FNBP1", "IGSF10", "IQCE", "KCNQ3", "TOX2")

pdf("plots/gene_sig_individual.pdf") 

for(i in 1:15 ){
  print(FeaturePlot_scCustom(sce.seurat, features = all_genes[i], pt.size=0.5,na_color = "lightgray",colors_use = pal))
}
dev.off()




##------------Assign labels----------------##

table <- sce.seurat.named@meta.data
write.xlsx(table,'plots/label_before_annotation.xlsx',colNames = TRUE)

#These are the labels AFTER cell annotation
labels <- read_excel('jansky_labels.xlsx',col_names = FALSE)
labels<-labels$...1

saveRDS(sce.seurat, file = 'jansky_clustered_unannotated')

sce.seurat.named <- sce.seurat

#moves the label from levels to 'ident'
sce.seurat.named[["ident"]] <- Idents(object = sce.seurat.named)

#renames label
#sce <- RenameIdents(object = sce, `Chromaffin` = "Chromaffin cells")


names(labels) <- levels(sce.seurat.named)
sce.seurat.named <- RenameIdents(sce.seurat, labels)

#UMAP plot with all celltypes annotated
# pdf("plots/annotated_all.pdf")
# plot <- DimPlot(sce.seurat.named, reduction = "UMAP", pt.size = 0.5, cols = palette)  +
#   theme(text=element_text(size=16,  family="Lato")) +
#   xlab('UMAP 1') +  ylab('UMAP 2')
# LabelClusters(plot = plot, id = 'ident', family="Lato")
# dev.off()

plot <- DimPlot(sce.seurat.named, reduction = "UMAP", pt.size = 0.5, cols = palette, raster=FALSE)  +
  theme(text=element_text(size=16,  family="Lato")) +
  xlab('UMAP 1') +  ylab('UMAP 2') + NoLegend()
plot <- LabelClusters(plot = plot, id = 'ident', family="Lato",size = 7) 




saveRDS(sce.seurat.named, file = 'jansky_annotated_seurat')

sce.named <- as.SingleCellExperiment(sce.seurat.named)

sce_label <- as.data.frame(colData(sce.named))
write.xlsx(sce_label,'plots/coldata3.xlsx',colNames = TRUE)




##------Plot cells with upreg and downreg pathway-------##


# upregulated<- WhichCells(object = sce.seurat.named, expression = pathway == TRUE)
# downregulated<- WhichCells(object = sce.seurat.named, expression = pathway_down == TRUE)
# cells <- list(upregulated  = upregulated, downregulated = downregulated)
# Cell_Highlight_Plot(seurat_object = sce.seurat.named, cells_highlight = cells, pt.size = 0.5) + ggtitle('All samples') 



pdf("plots/annotated_split_pathway.pdf")

upregulated<- WhichCells(object = sce.seurat.named, expression = pathway == TRUE)
downregulated<- WhichCells(object = sce.seurat.named, expression = pathway_down == TRUE)
cells <- list(upregulated  = upregulated, downregulated = downregulated)
Cell_Highlight_Plot(seurat_object = sce.seurat.named, cells_highlight = cells, pt.size = 0.5) + ggtitle('All samples') + NoLegend()

sce.nb <- (sce.seurat.named[,sce.seurat.named$batch == 'nb'])
#UMAP plot with only neuroblastoma samples
upregulated<- WhichCells(object = sce.nb, expression = pathway == TRUE)
downregulated<- WhichCells(object = sce.nb, expression = pathway_down == TRUE)
cells <- list(upregulated  = upregulated, downregulated = downregulated)
Cell_Highlight_Plot(seurat_object = sce.nb, cells_highlight = cells, pt.size = 0.5) + ggtitle('14 NB samples') + NoLegend()


sce.fetal <- (sce.seurat.named[,sce.seurat.named$batch == 'control'])
#UMAP plot with only fetal controls
upregulated<- WhichCells(object = sce.fetal, expression = pathway == TRUE)
downregulated<- WhichCells(object = sce.fetal, expression = pathway_down == TRUE)
cells <- list(upregulated  = upregulated, downregulated = downregulated)
Cell_Highlight_Plot(seurat_object = sce.fetal, cells_highlight = cells, pt.size = 0.5) + ggtitle('17 Fetal adrenal glands') + NoLegend()


dev.off()



pdf("plots/annotated_split_pathway2.pdf")

sce.old <- (sce.seurat.named[,sce.seurat.named$over_18_months == 'yes'])
#UMAP plot with only neuroblastoma samples
upregulated<- WhichCells(object = sce.old, expression = pathway == TRUE)
downregulated<- WhichCells(object = sce.old, expression = pathway_down == TRUE)
cells <- list(upregulated  = upregulated, downregulated = downregulated)
Cell_Highlight_Plot(seurat_object = sce.old, cells_highlight = cells, pt.size = 0.5) + ggtitle('over 18 months') 


sce.young <- (sce.seurat.named[,sce.seurat.named$over_18_months == 'no'])
#UMAP plot with only fetal controls
upregulated<- WhichCells(object = sce.young, expression = pathway == TRUE)
downregulated<- WhichCells(object = sce.young, expression = pathway_down == TRUE)
cells <- list(upregulated  = upregulated, downregulated = downregulated)
Cell_Highlight_Plot(seurat_object = sce.young, cells_highlight = cells, pt.size = 0.5) + ggtitle('under 18 months') 


dev.off()



##-----------Plot UMAP plots split by batches----------------##


sce_label <- as.data.frame(colData(sce))
sce_label <- sce_label[1:100,]

pdf("plots/annotated_split.pdf")

sce.nb <- (sce.seurat.named[,sce.seurat.named$batch == 'nb'])
#UMAP plot with only neuroblastoma samples
DimPlot(sce.nb, reduction = "UMAP", label = TRUE, repel=TRUE, pt.size = 0.5) + ggtitle('14 Neuroblastoma samples')

sce.control <- (sce.seurat.named[,sce.seurat.named$batch == 'control'])
#UMAP plot with only the controls
DimPlot(sce.control, reduction = "UMAP", label = TRUE, repel=TRUE, pt.size = 0.5)  + ggtitle('17 Fetal adrenal glands')

dev.off()




#another method for plotting the different batches i.e. control, nb

sce.seurat.named[["plot"]] <- Idents(object = sce.seurat.named)
table <- sce.seurat.named@meta.data


#i have only the cells of the neuroblastoma cells coloured, the rest are in grey?
#must turn metadata into regular table and then back again
table$plot <- as.character(table$plot)
table[table$batch == "control","plot"] <- "Controls"
sce.seurat.named@meta.data <- table
unique(table$plot) #see which cells are here

sce.seurat.named <- SetIdent(sce.seurat.named, value = "plot")

#I order so that the legend looks nicer (ordered by colour transition - makes legend look better)
levels(sce.seurat.named) <- c("Sympathoblasts","Chromaffin cells", "Adrenal cortex", "Mesenchyme",  "Macrophages", "Endothelium", "Schwann cell precursors", "Controls")
palette <- c("#F8766D", "#00B8E7", "#E68613", "#CD9600", "#00C19A", "#8494FF","#C77CFF",  "#D3D3D3")             
show_col(palette)

plot <- DimPlot(sce.seurat.named, reduction = "UMAP", pt.size = 0.5, cols = palette, raster=FALSE)  +
  theme(text=element_text(size=16,  family="Lato")) +
  xlab('UMAP 1') +  ylab('UMAP 2') + NoLegend()
plot <- LabelClusters(plot = plot, id = 'ident', family="Lato",size = 7) 
ggsave(plot, filename = "plots/annotated_all.pdf", device = cairo_pdf, 
       width = 8, height = 9, units = "in")




#Alternative with ggsave lets me save with custom font
sce.control <- (sce.seurat.named[,sce.seurat.named$batch == 'control'])
plot <- DimPlot(sce.nb, reduction = "UMAP", pt.size = 0.5, cols = palette, raster=FALSE)  +
  theme(text=element_text(size=16,  family="Lato")) +
  xlab('UMAP 1') +  ylab('UMAP 2') + NoLegend()
#plot <- LabelClusters(plot = plot, id = 'ident', family="Lato",size = 7) 
ggsave(plot, filename = "plots/annotated_all.pdf", device = cairo_pdf, 
       width = 8, height = 9, units = "in")





##-----------Plot UMAP plots by malignant cells----------------##


sce_malignant <- subset(x = sce.seurat.named, malignant == TRUE)

pdf("plots/malignantcells.pdf")
DimPlot(sce.malignant, reduction = "UMAP", label = TRUE, repel=TRUE, pt.size = 0.5)
dev.off()

sce_malignant_recluster <- as.SingleCellExperiment(sce_malignant)

#account for batch effects once again, after subsetting
gene_var <- modelGeneVar(sce_malignant_recluster, block=sce_malignant_recluster$patientID)
sce_malignant_recluster <- multiBatchNorm(sce_malignant_recluster, batch=sce_malignant_recluster$patientID)
sce_malignant.mnn <- fastMNN(sce_malignant_recluster, batch=sce_malignant_recluster$patientID, subset.row=getTopHVGs(gene_var))
reducedDim(sce_malignant_recluster, "corrected") <- reducedDim(sce_malignant.mnn, "corrected")




set.seed(100)
sce.neuro <- runUMAP(sce_malignant_recluster, 
                     dimred = "corrected",
                     n_neighbors = 50,
                     n_dimred = 10)


set.seed(100)
sce.neuro$leiden <- clusterCells(sce_malignant_recluster, 
                                 use.dimred = "corrected", 
                                 BLUSPARAM = SNNGraphParam(k = 20, 
                                                           cluster.fun = "leiden"))

table(sce_malignant_recluster$leiden)
table <- as.data.frame(colData(sce_malignant_recluster))


sce.malig.seurat <- as.Seurat(sce_malignant_recluster, counts = "counts", data = "logcounts")
#Set main Seurat label to the cluster numbers
sce.malig.seurat <- SetIdent(sce_malignant_recluster, value = "labels")


pdf("plots/malignant_cells_recluster.pdf")
DimPlot(sce_malignant_recluster, reduction = "UMAP", label = TRUE, repel = TRUE)
dev.off()


#Plotting cell types all together


n <- c(1,2,3,4,5,6,7,20,21,22)

pdf("plots/individual_malig.pdf")

for(i in n){
  
  celltype <- unlist(celltype_markers[[i]])
  celltype <- na.omit(celltype)
  sce_malignant_recluster <- AddModuleScore(sce_malignant_recluster,
                               features = list(celltype),
                               name= colnames(celltype_markers[i]))
  
  name <- (paste0(colnames(celltype_markers[i]),"1"))
  print(name)
  
  
  
  #must specify print() or the file will be corrupted
  print(FeaturePlot(sce_malignant_recluster,
                    features = name) +
          scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))))
  
  
}

dev.off()
