
#This script Z-transforms the data for the GSE62564 gene expression data 
#Then the expression data is concatenated with the clinical data of GSE49711 and GSE62564 from the SEQC cohort

library(openxlsx)
library(matrixStats)
library(tidyverse)



my_matrix <- read.delim("GSE62564_SEQC_NB_RNA-Seq_log2RPM.txt")
patient_data <- read.xlsx('patient_data.xlsx')


##---------------Calculating z-scores for the gene expression matrix-----------##

#The data is already log2 transformed, since the counts are in decimals and are not integers

#Calculating Z-score
#As discussed here at: https://www.biostars.org/p/153013/
#(Value of gene X for sample Y - Mean expression of gene X for all samples) / Standard deviation of gene X

matrix_zscore <- my_matrix
matrix_zscore <- data.frame(matrix_zscore[,-1], row.names = matrix_zscore[,1]) #change col 1 to be the index, so it doesn't affect the calculation of the mean and standard deviation

#Calculate z-score using the matrixStats package 
#As discussed here at:   https://stackoverflow.com/questions/34707527/improving-my-r-code-to-calculate-z-score-of-dataframe
matrix_zscore <- (matrix_zscore-rowMeans(matrix_zscore))/(rowSds(as.matrix(matrix_zscore)))[row(matrix_zscore)]



#------------Convert genes of interest to HGNC symbols--------------##

#Copy index into the first column and reset index to numbers
matrix_zscore <- cbind(RefSeqID = rownames(matrix_zscore), matrix_zscore)
rownames(matrix_zscore) <- 1:nrow(matrix_zscore)

#Obtain RefSeqID of pathway genes
pathway_genes <- c("NTRK1","PTPN6","TP53")
mart <- biomaRt::useMart(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl')

#searchFilters(mart = mart, pattern = "hgnc")  #use to look at filter options
gene_info <- biomaRt::getBM(attributes = c("refseq_mrna","hgnc_symbol"), 
                                 filters="hgnc_symbol", #change filter depending on output
                                 values = pathway_genes, 
                                 mart=mart)

#Using the annotables package, add additional information using dplyr package 
#gene_info <- left_join(refseq_mapping %>% dplyr::select(refseq_mrna, symbol = hgnc_symbol), grch38)

#There are many RefSeqIDs for the same gene name, so omit the RefSeqID columns with gaps
#For simplicity, retain only the first transcript variant
gene_info <-gene_info[!(gene_info$refseq_mrna==""),]
gene_info2 <- gene_info[!duplicated(gene_info$hgnc_symbol),]
rm(gene_info,refseq_mapping,mart)



#------------Using z-scores, determine in which patient samples the pathway is activated--------------##


#Subset the gene expression matrix so that only the genes you are interested in remain
matrix_subset <- matrix_zscore[1,]   #This line is necessary so that new lines can be added  (first line is overwritten)
for(i in 1:length(rownames(gene_info2)) ){
  val <- gene_info2$refseq_mrna[i]
  matrix_subset[i,] <- matrix_zscore[matrix_zscore$RefSeqID==val,]
}
rm(val)

#Add gene symbol name to data using tidyverse package
matrix_subset <- add_column(matrix_subset, Gene = gene_info2$hgnc_symbol, .after = "RefSeqID")


#Calculate where the pathway is activated for samples
#Pathway activated when  NTRK1 > 0 ; TP53 < 0 ; PTPN6 < 0  
matrix_subset[4,] <- "na" #Prepare new column
for(i in 3:length(colnames(matrix_subset)) ){
  if ((matrix_subset[1,i]>0) & (matrix_subset[2,i]<0) & (matrix_subset[3,i]<0)) {
    matrix_subset[4,i] <- "ACTIVATED"} else {
      matrix_subset[4,i] <- "NOT ACTIVATED"
    }
}

#Determine which samples have the pathway activated in long format
matrix_pathway <- t(matrix_subset[4,3:500])
colnames(matrix_pathway) <- "Pathway"


#Concatenate pathway activation information to patient data
patient_data <- cbind(patient_data,matrix_pathway)
#See how many patients have the pathway upregulated
table(patient_data$Pathway)

#Add column that has the patients stratified by age over and under 18 months 
months <- rep(NA, length(rownames(patient_data)))
for(i in 1:length(rownames(patient_data)) ){
  if (patient_data[i,3] > 540) {
    months[i] <- "over"
  } else {
    months[i] <- "under"
  }
}
patient_data <- add_column(patient_data, months = months, .after = "age")

#Write into Excel file
write.xlsx(patient_data,'patient_data_pathway.xlsx',colNames = TRUE)
write.table(matrix_zscore,"matrix_ztransform.txt",sep="\t")
