

#Load packages
library(openxlsx)
library(tidyverse)
library(annotables)
library(survival)
library(glmnet)



##-----Read in data------##

matrix_zscore <- read.delim("matrix_ztransform.txt")
patient_data <- read.xlsx('patient_data_pathway.xlsx')

#Read in differentially expressed genes (DEGs)
# Comparison 1: Genes that were differentially expressed in malignant only sympathoblasts with activated pathway vs inactivated pathway IN BOTH JANSKY AND KILDISUITE DATASETS
# Comparison 2: Genes that were differentially expressed in all sympathoblasts with activated pathway vs inactivated pathway IN BOTH JANSKY AND KILDISUITE DATASETS
degs <- read.table("Data/DEGs_Analysis1.txt")

#Log2 transform 'Fold change' column
log2FC_norm <- unlist(log(degs[14], 2))
degs <- add_column(degs, col = log2FC_norm, .after = "norm_foldChange")
names(degs)[names(degs) == "col"] <- "log2FC_norm"

#Filter by values with log2FC_norm of over magnitude 1.5 
genes <- rownames(degs[degs$pvalue.adj.FDR < 0.05 & (degs$log2FC_norm > 1.5 | degs$log2FC_norm < -1.5), ]) 

rm(degs,log2FC_norm)


##------------------Annotating the genes--------------------##


#Map hgnc symbols of the 'common genes' (identified from differential expression) to the RefSeqIDs

mart<- biomaRt::useMart(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl')

#searchFilters(mart = mart, pattern = "hgnc")  use to look at filter options
refseq_mapping <- biomaRt::getBM(attributes = c("refseq_mrna","hgnc_symbol"), 
                                 filters="hgnc_symbol", #change filter depending on output
                                 values = genes, 
                                 mart=mart)

#Using the annotables() package, add additional information
gene_info <- left_join(refseq_mapping %>% dplyr::select(refseq_mrna, symbol = hgnc_symbol), grch38)

#There are many RefSeqIDs for the same gene name, so we omit the RefSeqID columns with gaps and retain only the FIRST TRANSCRIPT VARIANT
gene_info <-gene_info[!(gene_info$refseq_mrna==""),]  #remove rows with missing RefSeqIDs
gene_info2 <- gene_info[!duplicated(gene_info$symbol),] #remove duplicates
all_genes <- (matrix_zscore$RefSeqID)
gene_info3 <- gene_info2[gene_info2$refseq_mrna %in% all_genes, ] #removes the genes that are not in the patient gene expression matrix data
rm(gene_info,gene_info2,all_genes,refseq_mapping,mart)


#Subset the gene expression matrix so that only the genes you are interested in remain
matrix_subset <- matrix_zscore[1,]   #This line is necessary so that new lines can be added (first line is overwritten)
#Subset original expression matrix for only the top differentially expressed genes
for(i in 1:length(rownames(gene_info3)) ){
  val <- gene_info3$refseq_mrna[i]
  matrix_subset[i,] <- matrix_zscore[matrix_zscore$RefSeqID==val,]
}
rm(val,i)


#Transpose the dataset so that the patient data can be joined onto it
matrix_subset <- matrix_subset %>% tidyr::gather(title, Expression,-RefSeqID) %>%
  tidyr::spread(RefSeqID, Expression)

data <- left_join(patient_data, matrix_subset)




##-----------------Univariate Coxpxh survival analysis--------------------##


#Can determine number of censored patients using the survival package (marked with a +) 
censor <- Surv(data$osday, data$b_os_all)  #393 cases censored, 105 not censored
count <- str_count(censor,"\\+")
table(count)

censor <- Surv(data$efsday, data$a_efs_all) #315 cases censored, 183 not censored
count <- str_count(censor,"\\+")
table(count)

rm(censor,count)





#Perform survival analysis by performing univariate Coxph regression in a loop, removing non-significant genes that are not associated with survival
results <- data.frame(gene=NA, pval=NA, hr=NA)
results <- data.frame(gene=character(),pval=numeric(),hr=numeric())
results[1:length(rownames(gene_info3)),] <- NA

for(i in 1:length(rownames(gene_info3)) ){
  gene <- gene_info3$refseq_mrna[i]
  gene_name <- gene_info3$symbol[i]
  
  fit <- coxph(Surv(efsday, a_efs_all)~eval(as.name(paste(gene))), data=data)
  summcph <- summary(fit)
  hr <- summcph$coefficients[2]
  result <- summcph$coefficients[5]
  
  results$hr[i] <- hr
  results$pval[i] <- result
  results$gene[i] <- gene_name
}
sig_genes <- results[results$pval < 0.05,] 
sig_genes <- sig_genes$gene
rm(genes,fit,result,results,summcph,gene,gene_name,hr,i)

#Running univariate Coxph on one gene
fit <- coxph(Surv(osday, b_os_all)~ NM_001099645, data=data)
summary(fit)

#Subset the matrix to only the significant genes
gene_info4 <- gene_info3[gene_info3$symbol  %in% sig_genes,]


#Subset the gene expression matrix so that only the genes you are interested in remain
matrix_subset <- matrix_zscore[1,]   #This line is necessary so that new lines can be added (first line is overwritten)
#Subset original expression matrix for only the top differentially expressed genes
for(i in 1:length(rownames(gene_info4)) ){
  val <- gene_info4$refseq_mrna[i]
  matrix_subset[i,] <- matrix_zscore[matrix_zscore$RefSeqID==val,]
}
rm(val,i)


#Transpose the dataset so that you can join the patient data onto it
matrix_subset <- matrix_subset %>% tidyr::gather(title, Expression,-RefSeqID) %>%
  tidyr::spread(RefSeqID, Expression)

data <- left_join(patient_data, matrix_subset)


#Order the column by increasing order 
gene_info4 <- gene_info4[order(gene_info4$refseq_mrna),]
#RefSeqIDs turned to row names
#gene_info4 <- data.frame(gene_info4[,-1], row.names = gene_info4[,1])


##------------------Multivariate Coxph survival analysis--------------------##


#Testing for interaction with the MYCN gene (unfavourable prognostic markers)
fit <- coxph(Surv(osday, b_os_all)~ mycn * NM_001099645, data=data)
summcph <- summary(fit)
as.data.frame(summcph$coefficients)
rm(fit,summcph)


##------------------Lasso regression--------------------##

#Lasso regression is used to construct prognostic risk scores

#Example from the glmnet package of how the input matrix should be formatted
#x and y must be in matrix form
# data(CoxExample)
# x <- CoxExample$x
# y <- CoxExample$y


#Prepare x
#Each row is a patient, each column is a gene
matrix_subset <- data.frame(matrix_subset[,-1], row.names = matrix_subset[,1])
x2 <- matrix_subset

#Prepare y
#Two columns, one with the time until event/censoring and the second column with the indicator of status (whether the event occurred or whether censoring happened)
y2 <- patient_data %>% dplyr::select(title,osday,b_os_all)
y2<- data.frame(y2[,-1], row.names = y2[,1])
str(y2)
y2$b_os_all <- as.factor(y2$b_os_all)
colnames(y2)<- c("time","status")

x2 <- data.matrix(x2)
y2 <- data.matrix(y2)


#First run cv.glmnet() to find the optimal lambda by running cross_validation on 10 folds (split by 10% validation and 90% training set)
lambda_seq <- 10^seq(2, -2, by = -.1)
cv_output <- cv.glmnet(x2, y2, alpha = 1, lambda = lambda_seq, nfolds = 10, family = "cox", type.measure = "C")
plot(cv_output)
c(cv_output$lambda.min,cv_output$lambda.1se)

#Use optimal lambda value determined from cross-validation 
#Possible to use either lambda min or lambda 1se
fit <- glmnet(x2, y2, family = "cox")
plot(fit)
coef(fit, s = cv_output$lambda.1se)

#Obtain the genes with non-zero coefficients
sig_genes2 <- (coef(fit, s = cv_output$lambda.1se))[coef(fit, s = cv_output$lambda.1se)[,1]!= 0,]
sig_genes2 <- as.data.frame(sig_genes2) 

#Genes associated with poor prognosis have a positive value and genes associated with good prognosis have negative values 


##---------Work out prognostic score----------##


#Subset to only the remaining genes
x2_2 <- x2[, unlist(rownames(sig_genes2))]  
#Keep y2 as it is


#For each sample, add together the gene expression value multiplied by the z-score value
#Transpose the dataset
x2_2 <- t(as.data.frame(x2_2))
#Obtain gene name instead of RefSeq ID
names <- gene_info4[gene_info4$refseq_mrna %in% rownames(sig_genes2),]
x2_2 <- data.frame(gene_name = names$symbol, x2_2)

#Create empty row to calculate risk score for each sample
x2_2[nrow(x2_2)+1, ] <- NA



#Calculate prognostic score
for(i in 2:length(x2_2)){  #number of patients  (omit first column with the gene names)
  total <- 1
  for(j in 1:nrow(x2_2)-1){  #number of genes(omit last row which is empty)
    coef <- sig_genes2[j,]  #lasso regression coefficient
    z_score <- x2_2[j,i]
    val <- coef*z_score
    total = sum(total, val)
    print(total)
  }
  x2_2[nrow(x2_2),i] <- total
}


#Create empty row to calculate risk group for each sample
x2_2[nrow(x2_2)+1, ] <- NA

values <- as.numeric(as.character(unlist(x2_2[11,2:length(x2_2)])))

low25 <- quantile(values,0.25) #determine the breakpoint
upper75 <- quantile(values,0.75) #determine the breakpoint
med <- median(values)

#Determine whether patients are high or low risk
for(i in 2:length(x2_2)) {  #number of patients 
  if (x2_2[11,i] < low25  ) {
    x2_2[12,i] <- "Q1"
  } else if (x2_2[11,i] > upper75) {
    x2_2[12,i] <- "Q4"
  } else if (x2_2[11,i] < upper75 & x2_2[11,i] > med) {
    x2_2[12,i] <- "Q3"
  } else {
    x2_2[12,i] <- "Q2"
  }}

rm(i,j,coef,val,total,z_score,low25,med,upper75)

rownames(x2_2)[11] <- "risk_score"
rownames(x2_2)[12] <- "risk_group"


#Turn row names into first column
patient_data <- data.frame(patient_data[,-1], row.names = patient_data[,1])

matrix_subset2 <- as.data.frame(t(x2_2))


