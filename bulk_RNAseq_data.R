
#This script uses the bulk RNA-seq datasets GSE49711 and GSE62564 from the SEQC cohort, both of which include the same sample data.
#Both datasets have been concatenated together to integrate the clinical information.

#Load packages
library(GEOquery)
library(tidyverse)
library(stringi)

##----------Obtaining the bulk RNA-seq data----------##

#Get GSE49711 data using GEOquery package
gse <- getGEO("GSE49711")[[1]]   

#Format into dataset and drop columns using tidyverse package
sample_info  <- pData(gse) %>% dplyr::select(title, contains("characteristics"), -characteristics_ch1, -characteristics_ch1.1)

#Rename columns using tidyverse library
sample_info2 <- sample_info  %>% 
  dplyr::rename(sex = characteristics_ch1.2,
                age = characteristics_ch1.3,
                mycn = characteristics_ch1.4,
                high_risk = characteristics_ch1.5,
                inss_stage = characteristics_ch1.6,
                class = characteristics_ch1.7,
                progression = characteristics_ch1.8,
                death_from_disease = characteristics_ch1.9)

sample_info2 <- sample_info2 %>% 
  mutate(sex = gsub("Sex: ", "", sex, fixed=TRUE),
         age = as.numeric(gsub("age at diagnosis: ","", age, fixed=TRUE)),
         mycn = as.factor(gsub("mycn status: ","", mycn, fixed=TRUE)),
         high_risk = as.factor(gsub("high risk: ","", high_risk, fixed=TRUE)),
         inss_stage = as.factor(gsub("inss stage: ","", inss_stage, fixed=TRUE)),
         class = as.factor(gsub("class label: ","", class, fixed=TRUE)),
         progression = as.factor(gsub("progression: ","", progression, fixed=TRUE)),
         death_from_disease = as.factor(gsub("death from disease: ","", death_from_disease, fixed=TRUE)))



#Get GSE62564 data using GEOquery package
gse <- getGEO("GSE62564")[[1]]

sample_info3 <- pData(gse) %>% 
  dplyr::select(title, contains("characteristics"), -characteristics_ch1,-characteristics_ch1.1,-characteristics_ch1.2, -characteristics_ch1.14)

#Gather data using the tidyverse package into long data format
sample_info3 <- sample_info3 %>% gather(key = "variable", value = "value", -title)
sample_info3$col <- 1 #create a new column filled with ones


#This for loop separates the text before the colon to make a new header e.g. Age, and then removes the text before the colon to isolate the value into the original column
#Strip first part of name using the stringi package to perform regex function
for (i in 1:nrow(sample_info3)) {
  var <- sample_info3$value[i]
  var <- stri_extract_first_regex(var, "^[^:]+")
  #print(var)
  name <- (paste0(var))
  sample_info3$col[i] <- name
  
  newvalue <- sample_info3$value[i]
  newvalue <- regmatches(newvalue,gregexpr("(?<=:).*",newvalue,perl=TRUE))
  name2 <- (paste0(newvalue))
  sample_info3$value[i] <- name2
}


#Since the words are separated by a space e.g. Age: 8, the space before the text AFTER the colon needs to be removed
#2 in the second parameter slot means correct by columns not rows
sample_info3 <- as.data.frame(
  apply(sample_info3,2, function(x) gsub("\\s+", "", x)))

#Remove all rows with missing values
sample_info3<-sample_info3[!(sample_info3$col=="NA"),]

#Convert back into wide format
sample_info3 <- sample_info3 %>%
  dplyr::select(title,value,col) %>%
  spread(key = col, value = value)

#Remove the [2] after the values in the title column
sample_info3 <- sample_info3 %>% mutate(title = gsub("[2]", "", title, fixed=TRUE))

#Remove unneeded columns
sample_info3 = subset(sample_info3, select = -c(tissue,efsbin,osbin) )

#Change the numeric columns from character to numeric
sample_info3 <- sample_info3 %>% 
  mutate_at(vars(efsday,osday,a_efs_all, b_os_all,c_sex_all,d_fav_all,e_efs_hr,f_os_hr), as.numeric)

#Check data
str(sample_info3)



#Merge the patient data from both datasets together
all_patientdata <- merge(sample_info2, sample_info3, on='title')

#Write into Excel file
write.xlsx(all_patientdata,'patient_data.xlsx',colNames = TRUE)

