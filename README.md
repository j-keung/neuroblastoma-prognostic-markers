# neuroblastoma-prognostic-markers

The following code was written for my first MSc project, under the supervision of Dr. Montano and Dr. Huntley. 

The following scripts analyse single-cell and bulk transcriptomic data to investigate the potential of putative neuroblastoma markers.

The single-cell RNA-seq datasets analysed include the Jansky dataset, which can be accessed here:
https://www.neuroblastomacellatlas.org/
https://www.science.org/doi/10.1126/sciadv.abd3311

The bulk RNA scripts analyse a publicly available bulk RNA-seq dataset which can be accessed using the following GEO codes: GSE49711 and GSE62564. GSE62564 is a re-analysis of the GSE49711 dataset, and contains additional clinical data on the same 498 patients in the SEQC cohort. The code concatenates information from the two datasets together.

Genes which are significantly associated with survival are filtered using univariate Cox proportional-hazards regression, and then a prognostic risk score is constructed using LASSO regression. The risk score is tested using multivariate Cox proportional-hazards regression and a log rank test on Kaplan Meier curves.

