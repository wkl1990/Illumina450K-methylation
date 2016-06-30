# illumina-450K-analysis
This R code uses minfi and wateRmelon preprocess the ILLUMINA 450K from idat to filtered and normalized beta-value.
The analysis includes reading 450K idat files, samples and probes filtered by detectionP and Nbead, type II probes normalized using BMIQ,   removing cross and snp probes using prepared list, and QC plots which includes PCA, PVCA and so on. 
