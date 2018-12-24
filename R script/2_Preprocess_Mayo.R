#
# Download Mayo data
#
rm(list=ls())
# source("http://bioconductor.org/biocLite.R")
# biocLite("ImmuneSpaceR")
library(ImmuneSpaceR)
sdy67 <- CreateConnection(study = "SDY67")
sdy67$listDatasets()
sdy67.titers <- sdy67$getDataset("hai")
sdy67.ge_info <- sdy67$getDataset("gene_expression_files")
sdy67.env <- sdy67$getGEMatrix(cohort = "healthy adults, 50-74 yo")
sdy67.exprs <- t(sdy67.env@assayData$exprs)
save(sdy67.exprs, file="sdy67.exprs.RData")
sdy67.feature <- sdy67.env@featureData@data
sdy67.pheno <- sdy67.env@phenoData@data

# Separating pheno data to three time points D0, D3, and D7
sdy67.pheno.ordered <- sdy67.pheno[order(sdy67.pheno$participant_id), ]
sdy67.mayo_pheno.d0 <- sdy67.pheno.ordered[which(sdy67.pheno.ordered$study_time_collected==0), ]
sdy67.mayo_pheno.d3 <- sdy67.pheno.ordered[which(sdy67.pheno.ordered$study_time_collected==3), ]
sdy67.mayo_pheno.d28 <- sdy67.pheno.ordered[which(sdy67.pheno.ordered$study_time_collected==28), ]

# ordering titer data
sdy67.titers.ordered <- sdy67.titers[order(sdy67.titers$participant_id), ]

# -------
sdy67.titer.Perth.d0 <- sdy67.titers.ordered[which(sdy67.titers.ordered$virus=="A/Perth/16/2009" & 
                                                                sdy67.titers.ordered$study_time_collected==0), ]
sdy67.titer.Perth.d28 <- sdy67.titers.ordered[which(sdy67.titers.ordered$virus=="A/Perth/16/2009" & 
                                                                 sdy67.titers.ordered$study_time_collected==28), ]
sdy67.titer.Perth <- data.frame(cbind(sdy67.titer.Perth.d0, sdy67.titer.Perth.d28))
sdy67.titer.Perth <- sdy67.titer.Perth[, c(1, 2, 3, 9, 18)]

sdy67.titer.Cal.d0 <- sdy67.titers.ordered[which(sdy67.titers.ordered$virus=="A/California/7/2009" & 
                                                                sdy67.titers.ordered$study_time_collected==0), ]
sdy67.titer.Cal.d28 <- sdy67.titers.ordered[which(sdy67.titers.ordered$virus=="A/California/7/2009" & 
                                                                 sdy67.titers.ordered$study_time_collected==28), ]
sdy67.titer.Cal <- data.frame(cbind(sdy67.titer.Cal.d0, sdy67.titer.Cal.d28))
sdy67.titer.Cal <- sdy67.titer.Cal[, c(1, 2, 3, 9, 18)]

# combine all titers
sdy67.mayo_titers <- cbind(sdy67.titer.Perth, sdy67.titer.Cal[, 4:5])
names(sdy67.mayo_titers)[4:7] <- c("(Perth) day0", "(Perth) day28","(Cal) day0", "(Cal) day28")
# save file
# write.csv(sdy67.mayo_titers, file = "GSE29619 mayo Flu Titers.csv")

# Separate gene expression data to three different time points
# we want to use day0 only, we might use day3 and day7 later.
sdy67.mayo_exprs.d0 <- sdy67.exprs[which(sdy67.mayo_pheno.d0$biosample_accession %in% rownames(sdy67.exprs)), ]
sdy67.mayo_exprs.d3 <- sdy67.exprs[which(sdy67.mayo_pheno.d3$biosample_accession %in% rownames(sdy67.exprs)), ]
sdy67.mayo_exprs.d28 <- sdy67.exprs[which(sdy67.mayo_pheno.d28$biosample_accession %in% rownames(sdy67.exprs)), ]
# checking feature ids of day0
identical(sdy67.feature$FeatureId, colnames(sdy67.mayo_exprs.d0))
# replace gene symbols with feature ids in gene expression matrix
colnames(sdy67.mayo_exprs.d0) <- sdy67.feature$gene_symbol
# remove NA columns in gene expression matrix
sdy67.mayo_exprs.d0 <- data.frame(sdy67.mayo_exprs.d0)
sdy67.mayo_exprs.d0 <- sdy67.mayo_exprs.d0[, grep("^(NA)", names(sdy67.mayo_exprs.d0), value = TRUE, invert = TRUE)]
# write.csv(sdy67.mayo_exprs.d0, file = "sdy67_mayo_exprs.d0.csv")

