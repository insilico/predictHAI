#
# Download Emory data 2007-2011
#
# Emory data 2007-2009, 63 samples
rm(list=ls())
# source("http://bioconductor.org/biocLite.R")
# biocLite("ImmuneSpaceR")
library(ImmuneSpaceR)
sdy269 <- CreateConnection(study = "SDY269")
sdy269$listDatasets()
sdy269.titers <- sdy269$getDataset("hai")
sdy269.ge_info <- sdy269$getDataset("gene_expression_files")
sdy269.env <- sdy269$getGEMatrix(cohort = "TIV Group 2008")
sdy269.exprs <- t(sdy269.env@assayData$exprs)
save(sdy269.exprs, file="sdy269.exprs.RData")
sdy269.feature <- sdy269.env@featureData@data
sdy269.pheno <- sdy269.env@phenoData@data

# Separating pheno data to three time points D0, D3, and D7
sdy269.pheno.ordered <- sdy269.pheno[order(sdy269.pheno$participant_id), ]
sdy269.emory_pheno.d0 <- sdy269.pheno.ordered[which(sdy269.pheno.ordered$study_time_collected==0), ]
sdy269.emory_pheno.d3 <- sdy269.pheno.ordered[which(sdy269.pheno.ordered$study_time_collected==3), ]
sdy269.emory_pheno.d7 <- sdy269.pheno.ordered[which(sdy269.pheno.ordered$study_time_collected==7), ]

# filtering titer data with TIV Group 2008
sdy269.titers_TIV.Group <- sdy269.titers[which(sdy269.titers$cohort == "TIV Group 2008"), ]
sdy269.titers_TIV.Group.ordered <- sdy269.titers_TIV.Group[order(sdy269.titers_TIV.Group$participant_id), ]

# -------
sdy269.titer.h1n1.d0 <- sdy269.titers_TIV.Group.ordered[which(sdy269.titers_TIV.Group.ordered$virus=="A/Brisbane/59/2007" & 
                                                         sdy269.titers_TIV.Group.ordered$study_time_collected==0), ]
sdy269.titer.h1n1.d28 <- sdy269.titers_TIV.Group.ordered[which(sdy269.titers_TIV.Group.ordered$virus=="A/Brisbane/59/2007" & 
                                                          sdy269.titers_TIV.Group.ordered$study_time_collected==28), ]
sdy269.titer.h1n1 <- data.frame(cbind(sdy269.titer.h1n1.d0, sdy269.titer.h1n1.d28))
sdy269.titer.h1n1 <- sdy269.titer.h1n1[, -c(4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17)]

sdy269.titer.h3n2.d0 <- sdy269.titers_TIV.Group.ordered[which(sdy269.titers_TIV.Group.ordered$virus=="A/Uruguay/716/2007  NYMC X-175C (H3N2)" & 
                                                         sdy269.titers_TIV.Group.ordered$study_time_collected==0), ]
sdy269.titer.h3n2.d28 <- sdy269.titers_TIV.Group.ordered[which(sdy269.titers_TIV.Group.ordered$virus=="A/Uruguay/716/2007  NYMC X-175C (H3N2)" & 
                                                          sdy269.titers_TIV.Group.ordered$study_time_collected==28), ]
sdy269.titer.h3n2 <- data.frame(cbind(sdy269.titer.h3n2.d0, sdy269.titer.h3n2.d28))
sdy269.titer.h3n2 <- sdy269.titer.h3n2[, -c(4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17)]

sdy269.titer.b.d0 <- sdy269.titers_TIV.Group.ordered[which(sdy269.titers_TIV.Group.ordered$virus=="B/Brisbane/3/2007" & 
                                                      sdy269.titers_TIV.Group.ordered$study_time_collected==0), ]
sdy269.titer.b.d28 <- sdy269.titers_TIV.Group.ordered[which(sdy269.titers_TIV.Group.ordered$virus=="B/Brisbane/3/2007" & 
                                                       sdy269.titers_TIV.Group.ordered$study_time_collected==28), ]
sdy269.titer.b <- data.frame(cbind(sdy269.titer.b.d0, sdy269.titer.b.d28))
sdy269.titer.b <- sdy269.titer.b[, -c(4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17)]

# combine all titers
sdy269.emory_titers <- cbind(sdy269.titer.h1n1, sdy269.titer.h3n2[, 4:5], sdy269.titer.b[4:5])
names(sdy269.emory_titers)[4:9] <- c("(H1N1) day0", "(H1N1) day28","(H3N2) day0", "(H3N2) day28", "(B) day0", "(B) day28")
# save file
# write.csv(sdy269.emory_titers, file = "GSE29619 Emory Flu Titers.csv")

# Separate gene expression data to three different time points
# we want to use day0 only, we might use day3 and day7 later.
sdy269.emory_exprs.d0 <- sdy269.exprs[which(sdy269.emory_pheno.d0$biosample_accession %in% rownames(sdy269.exprs)), ]
sdy269.emory_exprs.d3 <- sdy269.exprs[which(sdy269.emory_pheno.d3$biosample_accession %in% rownames(sdy269.exprs)), ]
sdy269.emory_exprs.d7 <- sdy269.exprs[which(sdy269.emory_pheno.d7$biosample_accession %in% rownames(sdy269.exprs)), ]
# checking feature ids of day0
identical(sdy269.feature$FeatureId, colnames(sdy269.emory_exprs.d0))
# replace gene symbols with feature ids in gene expression matrix
colnames(sdy269.emory_exprs.d0) <- sdy269.feature$gene_symbol
# remove NA columns in gene expression matrix
sdy269.emory_exprs.d0 <- data.frame(sdy269.emory_exprs.d0)
sdy269.emory_exprs.d0 <- sdy269.emory_exprs.d0[, grep("^(NA)", names(sdy269.emory_exprs.d0), value = TRUE, invert = TRUE)]
# write.csv(sdy269.emory_exprs.d0, file = "GSE29619_emory_exprs.d0.csv")

# boxplot(emory_exprs.d0)

### Emory data 2009-2011, 92 samples
sdy56 <- CreateConnection(study = "SDY56")
sdy56$listDatasets()
sdy56.titers <- sdy56$getDataset("hai")
sdy56.ge_info <- sdy56$getDataset("gene_expression_files")
sdy56.env_young <- sdy56$getGEMatrix("TIV_young")
sdy56.exprs_young <- t(sdy56.env_young@assayData$exprs)
sdy56.env_older <- sdy56$getGEMatrix("TIV_older")
sdy56.exprs_older <- t(sdy56.env_older@assayData$exprs)
sdy56.exprs <- rbind(sdy56.exprs_young, sdy56.exprs_older)
save(sdy56.exprs, file="sdy56.exprs.RData")
sdy56.feature <- sdy56.env_young@featureData@data

# Separating titer data to three different viruses
# H1N1
sdy56.titer.h1n1.d0_1 <- sdy56.titers[which(sdy56.titers$virus=="A/California/7/2009" & sdy56.titers$study_time_collected==0), ]
sdy56.titer.h1n1.d0_2 <- sdy56.titers[which(sdy56.titers$virus=="H1N1 A/California/07/09" & sdy56.titers$study_time_collected==0), ]
sdy56.titer.h1n1.d0 <- rbind(sdy56.titer.h1n1.d0_1, sdy56.titer.h1n1.d0_2)
sdy56.titer.h1n1.d30_1 <- sdy56.titers[which(sdy56.titers$virus=="A/California/7/2009" & sdy56.titers$study_time_collected==30), ]
sdy56.titer.h1n1.d30_2 <- sdy56.titers[which(sdy56.titers$virus=="H1N1 A/California/07/09" & sdy56.titers$study_time_collected==30), ]
sdy56.titer.h1n1.d30 <- rbind(sdy56.titer.h1n1.d30_1, sdy56.titer.h1n1.d30_2)
# common subjects
sdy56.titer.h1n1.d0_fltr <- sdy56.titer.h1n1.d0[!duplicated(sdy56.titer.h1n1.d0$participant_id), ]
sdy56.titer.h1n1.d0_ord <- sdy56.titer.h1n1.d0_fltr[order(sdy56.titer.h1n1.d0_fltr$participant_id), ]
sdy56.titer.h1n1.d30_ord <- sdy56.titer.h1n1.d30[order(sdy56.titer.h1n1.d30$participant_id), ]
sdy56.titer.h1n1 <- data.frame(cbind(sdy56.titer.h1n1.d0_ord, sdy56.titer.h1n1.d30_ord))
sdy56.titer.h1n1 <- sdy56.titer.h1n1[, -c(4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17)]

# H3N2
sdy56.titer.h3n2.d0_1 <- sdy56.titers[which(sdy56.titers$virus=="A/Perth/16/2009" & sdy56.titers$study_time_collected==0), ]
sdy56.titer.h3n2.d0_2 <- sdy56.titers[which(sdy56.titers$virus=="H3N2 A/Perth/16/09" & sdy56.titers$study_time_collected==0), ]
sdy56.titer.h3n2.d0 <- rbind(sdy56.titer.h3n2.d0_1, sdy56.titer.h3n2.d0_2)
sdy56.titer.h3n2.d30_1 <- sdy56.titers[which(sdy56.titers$virus=="A/Perth/16/2009" & sdy56.titers$study_time_collected==30), ]
sdy56.titer.h3n2.d30_2 <- sdy56.titers[which(sdy56.titers$virus=="H3N2 A/Perth/16/09" & sdy56.titers$study_time_collected==30), ]
sdy56.titer.h3n2.d30 <- rbind(sdy56.titer.h3n2.d30_1, sdy56.titer.h3n2.d30_2)
# common subjects
sdy56.titer.h3n2.d0_fltr <- sdy56.titer.h3n2.d0[!duplicated(sdy56.titer.h3n2.d0$participant_id), ]
sdy56.titer.h3n2.d0_ord <- sdy56.titer.h3n2.d0_fltr[order(sdy56.titer.h3n2.d0_fltr$participant_id), ]
sdy56.titer.h3n2.d30_ord <- sdy56.titer.h3n2.d30[order(sdy56.titer.h3n2.d30$participant_id), ]
sdy56.titer.h3n2 <- data.frame(cbind(sdy56.titer.h3n2.d0_ord, sdy56.titer.h3n2.d30_ord))
sdy56.titer.h3n2 <- sdy56.titer.h3n2[, -c(4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17)]

# B
sdy56.titer.b.d0_1 <- sdy56.titers[which(sdy56.titers$virus=="B/Brisbane/60/2008" & sdy56.titers$study_time_collected==0), ]
sdy56.titer.b.d0_2 <- sdy56.titers[which(sdy56.titers$virus=="B/Brisbane/60/08" & sdy56.titers$study_time_collected==0), ]
sdy56.titer.b.d0 <- rbind(sdy56.titer.b.d0_1, sdy56.titer.b.d0_2)
sdy56.titer.b.d30_1 <- sdy56.titers[which(sdy56.titers$virus=="B/Brisbane/60/2008" & sdy56.titers$study_time_collected==30), ]
sdy56.titer.b.d30_2 <- sdy56.titers[which(sdy56.titers$virus=="B/Brisbane/60/08" & sdy56.titers$study_time_collected==30), ]
sdy56.titer.b.d30 <- rbind(sdy56.titer.b.d30_1, sdy56.titer.b.d30_2)
# common subjects
sdy56.titer.b.d0_fltr <- sdy56.titer.b.d0[!duplicated(sdy56.titer.b.d0$participant_id), ]
sdy56.titer.b.d0_ord <- sdy56.titer.b.d0_fltr[order(sdy56.titer.b.d0_fltr$participant_id), ]
sdy56.titer.b.d30_ord <- sdy56.titer.b.d30[order(sdy56.titer.b.d30$participant_id), ]
sdy56.titer.b <- data.frame(cbind(sdy56.titer.b.d0_ord, sdy56.titer.b.d30_ord))
sdy56.titer.b <- sdy56.titer.b[, -c(4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17)]

# combine all titers
sdy56.emory_titers <- cbind(sdy56.titer.h1n1, sdy56.titer.h3n2[, 4:5], sdy56.titer.b[4:5])
names(sdy56.emory_titers)[4:9] <- c("(H1N1) day0", "(H1N1) day30","(H3N2) day0", "(H3N2) day30", "(B) day0", "(B) day30")

# intersect subjects
sdy56.emory_titers <- sdy56.emory_titers[which(sdy56.emory_titers$participant_id %in% unique(sdy56.ge_info$participant_id)), ]
# save titers as csv file
# write.csv(sdy56.emory_titers, file = "GSE74817 Emory Flu Titers.csv")

# Separating the gene expression data to five different time points
sdy56.ge_info.d0 <- sdy56.ge_info[which(sdy56.ge_info$study_time_collected==0), ]

# we want to use day0 gene expression only, we might use day1, day3, day7, and day14 later.
sdy56.emory_exprs.d0 <- sdy56.exprs[which(sdy56.ge_info$study_time_collected==0), ]
sdy56.emory_exprs.d1 <- sdy56.exprs[which(sdy56.ge_info$study_time_collected==1), ]
sdy56.emory_exprs.d3 <- sdy56.exprs[which(sdy56.ge_info$study_time_collected==3), ]
sdy56.emory_exprs.d7 <- sdy56.exprs[which(sdy56.ge_info$study_time_collected==7), ]
sdy56.emory_exprs.d14 <- sdy56.exprs[which(sdy56.ge_info$study_time_collected==14), ]
# checking feature ids of day0
identical(sdy56.feature$FeatureId, colnames(sdy56.emory_exprs.d0))
# replace gene symbols with feature ids in gene expression matrix
colnames(sdy56.emory_exprs.d0) <- sdy56.feature$gene_symbol
# remove NA columns in gene expression matrix
sdy56.emory_exprs.d0 <- data.frame(sdy56.emory_exprs.d0)
sdy56.emory_exprs.d0 <- sdy56.emory_exprs.d0[, grep("^(NA)", names(sdy56.emory_exprs.d0), value = TRUE, invert = TRUE)]
# intersect subjects
sdy56.emory_exprs.d0 <- sdy56.emory_exprs.d0[which(unique(sdy56.ge_info$participant_id) %in% sdy56.emory_titers$participant_id), ]
# save gene expression file
# write.csv(sdy56.emory_exprs.d0, file = "GSE74817_emory_exprs.d0.csv")

