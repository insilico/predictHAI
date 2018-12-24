# 
# Preprocessing of GSE48018(male) and GSE48023(female) RNA-Seq data
# and save for further analysis
#
# ----------------------------------------------------------------------------
# read data using read.ilmn function
# source("http://bioconductor.org/biocLite.R")
library(GEOquery)
library(Biobase)
library(preprocessCore)
# biocLite("cqn")
# library(cqn)
############ GSE48018 (male data) ##########
male.EListraw <- getGEO("GSE48018")
male.rawExpr <- exprs(male.EListraw[[1]])
# save(male.rawExpr, file = "male.rawExpr.RData")

male.normExpr <- normalize.quantiles(male.rawExpr)
# save(male.normExpr, file = "male.normExpr.RData")


# replace gene Symbol to ilumina_id
male.genes <- male.EListraw$GSE48018_series_matrix.txt.gz@featureData@data$Symbol

time_sbj <- pData(male.EListraw[[1]])[,c(12, 13)]
male.probeID <- rownames(male.rawExpr)
male.gene.name <- c("Subject", "Time",male.genes)
male.normExpr <- t(male.normExpr)
male.normExpr <- cbind(time_sbj, data.frame(male.normExpr))
colnames(male.normExpr) <- male.gene.name

############ GSE48023 (female data) ##########
female.EListraw <- getGEO("GSE48023")
female.rawExpr <- exprs(female.EListraw[[1]])
# save(female.rawExpr, file = "female.rawExpr.RData")

female.normExpr <- normalize.quantiles(female.rawExpr)
# save(female.normExpr, file = "female.normExpr.RData")

# replace gene Symbol to ilumina_id
male.genes <- female.EListraw$GSE48023_series_matrix.txt.gz@featureData@data$Symbol

time_sbj <- pData(female.EListraw[[1]])[,c(12, 13)]
female.probeID <- rownames(female.rawExpr)
female.gene.name <- c("Subject", "Time",male.genes)
female.normExpr <- t(female.normExpr)
female.normExpr <- cbind(time_sbj, data.frame(female.normExpr))
colnames(female.normExpr) <- female.gene.name

############ combine GSE48018 and GSE48023 ############
# extract four different data set, day0, day1, day3, day14
male_genes <- colnames(male.normExpr)
female_genes <- colnames(female.normExpr)

male.normExpr.matrix <- as.matrix(male.normExpr)
male.normExpr.commonGenes <- male.normExpr.matrix[ ,which(female_genes %in% male_genes)]

female.normExpr.matrix <- as.matrix(female.normExpr)
female.normExpr.commonGenes <- female.normExpr.matrix[ ,which(female_genes %in% male_genes)]

Baylor_flu_expr <- rbind(male.normExpr.commonGenes, female.normExpr.commonGenes)
rownames(Baylor_flu_expr) <- Baylor_flu_expr[,1]


flu_day0.expr <- Baylor_flu_expr[which(Baylor_flu_expr[,2]=="time: Day0"),]
flu_day1.expr <- Baylor_flu_expr[which(Baylor_flu_expr[,2]=="time: Day1"),]
flu_day3.expr <- Baylor_flu_expr[which(Baylor_flu_expr[,2]=="time: Day3"),]
flu_day14.expr <- Baylor_flu_expr[which(Baylor_flu_expr[,2]=="time: Day14"),]

# write.csv(flu_day0.expr[, 3:ncol(Baylor_flu_expr)],  file = "Flu_Baylor_Day0.csv")
# write.csv(flu_day1.expr[, 3:ncol(Baylor_flu_expr)],  file = "Flu_Baylor_Day1.csv")
# write.csv(flu_day3.expr[, 3:ncol(Baylor_flu_expr)],  file = "Flu_Baylor_Day3.csv")
# write.csv(flu_day14.expr[, 3:ncol(Baylor_flu_expr)],  file = "Flu_Baylor_Day14.csv")




