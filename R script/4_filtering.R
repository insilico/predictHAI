#
#
# Preprocess all data
#
# Developed by Saeid Parvandeh 12/25/2016
# 
#-------------------------------------------------
rm(list=ls())

#### ---- Baylor Day0 ---- ####
Baylor_d0.expr <- read.csv("flu_Baylor_Day0.csv", row.names=1, header=T)
# brush up row names
extr_id   <- function(x){substr(x, start=10, stop=15)} # function for manipulate row names
row_names <- rownames(Baylor_d0.expr); sbj_id  <- extr_id(row_names)
rownames(Baylor_d0.expr)<-sbj_id
# read titer data
library(readxl)
baylor_titer <- read_excel("Flu Baylor Titers.xlsx", sheet = 2)

# check OS to remove the 0's or NA's from data
baylor_titer <- baylor_titer[-which(is.na(baylor_titer[, 21])),]

# find common subject ids in both titers and expression data
baylor_titers <- baylor_titer[which(sbj_id %in% baylor_titer$IN_ID),]
save(baylor_titers, file = "baylor_titers.RData")
Baylor_expr <- Baylor_d0.expr[which(sbj_id %in% as.character(baylor_titer[,1])), ]
Baylor_expr <- Baylor_expr[, order(colnames(Baylor_expr))]
save(Baylor_expr, file = "Baylor_expr.RData")

#### ---- Emory Day0 Data ---- ####
# GSE29619
Emory_expr_1 <- read.csv("GSE29619_emory_exprs.d0.csv", row.names=1, header=T)
# GSE74817
Emory_expr_2 <- read.csv("GSE74817_emory_exprs.d0.csv", row.names=1, header=T)
identical(names(Emory_expr_1), names(Emory_expr_2))
Emory_expr <- rbind(Emory_expr_1, Emory_expr_2)
# save(Emory_expr, file = "Emory_expr.RData")
emory_titers <- read.csv("Emory Flu Titers.csv", header = TRUE, quote = "")
# save(emory_titers, file = "emory_titers.RData")

#### ---- Mayo Day0 Data ---- ####
Mayo_expr <- read.csv("sdy67_mayo_exprs.d0.csv", row.names=1, header=T)
Mayo_expr <- t(Mayo_expr)
extr_id   <- function(x){substr(x, start=3, stop=6)} # function for manipulate row names
row_names <- rownames(Mayo_expr); sbj_id  <- extr_id(row_names)
rownames(Mayo_expr)<-sbj_id
# save(Mayo_expr, file = "Mayo_expr.RData")

# read Mayo titers
mayo_titers <- read.table("GSE29619 mayo Flu Titers.csv", sep="\t", header=TRUE)
mayo_titers_idx <- as.character(mayo_titers[, 1])
mayo_titers <- mayo_titers[which(mayo_titers_idx %in% sbj_id), ]
# save(mayo_titers, file = "mayo_titers.RData")

# -------- 1. Coefficient of Variation (CoV) ---------
cv.filter <- function(dataMatrix, threshold) {
  # coefficient of variation filter
  mask <- apply(dataMatrix, 1, function(x) {(sd(x)/abs(mean(x))) < threshold})
  fdata <- dataMatrix[mask, ]
  # return the row mask and filtered data
  list(mask=mask, fdata=fdata)
}

# Mayo
Thresh <- .13396 # filtering 5000 genes by CoV
exprData1 <- t(Mayo_expr)
Mayo_expr.filter <- cv.filter(exprData1,Thresh)$fdata 

# Baylor
exprData2 <- t(Baylor_expr)
Bay_expr.filter <- cv.filter(exprData2,Thresh)$fdata

# Emory
exprData3 <- t(Emory_expr)
Emory_expr.filter <- cv.filter(exprData3,Thresh)$fdata 

# overlap genes
m_genes <- as.character(colnames(t(Mayo_expr.filter)))
b_genes <- as.character(colnames(t(Bay_expr.filter)))
e_genes <- as.character(colnames(t(Emory_expr.filter)))

common <- intersect(intersect(e_genes, m_genes), b_genes)# excluding Stanford

B_ex.cv.fltr <- t(Bay_expr.filter[which(b_genes %in% common), ])
B_ex.cv.fltr <- B_ex.cv.fltr[, order(colnames(B_ex.cv.fltr))]
dim(B_ex.cv.fltr)
# save(B_ex.cv.fltr, file = "B_ex.cv.fltr.RData")

M_ex.cv.fltr <- t(Mayo_expr.filter[which(m_genes %in% common), ])
M_ex.cv.fltr <- M_ex.cv.fltr[, order(colnames(M_ex.cv.fltr))]
dim(M_ex.cv.fltr)
# save(M_ex.cv.fltr, file = "M_ex.cv.fltr.RData")

E_ex.cv.fltr <- t(Emory_expr.filter[which(e_genes %in% common), ])
E_ex.cv.fltr <- E_ex.cv.fltr[, order(colnames(E_ex.cv.fltr))]
dim(E_ex.cv.fltr)
# save(E_ex.cv.fltr, file = "E_ex.cv.fltr.RData")


# ---------- Compute reGAIN with fc residuals as outcome ---------
# read Baylor titers
bay.d0 <- baylor_titers$`Matched Max day0`
bay.d28 <- baylor_titers$`Max day28`
bay.age <- baylor_titers$Age
bay.max.fc <- baylor_titers$`MAX FC`
bay.fc <- bay.max.fc

# read Emory titers
emory.d0 <- emory_titers$Matched.Max.day0
emory.d28 <- emory_titers$MAX.day28
emory.age <- emory_titers$age_reported
emory.max.fc <- emory_titers$MAX.FC
emory.fc <- emory.max.fc

# read Mayo titers
mayo.d0 <- mayo_titers$day0
mayo.d28 <- mayo_titers$day28
mayo.fc <- mayo.d28/mayo.d0

# compute fold change residuals
bay.df <- data.frame(fc = bay.fc, d0 = bay.d0)
bay.logfit <- lm(log2(fc)~log2(d0), data = bay.df)
bay.logfit.sum <- summary(bay.logfit)
cat("R-squared")
bay.logfit.sum$r.squared

a<-bay.logfit$coefficients[1]
b<-bay.logfit$coefficients[2]

# Baylor fold-change residuals
bay.fc.log.resid <- as.vector(bay.logfit$residuals)
bay.fc.log.resid <- log2(bay.fc) - (a + b*log2(bay.d0))
# save(bay.fc.log.resid, file="bay.fc.log.resid.RData")

# Baylor reGAIN
baylor_regain_resid <- Rinbix::regain(B_ex.cv.fltr, family = "gaussian")

