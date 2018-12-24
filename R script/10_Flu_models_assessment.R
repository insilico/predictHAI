# Assessing the flu vaccine prediction modeling
# 
#
library(bio3d)
# Baylor
bay_d0lmfit <- lm(bay.fc.log2~bay.d0.log2)
# Root Mean Squre Error test
rmse(bay.fc.log2, bay_d0lmfit$fitted.values) # d0
rmse(bay.log2.group$X1, bay.log2.group$obs) # d0+genes
# Root Mean Squre Deviation test
rmsd(bay.fc.log2, bay_d0lmfit$fitted.values) # d0
rmsd(bay.log2.group$X1, bay.log2.group$obs) # d0+genes

# plot the models
par(mfrow = c(3, 2))
plot(bay.fc.log2, bay_d0lmfit$fitted.values, main = "Baylor Day0", ylab = "Predicted (Fold-change)", xlab = "Observed (Fold-change)", add = TRUE) # d0
abline(lm(bay_d0lmfit$fitted.values~0+bay.fc.log2), col = "red") # 45 degree line
curve((x), bay_d0lmfit$fitted.values, from = min(bay.fc.log2), to = max(bay.fc.log2),col = "blue", add = TRUE) # fitted line

plot(bay.log2.group$obs, bay.log2.group$X1, main = "Baylor Day0+Genes", ylab = "Predicted (Fold-change)", xlab = "Observed (Fold-change)", add = TRUE) # d0+genes
abline(lm(bay.log2.group$X1~0+bay.log2.group$obs), col = "red") # 45 degree line
curve((x), bay.log2.group$X1, from = min(bay.log2.group$X1), to = max(bay.log2.group$X1),col = "blue", add = TRUE) # fitted line

# P-values
summary(bay_d0lmfit)$coefficients[2,4] # d0
summary(lm(bay.log2.group$X1~bay.log2.group$obs))$coefficients[2,4] # d0+genes
# F-statistics
summary(bay_d0lmfit)$fstatistic # d0
summary(lm(bay.log2.group$X1~bay.log2.group$obs))$fstatistic # d0+genes
# R-squared
summary(bay_d0lmfit)$r.squared # d0
summary(lm(bay.log2.group$X1~bay.log2.group$obs))$r.squared # d0+genes
# Intercept&Slope
summary(bay_d0lmfit)$coefficients[1:2] # d0
lm(bay.log2.group$X1~bay.log2.group$obs)$coefficients[1:2] # d0+genes


# Emory
emory_d0lmfit <- lm(emory.fc.log2~emory.d0.log2)
# Root Mean Squre Error test
rmse(emory.fc.log2, emory_d0lmfit$fitted.values) # d0
rmse(emory.log2.group$X1, emory.log2.group$obs) # d0+genes
# Root Mean Squre Deviation test
rmsd(emory.fc.log2, emory_d0lmfit$fitted.values) # d0
rmsd(emory.log2.group$X1, emory.log2.group$obs) # d0+genes

# plot the models
plot(emory.fc.log2, emory_d0lmfit$fitted.values, main = "Emory Day0", ylab = "Predicted (Fold-change)", xlab = "Observed (Fold-change)") # d0
abline(lm(emory_d0lmfit$fitted.values~0+emory.fc.log2), col = "red") # 45 degree line
curve((x), emory_d0lmfit$fitted.values, from = min(emory.fc.log2), to = max(emory.fc.log2),col = "blue", add = TRUE) # fitted line

plot(emory.log2.group$obs, emory.log2.group$X1, main = "Emory Day0+Genes", ylab = "Predicted (Fold-change)", xlab = "Observed (Fold-change)") # d0+genes
abline(lm(emory.log2.group$X1~0+emory.log2.group$obs), col = "red") # 45 degree line
curve((x), emory.log2.group$X1, from = min(emory.log2.group$X1), to = max(emory.log2.group$X1),col = "blue", add = TRUE) # fitted line

# P-values
summary(emory_d0lmfit)$coefficients[2,4] # d0
summary(lm(emory.log2.group$X1~emory.log2.group$obs))$coefficients[2,4] # d0+genes
# F-statistics
summary(emory_d0lmfit)$fstatistic # d0
summary(lm(emory.log2.group$X1~emory.log2.group$obs))$fstatistic # d0+genes
# R-squared
summary(emory_d0lmfit)$r.squared # d0
summary(lm(emory.log2.group$X1~emory.log2.group$obs))$r.squared # d0+genes
# Intercept&Slope
summary(emory_d0lmfit)$coefficients[1:2] # d0
summary(lm(emory.log2.group$X1~emory.log2.group$obs))$coefficients[1:2] # d0+genes


# Mayo
mayo_d0lmfit <- lm(mayo.fc.log2~mayo.d0.log2)
# Root Mean Squre Error test
rmse(mayo.fc.log2, mayo_d0lmfit$fitted.values) # d0
rmse(mayo.log2.group$X1, mayo.log2.group$obs) # d0+genes
# Root Mean Squre Deviation test
rmsd(mayo.fc.log2, mayo_d0lmfit$fitted.values) # d0
rmsd(mayo.log2.group$X1, mayo.log2.group$obs) # d0+genes

# plot the models
plot(mayo.fc.log2, mayo_d0lmfit$fitted.values, main = "Mayo Day0", ylab = "Predicted (Fold-change)", xlab = "Observed (Fold-change)", add = TRUE) # d0
abline(lm(mayo_d0lmfit$fitted.values~0+mayo.fc.log2), col = "red") # 45 degree line
curve((x), mayo_d0lmfit$fitted.values, from = min(mayo.fc.log2), to = max(mayo.fc.log2),col = "blue", add = TRUE) # fitted line

plot(mayo.log2.group$obs, mayo.log2.group$X1, main = "Mayo Day0+Genes", ylab = "Predicted (Fold-change)", xlab = "Observed (Fold-change)") # d0+genes
abline(lm(mayo.log2.group$X1~0+mayo.log2.group$obs), col = "red") # 45 degree line
curve((x), mayo.log2.group$X1, from = min(mayo.log2.group$X1), to = max(mayo.log2.group$X1),col = "blue", add = TRUE) # fitted line

# P-values
summary(mayo_d0lmfit)$coefficients[2,4] # d0
summary(lm(mayo.log2.group$X1~mayo.log2.group$obs))$coefficients[2,4] # d0+genes
# F-statistics
summary(lm(mayo_d0lmfit$fitted.values~mayo.fc.log2))$fstatistic # d0
summary(lm(mayo.log2.group$X1~mayo.log2.group$obs))$fstatistic # d0+genes
# R-squared
summary(mayo_d0lmfit)$r.squared # d0
summary(lm(mayo.log2.group$X1~mayo.log2.group$obs))$r.squared # d0+genes
# Intercept&Slope
summary(mayo_d0lmfit)$coefficients[1:2] # d0
summary(lm(mayo.log2.group$X1~mayo.log2.group$obs))$coefficients[1:2] # d0+genes

