##
## population characteristics

## distribution of the clinical measurements
pdf("data distribution.pdf", width = 10, height = 10)
par(mfrow = c(3,3))
for(i in c(clinical.valid, metabo.valid)){
  hist(data[,i], main = i)
}
dev.off()

other = colnames(data)[c(2, 6:11)]
other = c(other, "age")

group =  1:nrow(data)#which(data$Myocardial.scar==0)
rst = characteristics(data[group,clinical], data$clinically.Ischemia[group], 2, na.rm=T)
p=NULL
for(i in clinical){
  data$m = data[,i]
  p = c(p, wilcox.test(m~clinically.Ischemia, data[group,])$p.value)
}
rst = rbind(rst, p)
rownames(rst) = c("type0", "type1", "p")
rst.other = characteristics(data[group,other], data$clinically.Ischemia[group], 2, na.rm=T)
p.other = c(NA,
            sapply(data[group,other[-c(1,8)]], 
                   function(x) {
                     tmp = fisher.test(table(x, data$clinically.Ischemia[group]))
                     return(tmp$p.value)
                   }),
            age = wilcox.test(age ~ clinically.Ischemia, data, subset = group)$p.value
)
rst.other = rbind(rst.other, p.other)
rst = cbind(rst.other, rst)
write.csv(t(rst), "population characteristics_MyoScar1.csv")


sapply(data[which(data$Myocardial.scar==1),other], table, data$clinically.Ischemia[which(data$Myocardial.scar==1)])



## Correlation between the conventional markers and metabolites
cor.mtest <- function(mat, conf.level = 0.95, method = "pearson"){
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      tmp <- cor.test(mat[,i], mat[,j], conf.level = conf.level, method=method)
      p.mat[i,j] <- p.mat[j,i] <- tmp$p.value
      #lowCI.mat[i,j] <- lowCI.mat[j,i] <- tmp$conf.int[1]
      #uppCI.mat[i,j] <- uppCI.mat[j,i] <- tmp$conf.int[2]
    }
  }
  return(p.mat)
}

corr=cor(data[which(data$Myocardial.scar==1),c(clinical, metabo.valid)], use = "pairwise.complete.obs", method = "spearman")

corr.p = cor.mtest(data[,c(clinical, metabo.valid)], method = "spearman")
colnames(corr.p) = c(clinical, metabo.valid)
rownames(corr.p) = c(clinical, metabo.valid)

require(gplots)
pdf("differential correlation of the metabolites with other biomarkers_MyoScar.pdf", height=20, width=10)
heatmap.2(x=corr[metabo.valid,clinical],
          trace="none",keysize=0.5,
          main = "spearman correlations",
          col = greenred,
          margins = c(10, 5)
)
dev.off()

write.csv(corr[metabo.valid, clinical], file = "Correlation between metabolites and conventional markers.csv") 
write.csv(corr.p[metabo.valid, clinical], file = "Correlation (pvalue) between metabolites and conventional markers.csv")

pdf("distribution of the log transformed conventional markers.pdf", )
par(mfrow=c(3,3))
for(i in clinical){
  hist(log(data[,i]), main = i, xlab="value")
}
dev.off()


## Associations between conventionally measured biomarkers and clinical ischemia
rst = NULL;
for(i in clinical){
  data$m = scale(log(data[,i]))
  model = glm(clinically.Ischemia ~. 
              , data = data[, c("clinically.Ischemia", "m", other[-1])]
              , family = binomial
               , subset = data$Myocardial.scar==1
  )
  rst = rbind(rst, summary(model)$coef[2,])
}
rownames(rst) = clinical
write.table(rst, "results/association between cTnI and clinical ischemia.txt", sep = "\t", quote = F)
write.csv(rst, "association between the conventionally measured biomarkers and clinical ischemia.csv")

## associations between metabolites and cTnI
rst = NULL
for (i in metabo.valid){
  data$m = scale(log(data[,i]))
  model = lm(cTnI.pg.mL ~ m 
             + age + as.factor(sex)  
             + CVRF..aHT + CVRF..DMT2+ CVRF..Smoking + CVRF..Family
             + HDL.Cholesterin..mmol.L. + Cholesterin..mmol.L. 
             + Triglyceride..mmol.L.+Glucose..mg.dL.+HbA1c....
             + Insulin..mU.l.,
            subset = which(data$Myocardial.scar==1),
             data)
  coefs = summary(model)$coef
#   coefs = cbind(OR = exp(coefs[,1]), 
#                 lower = exp(coefs[,1] - 1.96*coefs[,2]),
#                 upper = exp(coefs[,1] + 1.96*coefs[,2]),
#                 P = coefs[,4])
  rst = rbind(rst, coefs[2,c(1,2,4)])
}
rownames(rst) = metabo.valid
write.table(rst, file = "results/associations between metabolites and cTnI.txt", sep = "\t", quote = F)


## association with clinical relevant ischemia
model.ref = glm(clinically.Ischemia~., data[,c("clinically.Ischemia",clinical, "sex","CVRF..aHT","CVRF..HLP","CVRF..DMT2","CVRF..Smoking","age", "Myocardial.scar")], family=binomial)

##############################################
##  adjustment:
##    unadjusted
##    crude: adjusted for age, sex
##    multivar: adjusted for aHT, DMT2, smoking, family history, HDL-cholesterol, total cholesterol, triglyceride, glucose, HbA1c, insulin
##    multivar2: adjusted for aHT, DMT2, smoking, family history, HDL-cholesterol, total cholesterol
##    multivar3: adjusted for aHT, DMT2, smoking, family history
rst = NULL
for(i in metabo.valid){
  data$m = scale(log(data[,i]))
  model = glm(clinically.Ischemia~ m * as.factor(Myocardial.scar) 
                 + age + as.factor(sex)  
                + CVRF..aHT + CVRF..DMT2+ CVRF..Smoking + CVRF..Family
                + HDL.Cholesterin..mmol.L. + Cholesterin..mmol.L. 
                + Triglyceride..mmol.L.+Glucose..mg.dL.+HbA1c....
                + Insulin..mU.l. 
                  , data = data, 
#                   , subset = which(data$Myocardial.scar==0)
                  , family=binomial
  )
  coefs = summary(model)$coef
  coefs = cbind(OR = exp(coefs[,1]), 
           lower = exp(coefs[,1] - 1.96*coefs[,2]),
           upper = exp(coefs[,1] + 1.96*coefs[,2]),
           P = coefs[,4])
  rst = rbind(rst, c(coefs[2,], coefs[nrow(coefs),]))
}
row.names(rst) = metabo.valid
write.csv(rst, file = "association between Ischemia and metabolites_multivar_interaction.csv")

## Association between clinical Ischemia and other clinical measurements
model = glm(clinically.Ischemia~ as.factor(Myocardial.scar) +
             age + as.factor(sex)  
            + CVRF..aHT + CVRF..DMT2+ CVRF..Smoking + CVRF..Family
            + HDL.Cholesterin..mmol.L. + Cholesterin..mmol.L. 
            + Triglyceride..mmol.L.+Glucose..mg.dL.+HbA1c....
            +Insulin..mU.l. #+ Chitinase.1
            , data = data, 
#           subset = which(data$Myocardial.scar==1),
            family=binomial
)

write.csv(summary(model)$coef, file = "association between the adjusted covariates with clinical ischemia_adjust MyoScar.csv")

###PLS
colnames(data)[5:10], 
require(pls)
pls.rst = plsr(clinically.Ischemia~., data = data.0impute[,c(metabo.valid, "clinically.Ischemia")], scale = T)
plot(pls.rst$scores[,1:2], col = c("green","red")[data$clinically.Ischemia+1])


### association with PCs of the data set
data.0impute = data
data.0impute[which(is.na(data), arr.ind = T)]=0

require(caret)
preProc = preProcess(data.0impute[,metabo.valid], 
           method = "pca",
           thresh = 0.95,
           )
pcs = scale(as.matrix(data.0impute[,metabo.valid])) %*% preProc$rotation

plot(100*diag(cov(pcs))/(sum(diag(cov(pcs))/0.95)), type = "h", 
     main = "Explained variance", ylab = "Percentage (%)", xlab = "PCs")
## select PCs with eigenvalue bigger than one
which(sqrt(diag(cov(pcs)))>1)
## total varianaces explained by the PCs
sum(100*diag(cov(pcs))[1:30]/(sum(diag(cov(pcs))/0.95)))
## 1. do varimax transformation of the rotation matrix;
## 2. Find the individual components in each of the factors 
rotation.varimax = varimax(preProc$rotation)
individual.components = apply(rotation.varimax$loadings[,1:30], 2,
       function(x) {
         metabolites= names(which(abs(x)>0.3))
         return(paste(metabolites, collapse = ","))
       }
)

data.pc = cbind(data.0impute[, c("clinically.Ischemia", clinical, other)], pcs)

## distribution of samples in each analysis batches
data.pc$batch = 2
data.pc$batch[which(data$ID_Blood.sample>=163)]=3
data.pc$batch[which(data$ID_Blood.sample<=81)]=1

fisher.test(table(data.pc$clinically.Ischemia, data.pc$batch))
fisher.test(table(data.pc$Myocardial.scar, data.pc$batch))

subset = which(data.pc$Myocardial.scar==1)
fisher.test(table(data.pc$clinically.Ischemia[subset], data.pc$batch[subset]))

## PCA analysis of potential batch effects
pdf("PCA for potential batch effects.pdf", width = 7, height = 7)
par(mfrow = c(2,2))
plot(data.pc$PC1, data.pc$PC2, col = c("red", "green", "blue")[data.pc$batch], pch = 19)
legend("topleft", legend = c("Apr.12", "Apr.18", "Apr.26"), pch = 19, col = c("red", "green", "blue"))
plot(data.pc$PC2, data.pc$PC3, col = c("red", "green", "blue")[data.pc$batch], pch = 19)
plot(data.pc$PC3, data.pc$PC4, col = c("red", "green", "blue")[data.pc$batch], pch = 19)
plot(data.pc$PC4, data.pc$PC5, col = c("red", "green", "blue")[data.pc$batch], pch = 19)
dev.off()


## Correlation between different PCs and the clinical variables
corr = cor( pcs, data.pc[, clinical], method = "spearman", use = "pair")
require(gplots)
pdf("Correlation of the PCs of the metabolites with other biomarkers.pdf", height=20, width=10)
heatmap.2(x=corr,
          trace="none",keysize=0.5,
          main = "spearman correlations",
          col = greenred,
          margins = c(10, 5)
)
dev.off()


rst = NULL
for(i in colnames(pcs)){
  data.pc$m = data.pc[,i]
  model = glm(clinically.Ischemia~ m #+ as.factor(Myocardial.scar) 
               + age + as.factor(sex)  
              + CVRF..aHT + CVRF..DMT2+ CVRF..Smoking + CVRF..Family
              + HDL.Cholesterin..mmol.L. + Cholesterin..mmol.L. 
              + Triglyceride..mmol.L.+Glucose..mg.dL.+HbA1c....
              + Insulin..mU.l. 
              , data = data.pc, 
              , subset = which(data$Myocardial.scar==0)
              , family=binomial
  )
  coefs = summary(model)$coef
  coefs = cbind(OR = exp(coefs[,1]), 
                lower = exp(coefs[,1] - 1.96*coefs[,2]),
                upper = exp(coefs[,1] + 1.96*coefs[,2]),
                P = coefs[,4])
  rst = rbind(rst, c(coefs[2,], coefs[nrow(coefs),]))
}
row.names(rst) = colnames(pcs)
write.csv(rst, "association between Ischemia and metabolite PCs_multivar_Myocardial.scar0.csv")

pdf("factor loadings for the associated factors.pdf", width = 7, height = 10)
par(mfrow = c(3,1))
plot(preProc$rotation[,3], type = "h", col = c("grey", "black")[as.numeric(abs(preProc$rotation[,3])>0.15)+1], main = "PC 3", ylab = "loading")
plot(preProc$rotation[,4], type = "h", col = c("grey", "black")[as.numeric(abs(preProc$rotation[,4])>0.15)+1], main = "PC 4", , ylab = "loading")
plot(preProc$rotation[,9], type = "h", col = c("grey", "black")[as.numeric(abs(preProc$rotation[,9])>0.15)+1], main = "PC 9", , ylab = "loading")
dev.off()


d = dist(scale(log(data[, c(metabo.valid)])))
hc = hclust(d)
hcd = as.dendrogram(hc)
labelColors = c( "black", "red")
colLab <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    labCol <- labelColors[data$clinically.Ischemia [which(data$ID_Blood.sample == a$label)]+1]
    attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
  }
  n
}
clusDendro = dendrapply(hcd, colLab)
plot(clusDendro,hang = -1)

pdf("clustering of the samples.pdf", height=20, width=10)
heatmap.2(x= scale(log(data[, c(metabo.valid)])),
          trace="none",keysize=0.5,
          col = greenred,
          margins = c(10, 5)
)
dev.off()


dendrapply(hcd, colLab.read)

colLab.read <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    return(a$nodePar)
  }
}