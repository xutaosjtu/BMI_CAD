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

rst = characteristics(data[which(data$Myocardial.scar==0),clinical], data$clinically.Ischemia[which(data$Myocardial.scar==0)], 2, na.rm=T)
p=NULL
for(i in clinical){
  data$m = data[,i]
  p = c(p, wilcox.test(m~clinically.Ischemia, data[which(data$Myocardial.scar==0),])$p.value)
}
rst = rbind(rst, p)
rownames(rst) = c("type0", "type1", "p")
write.csv(rst, "population characteristics_MyoScar0.csv")

rst = characteristics(data[which(data$Myocardial.scar==1),other], data$clinically.Ischemia[which(data$Myocardial.scar==1)], 2, na.rm=T)
sapply(data[which(data$Myocardial.scar==1),other], table, data$clinically.Ischemia[which(data$Myocardial.scar==1)])

sapply(data[which(data$Myocardial.scar==1),other[-c(1,8)]], 
       function(x) {
         (tmp = fisher.test(table(x, data$clinically.Ischemia[which(data$Myocardial.scar==1)])))
       }
)
wilcox.test(age ~ clinically.Ischemia, data, subset = data$Myocardial.scar==1)


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
  data$m = scale(data[,i])
  model = glm(clinically.Ischemia ~. 
              , data = data[, c("clinically.Ischemia", "m", other[-1])]
              , family = binomial
              , subset = data$Myocardial.scar==1
  )
  rst = rbind(rst, summary(model)$coef[2,])
}
rownames(rst) = clinical
write.csv(rst, "association between the conventionally measured biomarkers and clinical ischemia_Scar1.csv")


## association with clinical relevant ischemia
model.ref = glm(clinically.Ischemia~., data[,c("clinically.Ischemia",clinical, "sex","CVRF..aHT","CVRF..HLP","CVRF..DMT2","CVRF..Smoking","age", "Myocardial.scar")], family=binomial)


rst = NULL
for(i in metabo.valid){
  data$m = scale(log(data[,i]))
  model = glm(clinically.Ischemia~ m + as.factor(Myocardial.scar) 
#                + age + as.factor(sex)  
#               + CVRF..aHT + CVRF..DMT2+ CVRF..Smoking + CVRF..Family
#               + HDL.Cholesterin..mmol.L. + Cholesterin..mmol.L. 
#               + Triglyceride..mmol.L.+Glucose..mg.dL.+HbA1c....
#               +Insulin..mU.l. #+ Chitinase.1
                  , data = data, 
                  #, subset = which(data$Myocardial.scar==1)
                  , family=binomial
  )
  coefs = summary(model)$coef
  rst = rbind(rst,c(coefs["m",]))
}
row.names(rst) = metabo.valid
write.csv(rst, file = "association between Ischemia and metabolites_unadj_adj_MyoScar.csv")

## Association between clinical Ischemia and other clinical measurements
model = glm(clinically.Ischemia~ #as.factor(Myocardial.scar) 
             age + as.factor(sex)  
            + CVRF..aHT + CVRF..DMT2+ CVRF..Smoking + CVRF..Family
            + HDL.Cholesterin..mmol.L. + Cholesterin..mmol.L. 
            + Triglyceride..mmol.L.+Glucose..mg.dL.+HbA1c....
            +Insulin..mU.l. #+ Chitinase.1
            , data = data, 
            , subset = which(data$Myocardial.scar==0)
            , family=binomial
)

###PLS
colnames(data)[5:10], 
require(pls)
pls.rst = mvr(clinically.Ischemia~., data = data[,c(metabo.asso, "clinically.Ischemia")])
plot(pls.rst$scores[,1:2], col = c("green","red")[data$clinically.Ischemia+1])
