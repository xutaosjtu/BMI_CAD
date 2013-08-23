##
## population characteristics
rst = characteristics(data[,clinical], data$clinically.Ischemia, 2, na.rm=T)
p=NULL
for(i in clinical){
  data$m = data[,i]
  p = c(p, t.test(m~clinically.Ischemia, data)$p.value)
}

## Correlation between the conventional markers and metabolites
p.cor=cor(data[,c(clinical, metabo.valid)], use = "pairwise.complete.obs", method = "spearman")

require(gplots)
pdf("correlation of the metabolites with other biomarkers2.pdf", height=10, width=10)
heatmap.2(x=p.cor,
          trace="none",keysize=0.5,
          main = "spearman correlations",
          col = greenred
)
dev.off()

pdf("distribution of the log transformed conventional markers.pdf", )
par(mfrow=c(3,3))
for(i in clinical){
  hist((data[,i]), main = i, xlab="value")
}
dev.off()

## association with clinical relevant ischemia
model.ref = glm(clinically.Ischemia~., data[,c("clinically.Ischemia",clinical, "sex","CVRF..aHT","CVRF..HLP","CVRF..DMT2","CVRF..Smoking","age")], family=binomial)


rst = NULL
for(i in metabo.valid){
  data$m = scale(log(data[,i]))
  model = glm(clinically.Ischemia~ m + age + as.factor(sex)  
                + CVRF..aHT + CVRF..DMT2+ CVRF..Smoking + CVRF..Family
                +HDL.Cholesterin..mmol.L. + Cholesterin..mmol.L. + Triglyceride..mmol.L.+Glucose..mg.dL.+HbA1c....+Chitinase.1+Insulin..mU.l.,
                  data = data, #weights = data$weight,
              #subset = which(data$sex==1),
                  family=binomial
  )
  rst = rbind(rst,summary(model)$coef[2,])
}
row.names(rst) = metabo.valid


###PLS
colnames(data)[5:10], 
require(pls)
pls.rst = mvr(clinically.Ischemia~., data = data[,c(metabo.asso2, "clinically.Ischemia")])
plot(pls.rst$scores[,1:2], col = c("green","red")[data$clinically.Ischemia+1])
