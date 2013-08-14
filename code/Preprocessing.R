setwd("../../../Dropbox/BMI CAD/")
data = read.csv("data/BMI-CAD_Wang-Sattler.csv")

date.test = sapply(data$MRI..blood.draw, function(x) unlist(strsplit(as.character(x),split=".",fixed=T))[3])
date.birth = sapply(data$Birth.date, function(x) unlist(strsplit(as.character(x),split=".",fixed=T))[3])
date.birth = as.numeric(date.birth)
date.test = as.numeric(date.test)
data$age = date.test-date.birth

colnames(data)[10]="sex"
colnames(data)[11]="clinically.Ischemia"
data = subset(data, clinically.Ischemia!=2&!is.na(Hexose))

which(sapply(data,class)=="factor")

metabolites = colnames(data)[40:227]
clinical = colnames(data)[12:39]


index=apply(data[,metabolites],2, function(x) which(abs(x)>mean(x,na.rm=T)+4*sd(x,na.rm=T)|abs(x)<mean(x,na.rm=T)-4*sd(x,na.rm=T)))
for(i in names(index)){
  if(length(index[[i]])!=0) data[index[[i]],i]=NA
}


metabo.valid = metabolites[-which(sapply(data[,metabolites], function(x) sum(is.na(x)))>0.5*nrow(data))]
which(sapply(data[,clinical], function(x) sum(is.na(x)))>0.5*nrow(data))
metabo.valid = setdiff(metabo.valid,"PC.aa.C30.2")

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
  hist(log(data[,i]), main = i)
}
dev.off()