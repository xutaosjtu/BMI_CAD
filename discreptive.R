## box plot of the variables
pdf("results/distribution of the conventional biomarkers.pdf")
par(mfrow = c(3,3))
for(i in clinical){
  hist(data[,i], main = i)
  hist(data[which(data$Myocardial.scar==0),i], add = T, col = rgb(1,0,0,0.5))
  hist(data[which(data$Myocardial.scar==1),i], add = T, col = rgb(0,1,0,0.5))
  legend("topright", pch = 19, legend = c("MyoScar=0", "MyoScar=1"), col = c(rgb(1,0,0,0.5), rgb(0,1,0,0.5)))
}
dev.off()

metabo.na = sapply(data[,metabo.valid], function(x) sum(is.na(x)))
plot( metabo.na )
abline(h = 20, col = "red")
text(x = which(metabo.na>20), y = metabo.na[which(metabo.na>20)], labels = names(metabo.na[which(metabo.na>20)]))



