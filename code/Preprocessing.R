if (Sys.getenv("JAVA_HOME")!="") 
  Sys.setenv(JAVA_HOME="")

require(xlsx)

setwd("../../../Dropbox/BMI CAD/")
#data = read.csv("data/BMI-CAD_Wang-Sattler.csv")
data = read.csv(file="data//BMI-CAD_Wang-Sattler_04_2014_2.csv", stringsAsFactors = F)

## Calculate the age of the participants
date.test = sapply(data$MRI..blood.draw, function(x) unlist(strsplit(as.character(x),split="/",fixed=T))[3])
date.birth = sapply(data$Birth.date, function(x) unlist(strsplit(as.character(x),split="/",fixed=T))[3])
date.birth = as.numeric(date.birth)
date.test = as.numeric(date.test)
data$age = date.test-date.birth

data$age[193] = 69
data$CVRF..DMT2[123] = "1"
data$CVRF..DMT2 = as.numeric(data$CVRF..DMT2)

## rename some columns
colnames(data)[2] = "Myocardial.scar"
colnames(data)[11] = "sex"
colnames(data)[12]="clinically.Ischemia"

## subset of data which have metabolomics measurements and also have Clinically relevant Ischemia
data = subset(data, clinically.Ischemia!=2&!is.na(Hexose))

which(sapply(data,class)=="factor")

metabolites = colnames(data)[41:226]
clinical = colnames(data)[13:40]

### additional data of Stathmin(ng/ul), chitinase 1(ng/ul), chitinase 3(ng/ul) and EF-1a(AU)
data.add_1 = read.xlsx(file = "data/Finale Ergebnisse Rudolph August 2014.xlsx", sheetIndex = 2)
data.add_2 = read.xlsx(file = "data/Finale Ergebnisse Rudolph August 2014.xlsx", sheetIndex = 3)
data.add_3 = read.xlsx(file = "data/Finale Ergebnisse Rudolph August 2014.xlsx", sheetIndex = 4)

for(n in colnames(data.add_1)[-1]){
  boxplot(data.add_1[,n], data.add_2[,n], data.add_3[,n], main = n)
}

data.add = rbind(data.add_1, data.add_2, data.add_3)
data = merge(data, data.add, by.x = "ID_Blood.sample", by.y = "Probe.ID", all.x = T)
###
clinical = c(clinical, colnames(data.add)[c(3,5,7,9)])

## 
data$EF1a_2 = data$EF1a
data$EF1a_2[which(data$ID_Blood.sample %in% c(1:81))]=NA
data$EF.1a..AU._2 = data$EF.1a..AU.
data$EF.1a..AU._2[which(data$ID_Blood.sample %in% c(1:81))]=NA

## Preprocessing
data[,c(metabolites, clinical)] = sapply(data[,c(metabolites, clinical)], 
                                         function(x) {x[which(x==0)]=NA;return(x)})

index=apply(data[, c(metabolites,clinical)],2, 
            function(x) which(abs(x)>mean(x,na.rm=T)+4*sd(x,na.rm=T)|abs(x)<mean(x,na.rm=T)-4*sd(x,na.rm=T)))

for(i in names(index)){
  if(length(index[[i]])!=0) data[index[[i]],i]=NA
}

data$weight = rep(1, nrow(data))
data$weight[which(data$Birth.date %in% names(which(table(data$Birth.date)==2)))] = 0.5


metabo.valid = metabolites[-which(sapply(data[,metabolites], function(x) sum(is.na(x)))>0.5*nrow(data))]
clinical.valid = clinical[-which(sapply(data[,clinical], function(x) sum(is.na(x)))>0.3*nrow(data))]
which(sapply(data[,clinical], function(x) sum(is.na(x)))>0.5*nrow(data))
metabo.valid = setdiff(metabo.valid,"PC.aa.C30.2")

# data$total.AC = apply(data[, 82:121],1 ,sum, na.rm = T)
# data$total.PC = apply(data[, 122:197],1 ,sum, na.rm = T)
# data$total.lysoPC = apply(data[, 198:211],1 ,sum, na.rm = T)
# data$total.SM = apply(data[,212:226], 1, sum, na.rm = T)
# data$total.SMOH = apply(data[, 212:216],1 ,sum, na.rm = T)
# data$total.SMnonOH = apply(data[,217:226], 1, sum, na.rm = T)
# data$total.AA = apply(data[,40:60], 1, sum, na.rm = T)


# data$C2.C3.C0 = (data$C2+data$C3)/data$C0
# data$C2.C0 = data$C2/data$C0
# data$ADMA.Arg = data$ADMA/data$Arg
# data$totalDMA.Arg = data$total.DMA/data$Arg
# data$Cit.Arg = data$Cit/data$Arg
# data$Cit.Orn = data$Cit/data$Orn
# data$Kynurenine.Trp = data$Kynurenine/data$Trp
# data$Met.SO.Met = data$Met.SO/data$Met
# data$Orn.Arg = data$Orn/data$Arg
# data$Tyr.Phe = data$Tyr/data$Phe
# data$Putrescine.Orn = data$Putrescine/data$Orn
# data$SDMA.Arg = data$SDMA/data$Arg
# data$Serotonin.Trp = data$Serotonin/data$Trp
# data$Spermidine.Putrescine = data$Spermidine/data$Putrescine
# data$Spermine.Soermidine = data$Spermine/data$Spermidine
# data$total.AC.C0 = data$total.AC/data$C0
# 
# data$CPT.I.ratio = 
# data$essentialAA = 
# data$nonessentialAA =
# data$fisher.ratio
# data$MUFA.PC
# data$PUFA.PC
# data$SFA.PC
# data$total.PC.SM 







