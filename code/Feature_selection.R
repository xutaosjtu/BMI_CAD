## random forest
require(randomForest)
tmp = data[,c(metabo.asso,clinical.asso, "clinically.Ischemia")]
tmp = na.omit(tmp)
rf = randomForest(x = tmp[,1:(ncol(tmp)-1)], y = as.factor(tmp$clinically.Ischemia))
rfselect = rfcv(tmp[,1:(ncol(tmp)-1)], as.factor(tmp$clinically.Ischemia))

##variable selection in S2
tmp = S2[which(S2$subcoho==1|S2$inz_mi==1),c("ctalteru", "ccsex","ctbmi","my.diab","ctsysmm","my.cigreg","my.alkkon","cl_chola","cl_hdla","cl_crp",S2_valid_measures,"Arg.Trp","zz_nr","mi_time", "inz_mi","subcoho","prev_mi","ctantihy")]
tmp[,S2_valid_measures] = log(tmp[,S2_valid_measures])
clinical = c("ctalteru", "ccsex","ctbmi","my.diab","ctsysmm","my.cigreg","my.alkkon","cl_chola","cl_hdla","cl_crp")

require(penalized)
model.penal = penalized(
  clinically.Ischemia~., data = tmp[,c("clinically.Ischemia", metabo.asso, clinical.asso)], 
  standardize = T, ##centralize the covariates 
  #steps = "Park", trace = F, ## starting from the largest value of lambda1 untile the specificed lambda1 is reached, the function will return a list of penfit objects, and the chang of coefficients could be visualized by plotpath 
  #positive = T, ## positive restriction to all regression coefficients
  lambda1 = 2, lambda2= 0 #penalization parameter 
)

#### stepwise selection of cox regression
model = glm(clinically.Ischemia~., data=tmp[,c("clinically.Ischemia", feature.selected1)])
model.step = step(model, direction="backward")
fit = roc(model.step$y, predict(model.step))

model2 = glm(clinically.Ischemia~., data=tmp[,c("clinically.Ischemia", clinical.asso)])
model2.step = step(model2, direction="backward")
fit2 = roc(model2$y, predict(model2.step))

plot(fit)
plot(fit2,add=T)

# selectCox <- function(formula, data, rule = "aic") {
#   require("rms")
#   require("prodlim")
#   fit <- cph(formula, data, surv = TRUE)
#   bwfit <- fastbw(fit, rule = rule)
#   if (length(bwfit$names.kept) == 0) {
#     newform <- reformulate("1", formula[[2]])
#     newfit <- prodlim(newform, data = data)
#   } else{
#     newform <- reformulate(bwfit$names.kept, formula[[2]])
#     newfit <- cph(newform, data, surv = TRUE)
#   }
#   out <- list(fit = newfit,In = bwfit$names.kept)
#   out$call <- match.call()
#   class(out) <- "selectCox"
#   out
# }
# 
# model.selecCOx = selectCox(
#   formula = Surv(mi_time, inz_mi) ~  .,
#   data = data.frame(tmp[,c("inz_mi", "mi_time", metabo.asso2, clinical)]),
#   rule = "aic")