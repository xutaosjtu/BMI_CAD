## random forest
if(!require(randomForest)) install.packages("randomForest")
require(randomForest)
tmp = data[,c(metabo.valid,clinical, other, "clinically.Ischemia")]
tmp = na.omit(tmp)

table(data$clinically.Ischemia)
index.control = which(data$clinically.Ischemia==0)
index.case = which(data$clinically.Ischemia==1)
set.seed(12)
index.train = c(sample(index.control, 110), sample(index.case, 36))
index.test = c(1:nrow(data))[-index.train]

tmp = data[,c("clinically.Ischemia", clinical, other, metabo.valid)]
tmp[which(is.na(tmp), arr.ind = T)]=0

require(rpart)
fit.1 = rpart(as.factor(clinically.Ischemia) ~ ., 
              data = tmp, 
              subset = index.train, 
              method = "class"
              )
plot(fit.1); text(fit.1)
if(!require(ROCR)) install.packages("ROCR")
require(ROCR)
train.pred=predict(fit.1, type="prob", newdata=tmp[index.train,])
pred = prediction(train.pred[,2], tmp$clinically.Ischemia[index.train])
perf.train = performance(pred,  "tpr", "fpr")
auc.train = performance(pred, "auc")
plot(perf.train)

test.pred=predict(fit.1, type="prob", newdata=tmp[index.test,])
pred= prediction(test.pred[,2], tmp$clinically.Ischemia[index.test])
perf.test = performance(pred,  "tpr", "fpr")
auc.test = performance(pred, "auc")
plot(perf.test, add = T, col = "red")
abline(0,1, lty = 2)

require(randomForest)

fit.2 = randomForest(as.factor(clinically.Ischemia) ~ ., 
              data =  tmp,
              subset = index.train,
              na.action = na.omit
)
varImpPlot(fit.2)
train.pred = predict(fit.2, type="prob")
pred = prediction(train.pred[,2], fit.2$y)
perf.train = performance(pred,  "tpr", "fpr")
auc.train = performance(pred, "auc")
plot(perf.train)

test.pred=predict(fit.2, type="prob", newdata=tmp[index.test,])
pred= prediction(test.pred[,2], tmp$clinically.Ischemia[index.test])
perf.test = performance(pred,  "tpr", "fpr")
auc.test = performance(pred, "auc")
plot(perf.test, add = T, col = "red")
abline(0,1, lty = 2)

candidates = rownames(importance(fit.2))[order(importance(fit.2), decreasing = T)][1:20]
metabo.selected = intersect(candidates, metabo.valid)
clinical.selected = intersect(candidates, clinical)

##boosting
require(ada)
fit.3 = ada(as.factor(clinically.Ischemia) ~ ., 
    data = tmp, 
    subset = index.train, 
    bag.frac = 0.2, 
    nu=0.1,
    iter=100
    )
train.pred = predict(fit.3, newdata=tmp[index.train,], type="prob")
pred = prediction(train.pred[,2], tmp$clinically.Ischemia[index.train])
perf.train = performance(pred,  "tpr", "fpr")
auc.train = performance(pred, "auc")
plot(perf.train)

test.pred=predict(fit.3, type="prob", newdata=tmp[index.test,])
pred= prediction(test.pred[,2], tmp$clinically.Ischemia[index.test])
perf.test = performance(pred,  "tpr", "fpr")
auc.test = performance(pred, "auc")
plot(perf.test, add = T, col = "red")
abline(0,1, lty = 2)

## Caret
require(caret)
cl <- makePSOCKcluster(4)
registerDoParallel(cl)
cvCtrl <- trainControl(method = "repeatedcv", repeats = 3,
                       summaryFunction = twoClassSummary,
                       classProbs = TRUE)
model = train(as.factor(clinically.Ischemia) ~ .,
              data = tmp,
              subset = index.train,
              method = "ada",
#               trControl = cvCtrl,
#               metric = "ROC",
              tuneLength = 10)
stopCluster(cl)

train.pred=predict(model$finalModel, newdata = tmp[index.train,], type = "probs")
pred = prediction(train.pred[,2], tmp$clinically.Ischemia[index.train])
perf.train = performance(pred,  "tpr", "fpr")
auc.train = performance(pred, "auc")
plot(perf.train)

test.pred=predict(model$finalModel, newdata = tmp[index.test,], type = "probs")
pred= prediction(test.pred[,2], tmp$clinically.Ischemia[index.test])
perf.test = performance(pred,  "tpr", "fpr")
auc.test = performance(pred, "auc")
plot(perf.test, add = T, col = "red")
abline(0,1, lty = 2)

model = glm(as.factor(clinically.Ischemia)~., 
            data=tmp[,c("clinically.Ischemia", metabo.selected,other)],
            subset = index.train, 
            family = binomial,
            weights = data$weight
            )
test.pred = predict(model, newdata = data[index.test,])
fit = roc(data$clinically.Ischemia[index.test], test.pred)

model2 = glm(as.factor(clinically.Ischemia)~., 
             data=tmp[,c("clinically.Ischemia", clinical.selected, other)],
             subset = index.train,
             family = binomial,
             weights = data$weight
             )
test.pred = predict(model2, newdata = data[index.test,])
fit2 = roc(data$clinically.Ischemia[index.test], test.pred)

model3 = glm(as.factor(clinically.Ischemia)~., 
             data=tmp[,c("clinically.Ischemia", other)],
             subset = index.train,
             family = binomial,
             weights = data$weight
             )
test.pred = predict(model3, newdata = data[index.test,])
fit3 = roc(data$clinically.Ischemia[index.test], test.pred)

pdf("ROC curves_2.pdf")
plot(fit, col = "red")
plot(fit2,add=T, lty = 4)
plot(fit3, add = T, lty = 3)
dev.off()

##variable selection in S2
# tmp = S2[which(S2$subcoho==1|S2$inz_mi==1),c("ctalteru", "ccsex","ctbmi","my.diab","ctsysmm","my.cigreg","my.alkkon","cl_chola","cl_hdla","cl_crp",S2_valid_measures,"Arg.Trp","zz_nr","mi_time", "inz_mi","subcoho","prev_mi","ctantihy")]
# tmp[,S2_valid_measures] = log(tmp[,S2_valid_measures])
# clinical = c("ctalteru", "ccsex","ctbmi","my.diab","ctsysmm","my.cigreg","my.alkkon","cl_chola","cl_hdla","cl_crp")

require(penalized)
model.penal = penalized(
  clinically.Ischemia~., data = tmp[,c("clinically.Ischemia", metabo.asso, clinical.asso)], 
  standardize = T, ##centralize the covariates 
  #steps = "Park", trace = F, ## starting from the largest value of lambda1 untile the specificed lambda1 is reached, the function will return a list of penfit objects, and the chang of coefficients could be visualized by plotpath 
  #positive = T, ## positive restriction to all regression coefficients
  lambda1 = 2, lambda2= 0 #penalization parameter 
)

#### stepwise selection
model = glm(clinically.Ischemia~., data=tmp[,c("clinically.Ischemia", feature.selected1)])
model.step = step(model, direction="backward")
fit = roc(model.step$y, predict(model.step))
model = glm(clinically.Ischemia~., data=tmp[,c("clinically.Ischemia", feature.selected2, other)])
fit = roc(model$y, predict(model))

model2 = glm(clinically.Ischemia~., data=tmp[,c("clinically.Ischemia", clinical.asso)])
model2.step = step(model2, direction="backward")
fit2 = roc(model2$y, predict(model2.step))
clinical.selected = names(model2.step$coefficients)[-1]
model2 = glm(clinically.Ischemia~., data=tmp[,c("clinically.Ischemia", clinical.selected, other)])
fit2 = roc(model2$y, predict(model2))

model3 = glm(clinically.Ischemia~., data = tmp[,c("clinically.Ischemia", other)])
fit3= roc(model3$y, predict(model3))

pdf("ROC curves.pdf")
plot(fit, col = "red")
plot(fit2,add=T, lty = 4)
plot(fit3, add = T, lty = 3)
dev.off()

tmp = data[,c(metabo.asso,clinical.asso, other, "clinically.Ischemia")]
tmp = na.omit(tmp)
cases = which(tmp$clinically.Ischemia==1)
controls = which(tmp$clinically.Ischemia==0)
subset = c(sample(cases,30), sample(controls, 30))
rf = randomForest(x = tmp[subset,c(feature.selected.rf)], y = as.factor(tmp$clinically.Ischemia)[subset])
rfselect = rfcv(tmp[subset,c(metabo.asso,clinical.asso, other)], as.factor(tmp$clinically.Ischemia)[subset])
table(predict(rf, tmp[-subset, ]), tmp$clinically.Ischemia[-subset])

require(gbm)
model.boost = gbm(clinically.Ischemia~.,
    distribution = "bernoulli",
    data = tmp[,c("clinically.Ischemia", feature.selected2, other)],
    )

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