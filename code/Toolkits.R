# Functions could be used in the analysis
# 
# Author: tao.xu
###############################################################################

#function(data, subset = NULL, feature, disease, files){
#	if(is.null(subset)){
#		subset = dim(S4)[1]
#	}
#	for(i in 1:length(diseases)){
#		chars = apply(S4[subset,feature], 2, function(x) tapply(x, INDEX = as.factor(S4[subset, diseases[i]]), mean, na.rm= T))
#		write.csv(t(chars), file = paste(names(diseases)[i],".csv", sep = ""))
#	}
#}
#testfun = function(x, index, FUN){
# 	tapply(x, INDEX = index, FUN)
#}

preprocess=function(data,Metabolites)
{
	####log transform
	tmp=as.matrix(data[,Metabolites])
	tmp[which(tmp==0)]=NA
	data[,Metabolites]=tmp
	#data[,Metabolites]=log2(data[,Metabolites])
	####interpolation process
	tmp=data[,Metabolites];index1=NULL;index2=NULL;index=NULL
	for(i in 1:dim(tmp)[2]){
		tmp[which(is.na(tmp[,i])),i]=mean(tmp[,i],na.rm=T)
		index1=which(tmp[,i]>mean(tmp[,i])+3*sd(tmp[,i]))
		index2=which(tmp[,i]<mean(tmp[,i])-3*sd(tmp[,i]))
		index=unique(c(index,index1,index2))
	}
	data[,Metabolites]=tmp
	return(data)
}


characteristics = function(data , factor, d, ...)
{
	chars = sapply(
			data,  
			function(x) {
				m = tapply(x, INDEX = factor, mean, ...)
				s = tapply(x, INDEX = factor, sd, ...)
        print(m)
        print(s)
				m = round(m , digits = d)
				s = round(s, digits = d)
				return(paste(m, " (", s, ")", sep = ""))
			}
	)
	return(chars)
}

characteristics.disc <- function(data, factor)
{
	rst = sapply(data, function(x) tapply(x,INDEX = factor, function(a) table(a)/length(a)))
	return(rst)
}

## Correlation network analysis tools
## partial correlation calculation
cor.partial <- function(data, variables, ...)
{
	P = cor( data[,variables], ... )
	P.invers = solve( P )
	P.diag = diag( P.invers )
	Z = - P.invers / sqrt(P.diag %*% t(P.diag))
	diag(Z) = 1
	return(Z)
}

#correlation matrix to linkage pairs
cor2link = function(Zeta, threshold)
{
	links = NULL;
	for(i in 1:(dim(Zeta)[1]-1)){
		for(j in (i+1):dim(Zeta)[1]){
			if( abs(Zeta[i,j])  >= threshold){
				tmp = c( colnames( Zeta )[i], colnames( Zeta )[j], 
						Zeta[i,j])
				links = rbind(links, tmp)
			}
		}
	}
	return(links)
}

## convert number to coordinates in the matrix
## input 	num: the index number
##			ncol: number of columns
##			nrow: number of rows, by default equals to ncol
## output	vector with two elements indicating the (x,y) coordinate of the corresponding index

num2coord=function(num,ncol, nrow=ncol)
{
	y=ceiling(num/ncol)
	x=num%%ncol
	if(x==0){x=ncol}
	return(c(x,y))
}
#	Differential correxpression
#	The method was described in Cho, Sung, Jihun Kim, and Ju Kim. ?Identifying Set-wise Differential Co-expression in Gene Expression Microarray Data.BMC Bioinformatics 10, no. 1 (2009): 109.
#	Description of parameters:
#	cor1, cor2: two pearson correlation values
#	N1, N2: the number of samples used for the calculation of cor1 and cor2
diffcorr <- function(cor1, cor2, N1, N2)
{
	Z1 = 0.5*log((1+cor1)/(1-cor1))
	Z2 = 0.5*log((1+cor2)/(1-cor2))
	p = 1 - pnorm(abs((Z1-Z2) / sqrt(1/(N1-3) + 1/(N2-3))))
	return(p)
}

## coversion of parameters between S4 and F4
#function(parameters , optional = "l"){
#	substr(parameters, 1, 2) <- optional
#}

#############	residual calculation	##########
residue<-function(data,Metabolites,adj, control_group)
{
	data[,Metabolites]=scale(log(data[,Metabolites]))
  tmp=NULL;
	tmp=cbind(data[,adj])
	adj_f=paste(adj,collapse='+')
	for(i in 1:length(Metabolites)){
		model=lm(as.formula(paste(Metabolites[i],adj_f,sep="~")),data[control_group,])	
		data[,Metabolites[i]]=data[,Metabolites[i]]-predict(model,tmp)
	}
	return(data[,Metabolites])
}


# regression analysis
logisticRegression = function(meta , disease, valid_measures , feature.cont, feature.disc, metalog = TRUE, ...)
{
	rst = NULL 
	fdisc = NULL; fcont = NULL
	for (f in feature.disc){
		fdisc = cbind(fdisc, as.factor(meta[,f]))
	}
	for(m in valid_measures){
		if(metalog) { 
			metabolite = log(meta[,m])
			fcont = log(meta[ , feature.cont])
		}
		else { 
			metabolite = meta[,m]
			fcont = meta[ , feature.cont]
		}
		data = data.frame (disease, metabolite , fcont, fdisc)
		data = na.omit(data)
		model = glm(as.factor(disease) ~ ., data = data, family = binomial(link = "logit"), na.action = na.omit , ...)
		rst = rbind(rst , summary(model)$coefficients[2,])
	}
	rst = data.frame(rst, p.adjust(rst[,4], method = "BH"))
	rownames(rst) = valid_measures
	return(rst)
}

#########	longitudinal analysis using logistic regression	################
Comparison.prospective<- function(baseline, feature, metabo, adj, subset)
{
#	baseline	---		data frame of variables at baseline 
#	feature		---		list or matrix indicating the disease state in at different time points
#	metabo		---		names of metabolites
#
	rst = NULL
	
	for(i in 1:length(metabo)){
		
		data = data.frame(log(baseline[, metabo[i]]), baseline[,adj])
		
		model = glm(interaction(feature[,1], feature[,2]) ~ . , data = data, subset = subset,  family = binomial(link = "logit"))
		
		rst = rbind(rst, summary(model)$coefficients[2, ])
		
	}
	print(dim(rst)); print(i)
	rownames(rst) = metabo
	
	return (rst)
	
}

###############	cross validation for cox regression	################
theta.fit <- function(x, y, ...) 
{
	#d = data.frame(y, x)
	#print(dim(d))
	#print(colnames(d))
	coxph(y ~  ., data=x)
}
theta.predict <- function(fit, x)
{
	#if(is.null(dim(x))) x=t(x)
	#dim(x)
	#print(colnames(x))
	value=predict(fit,newdata= as.data.frame(x) , type="risk")
	return(value)
}

crossval.cox = function (x, y, theta.fit, theta.predict, ..., ngroup = n) 
{
	call <- match.call()
	#x <- as.matrix(x)
	n = dim(x)[1]
	ngroup <- trunc(ngroup)
	if (ngroup < 2) {
		stop("ngroup should be greater than or equal to 2")
	}
	if (ngroup > n) {
		stop("ngroup should be less than or equal to the number of observations")
	}
	if (ngroup == n) {
		groups <- 1:n
		leave.out <- 1
	}
	if (ngroup < n) {
		leave.out <- trunc(n/ngroup)
		o <- sample(1:n)
		groups <- vector("list", ngroup)
		for (j in 1:(ngroup - 1)) {
			jj <- (1 + (j - 1) * leave.out)
			groups[[j]] <- (o[jj:(jj + leave.out - 1)])
		}
		groups[[ngroup]] <- o[(1 + (ngroup - 1) * leave.out):n]
	}
	u <- NULL
	cv.fit <- rep(NA, n)
	for (j in 1:ngroup) {
		u <- theta.fit(x[-groups[[j]], ], y[-groups[[j]],])
		cv.fit[groups[[j]]] <- theta.predict(u, x[groups[[j]],])
	}
	if (leave.out == 1) 
		groups <- NULL
	return(list(cv.fit = cv.fit, ngroup = ngroup, leave.out = leave.out, 
					groups = groups, call = call))
}


## likelihood ratio test
likeli.test <- function(model1, model0)
{
	D = -2*logLik(model0)[1]+2*logLik(model1)[1]
	pvalue=pchisq(D, df = abs(model1$df.residual-model0$df.residual))
	return(data.frame(D,pvalue))
}

# Define the function
ggd.qqplot = function(pvector, main=NULL, ...) {
	o = -log10(sort(pvector,decreasing=F))
	e = -log10( 1:length(o)/length(o) )
	plot(e,o,pch=19,cex=1, main=main, ...,
			xlab=expression(Expected~~-log[10](italic(p))),
			ylab=expression(Observed~~-log[10](italic(p))),
			xlim=c(0,max(e)), ylim=c(0,max(o)))
	lines(e,e,col="red")
}

#pair-wise ratios among the data 
concen2ratio = function(data){
	rst = NULL;
	for(i in 1:dim(data)[2]){
		tmp =  t(apply(data, 1, function(x) x[i]/x[-i]))
		colnames(tmp) = paste(colnames(data)[i], colnames(data)[-i], sep = ".")
		rst = cbind(rst, tmp)
	}
	return(rst)
}

#detect outliers
outlier = function(x){
	Mean = mean(x, na.rm = T)
	SD = sd(x,na.rm =T)
	which(x>Mean+5*SD|x<Mean-5*SD)	
}

##Generate sampless for K fold cross validation
f_K_fold <- function(Nobs,K=5){
     rs <- runif(Nobs)
     id <- seq(Nobs)[order(rs)]
     k <- as.integer(Nobs*seq(1,K-1)/K)
     k <- matrix(c(0,rep(k,each=2),Nobs),ncol=2,byrow=TRUE)
     k[,1] <- k[,1]+1
     l <- lapply(seq.int(K),function(x,k,d) 
         list(train=d[!(seq(d) %in% seq(k[x,1],k[x,2]))],
              test=d[seq(k[x,1],k[x,2])]),k=k,d=id)
     return(l)
}

## Backward selection for cox model
selectCox <- function(formula, data, rule = "aic") {
  require("rms")
  require("prodlim")
  fit <- cph(formula, data, surv = TRUE)
  bwfit <- fastbw(fit, rule = rule)
  if (length(bwfit$names.kept) == 0) {
    newform <- reformulate("1", formula[[2]])
    newfit <- prodlim(newform, data = data)
  } else{
    newform <- reformulate(bwfit$names.kept, formula[[2]])
    newfit <- cph(newform, data, surv = TRUE)
  }
  out <- list(fit = newfit,In = bwfit$names.kept)
  out$call <- match.call()
  class(out) <- "selectCox"
  out
}

