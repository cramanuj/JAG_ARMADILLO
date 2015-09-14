lib.list=c("lme4","Rcpp","RcppArmadillo","plyr","rbenchmark")
for(i in 1:length(lib.list)){
	if(any(installed.packages()[,1]==lib.list[i])){
		library(lib.list[i],character.only=T)}else{
			source("http://bioconductor.org/biocLite.R");
			biocLite(lib.list[i]);
			library(lib.list[i],character.only=T)};
}

## Load the cpp programs
sourceCpp("jaguar.cpp")
sourceCpp("jaguar_ARMA.cpp")

## Initiate parameters
nobs = 150; k = 10;tau = 1;eps = 1;PVEg = 0;bta = 0;maf = 0.10;
miss_ind = 25; miss_k = 25;

## Generate datas
gamma = ((eps+tau)*PVEg) / (100-PVEg);
mu = rep(1,k); snp = rbinom(nobs,2,maf);
bkg = rnorm(k,0,gamma);
Y = do.call("rbind",lapply(1:nobs,function(i) mu + bta*snp[i] + rnorm(1,0,tau) + bkg*snp[i] + rnorm(k,0,eps)))

## Simulate missingNESS
miss_ind = (nobs*miss_ind)/100
if(miss_ind!=0){
	missOBS = sample(1:nobs,miss_ind,replace=F)
	missK = replicate(length(missOBS),sample(1:k,round(miss_k*k/100),replace=F),simplify=F)
	missV = cbind(rep(missOBS,each=round(miss_k*k/100)),as.vector(replicate(length(missOBS),sample(1:k,round(miss_k*k/100),replace=F),simplify=T)))
	Y[missV]<-NA
}

## Simulate data
snp = snp-mean(snp)
dataU = data.frame("IND"=as.factor(rep(1:nobs,each=k)),"Gene"=as.vector(t(Y)),"Geno" = rep(snp,each=k),"Tissue"=as.factor(rep(1:k,nobs)))
k_new = k - (apply(Y,1,function(x) sum(is.na(x))))
R = 1 - t(apply(is.na(Y),1,as.numeric));
YnewU = apply(Y,1,function(x) x-colMeans(Y,na.rm=T));
YnewU[is.na(YnewU)]<-0
Yhat = as.matrix(YnewU[YnewU!=0])
R_tmp = R * snp; R_tmp = split(R_tmp,1:NROW(R_tmp));
X = ldply(R_tmp,function(x) x*diag(length(x)),.id=NULL)
X = as.matrix(X[which(apply(X,1,sum)!=0),])

## Fit the model under the global null (H0: bta = gamma =0)
fit = lmer(Gene~0+Tissue+(1|IND),dataU,REML=F)

## Extract model components
eps = sigma(fit)^2; tau = VarCorr(fit)[[1]][1]

## Benchmarking
cols = c("test","replications","elapsed","relative")
benchmark(jaguarNEW(eps,tau,k_new,YnewU,snp,R),
		      jaguar(eps,tau,k_new,YnewU,snp,R), 
		      jaguar_ARMA(Yhat,k_new,snp,X,eps,tau),
		      jaguar_BALANCED(eps,tau,k,YnewU,snp),columns=cols,replications=1000)


