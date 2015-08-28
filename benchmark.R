library(lme4); library(Rcpp); library(RcppArmadillo); library(plyr); library(rbenchmark)

sourceCpp("jaguar.cpp")
sourceCpp("jaguar_ARMA.cpp")

nobs = 3; k = 3;tau = 1;eps = 1;PVEg = 0;bta = 0;maf = 0.10;
snp = c(1,1,2)	
snpC = snp-mean(snp)
gamma = ((eps+tau)*PVEg) / (100-PVEg); mu = rep(1,k); bkg = rnorm(k,0,gamma);
Y = do.call("rbind",lapply(1:nobs,function(i) mu + bta*snp[i] + rnorm(1,0,tau) + bkg*snp[i] + rnorm(k,0,eps)))
Y[2,1] = Y[3,3] <- NA
k_new = k - (apply(Y,1,function(x) sum(is.na(x))))
R = 1 - t(apply(is.na(Y),1,as.numeric));
data = data.frame("IND"=as.factor(rep(1:nobs,each=k)),"Gene"=as.vector(t(Y)),"Geno" = rep(snpC,each=k),"Tissue"=as.factor(rep(1:k,nobs)))
Ynew = apply(Y,1,function(x) x-colMeans(Y,na.rm=T));
Ynew[is.na(Ynew)]<-0
Yhat = as.matrix(Ynew[Ynew!=0])

## Fit the model under the global null (H0: bta = gamma =0)
fit_null = lmer(Gene~0+Tissue+(1|IND),data,REML=F)	
est.eps = sigma(fit_null)^2; est.tau = VarCorr(fit_null)[[1]][1]
X = getME(lmer(Gene~0+Tissue+(1|IND)+(0+Geno|Tissue),data,REML=F),"Z"); 
X = as.matrix(X[,c((nobs+1):(nobs+k))])

## Benchmarking
cols = c("test","replications","elapsed","relative")
benchmark(jaguar(est.eps,est.tau,k_new,Ynew,snp,R), jagARMA(Yhat,k_new,snpC,X,est.eps,est.tau),columns=cols,replications=10000)

