/*
 Author: Chaitanya Acharya
 Date: Sept 9, 2015
 ***RCPP code for JAGUAR using closed-form solution***
 
 Returns a p-value indicating the significance of association between a gene-SNP pair
 Our joint score test statistic is computed as --
 
 U_\psi = Y^t . V^{-1} . (0.5 a_\gamma XX^t + a_\beta GG^t) . V^{-1} . Y
 
 Arguments: 
 Eps -> \hat{\Epsilon} from the model
 Tau -> \hat{\Tau} from the model
 k -> vector indicating number of tissues per observation
 Y -> Residuals from the model
 snp -> vector of mean centered genotypes 
 R -> Matrix indicating the presence/absence of a sample
 
 Returns: P-value
 
 */

#include <Rcpp.h>
#include <stdint.h>
#include <vector>
#include <functional>
#include <algorithm>
#include <iostream>
#include <ctime>
using namespace Rcpp;
using namespace std;
using std::cerr;
using std::cout;
using std::endl;

// [[Rcpp::export]]
NumericMatrix addMat(NumericMatrix x, NumericMatrix y){
	int nr = x.nrow();
	int nc = x.ncol();
	Rcpp::NumericMatrix xplusy(nr, nc);
	for (int i = 0; i < nr; i++){
		for (int j = 0; j < nc; j++){
			xplusy[nr * j + i] = x[nr*j+i] + y[nr*j+i];
		}
	}
	return xplusy;
}


// [[Rcpp::export]]

SEXP jaguar(double Eps, double Tau, NumericVector k, NumericMatrix Y, NumericVector snp, NumericMatrix R){	
	
	double nobs = Y.ncol();
	double mG = mean(snp);
	NumericVector snpC = snp-mG;
	
	double kmax = R.ncol();
	NumericVector Yhat(kmax);
	NumericVector Ubta(nobs);
	NumericMatrix Ugam(kmax,nobs);
	NumericVector Ug(kmax);
	NumericVector mbeta(nobs);
	NumericVector mgamma(nobs);
	
	NumericVector V(kmax);
	NumericVector VG(kmax);
	NumericMatrix VGmat(kmax,kmax);
	NumericMatrix VGmat_fin(kmax,kmax);
	NumericMatrix initMAT(kmax,kmax); initMAT.fill(0);
	
	for(int i=0; i<nobs; i++){
		Yhat = Y(_,i);
		double k_ind = k[i];
		NumericVector R_vec = R(i,_);
		double Ysum = sum(Yhat);
		double G = snpC[i];
		double V1 = (Eps+(k_ind-1.0)*Tau) / ( (Eps*Eps)+k_ind*Eps*Tau );
		double V2 = - ((Tau)/(Eps*Eps+k_ind*Eps*Tau));
		std::vector<double> V(kmax-1);
		std::fill(V.begin(),V.begin()+(kmax-1),V2);
		V.insert(V.begin(),V1);
		for(int t=0; t<kmax; t++){
			Ug[t] = R_vec[t] * ( (V1*Yhat[t])+(V2*(Ysum-Yhat[t])) ) * G;
			for(int tt=0; tt<kmax; tt++){
				VG[tt] = G * G * (R_vec[t]*R_vec[tt]) * V[tt];
			}
			std::rotate(V.rbegin(),V.rbegin()+1.0,V.rend());
			VGmat(t,_)=VG;
		}
		VGmat_fin = addMat(initMAT,VGmat);
		initMAT = VGmat_fin;
		Ugam(_,i)=Ug;
		Ubta[i] = G * (V1 + (k_ind-1)*V2)*Ysum;
		mbeta[i] =  G * G * (V1+(k_ind-1)*V2) * k_ind;
		mgamma[i] = k_ind * G * G *V1;
	}	
	double U2beta = sum(Ubta)*sum(Ubta);
	double meanB = sum(mbeta);
	double varB = 2 * meanB * meanB;
	double meanG = 0.5 * sum(mgamma);
	NumericVector Ugamma_tmp(kmax);
	NumericVector VGsum(kmax);
	NumericVector covSum(kmax);
	for(int j=0; j<kmax; j++){
		Ugamma_tmp[j] = sum(Ugam(j,_))*sum(Ugam(j,_));
		VGsum[j] = sum(VGmat_fin(j,_)*VGmat_fin(j,_));
		covSum[j] = sum(VGmat_fin(j,_))*sum(VGmat_fin(j,_));
	}
	double varG = 0.5 * sum(VGsum);
	double cov = sum(covSum);
	double Ugamma = 0.5 * sum(Ugamma_tmp);
	double agam = (varB - cov)/(varB+varG-2*cov);
	double abeta = (varG - cov)/(varB+varG-2*cov);
	double Upsi = abeta*U2beta + agam*Ugamma;
	double meanPSI = (abeta*meanB)+ (agam*meanG);
	double varPSI = (abeta*varB*abeta) + (agam*agam*varG);
	double a1 = varPSI/(2*meanPSI);
	double a2 = (2*meanPSI*meanPSI)/varPSI;
	double scaled_upsi = Upsi/a1;
	double jag_pval = 1 - R::pchisq(scaled_upsi,a2,1,0);
	if(jag_pval==0) jag_pval = 2e-16;
	return wrap(jag_pval);
}

