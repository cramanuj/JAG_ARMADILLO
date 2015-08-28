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
NumericVector rowSums(NumericMatrix x) {
	
	int nrow = x.nrow(), ncol = x.ncol();
	NumericVector out(nrow);
	
	for (int i = 0; i < nrow; i++) {
		double total = 0;
		for (int j = 0; j < ncol; j++) {
			total += x(i, j);
		}
		out[i] = total;
	}
	return out;
}

// [[Rcpp::export]]
inline int randWrapper(const int n) { return floor(unif_rand()*n); }

// [[Rcpp::export]]
Rcpp::NumericVector randomShuffle(Rcpp::NumericVector a) {
    Rcpp::NumericVector b = Rcpp::clone(a);
    std::random_shuffle(b.begin(), b.end(), randWrapper);	
    return b;
}

// [[Rcpp::export]]

double jaguar(double Eps, double Tau, NumericVector k, NumericMatrix Y, NumericVector snp, NumericMatrix R){	
	
	double nobs = Y.ncol();
	double mG = mean(snp);
	double snp2 = sum((snp-mG)*(snp-mG));
	double kmax = R.ncol();
	NumericVector Yhat(kmax);
	NumericVector Ubta(nobs);
	NumericMatrix Ugam(kmax,nobs);
	NumericVector Ug(kmax);
	NumericVector mbeta(nobs);
	NumericVector mgamma(nobs);
	NumericVector V1(nobs); NumericVector V2(nobs);
	NumericVector covmat(nobs); NumericVector vgamma(nobs);
	NumericVector vgamma1(nobs);
	NumericVector covmat1(nobs); NumericVector covmat2(nobs);
	for(int i=0; i<nobs; i++){
		Yhat = Y(_,i);
		double k_ind = k[i];
		NumericVector R_vec = R(i,_);
		double Ysum = sum(Yhat);
		double G = snp[i]-mG;
		V1[i] = (Eps+(k_ind-1)*Tau) / ( (Eps*Eps)+k_ind*Eps*Tau );
		V2[i] = - (Tau)/(Eps*Eps+k_ind*Eps*Tau);
		double sbeta =  G * (V1[i] + (k_ind-1)*V2[i])*Ysum;
		for(int t=0; t<kmax; t++){
			Ug[t] = R_vec[t] * ( (V1[i]*Yhat[t])+(V2[i]*(Ysum-Yhat[t])) ) * G;
		}
		Ugam(_,i)=Ug;
		Ubta[i] = sbeta;
		mbeta[i] =  G * G * (V1[i]+(k_ind-1)*V2[i]) * k_ind;
		mgamma[i] = k_ind * G * G *V1[i];
		covmat[i] = sqrt(k_ind) * G * G * (V1[i]+(k_ind-1)*V2[i]);
		covmat1[i] = k_ind *G*G*(V1[i]+(k_ind-1)*V2[i]);
		covmat2[i] = k_ind *G*G*V1[i];
		
		vgamma[i] = k_ind * G * G * (V1[i]*V1[i]+(k_ind-1)*V2[i]*V2[i]);

	}
	
	NumericVector Ugamma_tmp(kmax);
	for(int j=0; j<kmax; j++){
		Ugamma_tmp[j] = sum(Ugam(j,_))*sum(Ugam(j,_));
	}

	double Ugamma = 0.5 * sum(Ugamma_tmp);
	double U2beta = sum(Ubta)*sum(Ubta);
	double meanB = sum(mbeta);
	double varB = 2 * meanB * meanB;
	double meanG = 0.5 * sum(mgamma);

//	double cov = (sum(covmat)*sum(covmat));
	double cov = sum(covmat2) * sum(covmat1);
	double varG = 0.5 * sum(vgamma) * snp2;
	
	Rcout<<"meanB: "<<meanB<<" varB: "<<varB<<" meanG: "<<meanG<<" varG: "<<varG<<" cov: "<<cov<<endl;
	double agam = (varB - cov)/(varB+varG-2*cov);
	double abeta = (varG - cov)/(varB+varG-2*cov);
	double Upsi = abeta*U2beta + agam*Ugamma;
	double meanPSI = (abeta*meanB)+ (agam*meanG);
	double varPSI = (abeta*varB*abeta) + (agam*agam*varG);
	double a1 = varPSI/(2*meanPSI);
	double a2 = (2*meanPSI*meanPSI)/varPSI;
	double scaled_upsi = Upsi/a1;
	return scaled_upsi;
}
