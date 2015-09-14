/*
 Author: Chaitanya Acharya
 Date: Aug 27, 2015
 ***RCPP ARMADILLO version of JAGUAR***

 Returns a p-value indicating the significance of association between a gene-SNP pair
 Our joint score test statistic is computed as --
 
 U_\psi = Y^t . V^{-1} . (0.5 a_\gamma XX^t + a_\beta GG^t) . V^{-1} . Y
 
 Arguments: Y -> Residuals from the model
			      k -> vector indicating number of tissues per observation
			      snp -> vector of mean centered genotypes 
			      X-> Stacked matrix of genotypes of dimension k*nobs by kmax where kmax is max(k)
			      Eps -> \hat{\Epsilon} from the model
			      Tau -> \hat{\Tau} from the model
 
 Returns: P-value
 
 UPDATED: Sept 14, 2015
*/

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <stdint.h>
#include <vector>
#include <functional>
#include <algorithm>
#include <iostream>

using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

SEXP jaguar_ARMA(const arma::mat Yhat, const arma::vec k, const arma::vec snp, const arma::mat X, double Eps, double Tau){
	
	// Initialize the parameters
	int ksum = sum(k); int kmax = max(k); int nobs = snp.size();
	arma::mat iV(ksum,ksum);				// V^{-1} variance-cov matrix of Y
	arma::mat G(ksum,1);					// G matrix (function of mean-centered genotypes)
	iV.fill(0);  G.fill(0);

	int row_st=0; int col_st=0; int row_end=0; int col_end=0; int nt = 0;
	int row_st_x=0; int row_end_x=0;
	
	//Loop over every observation
	for(int i=0; i<nobs; i++){
		double V1 = (Eps+(k[i]-1)*Tau) / ( (Eps*Eps)+k[i]*Eps*Tau );
		double V2 = - (Tau)/(Eps*Eps+k[i]*Eps*Tau);
		arma::mat id_iV(k[i],k[i]); id_iV.fill(V2);
		arma::vec V1_vec = rep(V1,k[i]);
		arma::vec iG_rep = rep(snp[i],k[i]);
		id_iV.diag() = V1_vec;
		nt+=k[i];
		row_end = col_end = nt;
		row_end_x = nt-1;
		G.submat(row_st_x,0,row_end_x,0) = iG_rep;
		row_st_x = row_end_x+1;
		iV.submat(row_st,col_st,row_end-1,col_end-1)=id_iV;
		row_st = col_st = col_end;
	}	
	// Compute variance-cov matrix and the optimal weights

	double meanB = trace(trans(G)*iV*G);
	double meanG = 0.5 * trace(trans(X)*iV*X);
	double varB = 2 * trace(trans(G)*iV*G*trans(G)*iV*G);
	double varG = 0.5 * trace(trans(X)*iV*X*trans(X)*iV*X);
	double cov = trace(trans(G)*iV*X*trans(X)*iV*G);
	double a_gam = (varB-cov)/(varB+varG-2*cov); 
	double a_beta = (varG-cov)/(varB+varG-2*cov);

	// Compute the score test statistic and the corresponding p-value using
	// Satterthwaite method
	
	const arma::mat U_psi = trans(Yhat)*iV*( (a_beta*G*trans(G)) + (0.5*a_gam*X*trans(X)) )*iV*Yhat;
	double meanPSI = (a_beta*meanB)+ (a_gam*meanG);
	double varPSI = (a_beta*varB*a_beta) + (a_gam*a_gam*varG);
	double a1 = varPSI/(2*meanPSI);
	double a2 = (2*meanPSI*meanPSI)/varPSI;
	double scaled_upsi = U_psi(0,0) / a1;
	return wrap(1 - R::pchisq(scaled_upsi,a2,1,0));
}


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

double jaguar_BALANCED(double Eps, double Tau, double k, NumericMatrix Y, NumericVector snp){	
	
	double nobs = Y.ncol();
	NumericVector Yhat(k);
	double V1 = (Eps+(k-1)*Tau) / ( (Eps*Eps)+k*Eps*Tau);
	double V2 = - (Tau)/(Eps*Eps+k*Eps*Tau);
	NumericVector Ug(k); 
	NumericVector Ubta(nobs);
	NumericMatrix Ugam(k,nobs);
	double mG = mean(snp);
	
	for(int i=0; i<nobs; i++){
		Rcpp::checkUserInterrupt();
		Yhat = Y(_,i);
		double G = snp[i]-mG;
		double Ysum = sum(Yhat);
		double sbeta =  (G * Ysum)/(Eps + k * Tau);
		for(int t=0; t<k; t++){
			Ug[t] = ( (V1*Yhat[t])+(V2*(Ysum-Yhat[t])) ) * G;
		}
		Ugam(_,i)=Ug;
		Ubta[i]=sbeta;
	}
	
	double Ugamma = 0.5 * sum(rowSums(Ugam)*rowSums(Ugam));
	double U2beta = sum(Ubta)*sum(Ubta);
	double snp2 = sum((snp-mG)*(snp-mG));
	double meanB = snp2 * (V1+(k-1)*V2) * k;
	double varB = 2 * meanB * meanB;
	double meanG = 0.5 * sum( (snp-mean(snp))*(snp-mean(snp))*V1*k );
	double varG = 0.5 * snp2 * snp2 * (V1*V1+(k-1)*V2*V2) * k;
	double cov = (snp2 * (V1+(k-1)*V2))*(snp2 * (V1+(k-1)*V2))*k;
	double agam = (varB - cov)/(varB+varG-2*cov);
	double abeta = (varG - cov)/(varB+varG-2*cov);
	double Upsi = abeta*U2beta + agam*Ugamma;
	double meanPSI = (abeta * meanB) + (agam * meanG); 
	double varPSI = (agam * agam * varG) + (abeta * abeta * varB);
	double a1 = varPSI/(2*meanPSI);
	double a2 = (2*meanPSI*meanPSI)/varPSI;
	double scaled_upsi = Upsi/a1;
	double pval = 1 - R::pchisq(scaled_upsi,a2,1,0);
	if(pval==0) pval = 2e-16;
	return pval;
}

// [[Rcpp::export]]

double jaguarNEW(const double Eps, const double Tau, const arma::vec k, const arma::mat& Y, const arma::vec snp, const arma::mat& R){
	
	double nobs = Y.n_cols;
	double kmax = R.n_cols;
	double Ub=0.0, meanB=0.0, meanG=0.0;
	
	arma::vec Yhat(kmax);
	arma::mat Ugam(kmax,nobs,arma::fill::zeros);
	arma::mat varGMAT(kmax,kmax,arma::fill::zeros);
	arma::rowvec Ysum = sum(Y,0);
	
	for(int i =0; i<nobs; i++){
		double G = snp[i]-mean(snp);
		double V1 = (Eps+(k[i]-1)*Tau) / ( (Eps*Eps)+k[i]*Eps*Tau );
		double V2 = - ( (Tau)/(Eps*Eps+k[i]*Eps*Tau) );
		arma::mat V(kmax,kmax); V.fill(V2);
		arma::vec V1_vec = rep(V1,kmax); V.diag() = V1_vec;
		arma::rowvec R_vec = R.row(i);
		arma::mat RR = trans(R.row(i)) * R.row(i);
		varGMAT+= G*G*(RR % V);
		arma::colvec Yhat = Y.col(i);
		Ugam.col(i) = (trans(R_vec)*G) % ( (V1*Yhat) + (V2*(Ysum[i]-Yhat)) );
		Ub+=G * (V1 + (k[i]-1)*V2)*Ysum[i];
		meanB+= G * G * (V1+(k[i]-1)*V2) * k[i];
		meanG+=k[i] * G * G * V1;
	}
	
	const double U2beta = Ub*Ub;
	const double Ugamma = 0.5 * accu( sum(Ugam,1) % sum(Ugam,1) );
	const double varB = 2 * meanB * meanB;
	meanG = 0.5 * meanG;
	const double varG = 0.5 * accu(varGMAT%varGMAT);
	const double cov = accu( sum(varGMAT,0)%sum(varGMAT,0) );
	const double agam = (varB - cov)/(varB+varG-2*cov);
	const double abeta = (varG - cov)/(varB+varG-2*cov);
	const double Upsi = abeta*U2beta + agam*Ugamma;
	const double meanPSI = (abeta*meanB)+ (agam*meanG);
	const double varPSI = (abeta*varB*abeta) + (agam*agam*varG);
	const double a1 = varPSI/(2*meanPSI);
	const double a2 = (2*meanPSI*meanPSI)/varPSI;
	const double scaled_upsi = Upsi/a1;
	double jag_pval = 1 - R::pchisq(scaled_upsi,a2,1,0);
	if(jag_pval==0) jag_pval = 2e-16;
	return jag_pval;
}