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
			X -> Stacked matrix of genotypes of dimension k*nobs by kmax where kmax is max(k)
			Eps -> \hat{\Epsilon} from the model
			Tau -> \hat{\Tau} from the model
 
 Returns: P-value
 
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

SEXP jagARMA(NumericMatrix Y, NumericVector k, NumericVector snp, NumericMatrix X, double Eps, double Tau){
	
	// Initialize the parameters
	int ksum = sum(k); int kmax = max(k); int nobs = snp.size();
	arma::mat bigV(ksum,ksum);		// V^{-1} variance-cov matrix of Y
	arma::mat bigG(ksum,1);			// G matrix (function of mean-centered genotypes)
	bigV.fill(0);  bigG.fill(0);

	int row_st=0; int col_st=0; int row_end=0; int col_end=0; int nt = 0;
	int row_st_x=0; int row_end_x=0;
	
	//Loop over every observation
	for(int i=0; i<nobs; i++){
		double V1 = (Eps+(k[i]-1)*Tau) / ( (Eps*Eps)+k[i]*Eps*Tau );
		double V2 = - (Tau)/(Eps*Eps+k[i]*Eps*Tau);
		arma::mat iV(k[i],k[i]); iV.fill(V2);
		arma::vec V1_vec = rep(V1,k[i]);
		arma::vec iG_rep = rep(snp[i],k[i]);
		iV.diag() = V1_vec;
		nt+=k[i];
		row_end = col_end = nt;
		row_end_x = nt-1;
		bigG.submat(row_st_x,0,row_end_x,0) = iG_rep;
		row_st_x = row_end_x+1;
		bigV.submat(row_st,col_st,row_end-1,col_end-1)=iV;
		row_st = col_st = col_end;
	}
	const arma::mat bigX = Rcpp::as<arma::mat>(X);								// X matrix (function of genotypes)
	const arma::mat GGt = bigG*trans(bigG);										// Transposing big matrices
	arma::mat XXt = bigX*trans(bigX);											// Transposing big matrices
	
	// Compute variance-cov matrix and the optimal weights
	double meanB = trace(bigV*GGt);
	double meanG = 0.5 * trace(bigV*XXt);
	double varB = 2 * trace(bigV*GGt*bigV*GGt);
	double varG = 0.5 * trace(bigV*XXt*bigV*XXt);
	double cov = trace(bigV*GGt*bigV*XXt);
	double a_gam = (varB-cov)/(varB+varG-2*cov); 
	double a_beta = (varG-cov)/(varB+varG-2*cov);

	// Compute the score test statistic and the corresponding p-value using
	// Satterthwaite method
	
	const arma::mat Yhat = Rcpp::as<arma::mat>(Y);
	
	arma::mat U_psi = trans(Yhat)*bigV*((a_beta*GGt) + (0.5*a_gam*XXt))*bigV*Yhat;
	double meanPSI = (a_beta*meanB)+ (a_gam*meanG);
	double varPSI = (a_beta*varB*a_beta) + (a_gam*a_gam*varG);
	double a1 = varPSI/(2*meanPSI);
	double a2 = (2*meanPSI*meanPSI)/varPSI;
	double scaled_upsi = U_psi(0,0) / a1;
	double jag_pval = 1 - R::pchisq(scaled_upsi,a2,1,0);	
	
	return wrap(jag_pval);

}
