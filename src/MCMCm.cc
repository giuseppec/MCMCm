// MCMCm.cc is C++ code to estimate a LVM with predictor
//
// Alexander Raach
// Ludwig-Maximilians-Universität München
//
//
// Copyright (C) 2005 Alexander Raach
// 

#include "matrix.h"
#include "distributions.h"
#include "stat.h"
#include "la.h"
#include "ide.h"
#include "smath.h"
#include "MCMCrng.h"
#include "MCMCfcds.h"
// #include "Windows.h"
#include "stdio.h"

#include<vector>
#include<iostream>
#include<algorithm>


#include <R.h>           // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts

using namespace SCYTHE;
using namespace std;

// Introduction of matrix functions which appear at the end of this code 
  template <class T>
  Matrix<T>
  inv_diag (const Matrix<T> &a);

  template <class T>
  Matrix<T>
  diag_mult_matrix (const Matrix<T> &a, const Matrix<T> &b);
  
  template <class T>
  Matrix<T>
  matrix_mult_diag (const Matrix<T> &a, const Matrix<T> &b);

  template <class T>
  Matrix<T>
  create_diag (const Matrix<T> &a);

	template <class T>
  Matrix<T>
  band_creation (const Matrix<T> &);

  template <class T>
  Matrix<T>
  band_reconversion (const Matrix<T> &a);
  
  template <class T>
  Matrix<T> 
  band_addition (const Matrix<T> &a, const Matrix<T> &b);

  template <class T>
  Matrix<T> 
  band_chol (const Matrix<T> &a);
  
  template <class T>
  Matrix<T> 
  band_solve (const Matrix<T> &chol, const Matrix<T> &vector);
  
  template <class T>
  Matrix<T> 
  band_solve_backward (const Matrix<T> &chol, const Matrix<T> &y);
  

extern "C"{

void MCMCm (double* sampledata, //const int* samplerow, const int* samplecol,
		  double* Ydata, const int* Yrow, const int* Ycol,											// Datenmatrizen Y, Omega, X
		 const double* Omegadata, const int* Omegacol,
		 const double* Offsdata,
		 const double* Xdata, const int* Xcol,
		 const int* MCMCparameter, const int* seedarray,											
		 const double* beta_sterndata, // const int* beta_sternrow, const int* beta_sterncol,	
		 const double* Beta_sterndata, // const int* Beta_sternrow, const int* Beta_sterncol,	
		 const double* Xdir_par_totaldata, const int* Xdir_par_totalcol,	
		 const double* Kdir_par_totaldata, const int* Kdir_par_totalrow, const int* Kdir_par_totalcol,	
		 const double* nonpardir_par_totaldata, const int* nonpardir_par_totalrow, const int* nonpardir_par_totalcol,
		 const double* X_par_totaldata, const int* X_par_totalcol,	
		 const double* K_par_totaldata, const int* K_par_totalrow, const int* K_par_totalcol,	
		 const double* nonpar_par_totaldata, const int* nonpar_par_totalrow, const int* nonpar_par_totalcol,	
		 const double* v_sterndata, const int* v_sternrow, const int* v_sterncol,	
		 const double* s_sterndata, const int* s_sternrow, const int* s_sterncol,
		 const double* zdata, const int* zspadata,    // abweichung, vorher: const int* zcol,																	// Startwerte für Parametervektoren
		 const double* etadata, const int* etacol,						// z, eta, gamma, beta, theta, nu 
		 const double* gammadata, const int* gammarow, const int* gammacol,
		 const double* betadata, const int* betarow, const int* betacol,		 
		 const double* thetadata, const int* thetarow, const int* thetacol,
		 const double* nudata, const int* nurow, const int* nucol,
		 const double* b0data, const int* b0row, const int* b0col,							// Zusätzlich Startwerte der einzelnen
		 const double* Adata, //const int* Arow, const int* Acol,									// Matrizen b0, A, B, Gamma
		 const double* Bdata, //const int* Brow, const int* Bcol,
		 const int*	ncatdata,																										// Sonstiges, z.B. Anzahl Kategorien
		 const double* tunedata,
		 int* accepts,                     
		 double* DICav, double* DICestimate,
		 const double* gamma_sterndata, const double* Gamma_sterndata,
		 const double* lambda_constraints
//		 const int* OtherConstants																									
		 ) {
	Rprintf("Wir sind im C++ - Code !!!\n");
	// SYSTEMTIME st;
	// SYSTEMTIME en;
        // GetSystemTime(&st);

	// Variablen und Matrizen initialisieren
	// Datenmatrizen Y, Omega, X
	const Matrix<double> Y = r2scythe(*Yrow, *Ycol, Ydata);	
  const Matrix<double> Omega = r2scythe(*Yrow, *Omegacol, Omegadata);
	const Matrix<double> X = r2scythe(*Yrow, *Xcol, Xdata);
	
	const Matrix<double> Xdir_par_total = r2scythe(*Yrow, *Xdir_par_totalcol, Xdir_par_totaldata);
	const Matrix<double> Kdir_par_total = r2scythe(*Kdir_par_totalrow, *Kdir_par_totalcol, Kdir_par_totaldata);
	const Matrix<double> nonpardir_par_total = r2scythe(*nonpardir_par_totalrow, *nonpardir_par_totalcol, nonpardir_par_totaldata);	
	
	const Matrix<double> X_par_total = r2scythe(*Yrow, *X_par_totalcol, X_par_totaldata);
	const Matrix<double> K_par_total = r2scythe(*K_par_totalrow, *K_par_totalcol, K_par_totaldata);
	const Matrix<double> nonpar_par_total = r2scythe(*nonpar_par_totalrow, *nonpar_par_totalcol, nonpar_par_totaldata);


	// MCMCparameter initialisieren
	const int* burnin				= MCMCparameter;																
	const int* mcmc					= MCMCparameter +  1;
	const int* thin					= MCMCparameter +  2;
	const int* lecuyer 			= MCMCparameter +  3;
	const int* lecuyerstream= MCMCparameter +  4;
	const int* verbose 			= MCMCparameter +  5;
	const int* storescores 	= MCMCparameter +  6;
	const int* storez 			= MCMCparameter +  7;
	const int* gm					 	= MCMCparameter +  8;
	const int* mh						= MCMCparameter +  9;
	const int Omega_exists  = *(MCMCparameter + 15);
	const int X_exists      = *(MCMCparameter + 16);
	const int X_par_exists  = *(MCMCparameter + 17);
	const int DIC_exists		=	*(MCMCparameter + 18);
	const int sim_nr				= *(MCMCparameter + 19);
	const int* ncatrow			= MCMCparameter + 20;
	const int* ncatcol			= MCMCparameter + 21;
  const int* zcol         = MCMCparameter + 24;
  const int Xdir_par_exists =  *(MCMCparameter + 25);
  const int* beta_sternrow  = MCMCparameter + 26;
  const int* beta_sterncol = MCMCparameter + 27;
  const int* Beta_sternrow = MCMCparameter + 28;
  const int* Beta_sterncol = MCMCparameter + 29;
  const int* Arow          = MCMCparameter + 30;
  const int* Acol          = MCMCparameter + 31;
  const int* Brow          = MCMCparameter + 32;
  const int* Bcol          = MCMCparameter + 33;
  const int* lambda_col = MCMCparameter + 38;
  const int* lambda_row = MCMCparameter + 39;
  const int* samplerow = MCMCparameter + 40;
  const int* samplecol = MCMCparameter + 41;
  	
	//Rprintf("burnin %d\nmcmc %d\nthin %d\nlecuyer %d\nlecuyerstream %d\nverbose %d\nstorescores %d\n",
	//				*burnin, *mcmc, *thin, *lecuyer, *lecuyerstream, *verbose, *storescores, * gm);

	// Prioris für Parametervektoren beta, gamma, theta initialisieren
	const Matrix<double> beta_stern = r2scythe(*beta_sternrow, *beta_sterncol, beta_sterndata);
	const Matrix<double> Beta_stern = r2scythe(*Beta_sternrow, *Beta_sterncol, Beta_sterndata);
	const Matrix<double> v_stern = r2scythe(*v_sternrow, *v_sterncol, v_sterndata);
	const Matrix<double> s_stern = r2scythe(*s_sternrow, *s_sterncol, s_sterndata);

	// Matrizen für Parametervektoren z, eta, gamma, beta, theta, nu initialisieren (inkl. Startwerten)
	Matrix<double> z = r2scythe(*Yrow, *zcol, zdata);
	Matrix<double> eta = r2scythe(*Yrow, *etacol, etadata);
	Matrix<double> gamma = r2scythe(*gammarow, *gammacol, gammadata);
	Matrix<double> beta = r2scythe(*betarow, *betacol, betadata);
	Matrix<double> theta = r2scythe(*thetarow, *thetacol, thetadata);
	Matrix<double> nu = r2scythe(*nurow, *nucol, nudata);
	
        const Matrix<double> lambdaConstr =  r2scythe(*lambda_row, *lambda_col, lambda_constraints);
	
	// Matrizen b0, A, B, Gamma initialisieren (inkl. Startwerten)
	Matrix<double> b0 = r2scythe(*b0row, *b0col, b0data);
	Matrix<double> A = r2scythe(*Arow, *Acol, Adata);
	Matrix<double> B = r2scythe(*Brow, *Bcol, Bdata);
	// für Gamma wird unten eine Matrix definiert
	// Sonstiges
	const Matrix<int> ncat = r2scythe(*ncatrow, *ncatcol, ncatdata);

  // initialize rng stream
  rng *stream = MCMCpack_get_rng(*lecuyer, seedarray, *lecuyerstream);
  
  // constants 
  const int N = Y.rows(); 					  // number of observations
  const int p = Y.cols();  						// number of response indicators y_j
  const int p1 = *(MCMCparameter + 10);			// number of ordinal indicators y_j
  const Matrix<int> zspa = r2scythe(2, p, zspadata);
  const int p3 = *(MCMCparameter + 22);		// number of metric indicators y_j
  const int p2 = *(MCMCparameter + 23);		// number of poisson-distributed indicators y_j
  const int p12 = (p1 + p2);
  Rprintf("p12 = %d\n",   p12);
  const int m = *(MCMCparameter + 11);		// number of latent factors eta
  const int d = *(MCMCparameter + 12);		// number of parametric direct effects
  const int sdir = *(MCMCparameter + 34);
  const int offset = *(MCMCparameter + 35);
  	Rprintf("Indikator offset= %d \n", offset);
  //const int offsetrow = *(MCMCparameter + 36);
  const int offsetcol = *(MCMCparameter + 37);
  const int q = *(MCMCparameter + 13);		// number of direct effects covariates
  //const int tot_ncut = *(MCMCparameter + 14); // total number of cutpoints
	const int s = *(MCMCparameter + 14);

  const int tot_iter = *burnin + *mcmc;		// total number of sampling iterations
  const int nsamp = *mcmc / *thin;		// number of stored samples
  Rprintf("tot_iter = %d\n",   tot_iter);
  Rprintf("nsamp = %d\n",   nsamp);
	const Matrix<double> tune = r2scythe(p1, 1, tunedata);
	
	const Matrix<double> Offs = r2scythe(*Yrow, offsetcol, Offsdata);
/*	
for (int i=0; i<  Offs.rows(); ++i){
	for(int j=0; j<  Offs.cols(); ++j){
		Rprintf("Matrix dim offs = %f\n", Offs(i,j));
}
}     	
      
		 return; // fehlersuche  
*/		
	
	// Define dynamic matrices for nonparametric direct covariates
	const Matrix<int> 	 	 pdir_nr = !nonpardir_par_total(0,0,0,sdir-1); 
	const Matrix<double> 	adir_par = !nonpardir_par_total(1,0,1,sdir-1);
	const Matrix<double> 	bdir_par = !nonpardir_par_total(2,0,2,sdir-1);
	const Matrix<double> rank_Kdir = !nonpardir_par_total(3,0,3,sdir-1);
	const Matrix<int> 		degrdir  = !nonpardir_par_total(4,0,4,sdir-1);
	const Matrix<int>	 var_cofdir  = !nonpardir_par_total(5,0,5,sdir-1);
	
	Rprintf("a priori (direct):\n");
  for (int i=0; i<adir_par.rows(); ++i) {
		for (int j=0; j<adir_par.cols(); ++j)
			Rprintf("%f  ", adir_par(i,j));
		Rprintf("\n");}Rprintf("\n");	
			
	Rprintf("b priori (direct):\n");
  for (int i=0; i<bdir_par.rows(); ++i) {
		for (int j=0; j<bdir_par.cols(); ++j)
			Rprintf("%f  ", bdir_par(i,j));
	Rprintf("\n");}Rprintf("\n");		
//	return;
	Matrix<double>* Kdir_par = new Matrix<double>[sdir];	
	Matrix<double>* Kdir_par_band = new Matrix<double>[sdir];	
	Matrix<double>* Xdir_par = new Matrix<double>[sdir];
	Matrix<double>* gammadir_par = new Matrix<double>[sdir];
	Matrix<double>* taudir_par = new Matrix<double>[sdir];
	
	for (int tdir=0, tempdir=0; tdir<sdir; ++tdir) {
		Kdir_par[tdir] = Kdir_par_total(0,tempdir, pdir_nr(tdir,0)-1, tempdir+pdir_nr(tdir,0)-1);
		Kdir_par_band[tdir] = band_creation(Kdir_par[tdir]);
		Xdir_par[tdir] = Xdir_par_total(0,tempdir, 				 N-1, tempdir+pdir_nr(tdir,0)-1);
		gammadir_par[tdir] = Matrix<double>(p, pdir_nr(tdir,0));
		taudir_par[tdir] = Matrix<double>(p,1,1, 0.5);				// Standardmäßig wird tau=0.5 gesetzt
		tempdir += pdir_nr(tdir,0);
	}

    Rprintf("zeilen   Xdir_par[0] rows = %d  spalten  = %d\n", Xdir_par[0].rows(), Xdir_par[0].cols());
/*  for (int i=0; i<  Xdir_par[0].rows(); ++i){
	for(int j=0; j<  Xdir_par[0].cols(); ++j){
		Rprintf("Matrix compvar = %f i= %d j= %d \n", Xdir_par[0](i,j),i,j);
	}
}     	
      
		return; // fehlersuche  */
		
		
		// Versuch: aufblähen von Xdir
int Xdir_pardim = p2 * sdir;
Rprintf("Xdir_pardim = %d \n", Xdir_pardim);
Matrix<double>* Xdir_parhelp = new Matrix<double>[Xdir_pardim];
Matrix<double>* Xdir_parhelptrans = new Matrix<double>[Xdir_pardim];

  // Zusatz Samplingstep4
	vector<int> dimpois(p2);
	for (int j=p1; j<p12; ++j){
		int helpdimpois = 0;
		for (int i=0; i<N; ++i){
			helpdimpois = helpdimpois + (int(Y(i,j))+1);
		}
		dimpois[(j-p1)] = helpdimpois;
	} // Ende Zusatz Samplingstep 4

//for (int test =0; test<p2; ++test)
//Rprintf("row = %d,  j= %d\n", dimpois[test], test+p1);
  
/*for (int tt=0; tt<sdir; ++tt){
Rprintf("cols = %d,  tt= %d\n", Xdir_par[tt].cols(), tt);
} */
  
for (int j=p1; j<(p1 + p2); ++j){
  for (int tempsdir=0; tempsdir<sdir; ++tempsdir){ 
  	Xdir_parhelp[(j-p1)*sdir+tempsdir] = Matrix<double>(dimpois[j-p1],Xdir_par[tempsdir].cols() );
  	Xdir_parhelptrans[(j-p1)*sdir+tempsdir] = Matrix<double>(Xdir_par[tempsdir].cols(),dimpois[j-p1] );
  }
}
  

int tdirhelp=0;
for (int j=p1; j<(p1 + p2); ++j){
for (int tempsdir=0; tempsdir<sdir; ++tempsdir){

int rowcount =0;
for (int i=0; i<N; ++i){
  for (int k=0; k<=int(Y(i,j)); ++k){
    for(int l=0; l<Xdir_par[tempsdir].cols(); ++l){
      Xdir_parhelp[tdirhelp](rowcount,l) = Xdir_par[tempsdir](i,l);
/*if(Xdir_par[tempsdir](i,l)>0.2){
            Rprintf("l = %d , tdirhelp =%d , k=%d , int(Y(i,j))=%d, rowcount = %d, Xdir_parhelp[tdirhelp](rowcount,l) = %f \n", l, tdirhelp, k, int(Y(i,j)),rowcount, Xdir_parhelp[tdirhelp](rowcount,l));
          } */
    }
    rowcount = rowcount +1;
  }
}

      tdirhelp = tdirhelp+1;
     }
   }    
   
for (int temptrans=0; temptrans<Xdir_pardim;  ++temptrans){
    Xdir_parhelptrans[temptrans] = !Xdir_parhelp[temptrans];
}



  /*    int testrow0 =  Xdir_parhelp[0].rows();
   int testcol0 =  Xdir_parhelp[0].cols();
   Rprintf("testrow0 = %d,   testcol0 = %d\n", testrow0, testcol0);
   int testrow1 =  Xdir_parhelp[1].rows();
   int testcol1 =  Xdir_parhelp[1].cols();
   Rprintf("testrow1 = %d,   testcol1 = %d\n", testrow1, testcol1); 
   
    for(int l=0; l<Xdir_parhelp[0].cols(); ++l){
   Rprintf("Xdir_parhelp[0]  = %f\n", Xdir_parhelp[0](3,l));  }
       for(int l=0; l<Xdir_parhelp[1].cols(); ++l){
   Rprintf("Xdir_parhelp[1]  = %f\n", Xdir_parhelp[1](2,l));  }
   Rprintf("Xdir_parhelp[1] rows = %d,   cols = %d\n", int(Xdir_parhelp[1].rows()), Xdir_parhelp[1].cols());
  	Rprintf("Xdir_par[0] rows = %d,   cols = %d\n", Xdir_par[0].rows(), Xdir_par[0].cols());  */
		

	
	// Define dynamic matrices for nonparametric indirect covariates
	const Matrix<int> 	 	 p_nr = !nonpar_par_total(0,0,0,s-1); 
	const Matrix<double> 	a_par = !nonpar_par_total(1,0,1,s-1);
	const Matrix<double> 	b_par = !nonpar_par_total(2,0,2,s-1);
	const Matrix<double> rank_K = !nonpar_par_total(3,0,3,s-1);
	const Matrix<int> 		degr  = !nonpar_par_total(4,0,4,s-1);
	const Matrix<int>	 var_cof  = !nonpar_par_total(5,0,5,s-1);
	Matrix<double> gamma_stern  = r2scythe(q, 1, gamma_sterndata);
	Matrix<double> Gamma_stern  = r2scythe(q, q, Gamma_sterndata);
	//Matrix<double> gamma_stern  = r2scythe(q+1, 1, gamma_sterndata);        //Änderung gamma_intercept
	//Matrix<double> Gamma_stern  = r2scythe(q+1, q+1, Gamma_sterndata);
	
	Rprintf("a priori (indirect):\n");
  for (int i=0; i<a_par.rows(); ++i) {
		for (int j=0; j<a_par.cols(); ++j)
			Rprintf("%f  ", a_par(i,j));
		Rprintf("\n");}Rprintf("\n");	
			
	Rprintf("b priori (indirect):\n");
  for (int i=0; i<b_par.rows(); ++i) {
		for (int j=0; j<b_par.cols(); ++j)
			Rprintf("%f  ", b_par(i,j));
	Rprintf("\n");}Rprintf("\n");		
//	return;
	Matrix<double>* K_par = new Matrix<double>[s];	
	Matrix<double>* K_par_band = new Matrix<double>[s];	
	Matrix<double>* X_par = new Matrix<double>[s];
	Matrix<double>* gamma_par = new Matrix<double>[s];
	Matrix<double>* tau_par = new Matrix<double>[s];
	
	for (int t=0, temp=0; t<s; ++t) {
		K_par[t] = K_par_total(0,temp, p_nr(t,0)-1, temp+p_nr(t,0)-1);
		K_par_band[t] = band_creation(K_par[t]);
		X_par[t] = X_par_total(0,temp, 				 N-1, temp+p_nr(t,0)-1);
		gamma_par[t] = Matrix<double>(m, p_nr(t,0));
		tau_par[t] = Matrix<double>(m,1,1, 0.5);				// Standardmäßig wird tau=0.5 gesetzt
		temp += p_nr(t,0);
	}
	
	Rprintf("Parameter Omega_exists: %d\n", Omega_exists);	
	Rprintf("Parameter X_exists: %d\n", X_exists);		
	Rprintf("Parameter Xdir_par_exists: %d\n", Xdir_par_exists);
	Rprintf("Parameter X_par_exists: %d\n", X_par_exists);		
	Rprintf("Parameter DIC_exists: %d\n", DIC_exists);
	Rprintf("Parameter *Yrow: %d\n", *Yrow);	
	Rprintf("Parameter *Ycol: %d\n", *Ycol);	
	Rprintf("Parameter *Omegacol: %d\n", *Omegacol);	
	Rprintf("Parameter *Xcol: %d\n", *Xcol);		
	Rprintf("Number of indicators:   p = %d,   p1 = %d,   p2 = %d,   p3 = %d\n", p, p1, p2, p3);
	Rprintf("Parameter m = %d\n", m);	
	Rprintf("Number of nonparametric direct covariates sdir = %d,   parametric direct covariates d = %d\n", sdir, d);
  Rprintf("Parameter *zcol = %d\n", z.cols());
	Rprintf("Number of nonparametric covariates s = %d,   parametric covariates q = %d\n", s, q);

	if (*mh) Rprintf("Durchführung Metropolis-Hastings-Algorithmus für Schwellenwerte: Ja !!!\n");
	else Rprintf("Durchführung Metropolis-Hastings-Algorithmus für Schwellenwerte: Nein !!!\n");
	if (*gm) Rprintf("Durchführung Grouped Move Step: Ja !!!\n");
	else Rprintf("Durchführung Grouped Move Step: Nein !!!\n");
	
	Matrix<double> Gamma;
	if (X_exists) {
		Gamma = Matrix<double>(m,q);
		for (int r=0; r<m; ++r) 
			for (int c=0; c<q; ++c)
				Gamma(r, c) = gamma(r*(q+1)+1+c,0);
	}		
	else { 
		Gamma = Matrix<double>(1,1);	
		Gamma(0,0) = -999;
	} 

  // storage matrices
  int sample_size = 0;
  Matrix<double> eta_store;
  if (*storescores==1) {
  	eta_store = Matrix<double>(nsamp, N*m);
  	sample_size += N*m;
  }
  Matrix<double> beta_store = Matrix<double>(nsamp, p*(m+d+1));
  sample_size += p*(m+d+1);
  
  
  // NUR ZUM TEST WIRD gamma-intercept gespeichert
  Matrix<double> gamma_intercept_store;
  if (X_exists || X_par_exists) {
  	gamma_intercept_store = Matrix<double>(nsamp, 1);
  	sample_size += 1;
  }  
  // TESTENDE
  
  
  // storage for nonparametric direct effects
    Matrix<double> gammadir_par_store, taudir_par_store;
  int pdir_nr_total = 0;
  if (Xdir_par_exists) {
  	for (int tdir=0; tdir<sdir; tdir++) 
  		pdir_nr_total += pdir_nr(tdir,0);
  	gammadir_par_store = Matrix<double>(nsamp, p*pdir_nr_total);
  	taudir_par_store   = Matrix<double>(nsamp, p*sdir);
  	sample_size += p*pdir_nr_total + p*sdir;
  } 
  
  // storage for indirect effects
  Matrix<double> gamma_store;
  if (X_exists) {
  	gamma_store = Matrix<double>(nsamp, m*q);
  	sample_size += m*q;
  }
  Matrix<double> gamma_par_store, tau_par_store;
  int p_nr_total = 0;
  if (X_par_exists) {
  	for (int t=0; t<s; t++) 
  		p_nr_total += p_nr(t,0);
  	gamma_par_store = Matrix<double>(nsamp, m*p_nr_total);
  	tau_par_store   = Matrix<double>(nsamp, m*s);
  	sample_size += m*p_nr_total + m*s;
  } 
  Matrix<double> theta_store;
  theta_store = Matrix<double>(nsamp, p);
  sample_size += p;
  
  Matrix<double> nu_store;
  int nr_nu = 0;

  if (p1 > 0 && nu.cols() > 3) {
  	for (int i=0; i<p; ++i) {
  		if (ncat(i,0) > 2) {
  			nr_nu += ncat(i,0) - 2;
  			sample_size += ncat(i,0) -2;
  		}
  	nu_store = Matrix<double>(nsamp, nr_nu);
  	}
  }

  Matrix<double> z_store;
  if ((*storez==1 && p1>0) || (*storez==1 && p2>0))  {
  	z_store = Matrix<double>(nsamp, (N*(*zcol - p3)) );
  	sample_size += (N*(*zcol - p3));
  }

  Matrix<double> DIC_store;
  if (DIC_exists) {
  	DIC_store = Matrix<double>(nsamp, 1);
  	sample_size += 1;
  }


  //////////////////////////////////////////////////////
  Rprintf("C++ - Size of sample: nrows = %d   ncols = %d\n", nsamp, sample_size);

 	//Rprintf("Vor Konstantendefinition...\n");
  /////////////////////
  // Constants, etc. //
  /////////////////////

  /// Allgemeine Konstanten/Variablen
  int count = 0;
  
  // allgemein benötigte Variablen

  int iter = 0;
  int i = 0;
  int j = 0;
  int k = 0;
  int l = 0;
  int r = 0;
  int t = 0;
  int w = 0;

  // samplingstep1 
  double sum_zwisch = 0.0; 
  double xi_exp = 0.0;
  double tempprop = 0.0;
  double compgleich = 0.0;
  double compsum = 0.0;
  int compind = 99;
  Matrix<double> nonpartemp;

  // samplingstep2
 Matrix<double> tempz;
  Matrix<double> tempb0;
  Matrix<double> tempB;
  Matrix<double> temptransB;
  Matrix<double> tempA;  

  //Matrix<double> temptheta;
  int rowcount = 0;
  Matrix<double> eta_var = Matrix<double>(m,m,1,0.0);
  Matrix<double> eta_var_chol = Matrix<double>(m,m,1,0.0);
  Matrix<double> tempvarwood = Matrix<double>(m,m,1,0.0);
  Matrix<double> Eta_mu = Matrix<double>(m,1,1,0.0);
  Matrix<double> invtemptheta;
  Matrix<double> tempvarwood2; // Woodbury's identity
  Matrix<double> tempinvthetaB;         

  // samplingstep3
  double average = 0.0;

  // samplingstep4            
  Matrix<double> tempzE;
  Matrix<double> tempU;
  Matrix<double> transtempU;
			//Matrix<double> tempW = Matrix<double>(dimpois[(j-p1)],1,1,0.0);
  Matrix<double> invtempW;
  Matrix<double> tempUWU = Matrix<double>((1+d+m),(1+d+m),1,0.0);
  int rowcountE = 0;
  int rowcountE2 = 0;

  // samplingstep5
  double dof = 0.0;
  double s_inf = 0.0;

  // samplingstep6
  Matrix<double> proposal;
  double prob_likelihood = 0.0;
  double prob_norm = 0.0;

  // Konstanten/Variablen/Vorberechnetes für Samplingschritt 1
  Matrix<double> temp1 = 0, temp2;
  
  Matrix<double> A_trans, B_trans, z_mu_vec, z_logmu_vec_pois, z_mu_vec_pois;
  // Zusatz Samplingstep 1
   // only for poisson distributed variables
	Matrix<double> IndVar = Matrix<double>(N,(z.cols()-(p1+p3)),1,-99); // Matrix für Varianzen, abhängig vom Komponentenindikator
        // Komponentenindikator  (Gewichte, Mittelwerte und Varianzen zur Approximation
        // der Dichte der log Exp(1)-Verteilung)
                   Matrix<double> compweight = Matrix<double>(10,1,1,0.0);           // Vektor der Gewichte
                   Matrix<double> compmean = Matrix<double>(10,1,1,0.0);              // Vektor der Mittelerte
                   Matrix<double> compvar = Matrix<double>(10,1,1,0.0);               // Vektor der Varianzen
                   Matrix<double> invcompvar = Matrix<double>(10,1,1,0.0);               // Vektor der inversen Varianzen
                   Matrix<double> sqrtcompvar = Matrix<double>(10,1,1,0.0);          // Vektor der Standardabweichungen

                   // Gewichte der 10 Komponenten
                   compweight(0,0) = 0.00397;
                   compweight(1,0) = 0.0396;
                   compweight(2,0) = 0.168;
                   compweight(3,0) = 0.147;
                   compweight(4,0) = 0.125;
                   compweight(5,0) = 0.101;
                   compweight(6,0) = 0.104;
                   compweight(7,0) = 0.116;
                   compweight(8,0) = 0.107;
                   compweight(9,0) = 0.088;

                   // Mittelwerte der Komponenten
                   compmean(0,0) = -5.09;
                   compmean(1,0) = -3.29;
                   compmean(2,0) = -1.82;
                   compmean(3,0) = -1.24;
                   compmean(4,0) = -0.764;
                   compmean(5,0) = -0.391;
                   compmean(6,0) = -0.0431;
                   compmean(7,0) = 0.306;
                   compmean(8,0) = 0.673;
                   compmean(9,0) = 1.06;

 		   // Varianzen der Komponenten
                   compvar(0,0) = 4.5;
                   compvar(1,0) = 2.02;
                   compvar(2,0) = 1.1;
                   compvar(3,0) = 0.422;
                   compvar(4,0) = 0.198;
                   compvar(5,0) = 0.107;
                   compvar(6,0) = 0.0778;
                   compvar(7,0) = 0.0766;
                   compvar(8,0) = 0.0947;
                   compvar(9,0) = 0.146;
                   
 		   // Inverse der Varianzen der Komponenten
                   invcompvar(0,0) = 1.0/4.5;
                   invcompvar(1,0) = 1.0/2.02;
                   invcompvar(2,0) = 1.0/1.1;
                   invcompvar(3,0) = 1.0/0.422;
                   invcompvar(4,0) = 1.0/0.198;
                   invcompvar(5,0) = 1.0/0.107;
                   invcompvar(6,0) = 1.0/0.0778;
                   invcompvar(7,0) = 1.0/0.0766;
                   invcompvar(8,0) = 1.0/0.0947;
                   invcompvar(9,0) = 1.0/0.146;
                   
      // Quadratwurzel der Varianzen
                   sqrtcompvar(0,0) = sqrt(4.5);
                   sqrtcompvar(1,0) = sqrt(2.02);
                   sqrtcompvar(2,0) = sqrt(1.1);
                   sqrtcompvar(3,0) = sqrt(0.422);
                   sqrtcompvar(4,0) = sqrt(0.198);
                   sqrtcompvar(5,0) = sqrt(0.107);
                   sqrtcompvar(6,0) = sqrt(0.0778);
                   sqrtcompvar(7,0) = sqrt(0.0766);
                   sqrtcompvar(8,0) = sqrt(0.0947);
                   sqrtcompvar(9,0) = sqrt(0.146);

        vector<double> compprop(10);
        vector<double> logcompprop(10);
        vector<double> explogcompprop(10);
        vector<double> comppropnorm(9);   // 9-dim Vektor der normierten diskreten Dichte
        vector<double> kumcomppropnorm(10);   // 10-dim Vektor der kumulierten
// Ende Zusatz Samplingstep1

// Konstanten/Variablen/Vorberechnetes für Samplingschritt 2
//  Matrix<double> eta_var;
  Matrix<double> eta_random;
	Matrix<double> temp3;
	Matrix<double> EINSER = Matrix<double>(N, 1, 1, 1);
	Matrix<double> tempnonpar2;
	Matrix<double> tempnonparind;
	Matrix<double> tempoffset;
        // Zusatz Samplingstep 2
        // berechnet Dimension von y_i* i=1,...,N
        vector<int> dim(N);
        for (int i=0; i<N; ++i){
            int poiszusatz=0;
            for (int j=p1; j<p12; ++j){
                poiszusatz = poiszusatz +(int(Y(i,j))+1);
            }
            dim[i] = p1 + poiszusatz + p3;
        } // Ende Zusatz Samplingstep2

// Konstanten/Variablen/Vorberechnetes für Samplingschritt 3
// 3.1. Parametric covariates - code is optimized and written only for diffuse priors ?wirklich
	Matrix<double> gamma_mu_r, gamma_r, praed_rest;
	Matrix<double> X_total, gamma_var_chol, gamma_mu_multiplicator;
	if (X_exists) {
		X_total = X;// X_total = cbind(EINSER, X);         //Änderung gamma_intercept
  		gamma_var_chol = cholesky (invpd(Gamma_stern + crossprod(X_total)));
 		gamma_mu_multiplicator = invpd(Gamma_stern + crossprod(X_total));
		// gamma_mu_multiplicator = invpd(Gamma_stern + crossprod(X_total)) * (!X_total);
 	}

 	// INTERCEPT TEST
 	double gamma_intercept = 0.0;
 	// 3.2. Nonparametric covariates
	Matrix<double> solution_vec, gamma_nonpar_mu_r, gamma_nonpar_chol, gamma_nonpar_rand;
	Matrix<double>* X_par_cross = new Matrix<double>[s];
	Matrix<double>* X_par_cross_band = new Matrix<double>[s];
	Matrix<double>* X_par_trans = new Matrix<double>[s];

	for (int t=0; t<s; ++t) {
		X_par_cross[t] = crossprod(X_par[t]);
		X_par_cross_band[t] = band_creation(X_par_cross[t]);
		X_par_trans[t] = !X_par[t];
	}

 	// 3.3. Smoothing parameters tau
	Matrix<double> b_tau;
	double a_tau;

 	// INTERCEPT TEST
 	//double gammadir_intercept = 0.0;
 	// 4.2. Nonparametric direct covariates
	Matrix<double> solutiondir_vec, gammadir_nonpar_mu_j, gammadir_nonpar_chol, gammadir_nonpar_rand, gammadir_nonpar_chol_help,gammadir_nonpar_chol_helpband;
	Matrix<double> tempeta, tempzE2, invtempW2, tempOmega;
	Matrix<double>* Xdir_par_cross = new Matrix<double>[sdir];
	Matrix<double>* Xdir_par_cross_band = new Matrix<double>[sdir];
	Matrix<double>* Xdir_par_trans = new Matrix<double>[sdir];

	for (int t=0; t<sdir; ++t) {
		Xdir_par_cross[t] = crossprod(Xdir_par[t]);
		Xdir_par_cross_band[t] = band_creation(Xdir_par_cross[t]);
		Xdir_par_trans[t] = !Xdir_par[t];
	}

 	//4.3. Smoothing parameters tau
	Matrix<double> bdir_tau;
	
	double adir_tau;	
	
  /// Konstanten/Variablen/Vorberechnetes für Samplingschritt 4
  Matrix<double> U_basis = EINSER;
  Matrix<double> U, beta_mu_j, beta_var_j, beta_j;
  Matrix<double> Beta_stern_j, beta_stern_j;
  Matrix<double> U_cross, U_trans, beta_var_j_chol;
  Matrix<double> praeddir_rest;
  Matrix<double> praeddir_restj;
  Matrix<double> tempdiff;
  
  double theta_inv_j;
  int dimension_beta_j = 1 + d + m;
  // Rprintf("zeilen  dimension_beta_j = %d  \n", dimension_beta_j);
  if (Omega_exists) {
  	U_basis = cbind(U_basis, Omega);
  }



  // Konstanten/Variablen/Vorberechnetes für Samplingschritt 5
  Matrix<double> temp4 = 0, temp5;
	double temp6;
  
  // Konstanten/Variablen/Vorberechnetes für Samplingschritt 6
  Matrix<double> nonpartempmh = Matrix<double>(N, 1, 1, 0.0); 	
  double LO, UP;
  
  // Konstanten/Variablen/Vorberechnetes für MH-Step 
  Matrix<double> A_row, B_row;
  double z_mu;

  // Konstanten/Variablen/Vorberechnetes für Grouped-Move-Step
  Matrix<double> temp7, temp8;
  // A_trans und B_trans sind bereits in "Konstanten für Samplingschritt 1" definiert
  Matrix<double> diff_temp;
  double a=0, b=0, g=0;
  
  // Konstanten/Variablen/Vorberechnetes für DIC-Berechnung
  double DIC = 0.0;
  
  //const double DIC_offset = - N*(p+m)*log(2.0*3.141592654) /2.0;  // von alter Methode (siehe ganz unten) benutzt
 //  double DIC_det = 0.0;
  Matrix<double> DIC_mu, DIC_mu1;
  Matrix<double> DIC_mu2;// = Matrix<double>(m, 1, 1, 0.0);
  Matrix<double> DIC_cov;
  Matrix<double> DIC_identity = eye<double>(m);
  Matrix<double> DIC_ysternmedium = Matrix<double>(N, p, 1, 0.0);
  Matrix<double> DIC_zmedium	= Matrix<double>(N, m, 1, 0.0);     // ZEILE FÜR DIC1 mit latent Variables
  Matrix<double> DIC_etamedium	= Matrix<double>(N, m, 1, 0.0);   // ZEILE FÜR DIC2 ohne latent Variables
  //Matrix<double> DIC_eta	= Matrix<double>(N, m, 1, 0.0);					// ZEILE FÜR DIC2 ohne latent Variables
  Matrix<double> DIC_b0medium	= Matrix<double>(p, 1, 1, 0.0);
	Matrix<double> DIC_Amedium	= Matrix<double>(p, d, 1, 0.0);
	Matrix<double> DIC_Bmedium	= Matrix<double>(p, m, 1, 0.0);
	Matrix<double> DIC_thetamedium = Matrix<double>(p, p, 1, 0.0);
	Matrix<double> DIC_numedium = Matrix<double>(*nurow, *nucol, 1, 0.0);
	
	
 
  ///////////////////
  // Gibbs Sampler //
  ///////////////////
	//Rprintf("Vor Gibbs Sampler....\n");
	
  for (iter=1; iter <= tot_iter; ++iter) {		// Start sampling loop
//	if (iter % 10 == 0)
//		Rprintf("Iteration Nr. %d\n", iter);
// 1. Sample  z_ij | eta_i, theta \ {gamma0, Gamma}, y_ij, omega_i
/*		if (p1>0) {		 // only for ordinal variables
			if (Omega_exists) A_trans = !A;
			B_trans = !B;
			for (j=0; j<p1; ++j) {
				if (Omega_exists)
					z_mu_vec = b0(j,0) + Omega * A_trans(_,j) + eta * B_trans(_,j);
				else
					z_mu_vec = b0(j,0) 												+ eta * B_trans(_,j);


				for (i=0; i<N; ++i)
					z(i,j) = stream->rtnorm_combo(z_mu_vec(i,0), 1.0, nu(j, static_cast<int>(Y(i,j))-1),
																								 		 				nu(j, static_cast<int>(Y(i,j))));
			}
   } // end of sampling step 1 ordinal (z)
 */
		if (p1>0) {		 // only for ordinal variables
			if (Omega_exists) A_trans = !A;
			B_trans = !B;
			for (j=0; j<p1; ++j) {
			 if(Xdir_par_exists){
          nonpartemp = Matrix<double>(N, 1, 1, 0.0); 	
					for (t=0; t<sdir; ++t) {
						nonpartemp = nonpartemp + Xdir_par[t] * !gammadir_par[t](j,_);
					}
				if (Omega_exists)
					z_mu_vec = b0(j,0) + nonpartemp +  Omega * A_trans(_,j) + eta * B_trans(_,j);
				else
					z_mu_vec = b0(j,0) + nonpartemp + 											 eta * B_trans(_,j);
			 }
       else{       // falls keine nonparametrischen direkten Effekte
        if (Omega_exists)
					z_mu_vec = b0(j,0) + Omega * A_trans(_,j) + eta * B_trans(_,j);
				else
					z_mu_vec = b0(j,0) 												+ eta * B_trans(_,j);
       }

       for (i=0; i<N; ++i)
					z(i,j) = stream->rtnorm_combo(z_mu_vec(i,0), 1.0, nu(j, static_cast<int>(Y(i,j))-1),
																								 		 				nu(j, static_cast<int>(Y(i,j))));

        
       }																					 		 				
   } // end of sampling step 1 ordinal (z)
      
  
                
////////////////////////////////////////////////////////////////////////////////
// Samplingstep 1 für poissonverteilte Indikatoren
// Matrix<double> IndVar = Matrix<double>(N,(z.cols()-p1-p3),1,-99); // Matrix für Varianzen, abhängig vom Komponentenindikator
//
                if (p2 > 0){                        // only for poisson distributed variables
                       for (j = p1; j < p12; ++j){        // für poissonverteilte Indikatoren
                       if (Omega_exists){                     // falls direkte Kovariaten vorhanden
                          A_trans = !A;
                       }
                       B_trans = !B;
                       if(Xdir_par_exists){            // nonparametrische direkte Effekte
                        nonpartemp = Matrix<double>(N, 1, 1, 0.0); 	
	                      for (t=0; t<sdir; ++t) {
						                nonpartemp = nonpartemp + Xdir_par[t] * !gammadir_par[t](j,_);
	                      }
                        if (Omega_exists){    // z_logmu_vec_pois: Vektor mit log(mu_ij), i =1,...,N
                            z_logmu_vec_pois = b0(j,0) + nonpartemp + Omega * A_trans(_,j) + eta * B_trans(_,j);
                        }
                        else {
                            z_logmu_vec_pois = b0(j,0) + nonpartemp + eta * B_trans(_,j);
                        }
                       }
                       else {  // falls keine nonparametrischen direkten Effekte
                        if (Omega_exists){    // z_logmu_vec_pois: Vektor mit log(mu_ij), i =1,...,N
                            z_logmu_vec_pois = b0(j,0) + Omega * A_trans(_,j) + eta * B_trans(_,j);
                        }
                        else {
                            z_logmu_vec_pois = b0(j,0) + eta * B_trans(_,j);
                        }
                       }
//Rprintf("vor dim: z_logmu_vec_pois_row = %d,   z_logmu_vec_pois_col  = %d\n",z_logmu_vec_pois.rows() , z_logmu_vec_pois.cols() ); 
//Rprintf("Matrix z_logmu_vec_pois = %f\n",z_logmu_vec_pois(10,0));
//Rprintf("Matrix offs = %f\n", #Offs(10,0));
                       if (offset) z_logmu_vec_pois = z_logmu_vec_pois + Offs;         // falls offset enthalten
//Rprintf("nach  dim: z_logmu_vec_pois_row = %d,   z_logmu_vec_pois_col  = %d\n",z_logmu_vec_pois.rows() , z_logmu_vec_pois.cols() ); 
//Rprintf("Matrix z_logmu_vec_pois = %f\n",z_logmu_vec_pois(10,0));
//return;
                       for (i=0; i<N; ++i){                // für alle Beobachtungen
                           vector<double> zwisch(int(Y(i,j))+1);   // (Y_ij+1)-dim. Vektor für Zwischenankunftszeiten
                           sum_zwisch = 0.0;            // summmierte Zwischenankunftszeiten
                           if (int(Y(i,j)) > 0){
                              vector<double> gleich(int(Y(i,j)));    // (Y_ij)-dim Vektor gleichverteilter Zufallszahlen
                              for (l = 0; l < int(Y(i,j)); ++l){
                                  gleich[l] = stream -> runif();   // ziehe auf (0,1)-gleichverteilte Zufallszahlen
                              }
                              sort(gleich.begin(),gleich.end());  // bildet Orderstatistik
                              zwisch[0] = gleich[0];              // erste Zwischenankunftszeit
                              //sum_zwisch = zwisch[0];             //alternative
                              for (l = 1; l < int(Y(i,j)); ++l){
                                  zwisch[l] = gleich[l] - gleich[l-1];   // berechnet alle Zwischenankunftszeiten
                                  //sum_zwisch = sum_zwisch + zwisch[l];   // summiert Zwischenankunftszeiten ??? anders lösbar?
                              }
                              sum_zwisch = gleich[int(Y(i,j))-1];       // bessere Alternative  ???
                           }  // Ende if (y(i,j) > 0)
                           z_mu_vec_pois(i,0) = exp(z_logmu_vec_pois(i,0));  // berechnet mu_ij
                           xi_exp = stream -> rexp(z_mu_vec_pois(i,0));      // zieht Exp(mu_ij)-verteilte Zufallszahl
                           zwisch[int(Y(i,j))] = 1 - sum_zwisch + xi_exp;         // abschließende Zwischenankunftszeit
                           for (l=0; l <= int(Y(i,j)); ++ l){
                               compsum=0.0;                   // summierte prop. diskrete Dichten
                               for (k = 0; k < 10; ++k){          // berechnet prop. diskrete Dichten
                                   tempprop =  (log(zwisch[l])  + z_logmu_vec_pois(i,0) - compmean(k,0))/ sqrtcompvar(k,0);
                                   logcompprop[k] = (-log(sqrtcompvar(k,0)) -0.5*(tempprop*tempprop)); // berechnet log f
                                   explogcompprop[k] = exp(logcompprop[k]); // berechnet aus log f durch exp die prop diskrete Dichte f
                                   compprop[k] = compweight(k,0) *  explogcompprop[k]; // berechnet f * weight
                                   compsum = compsum + compprop[k];   // summiert die k prop. diskreten Dichten
                               }
                               // Normierung der prop. diskreten Dichte
                               kumcomppropnorm[0] = 0.0;
                               for (k=0; k<9; ++k){
                                   comppropnorm[k] = compprop[k] / compsum;
                                   kumcomppropnorm[k+1] = kumcomppropnorm[k] + comppropnorm[k];
                               }
                               compgleich = stream -> runif();     // ziehe auf (0,1)-gleichverteilte Zahl
                               compind=99;				  // Komponentenindikator
                               // je nach gezogener Zahl, Zuweisung des Indikatorwertes
                               if(kumcomppropnorm[0] <= compgleich && compgleich < kumcomppropnorm[1]){
                                    compind = 0;
                               }
                               else if(compgleich < kumcomppropnorm[2]){
                                    compind = 1;
                               }
                               else if(compgleich < kumcomppropnorm[3]){
                                    compind = 2;
                               }
                               else if(compgleich < kumcomppropnorm[4]){
                                    compind = 3;
                               }
                               else if(compgleich < kumcomppropnorm[5]){
                                    compind = 4;
                               }
                               else if(compgleich < kumcomppropnorm[6]){
                                    compind = 5;
                               }
                               else if(compgleich < kumcomppropnorm[7]){
                                    compind = 6;
                               }
                               else if(compgleich < kumcomppropnorm[8]){
                                    compind = 7;
                               }
                               else if(compgleich < kumcomppropnorm[9]){
                                    compind = 8;
                               }
                               else if(compgleich <= 1){
                                    compind = 9;
                               }

                               // ziehe z_ijl
					                     if (j==0){ // Sonderbehandlung für j==0 nötig, da sonst in zspa(1,(j-1)) negative Spalte
						                    IndVar(i,l) = invcompvar(compind,0);
						                    z(i,l) = -log(zwisch[l]) + compmean(compind,0);
                               } // Ende j==0
					                     else {
                                IndVar(i,((zspa(1,(j-1)))+l-p1)) = invcompvar(compind,0);
                                z(i,((zspa(1,(j-1)))+l)) =  -log(zwisch[l]) + compmean(compind,0);
                               } // Ende else
                           } // Ende l-Schleife
                       } // Ende i-Schleife
                   } // Ende j-Schleife
                } // Ende if (p2 > 0)     end of step 1
                
////////////////////////////////////////////////////////////////////////////////
 //Rprintf("iter = %d   End of sampling step 1\n", iter);	
// 2. sample  eta_i | z_i, theta, x_i, omega_i
for (i=0; i<N; ++i){
    // individuellen Matrizen
    // tempA , tempB, tempb0, tempz, temptheta
    // Dimension individuell
    tempz = Matrix<double>(dim[i], 1, 1, 0.0);
    tempb0 = Matrix<double>(dim[i], 1, 1, 0.0);
    tempB = Matrix<double>(dim[i], m, 1, 0.0);
    temptransB = Matrix<double>(m,dim[i],1,0.0);
    tempA = Matrix<double>(dim[i], d, 1, 0.0);
    //temptheta = Matrix<double>(dim[i],1, 1, 0.0);
    invtemptheta = Matrix<double>(dim[i],1,1,0.0);
    tempnonpar2 = Matrix<double>(p,1,1,0.0);
    tempnonparind = Matrix<double>(dim[i],1,1,0.0);
    
    tempoffset = Matrix<double>(dim[i],1,1,0.0);      // nur für Poisson-verteilte Indikatoren Einträge ungleich 0
    
    rowcount = 0;
    if(Xdir_par_exists){
        for(t=0; t<sdir; ++t){
            tempnonpar2 = tempnonpar2 + gammadir_par[t]* !Xdir_par[t](i,_);
        }
    }
    for (j=0; j<p; ++j){
        if (j<p1){             // für ordinale Indikatoren
            tempz((rowcount),0) = z(i,j);
            tempb0((rowcount),0) = b0(j,0);
            for (k=0; k<m; ++k) {
                tempB((rowcount),k) = B(j,k);
            }
            if(Omega_exists){
              for (k=0; k<d; ++k) {
                tempA((rowcount),k) = A(j,k);
              }
            }
            invtemptheta((rowcount),0) = 1.0/theta(j,j);
            if(Xdir_par_exists){
            tempnonparind((rowcount),0) = tempnonpar2(j,0);
            }
            rowcount = (rowcount +1);
        } // Ende ordinal
        else if (j>=p1 && j<p12){            // für poissonverteilte Indikatoren
            for (l=0; l<(int(Y(i,j))+1); ++l){        // entsprechend dem Wert von Y_ij
                if (j == 0){   // Sonderbehandlung für j==0, da sonst negativer Spaltenindex in spa(1,(j-1))
                    tempz((rowcount),0) = z(i,l);
                    tempb0((rowcount),0) = b0(j,0);
                    for (k=0; k<m; ++k) {
                        tempB((rowcount),k) = B(j,k);
                    }
                    if(Omega_exists){
                      for (k=0; k<d; ++k) {
                        tempA((rowcount),k) = A(j,k);
                      }
                    }
                    tempoffset((rowcount),0) = Offs(i,0);
                    invtemptheta((rowcount),0) = IndVar(i,l);
                    if(Xdir_par_exists){
                        tempnonparind((rowcount),0) = tempnonpar2(j,0);
                    }
                    rowcount = (rowcount +1);
                }
                else {
                    tempz((rowcount),0) = z(i,((zspa(1,(j-1)))+l));
                    tempb0((rowcount),0) = b0(j,0);
                    for (k=0; k<m; ++k) {
                        tempB((rowcount),k) = B(j,k);
                    }
                    if(Omega_exists){
                      for (k=0; k<d; ++k) {
                        tempA((rowcount),k) = A(j,k);
                      }
                    }
                    tempoffset((rowcount),0)=Offs(i,0);
                    invtemptheta((rowcount),0) = IndVar(i,((zspa(1,(j-1)))+l-p1));
                    if(Xdir_par_exists){
                        tempnonparind((rowcount),0) = tempnonpar2(j,0);
                    }
                    rowcount = (rowcount +1);
                }
            }   // Ende for l
        } // Ende poisson
        else if (j>=p12 && j<p){        // metrische Indikatoren          
            tempz((rowcount),0) = z(i,(zspa(1,j)-1)); 
            tempb0((rowcount),0) = b0(j,0);
            for (k=0; k<m; ++k) {
                tempB((rowcount),k) = B(j,k);
            }
            if(Omega_exists){
              for (k=0; k<d; ++k) {
                tempA((rowcount),k) = A(j,k);
              }
            }
            invtemptheta((rowcount),0) = 1.0/theta(j,j);
            if(Xdir_par_exists){
              tempnonparind((rowcount),0) = tempnonpar2(j,0);
            }
            rowcount = (rowcount+1);
        }     // Ende metrisch
    } // Ende j

    // Hilfsmatrizen zur Berechnung von Eta_mu und eta_var
    // Matrix<double> eta_var = Matrix<double>(m,m,1,0.0);
    // Matrix<double> eta_var_chol = Matrix<double>(m,m,1,0.0);
   
    // Matrix<double> tempvarwood = Matrix<double>(m,m,1,0.0);
    tempvarwood2 = Matrix<double>(dim[i],dim[i],1,0.0); // Woodbury's identity
    tempinvthetaB = Matrix<double>(dim[i],m,1,0.0);
    temptransB = !tempB;
//    tempinvtheta = inv_diag(temptheta);
    tempinvthetaB = diag_mult_matrix(invtemptheta,tempB);
    tempvarwood = invpd(eye<double>(m) + temptransB * tempinvthetaB);
    tempvarwood2 = create_diag(invtemptheta) - tempinvthetaB * tempvarwood * matrix_mult_diag(temptransB,invtemptheta);
    temp3 = temptransB * tempvarwood2;
    eta_var = eye<double>(m) - temp3 * tempB;
    eta_var_chol = cholesky(eta_var);
    Eta_mu = Matrix<double>(m,1,1,0.0);
    //Matrix<double> Eta_mu = Matrix<double>(m,1,1,0.0);
    if(X_exists){
        Eta_mu = Gamma* !X(i,_);
    }
    if(X_par_exists){
        for(t=0; t<s; ++t){
            Eta_mu = Eta_mu + gamma_par[t]* !X_par[t](i,_);
        }
    }
    
// Version vor offset-Aufnahme:    
/*    if(Xdir_par_exists){        // falls nonparametrische direkte Effekte
      if(Omega_exists){
        Eta_mu = Eta_mu +  temp3 * (tempz - tempb0 - tempnonparind - tempA * !Omega(i,_) - tempB * Eta_mu);
      }
      else {
        Eta_mu = Eta_mu + temp3 * (tempz - tempb0 - tempnonparind - tempB * Eta_mu);
      }
    }
    else{        // falls keine nonparametrischen direkten Effekte
      if(Omega_exists){
        Eta_mu = Eta_mu +  temp3 * (tempz - tempb0 - tempA * !Omega(i,_) - tempB * Eta_mu);
      }
      else {
        Eta_mu = Eta_mu + temp3 * (tempz - tempb0 - tempB * Eta_mu);
      }
    }  */

    // Version mit offset:
    if(Xdir_par_exists){        // falls nonparametrische direkte Effekte
      if(Omega_exists){         // falls parametrische direkte Effekte
        if(offset){
          Eta_mu = Eta_mu +  temp3 * (tempz - tempb0 - tempoffset - tempnonparind - tempA * !Omega(i,_) - tempB * Eta_mu);
        }
        else{                    // kein offset
          Eta_mu = Eta_mu +  temp3 * (tempz - tempb0 - tempnonparind - tempA * !Omega(i,_) - tempB * Eta_mu);
        }
      }
      else {                  // keine paramatrischen direkten Effekte
        if(offset){
          Eta_mu = Eta_mu + temp3 * (tempz - tempb0 - tempoffset - tempnonparind - tempB * Eta_mu);
        }
        else{                 // kein offset
          Eta_mu = Eta_mu + temp3 * (tempz - tempb0 - tempnonparind - tempB * Eta_mu);
        }
      }
    }
    else{        // falls keine nonparametrischen direkten Effekte
      if(Omega_exists){
        if(offset){
          Eta_mu = Eta_mu +  temp3 * (tempz - tempb0 - tempoffset - tempA * !Omega(i,_) - tempB * Eta_mu);
        }
        else{      // kein offset
          Eta_mu = Eta_mu +  temp3 * (tempz - tempb0 - tempA * !Omega(i,_) - tempB * Eta_mu);
        }
      }
      else {          // keine parametrischen direkten Effekte
        if(offset){
          Eta_mu = Eta_mu + temp3 * (tempz - tempb0 - offset - tempB * Eta_mu);
        }
        else{            // kein offset
          Eta_mu = Eta_mu + temp3 * (tempz - tempb0 - tempB * Eta_mu);
        }
      }
    }
    
      eta_random = gaxpy(eta_var_chol, stream -> rnorm(m,1), Eta_mu);
      for (r=0; r<m; ++r){
        eta(i,r) = eta_random(r,0);
    }
} // Ende i-Schleife     // Ende Samplingstep 2
//Rprintf("End of sampling step 2 (eta):\n");
        
///////////////////////////////////////////////////////////////////////////////

// 3. sample  gamma | eta, x
// 3.1. Sampling of parametric covariates (gamma and X)
		if (X_exists) {
			for (r=0; r<m; ++r) {
				if (X_par_exists) {
					praed_rest = Matrix<double>(N, 1, 1, 0.0); 	
					for (t=0; t<s; ++t) {
						praed_rest = praed_rest + X_par[t] * !gamma_par[t](r,_);
					}
					gamma_mu_r = gamma_mu_multiplicator * (Gamma_stern * gamma_stern + !X_total * (eta(_, r) - praed_rest));
					//gamma_mu_r = gamma_mu_multiplicator * (eta(_, r) - praed_rest);
				}
				else gamma_mu_r = gamma_mu_multiplicator * (Gamma_stern * gamma_stern + !X_total * eta(_, r));
				// else gamma_mu_r = gamma_mu_multiplicator * eta(_, r);
				gamma_r = gaxpy(gamma_var_chol, stream->rnorm(q, 1), gamma_mu_r);	
        // gamma_r = gaxpy(gamma_var_chol, stream->rnorm(q+1, 1), gamma_mu_r);      //Änderung gamma_intercept
				
				// INTERCEPT-TEST
				//gamma_intercept = gamma_r(0,0); //Änderung gamma_intercept
				
				for (l=0; l<q; ++l)
					gamma(r*q+l,0) = gamma_r(l,0);
						//gamma(r*q+l,0) = gamma_r(l+1,0);  //Änderung gamma_intercept
			}
			
			// aus gamma wird Gamma gefüllt
			for (r=0; r<m; ++r) {
				for (l=0; l<q; ++l)
					Gamma(r, l) = gamma(r*q+l,0);
			}
		}
// Rprintf("iter = %d   End of sampling step 3.1\n", iter);
// 3.2. Sampling of nonparametric covariates (gamma_par and X_par and K_par)
		if (X_par_exists) {
			for (t=0; t<s; ++t) {
				for (r=0; r<m; ++r) {
					// Calculation of the rest of praedictor (for both old and new versions)
					if (s>1 || X_exists) {
						praed_rest = Matrix<double>(N, 1, 1, 0.0); 	
						if (s>1) {
							for (w=0; w<s; ++w)
								if (w != t) praed_rest = praed_rest + X_par[w] * !gamma_par[w](r,_);
						}

						if (X_exists)
							praed_rest = praed_rest + X * !Gamma(r,_); 
					}
//Rprintf("t = %d    r= % r   1. Grenze\n", t, r);
					//New version 
					if (s>1 || X_exists)
						solution_vec = X_par_trans[t] * (eta(_, r) - praed_rest);
					else 
						solution_vec = X_par_trans[t] * eta(_, r);			
						
					gamma_nonpar_chol = band_chol(band_addition(X_par_cross_band[t], K_par_band[t]/tau_par[t](r,0)));
					gamma_nonpar_rand = band_solve_backward(gamma_nonpar_chol, stream->rnorm(p_nr(t,0), 1));
					gamma_nonpar_mu_r = band_solve(gamma_nonpar_chol, solution_vec);
					gamma_nonpar_rand += gamma_nonpar_mu_r;
					
					// Datenübertragung und Mittelwertbildung
					average=0.0;
					for (l=0; l<p_nr(t,0); ++l) {
						gamma_par[t](r,l) = gamma_nonpar_rand(l,0);		// Datenübertragung (für alle nonpar. Effects gleich)
						average += gamma_par[t](r,l);
					}
				//Rprintf("Vorläufiger   average = %f\n", average/p_nr(t,0));
				// Ohne Zentrierung für P-Splines und Random walk schlechte Konfidenzintervalle
					if (var_cof[t] == 0) {
						switch(degr[t]) {
							case 1: average = average - 0.5 * gamma_par[t](r,0) - 0.5 * gamma_par[t](r,p_nr(t,0)-1);
											average /= p_nr(t,0) - 1;
											break;  // case 1
							case 2: average = average - (5.0/6.0) * gamma_par[t](r,0) - (1.0/6.0) * gamma_par[t](r,1) 
																				- (1.0/6.0) * gamma_par[t](r,p_nr(t,0)-2) - (5.0/6.0) * gamma_par[t](r,p_nr(t,0)-1);
											average /= p_nr(t,0) - 2;
											break;	// case 2
							case 3: average = average - (23.0/24.0) * gamma_par[t](r,0) - (12.0/24.0) * gamma_par[t](r,1)
																				- (1.0 /24.0) * gamma_par[t](r,2) - (1.0 /24.0) * gamma_par[t](r,p_nr(t,0)-3) 
																				- (12.0/24.0) * gamma_par[t](r,p_nr(t,0)-2) - (23.0/24.0) * gamma_par[t](r,p_nr(t,0)-1);
											average /= p_nr(t,0) - 3;
											break;	// case 3
							default: 		average /= p_nr(t,0);
						}

					//	Rprintf("Tatsächlicher average = %f\n", average);
						for (l=0; l<p_nr(t,0); ++l)
							gamma_par[t](r,l) -= average;
					}					
				}
			} 		
		}
//Rprintf("iter = %d   End of sampling step 3.2\n", iter);
		// 3.3. Sampling of smoothing parameters (tau_par)
		if (X_par_exists) {
			for (r=0; r<m; ++r) {
				for (t=0; t<s; ++t) {
					a_tau = a_par(t,0) + rank_K(t,0) / 2.0; 
					b_tau = b_par(t,0) + 0.5 * gamma_par[t](r,_) * K_par[t] * !gamma_par[t](r, _);
					
					tau_par[t](r,0) = stream->rigamma(a_tau, b_tau(0,0));
					//Rprintf("tau_par[t](r,0)  = %f    für t=%d  r=%d\n", tau_par[t](r,0), t, r);
				}
			}
		}		
		// End of sampling step 3 (gamma)
//Rprintf("iter = %d   End of sampling step 3.3\n", iter);  
////////////////////////////////////////////////////////////////////////////////

// 4.1 sample  beta | z, eta, theta, omega
	  U = cbind(U_basis, eta);
    U_cross = crossprod(U);
  	U_trans = !U;
  	for (j=0; j<p; ++j) { // für jeden Indikator ein Updateschritt
		if (j<p1){ // für ordinale Indikatoren
  			beta_stern_j = beta_stern(j*dimension_beta_j,0,(j+1)*dimension_beta_j-1,0);
  			Beta_stern_j = Beta_stern(j*dimension_beta_j,j*dimension_beta_j,(j+1)*dimension_beta_j-1,(j+1)*dimension_beta_j-1);
  			theta_inv_j = 1.0/theta(j,j);
			  beta_var_j = invpd(Beta_stern_j + theta_inv_j * U_cross);
			  if(Xdir_par_exists){
		      praeddir_rest = Matrix<double>(N, 1, 1, 0.0);
		      for (t=0; t<sdir; ++t) {
						    praeddir_rest = praeddir_rest + Xdir_par[t] * !gammadir_par[t](j,_);
		      }
          beta_mu_j = beta_var_j*(Beta_stern_j*beta_stern_j + theta_inv_j*U_trans*(z(_,j)-praeddir_rest));
        }
			  else{
			     beta_mu_j = beta_var_j*(Beta_stern_j*beta_stern_j + theta_inv_j*U_trans*z(_,j));
			  }
			  beta_var_j_chol = cholesky(beta_var_j);
		    beta_j = gaxpy(beta_var_j_chol, stream->rnorm(dimension_beta_j, 1), beta_mu_j);

        // aus beta_j wird beta, b0, B und A gefüllt
        b0(j,0) = beta_j(0,0);
        beta(j*dimension_beta_j,0) = b0(j,0);
        
        if (Omega_exists) {
				    for (i=0; i<d; ++i) {
					   A(j,i) = beta_j(1+i,0);
					   beta(j*dimension_beta_j+1+i,0) = A(j,i);
				    }
        }
        for (r=0; r<m; ++r) {
			    B(j,r) = lambdaConstr(j,r)*beta_j(1+d+r,0);
			    beta(j*dimension_beta_j+1+d+r,0) = B(j,r);
		    }
        //Rprintf("iter = %d ordinal = %d  End of sampling step 4 for ordinal \n", iter,j);

		} // Ende ordinale Indikatoren
		else if (j<p12){ //poissonverteilte Indikatoren
      // von j abhängige Matrizen tempzE, tempU, tempW, invtempW
			tempzE = Matrix<double>(dimpois[(j-p1)],1,1,0.0);
			tempU = Matrix<double>(dimpois[(j-p1)],(1+d+m),1,0.0);
			transtempU = Matrix<double>((1+d+m),dimpois[(j-p1)],1,0.0);
			//Matrix<double> tempW = Matrix<double>(dimpois[(j-p1)],1,1,0.0);
			invtempW = Matrix<double>(dimpois[(j-p1)],1,1,0.0);
			praeddir_restj = Matrix<double>(dimpois[(j-p1)],1,1,0.0);
			tempdiff =  Matrix<double>(dimpois[(j-p1)],1,1,0.0);
			//Matrix<double> tempUWU = Matrix<double>((1+d+m),(1+d+m),1,0.0);
			// Füllen der Matrizen
			rowcountE = 0;
			for (i=0; i<N; ++i){
				for (l=0; l<(int(Y(i,j))+1); ++l){
					if (j==0){ // falls p1==0 Sonderbehandlung da sonst in zspa(1,(j-i)) negativer Spaltenindex
						tempzE((rowcountE),0) = z(i,l);
						for (k=0; k<dimension_beta_j; ++k){
                tempU((rowcountE),k) = U(i,k);
            }
						invtempW((rowcountE),0) = IndVar(i,l);
						rowcountE = rowcountE + 1;
					}// Ende j==0
					else { // j>0
		        tempzE((rowcountE),0) = z(i,((zspa(1,(j-1)))+l));
            for (k=0; k<dimension_beta_j; ++k){
              tempU((rowcountE),k) = U(i,k);
            //Rprintf("Matrix U= %f i= %d k= %d \n", U(i,k),i,k);
            //Rprintf("Matrix tempU= %f (rowcountE+1)= %d k= %d \n", tempU((rowcountE+1),k),(rowcountE+1),k);
            }
						invtempW((rowcountE),0) = IndVar(i,((zspa(1,(j-1)))+l-p1));
            // Rprintf("Matrix tempW((rowcountE+1),(rowcountE+1))= %d\n",IndVar(i,((zspa(1,(j-1)))+l-p1)));
						rowcountE = rowcountE + 1;
					}
				} // Ende for l
			} // Ende for i
	    transtempU = !tempU;
			tempUWU = transtempU * diag_mult_matrix(invtempW,tempU);
      beta_stern_j = beta_stern(j*dimension_beta_j,0,(j+1)*dimension_beta_j-1,0);
      Beta_stern_j = Beta_stern(j*dimension_beta_j,j*dimension_beta_j,(j+1)*dimension_beta_j-1,(j+1)*dimension_beta_j-1);
     	beta_var_j = invpd(Beta_stern_j + tempUWU);
    	if(Xdir_par_exists){
     	    rowcountE2 = 0;
		      praeddir_rest = Matrix<double>(N, 1, 1, 0.0);
		      for (t=0; t<sdir; ++t) {
						    praeddir_rest = praeddir_rest + Xdir_par[t] * !gammadir_par[t](j,_);
		      }
          for (i=0; i<N; ++i){
            for(l=0; l<(int(Y(i,j))+1); ++l){
               praeddir_restj(rowcountE2,0) = praeddir_rest(i,0);
               rowcountE2 = rowcountE2 + 1;
             }
          }
          tempdiff = tempzE - praeddir_restj;
          beta_mu_j = beta_var_j * (Beta_stern_j * beta_stern_j + transtempU * diag_mult_matrix(invtempW,tempdiff));
      }
      else  {
        beta_mu_j = beta_var_j * (Beta_stern_j * beta_stern_j + transtempU * diag_mult_matrix(invtempW,tempzE));
     	}
      beta_var_j_chol = cholesky(beta_var_j);
			beta_j = gaxpy(beta_var_j_chol, stream->rnorm(dimension_beta_j, 1), beta_mu_j);
			// aus beta_j wird beta, b0, B und A gefüllt
			b0(j,0) = beta_j(0,0);
			beta(j*dimension_beta_j,0) = b0(j,0);

			if (Omega_exists) {
				for (i=0; i<d; ++i) {
					A(j,i) = beta_j(1+i,0);
					beta(j*dimension_beta_j+1+i,0) = A(j,i);
				}
			}
			for (r=0; r<m; ++r) {
					  B(j,r) = lambdaConstr(j,r)*beta_j(1+d+r,0);
					  beta(j*dimension_beta_j+1+d+r,0) = B(j,r);
			}
                //Rprintf("iter = %d poisson = %d  End of sampling step 4 for poisson \n", iter,j);

		} // Ende poissonverteilte Indikatoren
		else if (j<p){ // metrische Indikatoren
			beta_stern_j = beta_stern(j*dimension_beta_j,0,(j+1)*dimension_beta_j-1,0);
			Beta_stern_j = Beta_stern(j*dimension_beta_j,j*dimension_beta_j,(j+1)*dimension_beta_j-1,(j+1)*dimension_beta_j-1);
			theta_inv_j = 1.0/theta(j,j);
			beta_var_j = invpd(Beta_stern_j + theta_inv_j * U_cross);
			if(Xdir_par_exists){
		      praeddir_rest = Matrix<double>(N, 1, 1, 0.0);
		      for (t=0; t<sdir; ++t) {
						    praeddir_rest = praeddir_rest + Xdir_par[t] * !gammadir_par[t](j,_);
		      }
            beta_mu_j = beta_var_j*(Beta_stern_j*beta_stern_j + theta_inv_j*U_trans*(z(_,(zspa(1,j)-1))-praeddir_rest));

      }
      else{
          beta_mu_j = beta_var_j*(Beta_stern_j*beta_stern_j + theta_inv_j*U_trans*z(_,(zspa(1,j)-1)) );
      }
			beta_var_j_chol = cholesky(beta_var_j);
			beta_j = gaxpy(beta_var_j_chol, stream->rnorm(dimension_beta_j, 1), beta_mu_j);

			// aus beta_j wird beta, b0, B und A gefüllt
			b0(j,0) = beta_j(0,0);
			beta(j*dimension_beta_j,0) = b0(j,0);

			if (Omega_exists) {
				for (i=0; i<d; ++i) {
					A(j,i) = beta_j(1+i,0);
					beta(j*dimension_beta_j+1+i,0) = A(j,i);
				}
			}     

			for (r=0; r<m; ++r) {
					  B(j,r) = lambdaConstr(j,r)*beta_j(1+d+r,0);
					  beta(j*dimension_beta_j+1+d+r,0) = B(j,r);
			}
		} // Ende metrisch
  	} // End of sampling step 4.1 (beta)
//Rprintf("iter = %d   End of sampling step 4.1\n", iter);              
// 4.2. Sampling of nonparametric direct effects
		if (Xdir_par_exists) {
			for (t=0; t<sdir; ++t) {
				for (j=0; j<p; ++j) {
				if(j<p1){         // ordinale Indikatoren
					// Calculation of the rest of praedictor (for both old and new versions)
					if (sdir>1 || Omega_exists) {
						praeddir_rest = Matrix<double>(N, 1, 1, 0.0);
						praeddir_rest = b0(j,0) + eta * B_trans(_,j) + praeddir_rest;
						if (sdir>1) {
							for (w=0; w<sdir; ++w)
								if (w != t) praeddir_rest = praeddir_rest + Xdir_par[w] * !gammadir_par[w](j,_);
						}
/*Rprintf("praeddir_rest = %d   \n", praeddir_rest.cols());
Rprintf("praeddir_rest = %d   \n", praeddir_rest.rows());
Rprintf("Omega = %d   \n", Omega.cols());
Rprintf("Omega = %d   \n", Omega.rows());
Rprintf("A_trans(j,_) = %d   \n", A_trans(j,_).cols());
Rprintf("A_trans(j,_) = %d   \n", A_trans(j,_).rows()); */

            if (Omega_exists) //A_trans = !A;
							praeddir_rest = praeddir_rest + Omega * !A(j,_);
					}
        
//Rprintf("t = %d    r= % r   1. Grenze\n", t, r);
					//New version
/*Rprintf("Xdir_par_trans[t] = %d   \n", Xdir_par_trans[t].cols());
Rprintf("Xdir_par_trans[t] = %d   \n", Xdir_par_trans[t].rows());
Rprintf("z(_, j) = %d   \n", z(_, j).cols());
Rprintf("z(_, j) = %d   \n", z(_, j).rows());
Rprintf("praeddir_rest = %d   \n", praeddir_rest.cols());
Rprintf("praeddir_rest = %d   \n", praeddir_rest.rows());        */

					if (sdir>1 || Omega_exists)
						solutiondir_vec = Xdir_par_trans[t] * (z(_, j) - praeddir_rest);
					else
						solutiondir_vec = Xdir_par_trans[t] * z(_, j);
					
/*Rprintf("Xdir_par_cross_band[t].cols = %d   \n", Xdir_par_cross_band[t].cols());
Rprintf("Xdir_par_cross_band[t].rows = %d   \n", Xdir_par_cross_band[t].rows());
Rprintf("Kdir_par[0]cols = %d   \n", Kdir_par[0].cols());
Rprintf("Kdir_par[0] = %d   \n", Kdir_par[0].rows());
Rprintf("Kdir_par_band[0]cols = %d   \n", Kdir_par_band[0].cols());
Rprintf("Kdir_par_band[0] = %d   \n", Kdir_par_band[0].rows());
Rprintf("taudir_par[t]cols = %d   \n", taudir_par[t].cols());
Rprintf("taudir_par[t]] = %d   \n", taudir_par[t].rows());
Rprintf("Kdir_par_band[t].cols = %d   \n", Kdir_par_band[t].cols());
Rprintf("Kdir_par_band[t].rows = %d   \n", Kdir_par_band[t].rows());  */
					gammadir_nonpar_chol =  band_chol(band_addition(Xdir_par_cross_band[t], Kdir_par_band[t]/taudir_par[t](j,0)));
/*Rprintf("sdir = %d   \n", sdir);
Rprintf("Xdir_col = %d   \n", Xdir_par_trans[t].cols());
Rprintf("Xdir_rows = %d   \n", Xdir_par_trans[t].rows());
Rprintf("z.cols = %d   \n", z(_, j).cols());
Rprintf("z.rows = %d   \n", z(_, j).rows());
Rprintf("j= %d test:\n", j);  */

					gammadir_nonpar_rand = band_solve_backward(gammadir_nonpar_chol, stream->rnorm(pdir_nr(t,0), 1));
					gammadir_nonpar_mu_j = band_solve(gammadir_nonpar_chol, solutiondir_vec);
					gammadir_nonpar_rand += gammadir_nonpar_mu_j;
					//Rprintf("Ende ordinal nonparametric \n");   
         
					}
           else if (j<p12){ // Poisson-verteile Indikatoren
           tempzE2 = Matrix<double>(dimpois[(j-p1)],1,1,0.0);
			     invtempW2 = Matrix<double>(dimpois[(j-p1)],1,1,0.0);

                rowcountE = 0;
                
           			for (i=0; i<N; ++i){
				          for (l=0; l<(int(Y(i,j))+1); ++l){
					         if (j==0){ // falls p1==0 Sonderbehandlung da sonst in zspa(1,(j-i)) negativer Spaltenindex
						          tempzE2((rowcountE),0) = z(i,l);
						          for (k=0; k<eta.cols(); ++k){
                        tempeta((rowcountE),k) = eta(i,k);
                      }
                      if(Omega_exists){
                        for (l=0; l<Omega.cols(); ++l){
                          tempOmega(rowcountE,l) = Omega(i,l);
                        }
                      }
                      invtempW2((rowcountE),0) = IndVar(i,l);
						          rowcountE = rowcountE + 1;
			             }// Ende j==0
					         else { // j>0
		                  tempzE2((rowcountE),0) = z(i,((zspa(1,(j-1)))+l));
                      for (k=0; k<eta.cols(); ++k){
                        tempeta((rowcountE),k) = eta(i,k);
            //Rprintf("Matrix U= %f i= %d k= %d \n", U(i,k),i,k);
            //Rprintf("Matrix tempU= %f (rowcountE+1)= %d k= %d \n", tempU((rowcountE+1),k),(rowcountE+1),k);
                      }
                      if(Omega_exists){
                        for (l=0; l<Omega.cols(); ++l){
                          tempOmega(rowcountE,l) = Omega(i,l);
                        }
                      }
						          invtempW2((rowcountE),0) = IndVar(i,((zspa(1,(j-1)))+l-p1));
            // Rprintf("Matrix tempW((rowcountE+1),(rowcountE+1))= %d\n",IndVar(i,((zspa(1,(j-1)))+l-p1)));
						          rowcountE = rowcountE + 1;
					         }
				          } // Ende for l
	             } // Ende for i
					  if (sdir>1 || Omega_exists) {
						  praeddir_rest = Matrix<double>(tempzE.rows(), 1, 1, 0.0);
					   	praeddir_rest =  b0(j,0) + tempeta * B_trans(_,j) + praeddir_rest;
						if (sdir>1) {
							for (w=0; w<sdir; ++w)
								if (w != t) praeddir_rest = praeddir_rest + Xdir_parhelp[((j-p1)*sdir+w)] * !gammadir_par[w](j,_);
						}

            if (Omega_exists) 
							praeddir_rest = praeddir_rest + tempOmega * !A(j,_);
					}
//Rprintf("t = %d    r= % r   1. Grenze\n", t, r);
					//New version
					if (sdir>1 || Omega_exists)
						solutiondir_vec = !Xdir_parhelptrans[((j-p1)*sdir+t)] * diag_mult_matrix(invtempW2,(tempzE2(_, j) - praeddir_rest));
					else  {
/*Rprintf("invtempW_col = %d   \n", invtempW2.cols());
Rprintf("invtempW = %d   \n", invtempW2.rows());
Rprintf("tempzE(_, j).cols = %d   \n", tempzE2(_, j).cols());
Rprintf("tempzE(_, j).rows = %d   \n", tempzE2(_, j).rows());
Matrix<double> testmult = Matrix<double>(tempzE2.rows(),1,1,0.0);
testmult = diag_mult_matrix(invtempW2,tempzE2(_, j));
Rprintf("testmult.cols = %d   \n", testmult.cols());
Rprintf("testmult.rows = %d   \n",testmult.rows());
Rprintf("j = %d   \n",j);
Rprintf("p1 = %d   \n",p1);
Rprintf("sdir = %d   \n",sdir);
Rprintf("t = %d   \n",t);

Rprintf("((j-p1)*sdir+t) = %d   \n",((j-p1)*sdir+t));
Rprintf("Xdir_parhelp[((j-p1)*sdir+t)].cols = %d   \n", Xdir_parhelp[((j-p1)*sdir+t)].cols());
Rprintf("Xdir_parhelp[((j-p1)*sdir+t)].rows = %d   \n",Xdir_parhelp[((j-p1)*sdir+t)].rows());
Rprintf("Xdir_parhelptrans[((j-p1)*sdir+t)].cols = %d   \n", Xdir_parhelptrans[((j-p1)*sdir+t)].cols());
Rprintf("Xdir_parhelptrans[((j-p1)*sdir+t)].rows = %d   \n", Xdir_parhelptrans[((j-p1)*sdir+t)].rows());

Rprintf("Xdir_parhelp[0].cols = %d   \n", Xdir_parhelp[0].cols());
Rprintf("Xdir_parhelp[0].rows = %d   \n",Xdir_parhelp[0].rows());

Rprintf("Xdir_parhelp[1].cols = %d   \n", Xdir_parhelp[1].cols());
Rprintf("Xdir_parhelp[1].rows = %d   \n",Xdir_parhelp[1].rows()); */
            
						solutiondir_vec =  Xdir_parhelptrans[((j-p1)*sdir+t)] * diag_mult_matrix(invtempW2,tempzE2(_, j));
            
            }
          gammadir_nonpar_chol_help = Xdir_parhelptrans[((j-p1)*sdir+t)] *  diag_mult_matrix(invtempW2,Xdir_parhelp[((j-p1)*sdir+t)]);
          gammadir_nonpar_chol_helpband = band_creation(gammadir_nonpar_chol_help);

					gammadir_nonpar_chol = band_chol(band_addition(gammadir_nonpar_chol_helpband, Kdir_par_band[t]/taudir_par[t](j,0)));
					
					gammadir_nonpar_rand = band_solve_backward(gammadir_nonpar_chol, stream->rnorm(pdir_nr(t,0), 1));
					
					gammadir_nonpar_mu_j = band_solve(gammadir_nonpar_chol, solutiondir_vec);
					gammadir_nonpar_rand += gammadir_nonpar_mu_j;
					//Rprintf("Ende poisson nonparametric \n");
          
          }
          else {    //metrische Indikatoren
          theta_inv_j = 1.0/theta(j,j);
          					// Calculation of the rest of praedictor (for both old and new versions)
					if (sdir>1 || Omega_exists) {
						praeddir_rest = Matrix<double>(N, 1, 1, 0.0);
						praeddir_rest =  b0(j,0) + eta * B_trans(_,j) + praeddir_rest;
						if (sdir>1) {
							for (w=0; w<sdir; ++w)
								if (w != t) praeddir_rest = praeddir_rest + Xdir_par[w] * !gammadir_par[w](j,_);
						}

            if (Omega_exists) 
							praeddir_rest = praeddir_rest + Omega * !A(j,_);
					}
//Rprintf("t = %d    r= % r   1. Grenze\n", t, r);
					//New version
					if (sdir>1 || Omega_exists)
					{
						  solutiondir_vec = theta_inv_j * Xdir_par_trans[t] * (z(_,(zspa(1,j)-1)) - praeddir_rest);
					}
          else{
						  solutiondir_vec = theta_inv_j * Xdir_par_trans[t] * z(_,(zspa(1,j)-1));
          }  
					gammadir_nonpar_chol = band_chol(band_addition((Xdir_par_cross_band[t])/theta(j,j), Kdir_par_band[t]/taudir_par[t](j,0)));
					gammadir_nonpar_rand = band_solve_backward(gammadir_nonpar_chol, stream->rnorm(pdir_nr(t,0), 1));
					gammadir_nonpar_mu_j = band_solve(gammadir_nonpar_chol, solutiondir_vec);
					gammadir_nonpar_rand += gammadir_nonpar_mu_j;
					}

					// Datenübertragung und Mittelwertbildung
					average=0.0;
					for (l=0; l<pdir_nr(t,0); ++l) {
						gammadir_par[t](j,l) = gammadir_nonpar_rand(l,0);		// Datenübertragung (für alle nonpar. Effects gleich)
						average += gammadir_par[t](j,l);
					}


				//Rprintf("Vorläufiger   average = %f\n", average/p_nr(t,0));
				// Ohne Zentrierung für P-Splines und Random walk schlechte Konfidenzintervalle
					if (var_cofdir[t] == 0) {
						switch(degrdir[t]) {
							case 1: average = average - 0.5 * gammadir_par[t](j,0) - 0.5 * gammadir_par[t](j,pdir_nr(t,0)-1);
											average /= pdir_nr(t,0) - 1;
											break;  // case 1
							case 2: average = average - (5.0/6.0) * gammadir_par[t](j,0) - (1.0/6.0) * gamma_par[t](j,1)
																				- (1.0/6.0) * gammadir_par[t](j,pdir_nr(t,0)-2) - (5.0/6.0) * gammadir_par[t](j,pdir_nr(t,0)-1);
											average /= pdir_nr(t,0) - 2;
											break;	// case 2
							case 3: average = average - (23.0/24.0) * gammadir_par[t](j,0) - (12.0/24.0) * gammadir_par[t](j,1)
																				- (1.0 /24.0) * gammadir_par[t](j,2) - (1.0 /24.0) * gammadir_par[t](j,pdir_nr(t,0)-3)
																				- (12.0/24.0) * gammadir_par[t](j,pdir_nr(t,0)-2) - (23.0/24.0) * gammadir_par[t](j,pdir_nr(t,0)-1);
											average /= pdir_nr(t,0) - 3;
											break;	// case 3
							default: 		average /= pdir_nr(t,0);
						}

					//	Rprintf("Tatsächlicher average = %f\n", average);
						for (l=0; l<pdir_nr(t,0); ++l)
							gammadir_par[t](j,l) -= average;
					}
					}
					
/*for (int zeile=0; zeile< gammadir_par[0].rows(); ++zeile){
Rprintf("zeile = %d", zeile);
  for (int spalte=0; spalte< gammadir_par[0].cols(); ++spalte){
    Rprintf(" %f",   gammadir_par[0](zeile,spalte)) ;
    }
    Rprintf(" \n");
}  */

					
				}
		}
          
		// 4.3. Sampling of smoothing parameters (taudir_par)
		if (Xdir_par_exists) {
			for (j=0; j<p; ++j) {
				for (t=0; t<sdir; ++t) {
					adir_tau = adir_par(t,0) + rank_Kdir(t,0) / 2.0;
					bdir_tau = bdir_par(t,0) + 0.5 * gammadir_par[t](j,_) * Kdir_par[t] * !gammadir_par[t](j,_);
         
					taudir_par[t](j,0) = stream->rigamma(adir_tau, bdir_tau(0,0));
				//	Rprintf("taudir_par[t](j,0)  = %f    für t=%d  j=%d\n", taudir_par[t](j,0), t, j);
				}
			}
		}
		
/*     Rprintf("tau \n");
    for (int zeile=0; zeile<taudir_par[0].rows(); ++zeile) {
					for (int spalte=0; spalte<taudir_par[0].cols(); ++spalte) {
						Rprintf("%f ", taudir_par[0](zeile,spalte));    }
					Rprintf("\n");    
		 }*/
		// End of sampling step 4 (gammadir)
//Rprintf("iter = %d   End of sampling step 4.3\n", iter);

  
////////////////////////////////////////////////////////////////////////////////   
    	
//Rprintf("iter = %d   End of sampling step 4\n", iter);
			
////////////////////////////////////////////////////////////////////////////////

    // 5. sample  theta_j | z_j, eta, theta \ {theta_j}, omega
		if (p3>0) { 		// sampling only for metric indicators
			//double dof;
			//double s_inf;
			for (j=p12; j<p; ++j) {
				dof = N + v_stern(j-p12,0);
				s_inf = 0.0;
		    if(Xdir_par_exists){
          nonpartemp = Matrix<double>(N, 1, 1, 0.0); 	
					for (t=0; t<sdir; ++t) {
						nonpartemp = nonpartemp + Xdir_par[t] * !gammadir_par[t](j,_);
					}
				}
				for (i=0; i<N; ++i)
				{
					if (Omega_exists)
						temp4 = A(j,_) * (!Omega(i,_));
					temp5 = B(j,_) * (!eta(i,_));
					if (Xdir_par_exists) {
            temp6 = z(i,(zspa(1,j)-1)) - b0(j,0) - nonpartemp(i,0) - temp4(0,0) - temp5(0,0);
          }
          else {
					temp6 = z(i,(zspa(1,j)-1)) - b0(j,0) - temp4(0,0) - temp5(0,0);
					}
					s_inf += temp6*temp6;					
				}
				s_inf += v_stern(j-p12,0)*s_stern(j-p12,0)*s_stern(j-p12,0);
				s_inf = s_inf / dof;
				theta(j,j) = stream->rchisq(dof);
				theta(j,j) = dof * s_inf / theta(j,j);
			}
		}       // End of sampling step 5 (theta)
       //Rprintf("iter = %d   End of sampling step 5\n", iter);
		
////////////////////////////////////////////////////////////////////////////////
		
		// 6. sample  nu_jk | y_j, z_j, theta \ {nu_jk}
		//Matrix<double> proposal;
		//double prob_likelihood;
		//double prob_norm;
		
		if (p1 > 0) {
			if (*mh) {
				for (j=0; j<p1; ++j) {
					if (ncat(j,0) > 2) {
						proposal = nu(j, 0, j, ncat(j,0));			// proposal wird aus nu entnommen		
						for (k=2; k<ncat(j,0); ++k)					// proposal berechnen
							proposal(0, k) = stream->rtnorm_combo(proposal(0,k), tune(j,0)*tune(j,0), proposal(0, k-1), proposal(0, k+1));	
						prob_likelihood = 0.0;
						prob_norm 			= 0.0;
						if(Xdir_par_exists){
                for (t=0; t<sdir; ++t) {
						      nonpartempmh = nonpartempmh + Xdir_par[t] * !gammadir_par[t](j,_);
                }
            }

						for (i=0; i<N; ++i) {								// Calculation of the probability ratio of likelihood
							B_row = B(j,_);
              if(Xdir_par_exists){
							  if (Omega_exists) {						
								  A_row = A(j,_);
								  temp1 = A_row * (!Omega(i,_));
								  temp2 = B_row * (!eta(i,_));
								  z_mu = b0(j,0) + temp1(0,0) + temp2(0,0) + nonpartempmh(i,0);					
					      } else
						    {
								  temp2 = B_row * (!eta(i,_));
								  z_mu = b0(j,0) + temp2(0,0) + nonpartempmh(i,0);
							  }
              }
              else{     // falls keine nonparametrischen direkten Effekte
							if (Omega_exists) {						
								A_row = A(j,_);
								temp1 = A_row * (!Omega(i,_));
								temp2 = B_row * (!eta(i,_));
								z_mu = b0(j,0) + temp1(0,0) + temp2(0,0);					
							} else
							{
								temp2 = B_row * (!eta(i,_));
								z_mu = b0(j,0) + temp2(0,0);
							}
							}
							prob_likelihood = prob_likelihood
								+ log( pnorm(proposal(0, static_cast<int>(Y(i,j)  )) - z_mu)
									    -pnorm(proposal(0, static_cast<int>(Y(i,j)-1)) - z_mu)) 
								- log( pnorm(      nu(j, static_cast<int>(Y(i,j)  )) - z_mu)
									    -pnorm(      nu(j, static_cast<int>(Y(i,j)-1)) - z_mu));
						}	
						for (k=2; k<ncat(j,0); ++k)	{				// Calculation of normalization correction
							prob_norm = prob_norm
								+ log( pnorm(      nu(j,k+1), 			nu(j,k), tune(j,0))  
								      -pnorm(proposal(0,k-1), 			nu(j,k), tune(j,0)))	// ???proposal oder nu (1.)
										//	-pnorm(			nu(j,k-1), 			nu(j,k), tune(j,0)))
								- log( pnorm(proposal(0,k+1),proposal(0,k), tune(j,0)) 
								    	-pnorm(			 nu(j,k-1),proposal(0,k), tune(j,0)));	// ???nu oder proposal (1.)
								    //  -pnorm(proposal(0,k-1),proposal(0,k), tune(j,0)));	// ???nu oder proposal (1.)
						}
						if (stream->runif() <= exp(prob_likelihood + prob_norm)) {
							for (k=2; k < ncat(j,0); ++k) {
								nu(j,k)	= proposal(0,k);								
							}
							++accepts[0];
			}	}	}	}				
			else {																			// wenn kein MH-Algorithmus gewählt
				for (j=0; j<p1; ++j) {
					if (ncat(j,0) > 2) {										// for binary indicators this loop is not entered
						for (k=2; k<ncat(j,0); ++k) 
						{	
							LO  = -10000000;
							UP =  10000000;
							for (i=0; i<N; ++i) {
								if(static_cast<int>(Y(i,j)) == k)
									LO = SCYTHE::max(LO, z(i,j));
								else if(static_cast<int>(Y(i,j)) == (k+1))
									UP = SCYTHE::min(UP, z(i,j));
							}
							LO = SCYTHE::max(LO, nu(j,k-1));
							UP = SCYTHE::min(UP, nu(j,k+1));
							//if (LO > UP) Rprintf("LO > UP!\n");
							nu(j,k) = stream->runif()*(UP-LO) + LO ;
						} 
					}	
				}	
			}
		} // End of sampling step 6 (nu)	
//Rprintf("iter = %d   End of sampling step 6\n", iter);
   			
////////////////////////////////////////////////////////////////////////////////

		// Durchführung des Grouped-Move-Steps wenn gewünscht
//double testtest=0.0;
		if (*gm && p1>0) {
			if (Omega_exists) A_trans = !A;
			   B_trans = !B;
    	for (j=0; j<p1; ++j) {
				if (Omega_exists)
					diff_temp = z(_,j) - b0(j,0) - Omega * A_trans(_,j) - eta * B_trans(_,j);
				else 
					diff_temp = z(_,j) - b0(j,0) 												- eta * B_trans(_,j); 				
				diff_temp = crossprod(diff_temp);
      	a = (N+d+m+ncat(j,0)) / 2.0 ;			// Hier stand mal a = (N+d+m+ncat(j,0)+1) / 2.0 ;
				b = diff_temp(0,0) / 2.0;
				//Rprintf("Berechnetes a = %f für Var j=%d\n", a, j);
				//Rprintf("Berechnetes b = %f für Var j=%d\n", b, j);
				g = sqrt(stream->rgamma(a, b));
				//Rprintf("Zufallsvariable g = %f \n\n", g);
//if (j==0) testtest = g;
				//Update all relevant parameters

				// z_ij updaten
				for (i=0; i<N; ++i)
					z(i,j) *= g;
				// beta, b0, A und B updaten
				b0(j,0) *= g;
				beta(j*dimension_beta_j,0) = b0(j,0);
				if (Omega_exists)
					for (i=0; i<d; ++i) {
						A(j,i) *= g;
						beta(j*dimension_beta_j+1+i,0) = A(j,i);
					}
				for (r=0; r<m; ++r) {
					B(j,r) *= g;
					beta(j*dimension_beta_j+1+d+r,0) = B(j,r);
				}
				// nu updaten
				for (k=2; k<ncat(j,0); ++k)
					nu(j,k) *= g;

			} // End loop 0<=j<p1
		} // End grouped move step
		// print results to screen
    if (*verbose == 1 && iter % 250 == 0) {
      Rprintf("MCMCm iteration %i of %i \n\n", iter, tot_iter);
			if (X_exists) {
			 Rprintf("Gamma regressormatrix:\n");
      for (i=0; i<Gamma.rows(); ++i) {
				for (j=0; j<Gamma.cols(); ++j)
					Rprintf("%f  ", Gamma(i,j));
		    Rprintf("\n");}Rprintf("\n");}
			else Rprintf("No covariates X, therefore no regressormatrix Gamma.\n\n");
    	Rprintf("b0 intercept:\n");
      if (b0.cols() == 1 ) {
      	for (i=0; i<b0.rows(); ++i)
      		Rprintf("%f  ", b0(i,0));
      } else {
      	for (i=0; i<b0.rows(); ++i) {
					for (j=0; j<b0.cols(); ++j)
						Rprintf("%f  ", b0(i,j));
					Rprintf("\n");
				}
			}
			Rprintf("\n\n");
			Rprintf("B factor loadings:\n");
      if (B.cols() == 1) {
    		for (i=0; i<B.rows(); ++i) {
      		Rprintf("%15.5f ", B(i,0));
      	}
    	} else {
    		for (i=0; i<B.rows(); ++i) {
					for (j=0; j<B.cols(); ++j)
						Rprintf("%f ", B(i,j));
					Rprintf("\n");
				}	
			}
			
			
			/*

			Rprintf("\n\nDIC cov-matrix:\n");
      	for (int i=0; i<DIC_cov.rows(); ++i) {
					for (int j=0; j<DIC_cov.cols(); ++j)
						Rprintf("%f  ", DIC_cov(i,j));
				Rprintf("\n");}Rprintf("\n");

			Rprintf("\nDeterminante : %f\n\n", DIC_det);

						Rprintf("\n\nDIC mu2:\n");
      	for (int i=0; i<DIC_mu2.rows(); ++i) {
					for (int j=0; j<DIC_mu2.cols(); ++j)
						Rprintf("%f  ", DIC_mu2(i,j));
				Rprintf("\n");}Rprintf("\n");
			
			*/
			
			Rprintf("\n\n");
			if (Omega_exists) {
				Rprintf("A regressormatrix fixed effects:\n");
      	for (i=0; i<A.rows(); ++i) {
					for (j=0; j<A.cols(); ++j)
						Rprintf("%f  ", A(i,j));
				Rprintf("\n");}Rprintf("\n");}
			else Rprintf("No fixed effects Omega, therefore no regressormatrix A.\n\n");

      Rprintf("diag(theta) variances of e_ij:\n");
      for (i=0; i<theta.rows(); ++i) 
					Rprintf("%12.7f  ", theta(i,i));
			Rprintf("\n");Rprintf("\n");

      if (p1>0) {Rprintf("nu Schwellenwerte of z_ij:\n");
      for (i=0; i<nu.rows(); ++i) {
				for (j=0; j<nu.cols(); ++j)
					Rprintf("%f  ", nu(i,j));
			Rprintf("\n");}Rprintf("\n");}
			
			

			if (*mh)
				Rprintf("\nMetropolis-Hastings acceptance rate = %10.5f\n",
	      	static_cast<double>(*accepts)/(static_cast<double>(iter * p1)));

			Rprintf("MCMCm iteration %i of %i.   N = %d ,   MH-step = %d ,   GM-step = %d ", iter, tot_iter, N, *mh, *gm);
			if (sim_nr > 0) Rprintf("  , Sim-Nr. = %d .\n\n", sim_nr);
			  else Rprintf(" .\n\n");
    } // end print results to screen
  

		// store results
    if ((iter >= *burnin) && ((iter % *thin==0))) {

			// calculate DIC for this iteration
			if (DIC_exists) {
				DIC = 0.0;
				
				// eta (Prediktor der structural equation) berechnen   	 // ZEILE FÜR DIC2 ohne latent Variables
				Matrix<double> DIC_eta = Matrix<double>(N, m, 1, 0.0);	 // ZEILE FÜR DIC2 ohne latent Variables
				if (X_exists) DIC_eta = X * !Gamma;    									 // ZEILE FÜR DIC2 ohne latent Variables
				if (X_par_exists) {																			 // ZEILE FÜR DIC2 ohne latent Variables
					for (t=0; t<s; ++t)																 // ZEILE FÜR DIC2 ohne latent Variables
						DIC_eta += X_par[t] * !gamma_par[t];									 // ZEILE FÜR DIC2 ohne latent Variables
				}																												 // ZEILE FÜR DIC2 ohne latent Variables


				// Mittelwerte berechnen
				for (i=0; i<N; ++i) {
					if (Omega_exists)
						DIC_mu1 = b0 +  A * Omega(i,_) + B * (!DIC_eta(i,_));  // ZEILE FÜR DIC2 ohne latent Variables
						//DIC_mu1 = b0 +  A * Omega(i,_) + B * (!eta(i,_));  // ZEILE FÜR DIC1 mit latent Variables
					else
						DIC_mu1 = b0  								 + B * (!DIC_eta(i,_));	 // ZEILE FÜR DIC2 ohne latent Variables
						//DIC_mu1 = b0  								 + B * (!eta(i,_));	 // ZEILE FÜR DIC1 mit latent Variables
					
					for (j=0; j<p1; ++j) {   // Likelihood-Beitrag von ordinalen Indikatoren berechnen
						int obs_cat = static_cast<int>(Y(i,j));
						double prop = pnorm(nu(j,obs_cat)-DIC_mu1(j,0)) - pnorm(nu(j,obs_cat-1)-DIC_mu1(j,0));
						DIC = DIC + log(prop);
					}
					/*
					if (p1 < p) {
						for (int j=p1; j<p; ++j) { // Likelihood-Beitrag von metrischen Indikatoren berechnen

						}
					}
					*/
				}
				DIC = -2.0 * DIC;

				// Mittelwerte für die Berechnung von D(average(theta)) speichern
				//DIC_ysternmedium = DIC_ysternmedium + z;
				DIC_zmedium			 = DIC_zmedium + eta;  // ZEILE FÜR DIC1 mit latent Variables
				DIC_etamedium		 = DIC_etamedium + DIC_eta;  // ZEILE FÜR DIC2 ohne latent Variables
				DIC_b0medium	   = DIC_b0medium + b0;
				if (Omega_exists)
					DIC_Amedium		 = DIC_Amedium + A;
				DIC_Bmedium	  	 = DIC_Bmedium + B;
				//DIC_thetamedium	 = DIC_thetamedium + theta;
				DIC_numedium     = DIC_numedium + nu;

			} // end calculation of DIC


 			// store eta
 			if (*storescores==1)
 				for (r=0; r<m; ++r)
 					for (i=0; i<N; ++i) {
 						eta_store(count, r*N + i) = eta(i, r);
 					 	//Rprintf("eta_store i - r - value = %d - %d - %f\n", i, r, eta(i,r)); 
 					}
 		  // store beta
      for (i=0; i<p; ++i){	// Zuerst b0, dann b_1 speichern
      	for (j =0; j<(1+d+m); ++j){
      		if (j == 0)
						beta_store(count,i) = beta(i*(1+d+m)+j,0);
					else
						beta_store(count,p+i*(d+m)+j-1 ) = beta(i*(1+d+m)+j,0);
				}
			}
			// TEST SPEICHERE GAMMA-INTERCEPT
			if (X_exists || X_par_exists)
				gamma_intercept_store(count, 0) = gamma_intercept;

      // store gamma
      if (X_exists)
				for (i=0; i<m*q; ++i)
					gamma_store(count,i) = gamma(i,0);
			// store gamma_par
      if (X_par_exists) {
      	int temp = 0;
				for (r=0; r<m; ++r) // Reihenfolge: zuerst alle param für m=1 - t=0,...,s, dann m=2,......
					for (t=0; t<s; ++t)
						for (int prm=0; prm<p_nr(t,0); ++prm)
							{
								gamma_par_store(count, temp) = gamma_par[t](r,prm);
								++temp;
							}

				for (r=0; r<m; ++r) 	// Reihenfolge in tau_par_store: zuerst alle tau für m=1, dann m=2, ...
					for (t=0; t<s; ++t)
						//tau_par_store(count, r*m + t) = tau_par[t](r,0);   //Variante Alexander
						 tau_par_store(count, r*s + t) = tau_par[t](r,0);
			}
			
			// store gammadir_par
      if (Xdir_par_exists) {
      	int temp = 0;
				for (j=0; j<p; ++j) // Reihenfolge: zuerst alle param für j=1 - t=0,...,sdir, dann j=2,......
					for (t=0; t<sdir; ++t)
						for (int prm=0; prm<pdir_nr(t,0); ++prm)
							{
								gammadir_par_store(count, temp) = gammadir_par[t](j,prm);
								
								++temp;
							}


				for (j=0; j<p; ++j) {	// Reihenfolge in taudir_par_store: zuerst alle tau für p=1, dann p=2, ...
					for (t=0; t<sdir; ++t){
						taudir_par_store(count, j*sdir + t) = taudir_par[t](j,0);
					}
				}
			}

			// store theta
			for (j=0; j<p; ++j)   // ohne ordinale: int i=(p1+1); i<p
				{
					theta_store(count, j) = theta(j, j);
				}
			// store nu
  		if (p1 > 0 && nu.cols() > 3) {
				nr_nu = 0;
				for (i=0; i<p1; ++i)				// was p before
  				if (ncat(i,0) > 2) {
  						for (j=0; j<ncat(i,0)-2; ++j)
  							nu_store(count, nr_nu+j) = nu(i, j+2);
  						nr_nu += ncat(i,0) - 2;
  				}
     	}
     	// store z
 			if ((*storez==1 && p1>0) || (*storez==1 && p2>0))
 				for (j=0; j<(*zcol-p3); ++j) {
					for (i=0; i<N; ++i) {
 						z_store(count, j*N + i) = z(i, j);
                                        }
                                }
 			// store DIC
 			if (DIC_exists)
 				DIC_store(count, 0) = DIC;

      count++;
    } // End store results
        
		// allow user interrupts
    void R_CheckUserInterrupt(void);
  }		// end sampling-loop
   
  delete stream; // clean up random number stream

	// Berechnung von D(average(theta)) für DIC
	if (DIC_exists) { 
		DIC_zmedium			 = DIC_zmedium / (count*1.0);	   // ZEILE FÜR DIC1 mit latent Variables
		DIC_etamedium    = DIC_etamedium / (count*1.0);  // ZEILE FÜR DIC2 ohne latent Variables
		DIC_b0medium		 = DIC_b0medium / (count*1.0);
		if (Omega_exists)
			DIC_Amedium		 = DIC_Amedium / (count*1.0);
		DIC_Bmedium			 = DIC_Bmedium / (count*1.0);
		//DIC_thetamedium	 = DIC_thetamedium / (count*1.0);
		DIC_numedium	 	 = DIC_numedium / (count*1.0);
				
		DIC = 0.0;
				
		// Mittelwerte berechnen
				for (i=0; i<N; ++i) {
					if (Omega_exists)														
						DIC_mu1 = DIC_b0medium +  DIC_Amedium * Omega(i,_) + DIC_Bmedium * (!DIC_etamedium(i,_)); // ZEILE FÜR DIC2 ohne latent Variables
						//DIC_mu1 = DIC_b0medium +  DIC_Amedium * Omega(i,_) + DIC_Bmedium * (!DIC_zmedium(i,_)); // ZEILE FÜR DIC1 mit latent Variables
					else  
						DIC_mu1 = DIC_b0medium  								 					 + DIC_Bmedium * (!DIC_etamedium(i,_)); // ZEILE FÜR DIC2 ohne latent Variables				 	
						//DIC_mu1 = DIC_b0medium  								 					 + DIC_Bmedium * (!DIC_zmedium(i,_)); // ZEILE FÜR DIC1 mit latent Variables				 	
					
					for (j=0; j<p1; ++j) {   // Likelihood-Beitrag von ordinalen Indikatoren berechnen
						int obs_cat = static_cast<int>(Y(i,j));
						double prop = pnorm(DIC_numedium(j,obs_cat)-DIC_mu1(j,0)) - pnorm(DIC_numedium(j,obs_cat-1)-DIC_mu1(j,0));
						DIC = DIC + log(prop);
					}
					/*
					if (p1 < p) {
						for (int j=p1; j<p; ++j) { // Likelihood-Beitrag von metrischen Indikatoren berechnen
							
						}
					}
					*/
				}
				DIC = -2.0 * DIC;			
				
		DICestimate[0] = DIC;		
				
		DICav[0] = 0.0;
		for (i=0; i<count; i++)
			DICav[0] += DIC_store(i, 0);
		DICav[0] = DICav[0] / (count*1.0);
	}

	//Rprintf("\n\n D(average(theta)) = %f\n\n", DIC);
	
	    if (p1>0) {Rprintf("nu Schwellenwerte of z_ij:\n");
      for (i=0; i<DIC_numedium.rows(); ++i) {
				for (j=0; j<DIC_numedium.cols(); ++j)
					Rprintf("%f  ", DIC_numedium(i,j));
			Rprintf("\n");}Rprintf("\n");}
			
	
  // return output
  Matrix<double> output;

  if (*storescores==1) 
    output = cbind(eta_store, beta_store);
  else 
    output = beta_store;
  
  if (X_exists || X_par_exists)
  	output = cbind(output, gamma_intercept_store);
  
  if (X_exists) 
    output = cbind(output, gamma_store);

  if (X_par_exists) {
    output = cbind(output, gamma_par_store);
    output = cbind(output, tau_par_store);
  }
  
  if (Xdir_par_exists) {
    output = cbind(output, gammadir_par_store);
    output = cbind(output, taudir_par_store);
  }

 	//if (p3>0) 
 	output = cbind(output, theta_store);

	if (p1 > 0 && nu.cols() > 3) {
		output = cbind(output, nu_store);
	}
	
	if (*storez==1) {
  	output = cbind(output, z_store);
  }

	if (DIC_exists) {
  	output = cbind(output, DIC_store);
  }
  	
  const int size = *samplerow * *samplecol;
  for (int i=0; i<size; ++i)
    sampledata[i] = output[i];

  // Rprintf(" Startzeit: Hour:%d Min:%d Second:% d\n",st.wHour,st.wMinute,st.wSecond);
  // GetSystemTime(&en);
  // Rprintf(" Endzeit: Hour:%d Min:%d Second:% d\n",en.wHour,en.wMinute,en.wSecond);

	delete[] K_par;
	delete[] K_par_band;
	delete[] X_par;
	delete[] X_par_cross;
	delete[] X_par_cross_band;
	delete[] X_par_trans;
	delete[] gamma_par;
	delete[] tau_par;
	
	delete[] Kdir_par;
	delete[] Kdir_par_band;
	delete[] Xdir_par;
	delete[] Xdir_par_cross;
	delete[] Xdir_par_cross_band;
	delete[] Xdir_par_trans;
	delete[] gammadir_par;
	delete[] taudir_par;
	return;    
	
}

} // Ende extern c


  /* Band determination and creation of a symmetric nxn matrix */
  template <class T>
  Matrix<T>
  band_creation (const Matrix<T> &A){
  
    if (! A.isSquare() || A.rows()==1) {
      throw scythe_dimension_error(__FILE__, __PRETTY_FUNCTION__,
          __LINE__, "Matrix not square or 1x1 !");
    }
    
    // Determination of bandwidth
    int bandwidth = 0;
    for (int i = 1; i < A.rows(); ++i) {
    	for (int j = 0; j < i; ++j) {
    		if (A(i,j) != 0) {
    			if ( i-j > bandwidth) bandwidth = i-j;
    			break;
    		}
    	}
    }
		//Rprintf("Bandbreite = %d\n", bandwidth);
		
    // Initialization and construction of band matrix 
    Matrix<T> temp = Matrix<T> (A.rows(), bandwidth+1, true, 0.0);
    
    for (int i = 0; i < temp.rows(); ++i) {
    	for (int j = std::max(bandwidth-i,0) ; j < bandwidth+1; ++j) {
    		temp(i,j) = A(i, i-bandwidth+j);
    	}	
    }
  
    return temp;
  }
  
  // Reconversion from a band matrix to a normal matrix
  template <class T>
  Matrix<T>
  band_reconversion (const Matrix<T> &a)
  {
		Matrix<T> temp = Matrix<T>(a.rows(), a.rows(), true, 0.0);
 	 	int cols = a.cols();
  
 	 	// Reconversion      
 	  for (int i = 0; i < a.rows(); ++i) {
			for (int j = std::max(cols-i-1,0) ; j < cols; ++j) {
   			temp(i,i-cols+j+1) = a(i, j);
   		}	
 	 	}
  
		return temp;
  }

  //Band matrix addition
  template <class T>
  Matrix<T> band_addition (const Matrix<T> &a, const Matrix<T> &b)
	{
		if (a.rows() != b.rows()) {
     	throw scythe_dimension_error(__FILE__, __PRETTY_FUNCTION__,
       	  __LINE__, "Matrices don´t have equal number of rows !");
   	}
		
		//int max_cols;
		int diff_cols = a.cols() - b.cols();
		bool a_larger = false;
		Matrix<T> temp; 
		
		if (diff_cols >= 0) {
			//max_cols = a.cols();
			a_larger = true;
			temp = a;
		} 
		else {
			//max_cols = b.cols();
			temp = b;
		}
		
		if (a_larger == true) {
			for (register int i=0; i<a.rows(); ++i) {
				for (register int j=diff_cols; j<a.cols(); ++j) {
					temp(i,j) += b(i, j-diff_cols);
				}
			}
		} else {
			for (register int i=0; i<b.rows(); ++i) {
				for (register int j=-diff_cols; j<b.cols(); ++j) {
					temp(i,j) += a(i, j+diff_cols);
				}
			}
		}
		
		return temp;
	}
 /* invertiert in m x 1-dim Matrix gespeicherte Elemente einer
    Diagonalmatrix */
  template <class T>
  Matrix<T>
  inv_diag (const Matrix<T> &a){
	 Matrix<T> temp = Matrix<T>(a.rows(),1,true,0.0);
	 for(int i=0; i<a.rows(); ++i){
		temp(i,0) = 1/a(i,0);
	 }
	 return temp;
  }

/* multipliziert Diagonalmatrix mit Matrix, wobei elemente der
  Diegonalmatrix in m x 1 Matrix gespecihert sind */
  template <class T>
  Matrix<T>
  diag_mult_matrix (const Matrix<T> &a, const Matrix<T> &b){
	 if (a.rows() != b.rows()){
		throw scythe_dimension_error(__FILE__,__PRETTY_FUNCTION__,
			__LINE__, "Dimensionerror!");
	 }
	 Matrix<T> temp = Matrix<T>(a.rows(),b.cols(),true,0.0);
	 for (int j=0; j<a.rows(); ++j){
		for(int i=0; i<b.cols(); ++i){
			temp(j,i) = a(j,0)*b(j,i);
		}
	 }
   return temp;
  }
  
  
  /* multipliziert Matrix mit Diagonalmatrix , wobei Elemente der
  Diagonalmatrix in m x 1 Matrix gespeichert sind */
  template <class T>
  Matrix<T>
  matrix_mult_diag (const Matrix<T> &a, const Matrix<T> &b){
	 if (a.cols() != b.rows()){
		throw scythe_dimension_error(__FILE__,__PRETTY_FUNCTION__,
			__LINE__, "Dimensionerror!");
	 }
	 Matrix<T> temp = Matrix<T>(a.rows(),a.cols(),true,0.0);
	 for (int j=0; j<b.rows(); ++j){
		for(int i=0; i<a.rows(); ++i){
			temp(i,j) = a(i,j)*b(j,0);
		}
	 }
   return temp;
  }


  /* erzeugt aus in m x 1-dim Matrix abgespeicherten
	Elementen diagonalmatrix */ 
  template <class T>
  Matrix<T>
  create_diag (const Matrix<T> &a){
	 Matrix<T> temp = Matrix<T>(a.rows(),a.rows(),true,0.0);
	 for(int i=0; i<a.rows(); ++i){
		temp(i,i) = a(i,0);
	 }
	 return temp;
  }
  
	template <class T>
  Matrix<T> band_chol (const Matrix<T> &a)
  {
  	int b = a.cols()-1;
  	
  	if (b+1 >= a.rows()) {
     	throw scythe_dimension_error(__FILE__, __PRETTY_FUNCTION__,
       	  __LINE__, "Matrix is not a band matrix !");
   	}
   	  	
  	Matrix<T> temp = Matrix<T>(a.rows(), a.cols(), 1,0.0); 
  	if (b==0)  {
  		for (int i=0; i<a.rows(); ++i)
  			temp[i] = sqrt(a(i,0));                                   
  	}

  	register T h;
  	
  	temp(0,b) = sqrt(a(0,b));
  	for (int j=	b-1; j>=0; --j) {
  			temp(0+b-j,j) = a(0+b-j,j) / temp(0,b);
  	}
  	
  	for (int i=1; i<a.rows(); ++i) {
			// Calculation of l_ii or here l(i,b)
  		//temp(i,b) = sqrt(a(i,b));
			h = 0.0;

			for (int k=std::max(0,b-i); k<b; ++k)
				h += temp(i,k)*temp(i,k);

			temp(i,b) = sqrt(a(i,b)-h);
			
			//Calculation of l_ji für j>i

  		for (int j=b-1; j>=0; --j) {
				if (i+b-j >= a.rows()) break;
				h = 0.0;
  			for (int k=0; k<j; ++k) 
  				h += temp(i+b-j,j-k-1 ) * temp(i,b-k-1);

  			temp(i+b-j,j) = (a(i+b-j,j) - h) / temp(i,b);
  		}
  		
  	}  	
  	return temp;
  }
  
  template <class T>
  Matrix<T> 
  band_solve (const Matrix<T> &chol, const Matrix<T> &vector) {
  	if ((chol.rows() != vector.rows()) || (vector.cols() != 1)) {
   	 	throw scythe_dimension_error(__FILE__, __PRETTY_FUNCTION__,
       	  __LINE__, "Dimension error !");
  	}
  	int b = chol.cols()-1;
  	int N = chol.rows();
  	register T h;
  	Matrix<T> y = Matrix<T>(N,1, false);
  	Matrix<T> x = Matrix<T>(N,1, false);
  	
  	// Solve L*y = b  bzw. chol*y = vector
  	
  	y[0] = vector[0] / chol(0, b);
  	for (int i=1; i<N; ++i) {
  		h = 0.0;
  		for (int j=std::max(0,b-i); j<b; j++) {
  			h += chol(i,j) * y[i-b+j];
  		}
  		
  		y[i] = (vector[i] - h) / chol(i,b);
  	} 
  
  	//return y;
  	// Solve L'*x = y
  	x[N-1] = y[N-1] / chol(N-1, b);
  	
  	for (int i=N-2; i>-1; --i) {
  		h = 0.0;
  		for (int j=std::max(0,b-(N-1-i)); j<b; j++) {
  			h += chol(i+b-j,j) * x[i+b-j];
  		}
  		x[i] = (y[i] - h) / chol(i,b);
  	}   	
  	
  	return x;
  }
  
  template <class T>
  Matrix<T> 
  band_solve_backward (const Matrix<T> &chol, const Matrix<T> &y) {
  	if ((chol.rows() != y.rows()) || (y.cols() != 1)) {
   	 	throw scythe_dimension_error(__FILE__, __PRETTY_FUNCTION__,
       	  __LINE__, "Dimension error !");
  	}
  	int b = chol.cols()-1;
  	int N = chol.rows();
  	register T h;
  	Matrix<T> x = Matrix<T>(N,1, false);
  	
  	// Solve L'*x = y
  	x[N-1] = y[N-1] / chol(N-1, b);
  	
  	for (int i=N-2; i>-1; --i) {
  		h = 0.0;
  		for (int j=std::max(0,b-(N-1-i)); j<b; j++) {
  			h += chol(i+b-j,j) * x[i+b-j];
  		}
  		x[i] = (y[i] - h) / chol(i,b);
  	}   	
  	
  	return x;
  }



/*    Alte Version DIC berechnen
				DIC = DIC_offset;
				
				// Kovarianzmatrix zusammensetzen, Inverse und Determinante berechnen
				DIC_cov = (B*(!B))+theta;
				DIC_cov = cbind(DIC_cov, B);
				DIC_cov = rbind(DIC_cov, cbind(!B, DIC_identity ));
				
				DIC_det = ~DIC_cov;
				DIC_cov = invpd(DIC_cov);
				
				DIC = DIC - N*log(DIC_det)/2.0;   // Determinante zu DIC dazuzählen
				
				// Mittelwerte berechnen, zusammensetzen und Term zum DIC dazurechnen
				praed_rest = Matrix<double>(N, m, 1, 0.0);
				
				if (s>0) 
					for (int t=0; t<s; ++t) {
							praed_rest = praed_rest + X_par[t] * !gamma_par[t];
					}			
				if (X_exists)
							praed_rest = praed_rest + X * !Gamma; 
						
				
				
				
				
				for (int i=0; i<N; ++i) {
					DIC_mu2 = !praed_rest(i,_);
					
					if (Omega_exists)																		// in DIC_mu1 wird der Erwartungswert für y^*_i gespeichert
						DIC_mu1 = b0 +  A * Omega(i,_) + B * DIC_mu2;
					else  
						DIC_mu1 = b0  								 + B * DIC_mu2;					
					
					DIC_mu = rbind(DIC_mu1, DIC_mu2);   // Erwartungswert der Verteilung von y^* und z
					
					DIC_mu1 = !z(i,_);									// DIC_mu1 enthält die underlying variables y^*_i
					DIC_mu2 = !eta(i,_);						 		// DIC_mu2 enthält die latenten scores z_i
					DIC_mu2 = rbind(DIC_mu1, DIC_mu2);  // DIC_mu2 enthält oben y^*_i und unten z_i
					
					DIC_mu = DIC_mu2 - DIC_mu;
					DIC_mu2 = 0.5 * (!DIC_mu) * DIC_cov * DIC_mu;
					DIC -= DIC_mu2(0,0); 
						//DIC_mu2(r,0) = praed_rest;
				}
				DIC = -2.0 * DIC;			
				
				// Mittelwerte für die Berechnung von D(average(theta)) speichern
				DIC_ysternmedium = DIC_ysternmedium + z;
				DIC_zmedium			 = DIC_zmedium + eta;
				DIC_b0medium	   = DIC_b0medium + b0;
				if (Omega_exists) 
					DIC_Amedium		 = DIC_Amedium + A;
				DIC_Bmedium	  	 = DIC_Bmedium + B;					
				DIC_thetamedium	 = DIC_thetamedium + theta;		
				DIC_praed_rest   = DIC_praed_rest + praed_rest;		
				*/


