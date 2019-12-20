/*

References:

[1]https://matloff.wordpress.com/2015/01/16/openmp-tutorial-with-r-interface/
[2] Magwene, P. M. & Kim, J. "Estimating genomic coexpression networks using first-order conditional independance", Genome Biology, 5:R100, (2004).
[3] Zuo, Y. et. al. "Biological network inference using low order partial correlation", Methods, 69(3):266-273, (2014). doi:10.1016/j.ymeth.2014.06.010
[4] http://www.parallelr.com/r-and-openmp-boosting-compiled-code-on-multi-core-cpu-s/
[5] https://bisqwit.iki.fi/story/howto/openmp/

* build instructions:
1)
export R_LIBS_USER=/users/cmclean/R/x86_64-pc-linux-gnu-library/3.6/ # my R libraries
export PKG_LIBS="-lgomp"
export PKG_CXXFLAGS="-fopenmp -I/users/cmclean/R/x86_64-pc-linux-gnu-library/3.6/Rcpp/include"
R CMD SHLIB CIParCor.cpp

Or

2) 
PKG_CXXFLAGS="$(echo 'Rcpp:::CxxFlags()'| R --vanilla --slave) -fopenmp" R CMD SHLIB CIParCor.cpp


* to find the number of all installed cores/processors in linux: 
  nproc --all

* before start R, set number of threads, e.g. 

  export OMP_NUM_THREADS=2

in R:

> library(Rcpp)
> dyn.load("/users/cmclean/STUDIES/RCcpp/v5r1/CIParCor.so")
> GG = matrix(1,ncol=2,nrow=2)
> CC = GG
> .Call("firstOrder",GG,CC,NULL)
> dyn.unload("/users/cmclean/STUDIES/RCcpp/v5r1/CIParCorcpp.so")
*/


#include <Rcpp.h>
#include <iostream>
#include <math.h>
#include <limits.h>
#include <omp.h>

//namespace std;
//namespace Rcpp;

//--- RCpp Extensions
// Enable cpp11
// [[Rcpp::plugins("cpp11")]]

// Enable OpenMP (excludes macOS)
// [[Rcpp::plugins(openmp)]]
//--- 

//set omp env before running our subroutines
void setOMPenv(int print){

  //--- Setup OpenMP
  int max_threads = 0;
  int num_procs   = 0;
  max_threads     = omp_get_max_threads();
  num_procs       = omp_get_num_procs();
  if( print ){
    std::cout << "> OpenMP: " << " Set " << std::endl;
    std::cout << "> OpenMP: " << " Number of processors available " << num_procs << std::endl;
    std::cout << "> OpenMP: " << " Max number of threads " << max_threads << std::endl;
  }
  omp_set_num_threads(max_threads);
  //---

  
}


//create square matrix and fill with NA's
Rcpp::NumericMatrix na_matrix(int n){
  Rcpp::NumericMatrix m(n,n) ;
  std::fill( m.begin(), m.end(), Rcpp::NumericVector::get_na() ) ;
  return m;
}



void calSecondOrder(double pij,   double piq,   double pik,
		    double pjq,   double pjk,   double pqk,
		    double &pijk, double &piqk, double &pjqk,
		    double &pijq, double &pikq, double &pjkq ){

    
  double num, dem, pkq;

  //symmetric
  pkq  = pqk;
  
  pijk = piqk = pjqk = pijq = pikq = pjqk = std::numeric_limits<double>::quiet_NaN();
  
  //for pijk = [ pij, pik, pjk ]
  num = pij - pik*pjk;
  dem = sqrt( (1-pik*pik) ) * sqrt( (1-pjk*pjk) );
  if( (dem != 0) && !isnan(dem) ){
    pijk = num / dem;
  }

  //for piqk = [ piq, pik, pqk ]
  num = piq - pik*pqk;
  dem = sqrt( (1-pik*pik) ) * sqrt( (1-pqk*pqk) );
  if( (dem != 0) && !isnan(dem) ){
    piqk = num / dem;
  }

  //for pjqk = [ pjq, pjk, pqk ]
  num = pjq - pjk*pqk;
  dem = sqrt( (1-pjk*pjk) ) * sqrt( (1-pqk*pqk) );
  if( (dem != 0) && !isnan(dem) ){
    pjqk = num / dem;
  }

  //for pijq = [ pij, piq, pjq ]
  num = pij - piq*pjq;
  dem = sqrt( (1-piq*piq) ) * sqrt( (1-pjq*pjq) );
  if( (dem != 0) && !isnan(dem) ){
    pijq = num / dem;
  }


  //for pikq = [ pik, piq, kq ]
  num = pik - piq*pkq;
  dem = sqrt( (1-piq*piq) ) * sqrt( (1-pkq*pkq) );
  if( (dem != 0) && !isnan(dem) ){
    pikq = num / dem;
  }

  //for pjkq = [ pjk, pjq, pkq ]
  num = pjk - pjq*pkq;
  dem = sqrt( (1-pjq*pjq) ) * sqrt( (1-pkq*pkq) );
  if( (dem != 0) && !isnan(dem) ){
    pjkq = num / dem;
  }

  
}


// [[Rcpp::export]]
RcppExport SEXP secondOrder(SEXP GG, SEXP CC, SEXP N, SEXP nCPU, SEXP verbose )
{

  int i,j,P,k,K,q,Q,n,print;

  std::vector<double> cRES;
  std::vector<double> cPR;
  std::vector<double> cSD;
  std::vector<double> cK;
  std::vector<double> cQ;
  
  print = 0;
  if( !Rf_isNull(verbose) ){    
    Rcpp::NumericVector rVerbose(verbose);
    print = rVerbose[0];
    if( print != 0 ){ print = 1; }
  }
  
  if( !Rf_isNull(nCPU) ){
    int threads = 0;
    Rcpp::NumericVector rCPU(nCPU);
    threads = rCPU[0];
    omp_set_num_threads(threads);
    if( print ){
      std::cout << "> OpenMP: " << " number of threads " << threads << std::endl;
    }
  } else {  
    //--- default openmp setup
    setOMPenv( print );
  }

  if( !Rf_isNull(GG) && !Rf_isNull(CC) ){    
      
    Rcpp::NumericMatrix rGG(GG);
    Rcpp::NumericMatrix rCC(CC);
   
    //N == number of samples.
    double samples = std::numeric_limits<double>::quiet_NaN();
    if( !Rf_isNull(N) ){
      Rcpp::NumericVector rN(N);
      samples = rN[0];
    }    

    //get length of the square matrix
    n = rCC.nrow();

    //1) flatting R matrices to C++ 1d vectors
    std::vector<double> cGG = Rcpp::as<std::vector< double > >(rGG);
    std::vector<double> cCC = Rcpp::as<std::vector< double > >(rCC);

    cRES = std::vector<double>((n*n));
    cPR  = std::vector<double>((n*n));
    cSD  = std::vector<double>((n*n));
    cK   = std::vector<double>((n*n));
    cQ   = std::vector<double>((n*n));

    // initialise 
    #pragma omp parallel for schedule(static) default(none) firstprivate(n) private(k) shared(cRES,cPR,cSD,cK,cQ)
    for(k=0; k<(n*n); k++){
      cRES[k] = std::numeric_limits<double>::quiet_NaN();
      cPR[k]  = std::numeric_limits<double>::quiet_NaN();
      cSD[k]  = std::numeric_limits<double>::quiet_NaN();
      cK[k]   = std::numeric_limits<double>::quiet_NaN();
      cQ[k]   = std::numeric_limits<double>::quiet_NaN();      
    }
    
    //working code here
    #pragma omp parallel for schedule(static) default(none) firstprivate(n,samples) private(i,j,P,k,K,q,Q) shared(cGG,cCC,cRES,cPR,cSD,cK,cQ)
    for(P=0; P<(n*n); P++){

      i = floor(P/n);
      j = P % n;

      if( (i != j) && (!isnan(cGG[(i*n)+j])) ){
	//if( (i > j) && (!isnan(cGG[(i*n)+j])) ){
	
	double Rtemp, Rso, Rk, Rq;
	Rtemp = Rso = Rk = Rq = std::numeric_limits<double>::quiet_NaN();

	for( K=0; K<n; K++ ){
	  
	  k = K % n;
	  
	  for( Q=0; Q<n; Q++ ){
	      
	    q = Q % n;
		
	    if(  (i != k) && (j != k) &&
		 (i != q) && (j != q) &&
		 (k != q) &&		  
		 !isnan(cGG[(i*n)+j]) &&
		 !isnan(cGG[(i*n)+q]) &&
		 !isnan(cGG[(i*n)+k]) &&
		 !isnan(cGG[(j*n)+q]) &&
		 !isnan(cGG[(j*n)+k]) &&
		 !isnan(cGG[(q*n)+k]) ){
	      
	      double pijk, piqk, pjqk, pijq, pikq, pjkq;
		    
	      calSecondOrder( cCC[(i*n)+j], cCC[(i*n)+q], cCC[(i*n)+k],
			      cCC[(j*n)+q], cCC[(j*n)+k], cCC[(q*n)+k],
			      pijk, piqk, pjqk,
			      pijq, pikq, pjkq );
		    
	      if( !isnan(pijk) && !isnan(piqk) && !isnan(pjqk) &&
		  !isnan(pijq) && !isnan(pikq) && !isnan(pjkq) ){
		
		double num1    = pijk - (piqk*pjqk);
		double num1abs = abs(pijk) - (abs(piqk)*abs(pjqk));
		double dem1    = sqrt( (1-piqk*piqk) ) * sqrt( (1-pjqk*pjqk) );
		
		double num2    = pijq - (pikq*pjkq);
		double num2abs = abs(pijq) - (abs(pikq)*abs(pjkq));
		double dem2    = sqrt( (1-pikq*pikq) ) * sqrt( (1-pjkq*pjkq) );
						
		      
		if( (dem1 != 0) && (!isnan(dem1)) &&
		    (dem2 != 0) && (!isnan(dem2)) ){
		  
		  if( isnan(Rso) ){
		    //Rtemp = 0.5*( (num1abs/dem1) + (num2abs/dem2) );//Pij.kq and Pij.qk mutually exclusive?
		    Rtemp = (num1abs/dem1) * (num2abs/dem2);
		    Rso   = num1 / dem1;
		    Rk    = k;
		    Rq    = q;
		  } else {
		    //double r = 0.5*( (num1abs/dem1) + (num2abs/dem2));
		    double r = (num1abs/dem1) * (num2abs/dem2);
		    if( r < Rtemp ){ Rtemp = r; Rso = (num1/dem1);  Rk = k; Rq = q; }
		  }			
		}
	      }
		
	    }//i!=k & j!=k
	  }//Q
	}//K

	      
	cRES[(i*n)+j] = Rso;
	//cRES[(j*n)+i] = Rso;
	  
	cPR[(i*n)+j]  = Rtemp;
	//cPR[(j*n)+i]  = Rtemp;
	
	cK[(i*n)+j]   = Rk;
	//cK[(j*n)+i]   = Rk;

	cQ[(i*n)+j]   = Rq;
	//cQ[(j*n)+i]  = Rq;
	  
	if( !isnan(samples) ){
	  double norm = samples - 4;
	  if( norm >= 0 ){
	    cSD[(i*n)+j] = (1-Rso*Rso) / (sqrt(norm));
	    //cSD[(j*n)+i] = cSD[(i*n)+j];
	  }
	}
	  
      }//i!=j
	
    }//P      

  } else {

    n = 1;

    cRES = std::vector<double>((n*n));
    cPR  = std::vector<double>((n*n));
    cSD  = std::vector<double>((n*n));
    cK   = std::vector<double>((n*n));
    cQ   = std::vector<double>((n*n));
    
    for(k=0; k<(n*n); k++){
      cRES[k] = std::numeric_limits<double>::quiet_NaN();
      cPR[k]  = std::numeric_limits<double>::quiet_NaN();
      cSD[k]  = std::numeric_limits<double>::quiet_NaN();
      cK[k]   = std::numeric_limits<double>::quiet_NaN();
      cQ[k]   = std::numeric_limits<double>::quiet_NaN();      
    }
    
  }//if/else    

  
  Rcpp::NumericVector rRES = Rcpp::wrap(cRES);
  rRES.attr("dim")         = Rcpp::Dimension(n,n);
   
  Rcpp::NumericVector rPR  = Rcpp::wrap(cPR);
  rPR.attr("dim")          = Rcpp::Dimension(n,n);
   
  Rcpp::NumericVector rSD  = Rcpp::wrap(cSD);
  rSD.attr("dim")          = Rcpp::Dimension(n,n);
 
  Rcpp::NumericVector rK   = Rcpp::wrap(cK);
  rK.attr("dim")           = Rcpp::Dimension(n,n);
 
  Rcpp::NumericVector rQ   = Rcpp::wrap(cQ);
  rQ.attr("dim")           = Rcpp::Dimension(n,n);
   

  Rcpp::List OUT = Rcpp::List::create(Rcpp::Named("PC") =rRES,
				      Rcpp::Named("MPC")=rPR,
				      Rcpp::Named("SD") =rSD,
				      Rcpp::Named("k")  =rK,
				      Rcpp::Named("q")  =rQ);

  return OUT;

}
    
 

// SEXP is the internal data type for R objects, as seen from C/C++;
// here are input is an R matrix adjm, and the return value is another
// R matrix
// If you want the function to exist in an R session, included the
// following 'commented' line below and directly above the function def:
// [[Rcpp::export]]
// ----------------------------------
// "first order conditional independent partial correlation"
// ----------------------------------
// The partial correlation is what we get when we hold constant some third
// variable from two other variables. We know the correlation between i and
// j is .88. But k "accounts for" (or could account for) part of that. What would happen to the correlation if k were constant?
// We partial one variable out of a correlation.
RcppExport SEXP firstOrder(SEXP GG, SEXP CC, SEXP N, SEXP nCPU, SEXP verbose )
{

   int i,j,P,k,K,n,print;

   std::vector<double> cRES;
   std::vector<double> cPR;
   std::vector<double> cSD;
   std::vector<double> cK;
   
   print = 0;
   if( !Rf_isNull(verbose) ){    
     Rcpp::NumericVector rVerbose(verbose);
     print = rVerbose[0];
     if( print != 0 ){ print = 1; }
   }
   
   if( !Rf_isNull(nCPU) ){
     int threads = 0;
     Rcpp::NumericVector rCPU(nCPU);
     threads = rCPU[0];
     omp_set_num_threads(threads);
     if( print ){
       std::cout << "> OpenMP: " << " number of threads " << threads << std::endl;
     }
   } else {  
     //--- default openmp setup
     setOMPenv( print );
   }


    if( !Rf_isNull(GG) && !Rf_isNull(CC) ){    
      
    Rcpp::NumericMatrix rGG(GG);
    Rcpp::NumericMatrix rCC(CC);
   
    //N == number of samples.
    double samples = std::numeric_limits<double>::quiet_NaN();
    if( !Rf_isNull(N) ){
      Rcpp::NumericVector rN(N);
      samples = rN[0];
    }    

    //get length of the square matrix
    n = rCC.nrow();

    //1) flatting R matrices to C++ 1d vectors
    std::vector<double> cGG = Rcpp::as<std::vector< double > >(rGG);
    std::vector<double> cCC = Rcpp::as<std::vector< double > >(rCC);

    cRES = std::vector<double>((n*n));
    cPR  = std::vector<double>((n*n));
    cSD  = std::vector<double>((n*n));
    cK   = std::vector<double>((n*n));
   
    // initialise 
    #pragma omp parallel for schedule(static) default(none) firstprivate(n) private(k) shared(cRES,cPR,cSD,cK)
    for(k=0; k<(n*n); k++){
      cRES[k] = std::numeric_limits<double>::quiet_NaN();
      cPR[k]  = std::numeric_limits<double>::quiet_NaN();
      cSD[k]  = std::numeric_limits<double>::quiet_NaN();
      cK[k]   = std::numeric_limits<double>::quiet_NaN();
    }

    //working code here
     #pragma omp parallel for schedule(static) default(none) firstprivate(n,samples) private(i,j,P,k,K) shared(cGG,cCC,cRES,cPR,cSD,cK)
    for(P=0; P<(n*n); P++){

      i = floor(P/n);
      j = P % n;
      
      if( (i != j) && (!isnan(cGG[(i*n)+j])) ){
	//if( (i > j) && (!isnan(cGG[(i*n)+j])) ){

	double Rtemp, Rfo, Rk;
	Rtemp = Rfo = Rk = std::numeric_limits<double>::quiet_NaN();

	for( K=0; K<n; K++ ){

	  k = K % n;
	  
	  if( (i != k) && (j != k) &&
	      !isnan(cGG[(i*n)+k]) &&
	      !isnan(cGG[(j*n)+k]) ){
	      
	    double num    = cCC[(i*n)+j] - (cCC[(i*n)+k]*cCC[(j*n)+k]);
	    double numABS = abs(cCC[(i*n)+j]) - ( abs(cCC[(i*n)+k])*abs(cCC[(j*n)+k]) );
	    double dem    = sqrt( (1-cCC[(i*n)+k]*cCC[(i*n)+k]) ) *
	      sqrt( (1-cCC[(j*n)+k]*cCC[(j*n)+k]) );

		    if( (dem != 0) && (!isnan(dem)) ){

		      if( isnan(Rfo) ){
			Rtemp = numABS / dem;
			Rfo   = num / dem;
			Rk    = k;
		      } else {
			double r = numABS / dem;
			if( r < Rtemp ){ Rtemp = r; Rfo = (num/dem); Rk = k; }
		      }
		      
		    }
		  
		  
		
	  }//i!=k & j!=k
	      
	}//k

	
	cRES[(i*n)+j] = Rfo;
	//cRES[(j*n)+i] = Rfo;

	cPR[(i*n)+j]  = Rtemp;
	//cPR[(j*n)+i]  = Rtemp;

	cK[(i*n)+j]   = Rk;
	//rK[(j*n)+i]  = Rk;
	      
	if( !isnan(samples) ){
	  double norm = samples - 3;
	  if( norm >= 0 ){
	    cSD[(i*n)+j] = (1-Rfo*Rfo) / (sqrt(norm));
	    //cSD[(j*n)+i] = cSD[(i*n)+j];
	  }
	}
	
      }//i!=j
	
    }//P

    
    } else {

    n = 1;

    cRES = std::vector<double>((n*n));
    cPR  = std::vector<double>((n*n));
    cSD  = std::vector<double>((n*n));
    cK   = std::vector<double>((n*n));
        
    for(k=0; k<(n*n); k++){
      cRES[k] = std::numeric_limits<double>::quiet_NaN();
      cPR[k]  = std::numeric_limits<double>::quiet_NaN();
      cSD[k]  = std::numeric_limits<double>::quiet_NaN();
      cK[k]   = std::numeric_limits<double>::quiet_NaN();
    }
    
  }//if/else    

  
  Rcpp::NumericVector rRES = Rcpp::wrap(cRES);
  rRES.attr("dim")         = Rcpp::Dimension(n,n);
   
  Rcpp::NumericVector rPR  = Rcpp::wrap(cPR);
  rPR.attr("dim")          = Rcpp::Dimension(n,n);
   
  Rcpp::NumericVector rSD  = Rcpp::wrap(cSD);
  rSD.attr("dim")          = Rcpp::Dimension(n,n);
 
  Rcpp::NumericVector rK   = Rcpp::wrap(cK);
  rK.attr("dim")           = Rcpp::Dimension(n,n);
 
 
  Rcpp::List OUT = Rcpp::List::create(Rcpp::Named("PC") =rRES,
				      Rcpp::Named("MPC")=rPR,
				      Rcpp::Named("SD") =rSD,
				      Rcpp::Named("k")  =rK);

  return OUT;

}
    
 
