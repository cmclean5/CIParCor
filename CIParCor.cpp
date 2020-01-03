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
export PKG_CXXFLAGS="-w -fopenmp -I$R_LIBS_USER/Rcpp/include -I$R_LIBS_USER/RcppEigen/include"
R CMD SHLIB CIParCor.cpp geninv.cpp

option '-w' in PKG_CXXFLAGS to suppress all RcppEigen's warnings, see for more details: 
https://gcc.gnu.org/onlinedocs/gcc/Warning-Options.html

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
> dyn.unload("/users/cmclean/STUDIES/RCcpp/v5r1/CIParCor.so")
*/

#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <limits.h>
#include <omp.h>

//Moore-Penrose's Matrix inversion algorithm written in RcppEigen
#include "geninv.h"


//namespace std;
//namespace Rcpp;

//--- RCpp Extensions

// Enable cpp11
// [[Rcpp::plugins("cpp11")]]

// Enable OpenMP (excludes macOS)
// [[Rcpp::plugins(openmp)]]

//--- 


//---------------------------------------------------------------
// Keep code, but will try implementation in RcppEigen. 
// Moore–Penrose inverse in C++
// Reference:
// http://fractalytics.io/moore-penrose-matrix-optimization-cuda-c
// Used to find pesudo inverse the coefficients of a linear regression.
// if,                       Y = B * X + C
// then pesudo-inverse is, X^+ = (X^t * X)^-1 * X^t
// thus,                     B = X^+ * Y
//---------------------------------------------------------------
/*
void Transpose(float *matrix, float *t_matrix, int N){

  int i,j;
  
    for( i=0; i<N; i++ ){
      for( j=0; j<N; j++ ){
	t_matrix[(j*N)+i]=matrix[(i*N)+j]; 
      }
    }
    
    return;
}


void MatrixMult(float *matrix_1, float *matrix_2, float *matrix_product, int N){

  int i,j,k;
  for(i=0; i<N; i++){
    for(j=0; j<N; j++){             // not j<M
      matrix_product[(i*N)+j] = 0;
      for(k=0; k<N; k++){
	matrix_product[(i*N)+j] += matrix_1[(i*N)+k] *
	                           matrix_2[(k*N)+j];
      }
    }
  }
  
  return;
 
}
 

// Function to get cofactor 
void getCofactor(float *A, float *temp, int p, int q, int n, int N){

  int i,j,row,col;
  
  i = j = 0;
 
  // Looping for each element of the matrix
  for(row=0; row<n; row++){
    for(col=0; col<n; col++){
      // Copying into temporary matrix only those element
      // which are not in given row and column
      if (row != p && col != q){
	temp[(i*N)+j++] = A[(row*N)+col];
 
	// Row is filled, so increase row index and
	// reset col index
	if (j == (n - 1) ){
	  j = 0;
	  i++;
	}
      }
    }
  }

}


// Recursive function for finding determinant of matrix.
int determinant(float *A, int n, int N){

  int f, D, sign;
  
  D = 0; // Initialize result
 
  // Base case : if matrix contains single element
  if( n == 1 )
    return A[0];
 
  float temp[N*N]; // To store cofactors
 
  sign = 1; // To store sign multiplier
 
  // Iterate for each element of first row
  for(f=0; f<n; f++){
    
    // Getting Cofactor of A[0][f]
    getCofactor(A, temp, 0, f, n, N);
    D += sign * A[(0*N)+f] * determinant(temp, (n-1), N);
    
    // terms are to be added with alternate sign
    sign = -sign;

  }
 
  return D;

}

// Function to get adjoint
void adjoint(float *A, float *adj, int N){

  int i,j, sign;
  
  if (N == 1){
    adj[0] = 1;
    return;
  }
 
  // temp is used to store cofactors 
  sign = 1;
  float temp[N*N];
 
  for(i=0; i<N; i++){
    for(j=0; j<N; j++){

      // Get cofactor
      getCofactor(A, temp, i, j, N, N);
 
      // sign of adj positive if sum of row
      // and column indexes is even.
      sign = ((i+j)%2==0) ? 1 : -1;
 
      // Interchanging rows and columns to get the
      // transpose of the cofactor matrix
      adj[(j*N)+i] = (sign)*(determinant(temp, (N-1), N));
    }
  }

}

// Function to calculate and store inverse, returns false if
// matrix is singular
bool inverse(float *A, float *inverse, int N){

  int i,j;
  
  // Find determinant of A[][]
  int det = determinant(A, N, N); 
  if (det == 0){
    std::cout << "Singular matrix, can't find its inverse" << std::endl;
    return false;
  }
  
  // Find adjoint
  float adj[N*N];
  adjoint(A, adj, N);
 
  // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
  for(i=0; i<N; i++){
    for(j=0; j<N; j++){
      inverse[(i*N)+j] = adj[(i*N)+j]/float(det);
    }
  }
      
    return true;

}


// Generic function to display the matrix. We use it to display
// both adjoin and inverse. adjoin is integer matrix and inverse
// is a float.
template<class T>
void display(T *A, int N){

  int i,j;

  for(i=0; i<N; i++){
    for(j=0; j<N; j++)
      std::cout << A[(i*N)+j] << " ";
    std::cout << std::endl;
  }

}


// Driver program
void MRexample(){

  int i,j,N;

  N = 7;

  float A[N][N] = { {5, -2, 2, 7,9,8,0},
		    {1, 0, 0, 3,1,0,9},
		    {-3, 1, 5, 0,2,1,7},
		    {3, -1, -9, 4,6,5,2},
		    {1, 0, 4, 4,1,0,9},
		    {8, 0, 3, 8,6,5,2},
		    {5, 6, 4, 1,3,2,0}
                   };
 
    float *matrix         = new float[N*N];
    float *t_matrix       = new float[N*N];
    float *matrix_product = new float[N*N];
    float *pseudoinverse  = new float[N*N];
    float *adj            = new float[N*N]; // To store adjoint 
    float *inv            = new float[N*N]; // To store inverse 

    for(i=0; i<N; i++){
      for(j=0; j<N; j++){
	matrix[(i*N)+j] = A[i][j];
      }
    }
    
    Transpose(matrix, t_matrix, N);
    std::cout << "\nThe Transpose is :\n";
    display(t_matrix,N);
 
    std::cout << "The product of the matrix is: " << std::endl;
    MatrixMult(t_matrix, matrix, matrix_product, N);
    display(matrix_product,N);
 
    std::cout << "\nThe Inverse is :\n";
    if (inverse(matrix_product, inv, N))
      display(inv,N);
 
    MatrixMult(inv, t_matrix, pseudoinverse, N);
    std::cout << "\nThe Monroe-penrose inverse is :\n";
    display(pseudoinverse,N);
}

// END Moore–Penrose inverse in C++
//-----
*/


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

// Calculate partial correlation from inversion of the
// Covariance matrix, using Moore-Penrose algorithm.
// [[Rcpp::export]]
RcppExport SEXP parCov(SEXP CC, SEXP verbose )
{

  int n,print;

  Rcpp::NumericMatrix rCC;
  
  print = 0;
  if( !Rf_isNull(verbose) ){    
    Rcpp::NumericVector rVerbose(verbose);
    print = rVerbose[0];
    if( print != 0 ){ print = 1; }
  }

  // no omp used here for the moment
  int threads = 0;
  omp_set_num_threads(threads);
  if( print ){
    std::cout << "> OpenMP: " << " number of threads " << threads << std::endl;
  }

  if( !Rf_isNull(CC) ){

    if( print ){ std::cout << "Invert Covariance matrix using Moore-Penrose algorithm... "; }

    Rcpp::NumericMatrix temp = geninv(CC);

    if( print ){ std::cout << "...done. " << std::endl; }
       
    if( print ){ std::cout << "Calculate pairwise Partial Correlations... "; }
    
    rCC  = cov2cor(temp, -1.0);
    
    if( print ){ std::cout << "...done. " << std::endl; }
      
  } else {
  
    n = 1;

    rCC = Rcpp::NumericMatrix(n,n);
    
  }

  //Rcpp::NumericVector rRES = Rcpp::wrap(cRES);
  //rRES.attr("dim")         = Rcpp::Dimension(n,n);

  Rcpp::List OUT = Rcpp::List::create(Rcpp::Named("PC") =rCC);

  return OUT;
  
}

//----------  

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

    //struct to store min Conditional Independed Partial Correlation 
    struct minRso { double param1, param2, param3, param4; };

    //we'll need to define our own reduction clause
    #pragma omp declare reduction(storeMin : minRso : omp_in.param1 < omp_out.param1 ? omp_out = omp_in : omp_out)
    
    for(P=0; P<(n*n); P++){

      i = floor(P/n);
      j = P % n;

      //Correlation matrix symmetric... take the upper triangle
      if( (j > i) && (!isnan(cGG[(i*n)+j])) ){
	
	double Rtemp, Rso, Rk, Rq;
	Rtemp = Rso = Rk = Rq = std::numeric_limits<double>::quiet_NaN();

	minRso result{ std::numeric_limits<double>::infinity(),
		       std::numeric_limits<double>::infinity(),
		       std::numeric_limits<double>::infinity(),
		       std::numeric_limits<double>::infinity() };
	
  #pragma omp parallel for collapse(2) schedule(static) default(none) firstprivate(i,j,n,Rtemp,Rso,Rk,Rq) private(k,K,q,Q) \
  shared(cGG,cCC) reduction(storeMin:result)
	for( K=0; K<n; K++ ){
	  for( Q=0; Q<n; Q++ ){
	      
	    q = Q % n;
	    k = K % n;
	    
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
		    ////Rtemp  = (num1abs/dem1) * (num2abs/dem2);
		    Rtemp  = 0.5*( (num1abs/dem1) + (num2abs/dem2) );//Pij.kq and Pij.qk should be symmetric
		    Rtemp  = Rtemp * Rtemp;
		    Rso    = num1 / dem1;
		    Rk     = k;
		    Rq     = q;
		    result = minRso{Rtemp, Rso, Rk, Rq};
		  } else {
		    ////double param1 = (num1abs/dem1) * (num2abs/dem2);
		    double param1 = 0.5*( (num1abs/dem1) + (num2abs/dem2));
		    param1 = param1 * param1;
		    if( param1 < result.param1 ){
		      Rtemp  = param1; Rso = (num1/dem1);  Rk = k; Rq = q;
		      result = minRso{Rtemp, Rso, Rk, Rq};
		    }
		  }			
		}
	      }
		
	    }//i!=k & j!=k
	  }//Q
	}//K


	cPR[(i*n)+j]  = result.param1;//Rtemp
	cPR[(j*n)+i]  = result.param1;
	
	cRES[(i*n)+j] = result.param2;//Rso
	cRES[(j*n)+i] = result.param2;
	  
	cK[(i*n)+j]   = result.param3;//Rk
	cK[(j*n)+i]   = result.param3;

	cQ[(i*n)+j]   = result.param4;//Rq
	cQ[(j*n)+i]   = result.param4;
	  
	if( !isnan(samples) ){
	  double norm = samples - 4;
	  if( norm >= 0 ){
	    cSD[(i*n)+j] = (1-(result.param2*result.param2)) / (sqrt(norm));
	    cSD[(j*n)+i] = cSD[(i*n)+j];
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

    //struct to store min Conditional Independed Partial Correlation 
    struct minRfo { double param1, param2, param3; };

    //we'll need to define our own reduction clause
    #pragma omp declare reduction(storeMin : minRfo : omp_in.param1 < omp_out.param1 ? omp_out = omp_in : omp_out)

    
    //working code here

    for(P=0; P<(n*n); P++){

      i = floor(P/n);
      j = P % n;
      
      //Correlation matrix symmetric... take the upper triangle
      if( (j > i) && (!isnan(cGG[(i*n)+j])) ){

	double Rtemp, Rfo, Rk;
	Rtemp = Rfo = Rk = std::numeric_limits<double>::quiet_NaN();

	minRfo result{ std::numeric_limits<double>::infinity(),
		       std::numeric_limits<double>::infinity(),
		       std::numeric_limits<double>::infinity() };
		       

	#pragma omp parallel for schedule(static) default(none) firstprivate(i,j,n,Rtemp,Rfo,Rk) private(k,K) \
        shared(cGG,cCC) reduction(storeMin:result)
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
			Rtemp  = numABS / dem;
			Rtemp  = Rtemp * Rtemp;
			Rfo    = num / dem;
			Rk     = k;
			result = minRfo{Rtemp, Rfo, Rk};
		      } else {
			double param1 = numABS / dem;
			param1        = param1 * param1;
			if( param1 < result.param1 ){
			  Rtemp  = param1; Rfo = (num/dem); Rk = k;
			  result = minRfo{Rtemp, Rfo, Rk};
			}
		      }
		      
		    }
		  		  
		
	  }//i!=k & j!=k
	      
	}//k

	cPR[(i*n)+j]  = result.param1; //Rtemp
	cPR[(j*n)+i]  = result.param1; 
	
	cRES[(i*n)+j] = result.param2; //Rfo
	cRES[(j*n)+i] = result.param2; 

	cK[(i*n)+j]   = result.param3; //Rk
	cK[(j*n)+i]   = result.param3; 
	      
	if( !isnan(samples) ){
	  double norm = samples - 3;
	  if( norm >= 0 ){
	    cSD[(i*n)+j] = (1-(result.param2*result.param2)) / (sqrt(norm));
	    cSD[(j*n)+i] = cSD[(i*n)+j];
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
    
 
