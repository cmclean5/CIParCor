/*

References:

[1] Adapted from "OpenMP Tutorial, with R Interface", Matloff, 16/01/2015.
https://matloff.wordpress.com/2015/01/16/openmp-tutorial-with-r-interface/
[2] Magwene, P. M. & Kim, J. "Estimating genomic coexpression networks using first-order conditional independance", Genome Biology, 5:R100, (2004).
[3] Zuo, Y. et. al. "Biological network inference using low order partial correlation", Methods, 69(3):266-273, (2014). doi:10.1016/j.ymeth.2014.06.010

* build instructions:

export R_LIBS_USER=/users/cmclean/R/x86_64-pc-linux-gnu-library/3.6/ # my R libraries
export PKG_LIBS="-lgomp"
export PKG_CXXFLAGS="-fopenmp -I/users/cmclean/R/x86_64-pc-linux-gnu-library/3.6/Rcpp/include"
R CMD SHLIB CIPParCor.cpp

* to find the number of all installed cores/processors in linux: 
  nproc --all

* before start R, set number of threads, e.g. 

  export OMP_NUM_THREADS=2

in R:

> library(Rcpp)
> dyn.load("/users/cmclean/STUDIES/RCcpp/v4r1/CIParCor.so")
> GG = matrix(1,ncol=2,nrow=2)
> CC = GG
> .Call("firstOrder",GG,CC,NULL)
> dyn.unload("/users/cmclean/STUDIES/RCcpp/v4r1/CIParCorcpp.so")
*/


#include <Rcpp.h>
#include <math.h>
#include <limits.h>
#include <omp.h>

// Enable cpp11
// [[Rcpp::plugins("cpp11")]]

// finds the chunk of rows this thread will process
void findmyrange(int n, int nth, int me, int *myrange)
{  int chunksize = n / nth;
   myrange[0] = me * chunksize;
   if (me < nth-1) myrange[1] = (me+1) * chunksize - 1;
   else myrange[1] = n - 1;
}

//create square matrix and fill with NA's
Rcpp::NumericMatrix na_matrix(int n){
  Rcpp::NumericMatrix m(n,n) ;
  std::fill( m.begin(), m.end(), Rcpp::NumericVector::get_na() ) ;
  return m;
}

void calSecondOrder(double pij, double piq, double pik,
		    double pjq, double pjk, double pqk,
		    double &pijk, double &piqk, double &pjqk,
		    double &pijq, double &pikq, double &pjkq ){

  double num, dem, pkq;

  //symmetric
  pkq = pqk;
  
  pijk = std::numeric_limits<double>::quiet_NaN();
  piqk = std::numeric_limits<double>::quiet_NaN();
  pjqk = std::numeric_limits<double>::quiet_NaN();

  pijq = std::numeric_limits<double>::quiet_NaN();
  pikq = std::numeric_limits<double>::quiet_NaN();
  pjkq = std::numeric_limits<double>::quiet_NaN();
  
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
RcppExport SEXP secondOrder(SEXP GG, SEXP CC, SEXP N )
{

  int n;
  if( !Rf_isNull(GG) && !Rf_isNull(CC) ){    
    
    // make a C/C++ compatible view of the R matrix;
    // note: referencing is done with ( , ) not [ , ],
    // and indices start at 0
    Rcpp::NumericMatrix rGG(GG);
    Rcpp::NumericMatrix rCC(CC);

    //Problem here... probably need to check if null.
    //Rcpp::CharacterVector GN = rownames(rCC);

    //N == number of samples.
    double samples = std::numeric_limits<double>::quiet_NaN();
    if( !Rf_isNull(N) ){
      Rcpp::NumericVector rN(N);
      samples = rN[0];
    }    
    
    n = rCC.nrow();

    // create output matrices
    Rcpp::NumericMatrix rRES = na_matrix( n );
    Rcpp::NumericMatrix rPR  = na_matrix( n );
    Rcpp::NumericMatrix rSD  = na_matrix( n );

    #pragma omp parallel
    { int i,j,k,q;

      // each thread has an ID (starting at 0), so
      // determine the ID for this thread
      int me = omp_get_thread_num();
      // find total number of threads
      int nth = omp_get_num_threads();
      int myrows[2];
      // finds the chunk of rows this thread will process
      findmyrange(n,nth,me,myrows);
      
      //perform first order calculation with thread
      for( i=myrows[0]; i<=myrows[1]; i++ ){
	for( j=i; j<n; j++ ){

	  if( i != j ){

	    if( !Rcpp::NumericVector::is_na(rGG(i,j)) ){

	      double Rtemp = std::numeric_limits<double>::quiet_NaN();
	      double Rso   = std::numeric_limits<double>::quiet_NaN();
	      
	      for( k=0; k<n; k++ ){
		for( q=0; q<n; q++ ){
		
		  if( (i != k) && (j != k) &&
		      (i != q) && (j != q) && (q != k) ){

		    if( !Rcpp::NumericVector::is_na(rGG(i,k)) && 
			!Rcpp::NumericVector::is_na(rGG(j,k)) &&
			!Rcpp::NumericVector::is_na(rGG(i,q)) &&
			!Rcpp::NumericVector::is_na(rGG(j,q)) &&
			!Rcpp::NumericVector::is_na(rGG(q,k)) &&
			!Rcpp::NumericVector::is_na(rCC(i,j)) &&
			!Rcpp::NumericVector::is_na(rCC(i,q)) &&
			!Rcpp::NumericVector::is_na(rCC(i,k)) &&
			!Rcpp::NumericVector::is_na(rCC(j,q)) &&
			!Rcpp::NumericVector::is_na(rCC(j,k)) &&
			!Rcpp::NumericVector::is_na(rCC(q,k)) ){

		      double pijk, piqk, pjqk, pijq, pikq, pjkq;
		    
		      calSecondOrder( rCC(i,j), rCC(i,q), rCC(i,k),
				      rCC(j,q), rCC(j,k), rCC(q,k),
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
			    Rtemp = 0.5*( (num1abs/dem1) + (num2abs/dem2) );
			    Rso   = num1 / dem1;
			  } else {
			    double r = 0.5*( (num1abs/dem1) + (num2abs/dem2));
			    if( r < Rtemp ){ Rtemp = r; Rso = (num1/dem1); }
			  }			
			}
		      }
		    }		
		  }//i!=k & j!=k
		}//q
	      }//k

	      rRES(i,j) = Rso;
	      rRES(j,i) = Rso;

	      rPR(i,j)  = Rtemp;
	      rPR(j,i)  = Rtemp;

	      if( !isnan(samples) ){
		double norm = samples - 4;
		if( norm >= 0 ){
		  rSD(i,j) = (1-Rso*Rso) / (sqrt(norm));
		  rSD(j,i) = rSD(i,j);
		}
	      }
	      
	      
	    }
	  
	  }//i!=j
	
	}//j
      }//i

    }//omp
      
    Rcpp::List OUT = Rcpp::List::create(Rcpp::Named("PC") =rRES,
					Rcpp::Named("MPC")=rPR,
					Rcpp::Named("SD") =rSD);

    return OUT;

  }

  Rcpp::List OUT = Rcpp::List::create(Rcpp::Named("PC") =na_matrix(1),
				      Rcpp::Named("MPC")=na_matrix(1),
				      Rcpp::Named("SD") =na_matrix(1));

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
RcppExport SEXP firstOrder(SEXP GG, SEXP CC, SEXP N )
{

  int n;
  if( !Rf_isNull(GG) && !Rf_isNull(CC) ){    
    
    // make a C/C++ compatible view of the R matrix;
    // note: referencing is done with ( , ) not [ , ],
    // and indices start at 0
    Rcpp::NumericMatrix rGG(GG);
    Rcpp::NumericMatrix rCC(CC);

    //Problem here... probably need to check if null.
    //Rcpp::CharacterVector GN = rownames(rCC);

    //N == number of samples.
    double samples = std::numeric_limits<double>::quiet_NaN();
    if( !Rf_isNull(N) ){
      Rcpp::NumericVector rN(N);
      samples = rN[0];
    }
    
    n = rCC.nrow();

    // create output matrices
    Rcpp::NumericMatrix rRES = na_matrix( n ); //first order partial correlation
    Rcpp::NumericMatrix rPR  = na_matrix( n ); //first order modified partial correlation
    Rcpp::NumericMatrix rSD  = na_matrix( n ); // standard error on first order partial correlation

    #pragma omp parallel
    { int i,j,k;

      // each thread has an ID (starting at 0), so
      // determine the ID for this thread
      int me = omp_get_thread_num();
      // find total number of threads
      int nth = omp_get_num_threads();
      int myrows[2];
      // finds the chunk of rows this thread will process
      findmyrange(n,nth,me,myrows);
      
      //perform first order calculation with thread
      for( i=myrows[0]; i<=myrows[1]; i++ ){
	for( j=i; j<n; j++ ){

	  if( i != j ){

	    if( !Rcpp::NumericVector::is_na(rGG(i,j)) ){

	      double Rtemp = std::numeric_limits<double>::quiet_NaN();
	      double Rfo   = std::numeric_limits<double>::quiet_NaN();
	      
	      for( k=0; k<n; k++ ){

		if( (i != k) && (j != k) ){

		  if( !Rcpp::NumericVector::is_na(rGG(i,k)) && 
		      !Rcpp::NumericVector::is_na(rGG(j,k)) ){

		    double num    = rCC(i,j) - (rCC(i,k)*rCC(j,k));
		    double numABS = abs(rCC(i,j)) - ( abs(rCC(i,k))*abs(rCC(j,k)) );
		    double dem    = sqrt( (1-rCC(i,k)*rCC(i,k)) ) *
		      sqrt( (1-rCC(j,k)*rCC(j,k)) );

		    if( (dem != 0) && (!isnan(dem)) ){

		      if( isnan(Rfo) ){
			Rtemp = numABS / dem;
			Rfo   = num / dem;
		      } else {
			double r = numABS / dem;
			if( r < Rtemp ){ Rtemp = r; Rfo = (num/dem); }
		      }
		      
		    }
		  
		  }
		
		}//i!=k & j!=k
	      
	      }//k

	      rRES(i,j) = Rfo;
	      rRES(j,i) = Rfo;

	      rPR(i,j)  = Rtemp;
	      rPR(j,i)  = Rtemp;

	      if( !isnan(samples) ){
		double norm = samples - 3;
		if( norm >= 0 ){
		  rSD(i,j) = (1-Rfo*Rfo) / (sqrt(norm));
		  rSD(j,i) = rSD(i,j);
		}
	      }
	      
	  }
	  
	  }//i!=j
	
	}//j
      }//i

    }//omp
      
    Rcpp::List OUT = Rcpp::List::create(Rcpp::Named("PC") =rRES,
					Rcpp::Named("MPC")=rPR,
					Rcpp::Named("SD") =rSD);

    return OUT;

  }

  Rcpp::List OUT = Rcpp::List::create(Rcpp::Named("PC") =na_matrix(1),
				      Rcpp::Named("MPC")=na_matrix(1),
				      Rcpp::Named("SD") =na_matrix(1));

  return OUT;

}
    
 
