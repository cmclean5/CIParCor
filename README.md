# CIParCor
Conditional Independent first- and second-order Partial Correlation 

Subroutines (in R) to calculate the first and second order conditional independent partial correlation. 

Subroutines are given as inputed an NxN zeorth order correlation matrix called "C", possible NxN 'filter' matrix called "G", and possible number of samples called "N". 

Subroutines output the NxN first or second order partial correlation matrix called "PC", the corresponding 
NxN conditional independent test score called "MPC", and NxN standard deviation for "PC" given "N".


"first order conditional independent partial correlation"


The first order partial correlation is what we get when we hold constant some third variable from two other variables. For example,  we may have a set of three variables [i,j,k] and know the zeroth order correlation between i and j is .88. But k accounts for (or could account for) part of that. What would happen to the correlation if k were constant? In other words, 
we partial one variable out of a correlation.

References:

[1] Adapted from "OpenMP Tutorial, with R Interface", Matloff, 16/01/2015.
https://matloff.wordpress.com/2015/01/16/openmp-tutorial-with-r-interface/

[2] Magwene, P. M. & Kim, J. "Estimating genomic coexpression networks using first-order conditional independance", Genome Biology, 5:R100, (2004).

[3] Zuo, Y. et. al. "Biological network inference using low order partial correlation", Methods, 69(3):266-273, (2014). doi:10.1016/j.ymeth.2014.06.010


Build instructions:

export R_LIBS_USER=/users/cmclean/R/x86_64-pc-linux-gnu-library/3.6/ # your R libraries

export PKG_LIBS="-lgomp"

export PKG_CXXFLAGS="-fopenmp -I/users/cmclean/R/x86_64-pc-linux-gnu-library/3.6/Rcpp/include"

R CMD SHLIB CIParCor.cpp

* to find the number of all installed cores/processors in linux: 

nproc --all

* before start R, set number of threads, e.g. 

export OMP_NUM_THREADS=2


Test in R:

> library(Rcpp)

> dyn.load("/users/cmclean/STUDIES/RCcpp/v4r1/CIParCor.so")

> GG = matrix(1,ncol=2,nrow=2)

> CC = GG

> .Call("firstOrder",GG,CC,NULL)

> dyn.unload("/users/cmclean/STUDIES/RCcpp/v4r1/CIParCorcpp.so")

