# CIParCor
Conditional Independent first- and second-order Partial Correlation 

Subroutines (in R) to calculate the first and second order conditional independent partial correlation. 

Subroutines are given as inputed an NxN zeorth order correlation matrix called "C", possible NxN 'filter' matrix called "G", and possible number of samples called "N". 

Subroutines output the NxN first or second order partial correlation matrix called "PC", the corresponding 
NxN conditional independent test score called "MPC", and NxN standard deviation for "PC" given "N".

Subrountine to calculate the pair-wise partial correlation, given the full set, from the covariance matrix using the Moore-Penrose matrix inversion method.

"first order conditional independent partial correlation"


The first order partial correlation is what we get when we hold constant some third variable from two other variables. For example,  we may have a set of three variables [i,j,k] and know the zeroth order correlation between i and j is .88. But k accounts for (or could account for) part of that. What would happen to the correlation if k were constant? In other words, 
we partial one variable out of a correlation.

### Test

$ Rscript test.R 200 2 2

<table class="tg">
  <tr>
    <th class="tg-yw4l"><b>NxN</b></th>
    <th class="tg-yw4l"><b>CPUs</b></th>
    <th class="tg-yw4l"><b>firstOrder (s)</b></th>
     <th class="tg-yw41"><b>secondOrder (s)</b></th>
  </tr>
  <tr>
    <td class="tg-yw4l">200x200</td>
    <td class="tg-yw4l">1</td>
    <td class="tg-yw4l">0.37</td>
     <td class="tg-yw41">205.3</td>
</tr>
  <tr>
    <td class="tg-yw4l">200x200</td>
    <td class="tg-yw4l">2</td>
    <td class="tg-yw4l">0.32</td>
     <td class="tg-yw41">103.1</td>
</tr>
  <tr>
    <td class="tg-yw4l">200x200</td>
    <td class="tg-yw4l">4</td>
    <td class="tg-yw4l">0.30</td>
     <td class="tg-yw41">52.1</td>
</tr>
  <tr>
    <td class="tg-yw4l">200x200</td>
    <td class="tg-yw4l">8</td>
    <td class="tg-yw4l">0.32</td>
     <td class="tg-yw41">26.0</td>
</tr>
</table>

---

References:

[1] https://matloff.wordpress.com/2015/01/16/openmp-tutorial-with-r-interface/

[2] Magwene, P. M. & Kim, J. "Estimating genomic coexpression networks using first-order conditional independance", Genome Biology, 5:R100, (2004).

[3] Zuo, Y. et. al. "Biological network inference using low order partial correlation", Methods, 69(3):266-273, (2014). doi:10.1016/j.ymeth.2014.06.010

[4] http://www.parallelr.com/r-and-openmp-boosting-compiled-code-on-multi-core-cpu-s/

[5] https://bisqwit.iki.fi/story/howto/openmp/

---

Build instructions:

export R_LIBS_USER=/users/cmclean/R/x86_64-pc-linux-gnu-library/3.6/ # your R libraries

export PKG_LIBS="-lgomp"

export PKG_CXXFLAGS="-fopenmp -I$R_LIBS_USER/Rcpp/include -I$R_LIBS_USER/RcppEigen/include"

R CMD SHLIB CIParCor.cpp geninv.cpp

OR:

PKG_CXXFLAGS="$(echo 'Rcpp:::CxxFlags()'| R --vanilla --slave) -fopenmp" R CMD SHLIB CIParCor.cpp

* to find the number of all installed cores/processors in linux: 

nproc --all

* Test in R:

---

