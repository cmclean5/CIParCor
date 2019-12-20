# CIParCor
Conditional Independent first- and second-order Partial Correlation 

Subroutines (in R) to calculate the first and second order conditional independent partial correlation. 

Subroutines are given as inputed an NxN zeorth order correlation matrix called "C", possible NxN 'filter' matrix called "G", and possible number of samples called "N". 

Subroutines output the NxN first or second order partial correlation matrix called "PC", the corresponding 
NxN conditional independent test score called "MPC", and NxN standard deviation for "PC" given "N".


"first order conditional independent partial correlation"


The first order partial correlation is what we get when we hold constant some third variable from two other variables. For example,  we may have a set of three variables [i,j,k] and know the zeroth order correlation between i and j is .88. But k accounts for (or could account for) part of that. What would happen to the correlation if k were constant? In other words, 
we partial one variable out of a correlation.

---

### Test

<table class="tg">
  <tr>
    <th class="tg-yw4l"><b>NxN</b></th>
    <th class="tg-yw4l"><b>CPUs</b></th>
    <th class="tg-yw4l"><b>firstOrder (s)</b></th>
     <th class="tg-yw41"><b>secondOrder (s)</b></th>
  </tr>
  <tr>
    <td class="tg-yw4l">100x100</td>
    <td class="tg-yw4l">1</td>
    <td class="tg-yw4l">0.23</td>
     <td class="tg-yw41">23.7</td>
</tr>
  <tr>
    <td class="tg-yw4l">100x100</td>
    <td class="tg-yw4l">2</td>
    <td class="tg-yw4l">0.17</td>
     <td class="tg-yw41">4.58</td>
</tr>
  <tr>
    <td class="tg-yw4l">100x100</td>
    <td class="tg-yw4l">4</td>
    <td class="tg-yw4l">0.2</td>
     <td class="tg-yw41">5.9</td>
</tr>
  <tr>
    <td class="tg-yw4l">100x100</td>
    <td class="tg-yw4l">8</td>
    <td class="tg-yw4l">0.15</td>
     <td class="tg-yw41">2.95</td>
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

export PKG_CXXFLAGS="-fopenmp -I/users/cmclean/R/x86_64-pc-linux-gnu-library/3.6/Rcpp/include"

R CMD SHLIB CIParCor.cpp

OR:

PKG_CXXFLAGS="$(echo 'Rcpp:::CxxFlags()'| R --vanilla --slave) -fopenmp" R CMD SHLIB CIParCor.cpp

* to find the number of all installed cores/processors in linux: 

nproc --all

---

* Test in R:

$ Rscript test.R 100 0\

args...\
N    = 100\
CPUs = 0\

run firstOrder...\ 
OpenMP:  number of threads 0\

user  system elapsed\    
0.208   0.000   0.208\   
...done.

> run secondOrder...\ 
> OpenMP:  number of threads 0\

>   user  system elapsed\    
> 23.742   0.000  23.762\ 
> ...done.\

> print mean.

>         TT       fO.PC       sO.PC 
> 0.001865987 0.014508042 0.161405167 

> print sd.

>       TT     fO.PC     sO.PC 
> 0.2471588 0.5556790 4.3702410 


$ Rscript test.R 100 4

> args...

> N    = 100 

> CPUs = 4 

> run firstOrder... 

> OpenMP:  number of threads 4

>   user  system elapsed 
  
>  0.204   0.007   0.171 
  
> ...done.

> run secondOrder... 

> OpenMP:  number of threads 4

>   user  system elapsed 
 
> 23.605   0.000   5.912 
 
> ...done.

> print mean.

>          TT        fO.PC        sO.PC 
> 0.0001579991 0.0096378304 0.2208161727 

> print sd.

>       TT     fO.PC     sO.PC 
> 0.2493581 0.7044041 7.5740866 
