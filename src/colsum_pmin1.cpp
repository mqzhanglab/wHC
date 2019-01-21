#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericMatrix colsum_pmin1(NumericMatrix pm, NumericMatrix wm) {
  int l = pm.nrow();
  double invl = 1.0 / l;
  int m = pm.ncol();
  int n = wm.ncol();
  NumericMatrix res(l, m);
  for (int j1 = 0; j1 != m; j1++) {
    for (int i = 0; i != l; i++) {
      for (int j2 = 0; j2 != n; j2++) {
        res(i, j1) += std::min(1.0, pm(i, j1) * wm(i, j2)) * invl;
      }
    }
  }
  return res;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
mat1 = matrix(1:9, nrow=3) / 20
wm  = matrix(11:19, nrow=3) / 4
l = nrow(mat1)
fun = function(y) {apply(y,2,function(x){return(rowSums(pmin(x*wm,1))/l)})}
fun(mat1)
colsum_pmin1(mat1, wm)
microbenchmark::microbenchmark(fun(mat1), colsum_pmin1(mat1, wm))
*/
