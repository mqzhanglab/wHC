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
NumericVector get_2x2xk(IntegerVector dis, IntegerVector cur_g, IntegerVector strata_int, int k) {
  int n = dis.size();
  NumericMatrix mat_2x2xk(4, k);
  for (int i = 0; i != n; i++) {
    if (dis[i] == 0) {
      if (cur_g[i] == 0) {
        ++(mat_2x2xk(0, strata_int[i]));
      } else {
        ++(mat_2x2xk(2, strata_int[i]));
      }
    } else {
      if (cur_g[i] == 0) {
        ++(mat_2x2xk(1, strata_int[i]));
      } else {
        ++(mat_2x2xk(3, strata_int[i]));
      }
    }
  }
  IntegerVector dim(3);
  dim[0] = 2;
  dim[1] = 2;
  dim[2] = k;
  mat_2x2xk.attr("dim") = dim;
  return mat_2x2xk;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
g=  c(1,0,0,1,0,0,1)
dis=c(1,0,1,1,1,0,0)
strata=c(0,0,0,0,0,0,0)
k=1
get_2x2xk(dis,g,strata,k)
*/
