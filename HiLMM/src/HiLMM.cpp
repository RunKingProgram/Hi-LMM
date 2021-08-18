// [[Rcpp::depends(RcppArmadillo)]]
#include <fstream>
#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace std;
using namespace arma;


class Params
{
public:
   arma::colvec beta    ;
   arma::colvec err ;
};

Params sfastLm(const arma::mat& X, const arma::colvec& y) {
   int n = X.n_rows, k = X.n_cols;
   arma::colvec coef = arma::solve(X, y);    // fit model y ~ X
   arma::colvec res  = y - X * coef;           // residuals
   // std.errors of coefficients
   double s2 = std::inner_product(res.begin(), res.end(), res.begin(), 0.0)/(n - k);
   arma::colvec err = s2 * arma::diagvec(arma::pinv(arma::trans(X)*X));
   Params results;
   results.beta = coef;
   results.err = err;
   return results;
}
 // [[Rcpp::export]]
void Hilmm_multithreads(arma::mat y,arma::mat b0,SEXP bedfile_in,
 SEXP outfile_in,int thread,int nmar,int threadnumb) {
   string bedfile = Rcpp::as<string>(bedfile_in);
   string outfile = Rcpp::as<string>(outfile_in);
   size_t n = y.size();
   arma::mat ty(n, 1);
   size_t ncount;
   string *tmpout = new string[n];
   char *cp;
   int st,end,range;
   ifstream readbedfile (bedfile.c_str(), ios::binary);
   ofstream writefile (outfile.c_str(), ofstream::out);
   int nblocks = (n + 3) / 4, pos;
   unsigned char magic[3], temp[2], buffer[nblocks];
   readbedfile.read((char *)magic, 3);
   arma::mat  cholK(n, n);
   arma::mat  cholK2(n,n);
    FILE *fp;
    fp=fopen("./output/cholK.txt","r");
    if(fp==NULL){
      printf("No cholK file!");
    }
    char* line = NULL;
    size_t len = 0;
    for(int count = 0;count < n;count++) {
      getline(&line, &len, fp);
      char *ptr = NULL;
      ptr = strtok(line, "\t");
      for(int i = 0;i < n;i++) {
        cholK(count,i) = atof(ptr);
        cholK2(count,i)= 2 * cholK(count,i);
        ptr = strtok(NULL, "\t");
      }
      ty(count,0) = dot(cholK.row(count),y);
    }
   arma::mat G(n, 2);
   range =nmar / thread;
   st = range * ( threadnumb - 1 );
   end = range * (threadnumb);
   if(threadnumb == thread) {end = nmar;}
   cout << "start=" << st << endl;
   cout << "end=" << end << endl;

   for(size_t i=st; i<end; ++i) {
     G.zeros(n,2);
     G.col(1) = b0.col(0);
     stringstream writeout;
     readbedfile.seekg(i * nblocks+3);
     ncount = 0;
     readbedfile.read((char *)buffer, nblocks);
     for(size_t j = 0; j < nblocks; ++j) {
       pos = 0;
       for(size_t k = 0; k < 4; ++k) {
         if(ncount == n && j == nblocks - 1) {break;}
         for(size_t l = 0; l < 2; ++l) {
           temp[l] = (buffer[j] >> pos) & 1;
           pos++;
         }
         if(temp[0] == 0 && temp[1] == 0){
           G.col(0) += cholK2.col(ncount);
           ncount++;
           continue;
         } else if (temp[0] == 1 && temp[1] == 1){
           ncount++;
           continue;
         }
         G.col(0) += cholK.col(ncount);
         ncount++;
       }
     }
    Params rgres;
    rgres = sfastLm(G, ty);
    writefile << rgres.beta(0) << "\t"<<rgres.beta(0)*rgres.beta(0)/rgres.err(0) << "\n";
   }

   return;
 }

