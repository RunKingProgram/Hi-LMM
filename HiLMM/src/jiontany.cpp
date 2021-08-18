// [[Rcpp::depends(RcppArmadillo)]]
#include <fstream>
#include <RcppArmadillo.h>
#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>

using namespace Rcpp;
using namespace std;
using namespace arma;

// [[Rcpp::export]]
arma::mat Cposi_Choice(arma::mat b, float ws){
  if(b.n_rows<3){return b;}
  int nmar = b.n_rows;
  int pos,pos2,indexpos;
  arma::vec b1(nmar);
  arma::vec p1(nmar);
  arma::mat result(nmar,2);
  result.zeros(nmar,2);
  b1=b.col(0);
  p1=b.col(1);
  int count =0;
  for(int i=0;i<nmar-1;){
      pos = b1(i);
      int j=i+1;
      pos2=b1(j);
      if(pos2-pos>=ws){
        i = j;
        continue;
      }
      while(pos2-pos<ws & j<nmar ){
        pos2=b1(j);
        j++;
        }
      indexpos = p1.subvec(i,j-1).index_min();
      result(count,0) = b1(indexpos+i);
      result(count,1) = p1.subvec(i,j-1).min();
      count++;
      i = j;
  }
  return(result.rows(0,count-1));
}
arma::colvec fastLm(const arma::mat& X, const arma::colvec& y) {
    int n = X.n_rows, k = X.n_cols;    
    arma::colvec coef = arma::solve(X, y);    // fit model y ~ X
    arma::colvec res  = y - X*coef; // residuals  
    // std.errors of coefficients
    double s2 = std::inner_product(res.begin(), res.end(), res.begin(), 0.0)/(n - k);                                                     
    arma::colvec err = s2 * arma::diagvec(arma::pinv(arma::trans(X)*X));  
    
    return (coef%coef)/err;
}
int FindMin(arma::colvec a, int n, int *pMinPos)
{
    int i, min;
    min = a(1);              
    *pMinPos = 0;           
     
    for (i=1; i<n; i++)
    {
        if (a(i) < min)
        {
            min = a[i];
            *pMinPos = i;  
        }
    }
    return min ;
}


 // [[Rcpp::export]]
Rcpp::List fast_jiont(arma::mat y,arma::mat b0,arma::mat position,float thrd, SEXP bedfile_in) {
   string bedfile = Rcpp::as<string>(bedfile_in);
   size_t n = y.size();
   size_t nmar = position.size();
   size_t ncount;
   arma::mat ty(n, 1);
   char *cp;
   int geno;
   float minValue;
   ifstream readbedfile (bedfile.c_str(), ios::binary);
   int nblocks = (n + 3) / 4, pos;
   unsigned char magic[3], temp[2], buffer[nblocks];
   readbedfile.read((char *)magic, 3);
   arma::mat cholK(n, n);
   
    FILE *fp;
    fp=fopen("./output/cholK.txt","r");
    if(fp==NULL){
      printf("No cholK file!");
    }   
    char* line = NULL;
    size_t len = 0;   
    for(int count = 0;count<n;count++) { 
      getline(&line, &len, fp);
      char *ptr = NULL;
      ptr = strtok(line, "\t");
      for(int i = 0;i < n;i++) {
        cholK(count,i) = atof(ptr);
        ptr = strtok(NULL, "\t");
      }
      ty(count,0)=dot(cholK.row(count),y);
    }
   arma::mat Gall(n, nmar+1);
   Gall.zeros(n,nmar+1);
   for(size_t i=0; i<nmar; ++i) {
     stringstream writeout;
     readbedfile.seekg(position(i)*nblocks+3);
     ncount = 0;
     readbedfile.read((char *)buffer, nblocks);
     for(size_t j=0; j<nblocks; ++j) {
       pos = 0;
       for(size_t k=0; k<4; ++k) {
         geno=1;
         if(ncount == n && j == nblocks - 1) {break;}
         for(size_t l=0; l<2; ++l) {
           temp[l] = (buffer[j] >> pos) & 1;
           pos++;
         }
         if(temp[0] ==0 && temp[1] ==0){
           geno = 2;
         } else if (temp[0] ==1 && temp[1] ==1){
           geno = 0;
         } 
         if (i==0){
          Gall(ncount,0) = b0(ncount,0);
          }
         if(geno != 0){
            Gall.col(i+1) += cholK.col(ncount)*geno;
         }
         ncount++;
       }
     }
     }
    
     arma::colvec chi;
     chi = fastLm(Gall, y);
     int nmarnow = nmar;
     for (int i = 0; i < nmar; i++){
      int minPos;
      minValue = FindMin(chi, nmarnow-i, &minPos); 
      if(minValue>thrd){
      return Rcpp::List::create(Rcpp::Named("chi") = chi,
                            Rcpp::Named("position") = position);
      }
      Gall.shed_col(minPos);
      position.shed_row(minPos);
      chi = fastLm(Gall, y);
     }
   return Rcpp::List::create(Rcpp::Named("chi") = chi,
                            Rcpp::Named("position") = position);
 }






#include <Rcpp.h>
// fast_jiont
Rcpp::List fast_jiont(arma::mat y, arma::mat b0, arma::mat position, float thrd, SEXP bedfile_in);
RcppExport SEXP sourceCpp_19_fast_jiont(SEXP ySEXP, SEXP b0SEXP, SEXP positionSEXP, SEXP thrdSEXP, SEXP bedfile_inSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type b0(b0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type position(positionSEXP);
    Rcpp::traits::input_parameter< float >::type thrd(thrdSEXP);
    Rcpp::traits::input_parameter< SEXP >::type bedfile_in(bedfile_inSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_jiont(y, b0, position, thrd, bedfile_in));
    return rcpp_result_gen;
END_RCPP
}
