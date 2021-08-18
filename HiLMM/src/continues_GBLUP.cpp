// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <fstream>


using namespace std;

// [[Rcpp::export]]
Eigen::VectorXf gblup_ebv(Eigen::VectorXf y, const float h2) {
  const int nobs = y.size();
  Eigen::MatrixXf X(nobs+1,nobs+1);
  Eigen::VectorXf Y(nobs+1);
  Eigen::VectorXf ebv(nobs+1);
  Eigen::VectorXf tempv(nobs); 
  X.setOnes(nobs+1,nobs+1);
  FILE *fp;
  fp=fopen("./output/kinship.sXX.txt","r");
  if(fp==NULL){
    printf("No eigenU file!");
    return y;
  }   
  char* line = NULL;
  size_t len = 0;
  for(int count=0;count<nobs;count++) { 
      getline(&line, &len, fp);
      char *ptr=NULL;
      ptr = strtok(line, "\t");
      for(int i=0;i<nobs;i++) {
         tempv(i)=atof(ptr);
          if(i == count){ tempv(i)+=0.001;}
          ptr = strtok(NULL, "\t");
      }
      X.row(count)=tempv;
      X(count,nobs)=tempv.sum();
      Y(count)=y.dot(tempv);
      X(count,count)=X(count,count)+(1-h2)/h2;
  }
  free(line);
  Y(nobs)=y.sum();
  X(nobs,nobs)=nobs;
  ebv = X.lu().solve(Y);
  return ebv.head(nobs);
}

