// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <fstream>
using namespace std;
using Eigen::SelfAdjointEigenSolver;
float cal_nobslogVe(Eigen::VectorXf x){
    int nobs = x.size();
    float Ve = x.array().square().sum()/(nobs-2);
    return nobs*log(Ve);
}
// [[Rcpp::export]]
float Loglike(Eigen::VectorXf sg,Eigen::VectorXf ynew,
              Eigen::VectorXf gnew,float h2){
    int nobs =ynew.size();
    Eigen::VectorXf V050(nobs);
    V050 =sg.array()*h2/(1-h2)+1;
    float logV0=V050.array().log().sum();
    V050=V050.array().sqrt().inverse();
    ynew=ynew.array()*V050.array();
    gnew=gnew.array()*V050.array();
    float coef=gnew.dot(ynew)/gnew.dot(gnew);
    float L=logV0+cal_nobslogVe(ynew.array()-gnew.array()*coef);
    return L;
}
// [[Rcpp::export]]
Rcpp::List tugxy (Eigen::VectorXf y,Eigen::VectorXf tugx) {
  const int nobs = y.size();
  Eigen::VectorXf tempv(nobs);
  Eigen::MatrixXf G(nobs,nobs);
  Eigen::MatrixXf H(nobs,nobs);
  Eigen::MatrixXf ug(nobs,nobs);
  Eigen::VectorXf tugy(nobs);
  Eigen::VectorXf gnewb0(nobs);
  Eigen::VectorXf sg(nobs);
  Eigen::VectorXf sginv(nobs);
  sg.setOnes(nobs);
  gnewb0.setOnes(nobs);
  FILE *fp;
  fp=fopen("./output/kinship.sXX.txt","r");
  if(fp==NULL){
    printf("No eigenU file!");
    return 0;
  }
  char* line = NULL;
  size_t len = 0;

  for(int count = 0;count<nobs;count++) {
    getline(&line, &len, fp);
    char *ptr = NULL;
    ptr = strtok(line, "\t");
    for(int i = 0;i < nobs;i++) {
      tempv(i) = atof(ptr);
      ptr = strtok(NULL, "\t");
    }
    G.row(count) = tempv;
    G(count,count) += 0.0001;
  }
  cout<<"Starting spectral decomposition"<<endl;
  SelfAdjointEigenSolver<Eigen::MatrixXf> eig(G);
  cout<<"Spectral decomposition finished"<<endl;
  for(int count = 0;count<nobs;count++) {
    sg(count)=eig.eigenvalues()(nobs-count-1);
    ug.col(count)=(eig.eigenvectors()).col(nobs-count-1);
  }
  sginv = sg.array().sqrt();
  sginv =1/sginv.array();
  tugx = ug.transpose()*tugx;
  tugy = ug.transpose()*y;
  for(int count = 0;count<nobs;count++) {
    H.col(count) = sginv.array() *
      (ug.transpose()).col(count).array();
  }
  gnewb0 = H*gnewb0;
  fclose(fp);
  fp = fopen("./output/cholK.txt","w");
  for(int i = 0;i < nobs;i++){
    for(int j= 0;j < nobs;j++){
      fprintf(fp,"%f",H(i,j));

       if(j < nobs-1){fprintf(fp,"%s","\t");}
      }
    fprintf(fp,"%s","\n");
  }
  fclose(fp);

  return Rcpp::List::create(Rcpp::Named("tugx") = tugx,
                            Rcpp::Named("tugy") = tugy,
                            Rcpp::Named("sg") = sg,
                            Rcpp::Named("gnewb0") = gnewb0);
}

// [[Rcpp::export]]
Rcpp::List tugxy2 (Eigen::VectorXf y) {
  const int nobs = y.size();
  Eigen::VectorXf tempv(nobs);
  Eigen::MatrixXf G(nobs,nobs);
  Eigen::MatrixXf H(nobs,nobs);
  Eigen::MatrixXf ug(nobs,nobs);
  Eigen::VectorXf tugy(nobs);
  Eigen::VectorXf gnewb0(nobs);
  Eigen::VectorXf sg(nobs);
  Eigen::VectorXf sginv(nobs);
  Eigen::VectorXf tugx(nobs);
  sg.setOnes(nobs);
  gnewb0.setOnes(nobs);
  FILE *fp;
  fp=fopen("./output/kinship.sXX.txt","r");
  if(fp==NULL){
    printf("No eigenU file!");
    return 0;
  }
  char* line = NULL;
  size_t len = 0;

  for(int count = 0;count<nobs;count++) {
    getline(&line, &len, fp);
    char *ptr = NULL;
    ptr = strtok(line, "\t");
    for(int i = 0;i < nobs;i++) {
      tempv(i) = atof(ptr);
      ptr = strtok(NULL, "\t");
    }
    G.row(count) = tempv;
    G(count,count) += 0.0001;
  }
  cout<<"Starting spectral decomposition"<<endl;
  SelfAdjointEigenSolver<Eigen::MatrixXf> eig(G);
  cout<<"Spectral decomposition finished"<<endl;
  for(int count = 0;count<nobs;count++) {
    sg(count) = eig.eigenvalues()(nobs-count-1);
    ug.row(count) = (eig.eigenvectors()).col(nobs-count-1);
    tugx(count) = ug.row(count).sum();
    tugy(count) = y.dot(ug.row(count));
  }
  sginv = sg.array().sqrt();
  sginv =1 / sginv.array();

  for(int count = 0;count<nobs;count++) {
    H.col(count) = sginv.array() *
      ug.col(count).array();
  }
  gnewb0 = H * gnewb0;
  fclose(fp);
  fp = fopen("./output/cholK.txt","w");
  for(int i = 0;i < nobs;i++){
    for(int j= 0;j < nobs;j++){
      fprintf(fp,"%f",H(i,j));
      if(j < nobs-1){fprintf(fp,"%s","\t");}
    }
    fprintf(fp,"%s","\n");
  }
  fclose(fp);
  fp = fopen("./output/tug.txt","w");
  for(int i = 0;i < nobs;i++){
    for(int j= 0;j < nobs;j++){
      fprintf(fp,"%f",ug(i,j));
      if(j < nobs-1){fprintf(fp,"%s","\t");}
    }
    fprintf(fp,"%s","\n");
  }
  fclose(fp);
  return Rcpp::List::create(Rcpp::Named("tugx") = tugx,
                            Rcpp::Named("sg") = sg,
                            Rcpp::Named("gnewb0") = gnewb0);
}
