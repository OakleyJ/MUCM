#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;

using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]

//' @param delta [q] is a vector the log of the correlation length
//' @param matX [p x q] is the input data matrix with p observations and q variables
//' @param matH [q x ?] is the H matrix, the number of columns depends on the prior model for the mean function
//' @param matY [p x r] is the output data matrix with p observations and r variables

//' @export
// [[Rcpp::export]]
double GetNllDelta (const NumericVector & delta,
                  const NumericMatrix & matX,
                  const Eigen::Map<Eigen::MatrixXd> & matY,
                  const Eigen::Map<Eigen::MatrixXd> & matH){
  

  
  unsigned int p = matX.nrow();
  unsigned int q = matH.cols();
  unsigned int r = matY.cols();
  
  const NumericVector psi2 = exp(delta)*exp(delta);
  const double NLL_CAP = 1000000.0;
  
  // Get A
  NumericMatrix matA(p, p);
  std::fill(matA.begin(), matA.end(), 1);
  unsigned int i=0, j=0;
  double aij; 
  
  for (i = 0; i < p - 1; i++){
    NumericVector vi = matX.row(i);
    for (j = i + 1; j < p ; j ++){
      aij = exp(-sum((vi-matX.row(j))*(vi-matX.row(j))/psi2));
      matA(j,i)=aij;
      matA(i,j)=aij;
    }
  }
  
  // Get matL := chol(matA)
  Map<MatrixXd> eigenMatA = as<Map<MatrixXd> >(matA);
  LLT<MatrixXd> lltA(eigenMatA);
  if(lltA.info() == NumericalIssue){
    // std::cout<<"Non-PSD matrix A!"<<std::endl;
    return(NLL_CAP);
  }
  MatrixXd matL = lltA.matrixL();
  
  // Get matQ := t(matH) %*% solve(matA) %*% matH
  MatrixXd mat1 = lltA.matrixL().solve(matH);
  MatrixXd matQ = MatrixXd(q,q).setZero().selfadjointView<Lower>()
                               .rankUpdate(mat1.adjoint());
  
  // Get matW := solve(A) %*% matY
  MatrixXd matW = lltA.solve(matY);
  
  // Get matBeta := solve(matQ) %*% t(matH) %*% matW
  LLT<MatrixXd> lltQ(matQ);
  if(lltQ.info() == NumericalIssue){
    // std::cout<<"Non-PSD matrix Q!"<<std::endl;
    return(NLL_CAP);
  }
  MatrixXd matK = lltQ.matrixL();
  MatrixXd mat2 = lltQ.matrixL().solve(matH.transpose());
  MatrixXd mat3 = mat2 * matW;
  MatrixXd matBeta = lltQ.matrixU().solve(mat3);

  // Get vecU := solve(matL) %*% matY - matH %*% matBeta
  MatrixXd mat4 = matY - matH*matBeta;
  MatrixXd matU = lltA.matrixL().solve(mat4);
  
  // Get nll
  MatrixXd sigHat2 = MatrixXd(r,r).setZero().selfadjointView<Lower>()
                                  .rankUpdate(matU.adjoint());

  MatrixXd sigHat2Chol = sigHat2.llt().matrixL();
  double comp3 = sigHat2Chol.diagonal().array().log().sum()*(p-q);
  double comp1 = matL.diagonal().array().log().sum()*r;
  double comp2 = matK.diagonal().array().log().sum()*r;
  double nll = comp1 + comp2 + comp3;
  
  return(nll);
}