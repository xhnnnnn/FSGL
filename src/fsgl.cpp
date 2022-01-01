#include <iostream>
#include <math.h>
#include <R.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
//' @import Rcpp
using namespace arma;
using namespace std;
using namespace Rcpp;

double l2norm(arma::mat x){
    int n = x.n_rows;
    double norm = 0;

    for(int i = 0; i < n; i++){
        norm += accu(pow(x.row(i),2));
    }
    norm = sqrt(norm);
    return norm;
}

double mutilply_l2norm(arma::mat X, arma::mat Y){
    arma::mat K(X.n_cols, Y.n_cols);
    K = X.t()*Y;

    double norm;
    norm = l2norm(K);
    return norm;
}

void PreCompute_K_XX(arma::mat X,  arma::mat& K, arma::field<arma::mat>& XX,int numGroup, arma::vec rangeGroupInd, arma::vec groupLen){
    for(int i = 0; i < numGroup; i++){
        for(int j = 0; j < numGroup; j++){
            if(i == j){
                XX(i) = X.cols(rangeGroupInd(i), rangeGroupInd(i)+groupLen(i)-1).t() * X.cols(rangeGroupInd(j), rangeGroupInd(j)+groupLen(j)-1);
                K(i, i) = l2norm(XX(i));
            }else{
                K(i, j) = mutilply_l2norm(X.cols(rangeGroupInd(i), rangeGroupInd(i) +groupLen(i)-1), X.cols(rangeGroupInd(j), rangeGroupInd(j)+groupLen(j)-1));
            }
        }
    }
}
void PreCompute_XY(arma::mat X,  arma::vec Y, arma::field<arma::mat>& XY,int numGroup, arma::vec rangeGroupInd, arma::vec groupLen){
    for(int i = 0; i< numGroup; i++){
        XY(i) = X.cols(rangeGroupInd(i), rangeGroupInd(i)+groupLen(i)-1).t() * Y;
    }
}


//subset function
arma::vec subsetVector(arma::vec v, int i,int j){
    arma::vec result(j-i);
    for(int k = 0; k < j-i; k++){
        result(k) = v(k+i);
    }

    return result;
}

//calculate loglikelihood
double linNegLogLikelihoodCalc(int nrow,  arma::vec eta,  arma::vec y){
    double squareSum = 0;
    squareSum = accu(pow(eta - y, 2));
    return squareSum/nrow/2;
}



//main algorithm
//' @useDynLib FSGL
void linear_solve( arma::mat X,  arma::vec y,  arma::vec &beta,int numGroup,  arma::vec rangeGroupInd,
                       arma::vec groupLen, double lambda1,  double lambda2,  int innerIter,  double thresh,
                       double gamma, arma::vec & betaIsZero, int & groupChange, arma::vec& isActive, arma::vec & useGroup,
                       double step,  int reset, arma::mat K, arma::field <arma::mat> XX, arma::field<arma::mat> XY ,
                       arma::vec &eta, arma::vec IsCandidate, int UseUpperBound, arma::vec & U, arma::vec delta, arma::vec betaNull)
{
    int nrow = X.n_rows;
    int ncol = X.n_cols;
    arma::vec theta(ncol);
    int startInd = 0;
    double zeroCheck = 0;
    double check = 0;
    int count = 0;
    double t = step;
    double diff = 1;
    double norm = 0;
    double u0p = 0;
    double Lnew = 0;
    double Lold = 0;
    double sqNormG = 0;
    double iProd = 0;
    arma::vec etaNew(nrow);
    arma::vec etaNull(nrow);
    arma::vec grad;


    for(int i = 0; i < numGroup; i++){
        if(useGroup(i)==1 && IsCandidate(i)==1){
            startInd = rangeGroupInd(i);
            arma::vec tempGroupBeta(groupLen(i));
            arma::mat tempGroupX;
            tempGroupX = X.cols(startInd, startInd+groupLen(i)-1);
            tempGroupBeta = subsetVector(beta, startInd, startInd+groupLen(i));
            grad = (XX(i) * tempGroupBeta - XY(i))/nrow;
            //Use Upper Bound
            if(UseUpperBound==1){
                if(U(i) < sqrt(groupLen(i))*lambda2){
                    if(betaIsZero(i) == 0){
                        eta = eta - tempGroupX * tempGroupBeta;
                    }
                    betaIsZero(i) = 1;
                    for(int j = 0; j < groupLen(i); j++){
                        beta(j+startInd) = 0;
                    }
                    U(i) = U(i) - 2*delta(i)  + 2*K(i,i)*sqrt(accu(pow(beta-betaNull,2)));
                    delta(i) = K(i,i)*sqrt(accu(pow(beta-betaNull,2)));
                    continue;
                }
            }

            //soft max calculas
            for(int j = 0; j < groupLen(i); j++){
                if((grad(j) < lambda1)&& (grad(j) > -lambda1)){
                    grad(j) = 0;
                }else if(grad(j) > lambda1){
                    grad(j) = grad(j) -lambda1;
                }else if(grad(j) < -lambda1){
                    grad(j) = grad(j) + lambda1;
                }else if(pow(grad(j), 2) == pow(lambda1,2)){
                    grad(j) = 0;
                }
            }

            zeroCheck = accu(pow(grad, 2));
            //group level decision
            if(zeroCheck <= pow(lambda2, 2)*groupLen(i)){
                if(betaIsZero(i) == 0){
                    eta = eta - tempGroupX * tempGroupBeta;
                }
                betaIsZero(i) = 1;
                for(int j = 0; j < groupLen(i); j++){
                    beta(j+startInd) = 0;
                }
            }else{
                if(isActive(i) == 0){
                    groupChange = 1;
                }
                isActive(i) = 1;
                theta = beta;
                arma::vec z(groupLen(i));
                arma::vec U(groupLen(i));
                arma::vec G(groupLen(i));
                //inner loop
                count = 0;
                check = 100000;
                while(count <= innerIter && check > thresh){
                    count++;
                    tempGroupBeta = subsetVector(beta, startInd, startInd+groupLen(i));
                    grad =  (XX(i) * tempGroupBeta - XY(i))/nrow;
                    diff = -1;
                    Lold = linNegLogLikelihoodCalc(nrow, eta, y);
                    //check the QM condition
                    while(diff < 0){
                        for(int j = 0; j < groupLen(i); j++)
                        {
                            z(j) = beta(j + startInd) - t * grad(j);
                            if(z(j) < lambda1 * t && z(j) > -lambda1 * t)
                            {
                                z(j) = 0;
                            }
                            if(z(j) > lambda1 * t)
                            {
                                z(j) = z(j) - lambda1 * t;
                            }
                            if(z(j) < -lambda1 * t)
                            {
                                z(j) = z(j) + lambda1 * t;
                            }
                        }
                        norm = sqrt(accu(pow(z,2)));
                        if(norm!=0){
                            u0p = (1-lambda2*sqrt(double(groupLen(i)))*t/norm);
                        }else{u0p = 0;}
                        if(u0p < 0){u0p=0;}
                        //update U G
                        U = u0p * z;
                        G = 1/t*(tempGroupBeta - U);
                        //Setting up betaNew and etaNew in direction of Grad for descent step
                        etaNew = eta-t*tempGroupX*G;
                        Lnew = linNegLogLikelihoodCalc(nrow, etaNew, y);

                        sqNormG = accu(pow(G, 2));
                        iProd = accu(dot(grad, G));
                        diff = Lold - Lnew - t * iProd + t/2 * sqNormG;
                        t = t * gamma;
                    }
                    t = t/gamma;

                    //Nesterov step to update theta
                    check= accu(abs(subsetVector(theta, startInd, startInd+groupLen(i)) - U));
                    eta = eta - tempGroupX * tempGroupBeta;
                    for(int j = 0; j < groupLen(i); j++){
                        beta(j+startInd) = U(j) + count%reset/(count%reset+3) * (U(j) - theta(j+startInd));
                        theta(j+startInd) = U(j);
                    }

                    eta = eta + tempGroupX * subsetVector(beta, startInd, startInd+groupLen(i));
                }
            }
            //update U(g)
            if(UseUpperBound==1){
                U(i) = U(i) - 2*delta(i)  + 2*K(i,i)*sqrt(accu(pow(beta-betaNull,2)));
                delta(i) = K(i,i)*sqrt(accu(pow(beta-betaNull,2)));
            }
        }
    }
}


//main  algorithms of Fast Sparse Group Lasso
//' @useDynLib FSGL
// [[Rcpp::export()]]
Rcpp::List FastSGL   ( arma::mat X,  arma::vec y, int numGroup,  arma::vec rangeGroupInd,
                      arma::vec groupLen, arma::vec lambdas,  double alpha,int innerIter, int outerIter,
                      double thresh, double outerThresh,double gamma,  double step, int reset, int UseUPB)
{
    //parameter initial
    int outermostCounter;
    int n = X.n_rows;
    int p = X.n_cols;
    int nlam = lambdas.n_elem;
    int groupChange;
    double lambda1;
    double lambda2;
    int UseUpperBound;
    arma::mat beta(p, nlam);
    arma::mat betaIsZero(p, nlam);

    arma::vec isActive(numGroup);
    arma::vec useGroup(numGroup);
    arma::vec tempIsActive(numGroup);
    arma::vec tempBeta(p);
    arma::vec tempBetaIsZero(p);
    arma::vec betaNull(p);
    arma::vec delta(numGroup); // delta(g,g) for group g

    //outer loop paramter
    double outermostCheck;
    double C;
    arma::vec outerOldBeta;
    arma::vec eta(n);
    arma::vec IsCandidate(numGroup);
    arma::vec U(numGroup);

    //precompute
    arma::mat K(numGroup, numGroup);
    arma::field<arma::mat> XX(numGroup);
    arma::field<arma::mat> XY(numGroup);
    PreCompute_K_XX(X, K, XX, numGroup, rangeGroupInd, groupLen);
    PreCompute_XY(X, y, XY, numGroup, rangeGroupInd, groupLen);

    for(int i = 0; i < nlam; i++){
        UseUpperBound = 0;
        outermostCounter = 0;
        outermostCheck = 100000;
        groupChange = 1;
        isActive.fill(0);
        useGroup.fill(1);
        U.fill(0);
        betaNull.fill(0);
        tempBeta.fill(0);
        tempBetaIsZero.fill(1);
        delta.fill(0);
        eta.fill(0);
        lambda1 = lambdas(i)*alpha;
        lambda2 = lambdas(i)*(1-alpha);
        //candidate set
        for(int j = 0; j < numGroup; j++){
            C = sqrt(accu(pow(XY(j),2))) - lambda1 * sqrt(groupLen(j)/2);
            if(C > sqrt(groupLen(j)) * lambda2){
                IsCandidate(j) = 1;
            }else{
                IsCandidate(j) = 0;
            }
        }
        //candidate main loop
        while(groupChange==1){
            groupChange = 0;
            linear_solve(X, y, tempBeta, numGroup,rangeGroupInd, groupLen, lambda1, lambda2, innerIter, thresh, gamma, tempBetaIsZero, groupChange, isActive, useGroup, step, reset, K, XX, XY, eta, IsCandidate, UseUpperBound, U, delta, betaNull);
            while(outermostCounter < outerIter && outermostCheck > outerThresh ){
                outermostCounter++;
                outerOldBeta = tempBeta;
                tempIsActive = isActive;
                linear_solve(X, y, tempBeta, numGroup,rangeGroupInd, groupLen, lambda1, lambda2, innerIter, thresh, gamma, tempBetaIsZero, groupChange, isActive, tempIsActive, step, reset, K, XX, XY, eta, IsCandidate, UseUpperBound, U, delta, betaNull);
                outermostCheck = accu(abs(outerOldBeta-tempBeta));
            }
        }
        //reset parameter
        betaNull=tempBeta;
        if(UseUPB==1){
          UseUpperBound = 1;
          //compute U(g)
          for(int k = 0; k< numGroup; k++){
            delta(k) = 0;
            U(k) = sqrt(accu(pow((XY(k)-XX(k)*subsetVector(betaNull, rangeGroupInd(k), rangeGroupInd(k)+groupLen(k)))/n, 2))) + delta(k);
          }
        }else{
          UseUpperBound = 0;
        }

        outermostCounter = 0;
        outermostCheck = 100000;
        groupChange = 1;
        isActive.fill(0);
        useGroup.fill(1);
        IsCandidate.fill(1);

        //main loop of upperErrorBound
        while(groupChange==1){
            groupChange = 0;
            linear_solve(X, y, tempBeta, numGroup,rangeGroupInd, groupLen, lambda1, lambda2, innerIter, thresh, gamma, tempBetaIsZero, groupChange, isActive, useGroup, step, reset, K, XX, XY, eta, IsCandidate, UseUpperBound, U , delta, betaNull);
            while(outermostCounter < outerIter && outermostCheck > outerThresh ){
                outermostCounter++;
                outerOldBeta = tempBeta;
                tempIsActive = isActive;
                linear_solve(X, y, tempBeta, numGroup,rangeGroupInd, groupLen, lambda1, lambda2, innerIter, thresh, gamma, tempBetaIsZero, groupChange, isActive, tempIsActive, step, reset, K, XX, XY, eta, IsCandidate, UseUpperBound, U, delta, betaNull);
                outermostCheck = accu(abs(outerOldBeta-tempBeta));
            }
        }

        beta.col(i) = tempBeta;
        betaIsZero.col(i) = tempBetaIsZero;
    }

    //输出结果
    return Rcpp::List::create(
        Rcpp::Named("beta") = beta,
        // Rcpp::Named("betaIsZero") = betaIsZero
        //    Rcpp::Named("eta") = eta,
        Rcpp::Named("outerIter") = outermostCounter
    );
}



