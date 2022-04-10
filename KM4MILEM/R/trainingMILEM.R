#' trainingMILEM
#'
#' This is training function for multiple instance learning using EM algorithm
#' MILdat,tol,maxIter,kerType,sigma,lambda
#' @param MILdat multiple instance learning data
#' @param p dimensions of instance covariance
#' @param ki number of instances in each bags
#' @param tol error tolerance for EM algorithm
#' @param kerType type of kerenl function, "RBF"(if choose) or "linear" (else choose, automatically)
#' @param lambda tuning parameter
#' @param sigma parameter for RBF kernel, also act as tuning parameter
#' @param testPurpose test purpose, for crossivalidation, the loss function is chosen as the phi loss;
#'                                  for test (default value) , the loss function is chosen as the classification loss.
#' @return trainRes
#'
#' @importFrom kernlab rbfdot vanilladot kernelMatrix
#' @importFrom rootSolve multiroot
#' @export
#'
trainingMILEM=function(MILdat,p,ki,tol,kerType,sigma,lambda,testPurpose)
{
  X=as.matrix(MILdat[,1:(ki*p)])
  Y=MILdat[,(ki*p)+1]
  R=MILdat[,(ki*p)+2]
  sampleSize=length(Y)
  index1=c(1:p)
  index2=c((p+1):(ki*p))
  X_inst=t(matrix(t(X),p))
  n_inst=dim(X_inst)[1]
  # set kernel function
  if(kerType=="RBF"){
    kerFun=rbfdot(sigma)# The Gaussian RBF kernel k(x,x') = exp(-sigma||x - x'||^2)
  }else{kerFun=vanilladot()} # RBF kernel or linear kernel
  
  ## get kernel matrix
  K=kernelMatrix(kerFun,X_inst)
  ##logit function
  logistic<-function(x){return(1/(1+exp(-x)))} # define lositic function
  
  # function g(X_i, beta^t) for the conditional distribution
  g=function(Ker_X,alpha_t){
    logistic1=logistic(sum(Ker_X[,1]*alpha_t))
    logistic2=logistic(sum(Ker_X[,2]*alpha_t))
    g_x_beta=(1-logistic1)*logistic2+(1-logistic2)*logistic1+logistic1*logistic2
    return(g_x_beta)
  }
  
  ##omega function
  wFun=function(Ker_X,alpha_t,group){
    logistic1=logistic(sum(Ker_X[,1]*alpha_t))
    logistic2=logistic(sum(Ker_X[,2]*alpha_t))
    return(switch(group, # 1: biv normal; 2: biv Bernoulli
                  "1" = {(1-logistic1)*logistic2/g(Ker_X,alpha_t) },
                  "2" = {(1-logistic2)*logistic1/g(Ker_X,alpha_t)} ) # end of switch
    )# end of return
  }# end of w fucntion
  
  ## Two element in derivative Q
  derivativeQ12=function(r,Ker_X,alpha,alpha_t,group){
    if(group==1){
      ##first part of derivativeQ
      return(as.numeric((1-r)*(1-wFun(Ker_X,alpha_t,1)))-as.numeric(logistic(sum(Ker_X[,1]*alpha))))}else{
        ##second part of derivativeQ
      return(as.numeric((1-r)*(1-wFun(Ker_X,alpha_t,2)))-as.numeric(logistic(sum(Ker_X[,2]*alpha))))}
  }# end of derivativeQ12 function
  
  nonLinVec=rep(0,n_inst)
  derivativeQ=function(alpha,alpha_t){
    for(iSample in 1:sampleSize){
      startLabel=(iSample-1)*ki+1
      endLabel=iSample*ki
      nonLinVec[startLabel]=derivativeQ12(R[iSample],K[,startLabel:endLabel],alpha,alpha_t,1)
      nonLinVec[endLabel]=derivativeQ12(R[iSample],K[,startLabel:endLabel],alpha,alpha_t,2)}# end of recycling
    derivative=1/(n_inst)*nonLinVec-2*lambda*alpha # the estimating equation for alpha
    return(derivative)
  } # end of function
  
  ## EM algorithm
  alpha_Previous=rep(-1,n_inst) # dummy start
  alpha_Current=rep(0,n_inst)
  maxIter=100
  alphaMatrix=matrix(0,maxIter,n_inst) # matrix for keep beta
  numIter=0
  
  ## E-M algorithm
  while(sum(abs(alpha_Current-alpha_Previous))>tol & numIter<maxIter)
  {
    numIter=numIter+1
    alpha_Previous=alpha_Current
    ## E step
    solveFun=function(alpha){
      return(derivativeQ(alpha,alpha_Current))} # end of solveFun
    ## M step
    Res=multiroot(solveFun,rep(0,n_inst),maxiter = 200) # M-setp
    alphaMatrix[numIter,]=Res$root
    alpha_Current=Res$root
  }# end of while
  alphaEM=alphaMatrix[which(rowSums(alphaMatrix)!=0),]
  iterConvergNumKern=dim(alphaEM)[1]
  alphaHat=alphaEM[iterConvergNumKern,]
  trainRes=list(p=p,ki=ki,alphaHat=alphaHat,kerType=kerType,sigma=sigma,
                lambda=lambda,testPurpose=testPurpose,Xtrain_inst=X_inst)
  return(trainRes)
}# end of function MIEMkernel