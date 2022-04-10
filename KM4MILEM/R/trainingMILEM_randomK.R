#' trainingMILEM_randomK
#'
#' This is training function for multiple instance learning using EM algorithm with multiple instances case
#' MILdat,tol,maxIter,kerType,sigma,lambda
#' @param MILdat multiple instance learning data
#' @param p dimensions of instance covariance
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
trainingMILEM_randomK=function(MILdat,p,tol,kerType,sigma,lambda,testPurpose)
{
  trainBags=MILdat$Y1
  nPB=sum(trainBags==1)
  nNB=sum(trainBags==-1)
  sampleSize=nPB+nNB # number of bags
  kVec=MILdat$numInst
  R=MILdat$R
  MILdat_inst=matrix(0,sum(kVec),p)
  MILdat_inst[1:kVec[1],]=t(matrix(t(as.matrix(MILdat[1,1:(p*kVec[1])])),p))
  for(i in 2:sampleSize){
    startLabel=sum(kVec[1:(i-1)])+1
    endLabel=sum(kVec[1:i])
    MILdat_inst[startLabel:endLabel,]=t(matrix(t(as.matrix(MILdat[i,1:(p*kVec[i])])),p))
  }# end of for
  
  n_inst=dim(MILdat_inst)[1]
  
  # set kernel function
  if(kerType=="RBF"){
    kerFun=rbfdot(sigma)# The Gaussian RBF kernel k(x,x') = exp(-sigma||x - x'||^2)
  }else{kerFun=vanilladot()} # RBF kernel or linear kernel
  
  ## get kernel matrix
  K=kernelMatrix(kerFun,MILdat_inst)
  
  ##logit function
  logistic<-function(x){return(1/(1+exp(-x)))} # define lositic function
  
  # ##Traverse matrix
  # bin01=as.data.frame(matrix(c(rep(0, ki),rep(1, ki)),ncol=ki, byrow = TRUE))
  # binMatrix= expand.grid(bin01)[-1,]
  # caseNum=dim(binMatrix)[1]
  
  
  ## product function
  product <- function(vec){
    out <- 1
    for(i in 1:length(vec)){out <- out*vec[i]}
    return(out) }# end of product function
  
  
  ## summmation of weights for each instance
  weiSum=function(x,alpha_t,ki){
    logiProb=logistic(alpha_t%*%x) # logit probability of instances for each bag
    ##Traverse matrix
    bin01=as.data.frame(matrix(c(rep(0, ki),rep(1, ki)),ncol=ki, byrow = TRUE))
    binMatrix= expand.grid(bin01)[-1,]
    caseNum=dim(binMatrix)[1]
    posLogit=matrix(logiProb,nrow=caseNum,ncol=ki, byrow = TRUE)
    weightMatrix=binMatrix*posLogit+(1-binMatrix)*(1-posLogit)
    weightVec=apply(weightMatrix,1,product)/(1-product(1-logiProb)) # ogema
    weightSum=apply(weightVec*binMatrix,2,sum) # summmation of weights for each instance
    return(weightSum)
  }
  
  nonLinVec1=rep(0,n_inst)
  derivativeQ1=function(alpha_t){
    if(R[1]==0){nonLinVec1[1:kVec[1]]=weiSum(K[,1:kVec[1]],alpha_t,kVec[1])} 
    
    for(iSample in 2:sampleSize){
      startLabel=sum(kVec[1:(iSample-1)])+1
      endLabel=sum(kVec[1:iSample])
      if(R[iSample]==0){
        nonLinVec1[startLabel:endLabel]=weiSum(K[,startLabel:endLabel],alpha_t,kVec[iSample])} #end of else.
    } # end of for
    derivative1=1/(n_inst)*nonLinVec1 # the estimating equation for alpha
    return(derivative1)
  } # end of function
  
  nonLinVec2=rep(0,n_inst)
  derivativeQ2=function(alpha){
    logit_kernelVec=logistic(alpha%*%K) ## exp(alpha^{t}K)/(1+exp(alpha^{t}K))
    nonLinVec2[1:kVec[1]]=-logit_kernelVec[1:kVec[1]]
    for(iSample in 2:sampleSize){
      startLabel=sum(kVec[1:(iSample-1)])+1
      endLabel=sum(kVec[1:iSample])
      nonLinVec2[startLabel:endLabel]=-logit_kernelVec[startLabel:endLabel]
    } # end of for
    derivative2=1/(n_inst)*nonLinVec2-2*lambda*alpha # the estimating equation for alpha
    return(derivative2)
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
    Q1=derivativeQ1(alpha_Current)
    ## E step
    solveFun=function(alpha){
      return(Q1+derivativeQ2(alpha))} # end of solveFun
    ## M step
    Res=multiroot(solveFun,rep(0,n_inst),maxiter = 200) # M-setp
    alphaMatrix[numIter,]=Res$root
    alpha_Current=Res$root
  }# end of while
  if(max(which(rowSums(alphaMatrix)!=0))>1){
    alphaEM=alphaMatrix[which(rowSums(alphaMatrix)!=0),]
    iterConvergNumKern=dim(alphaEM)[1]
    alphaHat=alphaEM[iterConvergNumKern,]}else{
      alphaEM=alphaMatrix[which(rowSums(alphaMatrix)!=0),]
      alphaHat=alphaEM}#end of if else
  
  trainRes=list(p=p,kVec=kVec,alphaHat=alphaHat,kerType=kerType,sigma=sigma,
                lambda=lambda,testPurpose=testPurpose,Xtrain_inst=MILdat_inst)
  return(trainRes)
}

