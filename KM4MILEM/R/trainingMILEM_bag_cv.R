#' trainingMILEM_GK_multik
#'
#' This is training function for multiple instance learning using EM algorithm with multiple instances case
#' MILdat,tol,maxIter,kerType,sigma,lambda
#' @param MILdat multiple instance learning data
#' @param p dimensions of instance covariance
#' @param seperateNum number of seperating positive bags and negative bags
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
trainingMILEM_bag_cv=function(MILdat,p,seperateNum,tol,kerType,sigma,lambda,testPurpose)
{
  MILdat=MILdat[order(MILdat$bagName),]
  trainBags=MILdat$bagName
  nPB=sum(trainBags<seperateNum)
  nNB=sum(trainBags>=seperateNum)
  sampleSize=nPB+nNB # number of bags
  kVec=MILdat$numInst
  MILdat_inst=matrix(0,sum(kVec),p)
  MILdat_inst[1:kVec[1],]=t(matrix(t(as.matrix(MILdat[1,3:(p*kVec[1]+2)])),p))
  for(i in 2:sampleSize){
    startLabel=sum(kVec[1:(i-1)])+1
    endLabel=sum(kVec[1:i])
    MILdat_inst[startLabel:endLabel,]=t(matrix(t(as.matrix(MILdat[i,3:(p*kVec[i]+2)])),p))
  }# end of for
  
  n_inst=dim(MILdat_inst)[1]
  
  nPB_ins=sum(kVec[1:nPB])
  nNB_ins=sum(kVec[(nPB+1):(nPB+nNB)])
  X_PBins=MILdat_inst[1:nPB_ins,]
  
  X_NBins=MILdat_inst[(nPB_ins+1):n_inst,]
  Y_NBins=rep(-1,nNB_ins)
  X_inst=rbind(X_PBins,X_NBins) # covariates of all instances
  
  
  Y=c(rep(1,nPB),rep(-1,nNB)) #Labels of bags
  R=c(rep(0,nPB),rep(1,nNB)) #Labels of bags
  
  # set kernel function
  if(kerType=="RBF"){
    kerFun=rbfdot(sigma)# The Gaussian RBF kernel k(x,x') = exp(-sigma||x - x'||^2)
  }else{kerFun=vanilladot()} # RBF kernel or linear kernel
  
  ## get kernel matrix
  K=kernelMatrix(kerFun,X_inst)
  
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
                lambda=lambda,testPurpose=testPurpose,Xtrain_inst=X_inst,seperateNum=seperateNum)
  return(trainRes)
}



