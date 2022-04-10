#' testingMILEM_GK_multik
#'
#' Compute the empirical risk of cv or testing for MIL with multiple instance
#' @param MILdatTst testing data set
#' @param trainRes  Results results from trainging results
#' @return YHat the empirical risk of the test
#'
#' @importFrom kernlab rbfdot vanilladot
#' @import stats
#' @export
#'
#'
testingMILEM_bag_cv=function(MILdatTst,trainRes)
{
  ## results obtained from thre training function
  p=trainRes$p
  alphaHat=trainRes$alphaHat
  kerType=trainRes$kerType
  sigma=trainRes$sigma
  lambda=trainRes$lambda
  testPurpose=trainRes$testPurpose
  seperateNum=trainRes$seperateNum
  Xtrain_inst=trainRes$Xtrain_inst
  
  
  MILdatTst=MILdatTst[order(MILdatTst$bagName),]
  testBags=MILdatTst$bagName
  nPB_test=sum(testBags<seperateNum)
  nNB_test=sum(testBags>=seperateNum)
  nTst=nPB_test+nNB_test
  kVec_test=MILdatTst$numInst
  
  YTst=c(rep(1,nPB_test),rep(-1,nNB_test))
  RTst=c(rep(0,nPB_test),rep(1,nNB_test))
  
  MILdatTst_inst=matrix(0,sum(kVec_test),p)
  MILdatTst_inst[1:kVec_test[1],]=t(matrix(t(as.matrix(MILdatTst[1,3:(p*kVec_test[1]+2)])),p))
  for(i in 2:nTst){
    startLabel=sum(kVec_test[1:(i-1)])+1
    endLabel=sum(kVec_test[1:i])
    MILdatTst_inst[startLabel:endLabel,]=t(matrix(t(as.matrix(MILdatTst[i,3:(p*kVec_test[i]+2)])),p))
  }# end of for
  XTst_inst=as.matrix(MILdatTst_inst)
  
  
  if(kerType=="RBF"){
    kerFun=rbfdot(sigma)
  }else{kerFun=vanilladot()}
  ## get kernel vector (K(X_inst[1,],x),K(X_inst[2,],x),.....K(X_inst[n*ki,],x))
  KTst=kernelMatrix(kerFun,Xtrain_inst,XTst_inst)
  fHat=alphaHat%*%KTst
  Y_instHat=c(sign(fHat))
  
  Y_Hat=rep(0,nTst)
  Y_Hat[1]=max(Y_instHat[1:kVec_test[1]])
  for(iB in 2:nTst){
    labelStart=sum(kVec_test[1:(iB-1)])+1
    labelEnd=sum(kVec_test[1:iB])
    Y_Hat[iB]=max(Y_instHat[labelStart:labelEnd])
  }# estimated bag label
  
  switch(testPurpose,
         "cv"= {empiricalError=1-mean(Y_Hat==YTst)
         Y_labelHat=c(empiricalError,rep(0,(nTst-1)))
         },
         "testing"={Y_labelHat=Y_Hat
         }
  )# end of switch
  return(Y_labelHat)
}# enf of testingMILEM function
