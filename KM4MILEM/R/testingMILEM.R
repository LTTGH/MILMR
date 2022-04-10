#' testingMILEM
#' 
#' Compute the empirical risk of cv or testing for kernel machies with missing covariates
#'
#' @param trainRes  Results results from trainging results
#' @param testData testing data set
#' @return YHat the empirical risk of the test
#'
#' @importFrom kernlab rbfdot vanilladot
#' @import stats
#' @export
#'
#'
testingMILEM=function(MILdatTst,trainRes)
{
  ## results obtained from thre training function
  p=trainRes$p
  ki=trainRes$ki
  alphaHat=trainRes$alphaHat
  kerType=trainRes$kerType
  sigma=trainRes$sigma
  lambda=trainRes$lambda
  testPurpose=trainRes$testPurpose
  Xtrain_inst=trainRes$Xtrain_inst
  
  ##
  XTst=as.matrix(MILdatTst[,1:(ki*p)])
  nTst=dim(XTst)[1]
  YTst=MILdatTst[,(ki*p)+1]
  RTst=MILdatTst[,(ki*p)+2]
  XTst_inst=t(matrix(t(XTst),p))
  
  if(kerType=="RBF"){
    kerFun=rbfdot(sigma)
  }else{kerFun=vanilladot()} 
  ## get kernel vector (K(X_inst[1,],x),K(X_inst[2,],x),.....K(X_inst[n*ki,],x))
  KTst=kernelMatrix(kerFun,Xtrain_inst,XTst_inst)
  fHat=alphaHat%*%KTst
  Y_instHat=c(sign(fHat))
  Y_instHatMatrix=matrix(Y_instHat,nTst,ki,byrow = TRUE)
  Y_Hat=apply(Y_instHatMatrix,1,max)
  switch(testPurpose, 
         "cv"= {empiricalError=1-mean(Y_Hat==YTst)
         Y_labelHat=c(empiricalError,rep(0,(nTst-1)))
         },
         "testing"={Y_labelHat=Y_Hat
         }
         )# end of switch
  return(Y_labelHat)
}# enf of testingMILEM function
