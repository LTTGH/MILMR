Sys.setenv(RETICULATE_PYTHON = "C:/Users/liutiantian/anaconda3") 
library(reticulate)
source_python("C:/Users/liutiantian/Desktop/50/models_S1k250R.py")
source_python("C:/Users/liutiantian/Desktop/50/DataLoader_S1k250R.py")
source_python("C:/Users/liutiantian/Desktop/50/train_test_script_S1k250R.py")
source_python("C:/Users/liutiantian/Desktop/50/main_S1k250R.py")

source_python("C:/Users/liutiantian/Desktop/MILpy/S1/k2/100/models_S1k2100R.py")
source_python("C:/Users/liutiantian/Desktop/MILpy/S1/k2/100/DataLoader_S1k2100R.py ")
source_python("C:/Users/liutiantian/Desktop/MILpy/S1/k2/100/train_test_script_S1k2100R.py")
source_python("C:/Users/liutiantian/Desktop/MILpy/S1/k2/100/main_S1k2100R.py")


source_python("C:/Users/liutiantian/Desktop/MILpy/S1/k2/200/models_S1k2200R.py")
source_python("C:/Users/liutiantian/Desktop/MILpy/S1/k2/200/DataLoader_S1k2200R.py ")
source_python("C:/Users/liutiantian/Desktop/MILpy/S1/k2/200/train_test_script_S1k2200R.py")
source_python("C:/Users/liutiantian/Desktop/MILpy/S1/k2/200/main_S1k2200R.py")



source_python("C:/Users/liutiantian/Desktop/MILpy/S1/k4/50/models_S1k450R.py")
source_python("C:/Users/liutiantian/Desktop/MILpy/S1/k4/50/DataLoader_S1k450R.py ")
source_python("C:/Users/liutiantian/Desktop/MILpy/S1/k4/50/train_test_script_S1k450R.py")
source_python("C:/Users/liutiantian/Desktop/MILpy/S1/k4/50/main_S1k450R.py")


source_python("C:/Users/liutiantian/Desktop/MILpy/S1/k4/200/models_S1k4200R.py")
source_python("C:/Users/liutiantian/Desktop/MILpy/S1/k4/200/DataLoader_S1k4200R.py ")
source_python("C:/Users/liutiantian/Desktop/MILpy/S1/k4/200/train_test_script_S1k4200R.py")
source_python("C:/Users/liutiantian/Desktop/MILpy/S1/k4/200/main_S1k4200R.py")



## Check the result using the data read from csv
library(MASS)
library(e1071)
library(rootSolve)
p=5
ki=2
maxIter=100
tol=0.001
## SVM based on maximum pattern margin
source('C:/Users/liutiantian/Desktop/MILcode/1711/computer2/MILcode/miSVM.R')
source('C:/Users/liutiantian/Desktop/MILcode/1711/computer2/MILcode/miSVM_method.R')
## SVM based on maximum bag margin
source('C:/Users/liutiantian/Desktop/MILcode/1711/computer2/MILcode/MI_SVM.R')
source('C:/Users/liutiantian/Desktop/MILcode/1711/computer2/MILcode/MI_SVM_method.R')

MILtrain=read.csv('C:/Users/liutiantian/Desktop/Setting1DATA/k2/200/s1k2_200_5.csv')[,-1]
MILtest=read.csv('C:/Users/liutiantian/Desktop/Setting1DATA/k2/200/s1k2_200_Test.csv')[,-1]
set.seed(5)
MISVM_method(MILtrain,MILtest,p,ki,maxIter)
miSVM_method(MILtrain,MILtest,p,ki,maxIter)
