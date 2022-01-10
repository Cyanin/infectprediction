data_comp=data

########分层抽样，抽取20%为测试集，80%为训练集########
library(sampling)
#n=round(0.8*nrow(data_comp))
train_indel<-strata(data_comp,stratanames = ("XUEPY_RESULT"),size=c(350,350),method='srswor')
#write.csv(train,"train_sampling.csv")
indel<-sort(train_indel$ID_unit)
train<-data_comp[indel,]
varify_test<-data_comp[-indel,]

######################################
library(xgboost)
library(Matrix)
set.seed(100)

###参数调试
traindata2<-as(as.matrix(train[,c(2:ncol(train))]),"dgCMatrix") 
traindata3<-as.numeric(as.character(train[,1]))
traindata4<-list(data=traindata2,label=traindata3)
dtrain<-xgb.DMatrix(data=traindata4$data,label=traindata4$label)

varifyset2<-as(as.matrix(varify_test[,c(2:ncol(train))]),"dgCMatrix")
varifyset3<-as.numeric(as.character(varify_test[,1]))
varifyset4<-list(data=varifyset2,label=varifyset3)
dvarifyset<-xgb.DMatrix(data=varifyset4$data,label=varifyset4$label)

#####################网格搜索寻找最优值#######################
max_depth.seq<-seq(3,10,1)
subsample.seq<-seq(0.5,0.95,0.05)
colsample_bytree.seq<-seq(0.5,0.95,0.05)
min_child_weight.seq<-seq(1,6,1)
alpha.seq<-c(1e-5,1e-2,0.1,1,100)
gamma.seq<-seq(0,0.5,0.1)

n1=length(max_depth.seq)  #max_depth的备选数
n2=length(subsample.seq)
n3=length(colsample_bytree.seq)
n4=length(min_child_weight.seq)
n5=length(alpha.seq)
n6=length(gamma.seq)
r=array(0,c(n1,n2,n3,n4,n5,n6))

n_total=n1*n2*n3*n4*n5*n6
n_total=as.character(n_total)
count=0

system.time({
  for(ind.max_depth in 1:n1){
    max_depth<-max_depth.seq[ind.max_depth]
    for(ind.subample in 1:n2){
      subsample<-subsample.seq[ind.subample]
      for(ind.colsample_bytree in 1:n3){
        colsample_bytree<-colsample_bytree.seq[ind.subample]
        for(ind.min_child_weight in 1:n4){
          min_child_weight<-min_child_weight.seq[ind.min_child_weight]
          for(ind.alpha in 1:n5){
            reg_alpha=alpha.seq[ind.alpha]
            for(ind.gamma in 1:n6){
              reg_gamma=gamma.seq[ind.gamma]
              
              param=list(
                objective="binary:logistic",
                booster='gbtree',
                max_depth=max_depth,
                subsample=subsample,
                colsample_bytree=colsample_bytree,
                min_child_weight=min_child_weight,
                alpha=reg_alpha,
                gamma=reg_gamma
              )
              fit=xgb.train(
                params = param,
                data=dtrain,
                nrounds=10,
                #early_stopping_rounds = 50,
                watchlist=list(validation1=dvarifyset)
              )
              
              pred=predict(fit,varifyset2)
              pred=round(pred)
              
              acc=mean(varifyset3-pred==0)
              r[ind.max_depth,ind.subample,ind.colsample_bytree,ind.min_child_weight,ind.alpha,ind.gamma]=acc
              
              count=count+1
              n_now=formatC(count,width=nchar(n_total),flag='0')
              cat('\n','????:',n_now,'/',n_total,'\n\n',seq='')
            }
          }
        }
      }
    }
  }
})


rr=apply(r,1:6,mean)
max_acc=max(rr)
ind=which(rr==max_acc,arr.ind = TRUE)[1,]
max_depth=max_depth.seq[ind[1]]   
subsample=subsample.seq[ind[2]]   
colsample_bytree=colsample_bytree.seq[ind[3]]  
min_child_weight=min_child_weight.seq[ind[4]]  
reg_alpha=alpha.seq[ind[5]]  
reg_gamma=gamma.seq[ind[6]]  
max_depth;subsample;colsample_bytree;min_child_weight;reg_alpha;reg_gamma







####################################################################
###################################################################
######################对top变量进行调参工作########################

data_comp=data_top10

########分层抽样，抽取20%为测试集，80%为训练集########
library(sampling)
#n=round(0.8*nrow(data_comp))
train_indel<-strata(data_comp,stratanames = ("XUEPY_RESULT"),size=c(350,350),method='srswor')
#write.csv(train,"train_sampling.csv")
indel<-sort(train_indel$ID_unit)
train<-data_comp[indel,]
varify_test<-data_comp[-indel,]

######################################
library(xgboost)
library(Matrix)
set.seed(100)

###参数调试
traindata2<-as(as.matrix(train[,c(2:ncol(train))]),"dgCMatrix") 
traindata3<-as.numeric(as.character(train[,1]))
traindata4<-list(data=traindata2,label=traindata3)
dtrain<-xgb.DMatrix(data=traindata4$data,label=traindata4$label)

varifyset2<-as(as.matrix(varify_test[,c(2:ncol(train))]),"dgCMatrix")
varifyset3<-as.numeric(as.character(varify_test[,1]))
varifyset4<-list(data=varifyset2,label=varifyset3)
dvarifyset<-xgb.DMatrix(data=varifyset4$data,label=varifyset4$label)

#####################网格搜索寻找最优值#######################
max_depth.seq<-seq(3,10,1)
subsample.seq<-seq(0.5,0.95,0.05)
colsample_bytree.seq<-seq(0.5,0.95,0.05)
min_child_weight.seq<-seq(1,6,1)
alpha.seq<-c(1e-5,1e-2,0.1,1,100)
gamma.seq<-seq(0,0.5,0.1)

n1=length(max_depth.seq)  #max_depth的备选数
n2=length(subsample.seq)
n3=length(colsample_bytree.seq)
n4=length(min_child_weight.seq)
n5=length(alpha.seq)
n6=length(gamma.seq)
r=array(0,c(n1,n2,n3,n4,n5,n6))

n_total=n1*n2*n3*n4*n5*n6
n_total=as.character(n_total)
count=0

system.time({
  for(ind.max_depth in 1:n1){
    max_depth<-max_depth.seq[ind.max_depth]
    for(ind.subample in 1:n2){
      subsample<-subsample.seq[ind.subample]
      for(ind.colsample_bytree in 1:n3){
        colsample_bytree<-colsample_bytree.seq[ind.subample]
        for(ind.min_child_weight in 1:n4){
          min_child_weight<-min_child_weight.seq[ind.min_child_weight]
          for(ind.alpha in 1:n5){
            reg_alpha=alpha.seq[ind.alpha]
            for(ind.gamma in 1:n6){
              reg_gamma=gamma.seq[ind.gamma]
              
              param=list(
                objective="binary:logistic",
                booster='gbtree',
                max_depth=max_depth,
                subsample=subsample,
                colsample_bytree=colsample_bytree,
                min_child_weight=min_child_weight,
                alpha=reg_alpha,
                gamma=reg_gamma
              )
              fit=xgb.train(
                params = param,
                data=dtrain,
                nrounds=10,
                #early_stopping_rounds = 50,
                watchlist=list(validation1=dvarifyset)
              )
              
              pred=predict(fit,varifyset2)
              pred=round(pred)
              
              acc=mean(varifyset3-pred==0)
              r[ind.max_depth,ind.subample,ind.colsample_bytree,ind.min_child_weight,ind.alpha,ind.gamma]=acc
              
              count=count+1
              n_now=formatC(count,width=nchar(n_total),flag='0')
              cat('\n','????:',n_now,'/',n_total,'\n\n',seq='')
            }
          }
        }
      }
    }
  }
})


rr=apply(r,1:6,mean)
max_acc=max(rr)
ind=which(rr==max_acc,arr.ind = TRUE)[1,]
max_depth=max_depth.seq[ind[1]]   
subsample=subsample.seq[ind[2]]   
colsample_bytree=colsample_bytree.seq[ind[3]]  
min_child_weight=min_child_weight.seq[ind[4]]  
reg_alpha=alpha.seq[ind[5]]  
reg_gamma=gamma.seq[ind[6]]  
max_depth;subsample;colsample_bytree;min_child_weight;reg_alpha;reg_gamma