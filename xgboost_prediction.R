
########################################################################################################
######导入数据dataset，从0和1中随机分别抽取220例样本，并与result=2的样本组合为新的样本集################
########################################################################################################

dataset<-read.csv("F:/侵袭性真菌病项目—罗燕萍主任/2-第二次提取数据/3-数据分析/1-数据预处理/20190918-xgboos自动补值/dataset.csv")
dataset<-dataset[,-3]

##随机抽取对照组440例
indel_0<-sample(which(dataset$XUEPY_RESULT==0),220)
indel_1<-sample(which(dataset$XUEPY_RESULT==1),220)
indel<-c(indel_0,indel_1,which(dataset$XUEPY_RESULT==2))
dataset_new<-dataset[indel,]

#乱序排列
library(plyr)
dataset_new<-arrange(dataset_new,NUM)

#将result=0或1改为类别为0，result=2改为类别为1
dataset_new[dataset_new$XUEPY_RESULT==1,2]=0
dataset_new[dataset_new$XUEPY_RESULT==2,2]=1

#删掉序列号
dataset_new<-dataset_new[,-1]
write.csv(dataset_new,file='F:/侵袭性真菌病项目—罗燕萍主任/2-第二次提取数据/3-数据分析/1-数据预处理/20190918-xgboos自动补值/dataset_new.csv')


######################################################################################
############################数据分析建模##############################################
######################################################################################
setwd('F:/侵袭性真菌病项目—罗燕萍主任/2-第二次提取数据/3-数据分析/2-分析结果')
dataset_new<-read.csv('F:/侵袭性真菌病项目—罗燕萍主任/2-第二次提取数据/3-数据分析/1-数据预处理/20190918-xgboos自动补值/dataset_new.csv')[,-1]

data<-dataset_new[,-c(2,3,4,6)]

#转换格式
indel<-c(1,3,4,6:15,28:40)
data[,indel]<-lapply(data[,indel],as.factor)

##利用mice函数进行缺失值填补
library(lattice)
library(mice)
##填入数据
imp<-mice(data,m=5,seed=1234,maxit=10)
print(imp)

#############数据评价##################
#完整的数据集
data_comp1<-complete(imp)   #第一个完成数据集，默认有五个
data_comp2<-complete(imp,2)
data_comp3<-complete(imp,3)
data_comp4<-complete(imp,4)
data_comp5<-complete(imp,5)
for(i in 1:5){
  write.csv(complete(imp,m=i),file=paste0("data_comp",i,".csv"))
}
#检查原始数据集与插入数据的分布
stripplot(imp,pch=20,cex=1.2)  #散点图
densityplot(imp,col=c("orange","red","purple","blue","green"))  #密度函数图



##################################################################
##################################################################
##利用填补的数据3进行分析
#data<-data_comp3
data<-read.csv("F:/侵袭性真菌病项目—罗燕萍主任/2-第二次提取数据/3-数据分析/2-分析结果/data_comp3.csv")[,-c(1,30)]
#转换格式
indel<-c(1,3,4,6:15,28:39)
data[,indel]<-lapply(data[,indel],as.factor)

#######################十折交叉验证############
library(lattice)
library(ggplot2)
library(caret)
mm<-createFolds(y=data[,1],k=10)

###########################################################
#################logisit回归###############################
###########################################################
library(e1071)
library(stats)
library(ROCR)
library(pROC)
set.seed(100)

#利用全数据集得到最后的模型
lgst_model<-glm(data$XUEPY_RESULT ~ .,data,family=binomial(link = "logit"),control = list(maxit=100))
summary(lgst_model)

lgst_model<-glm(data$XUEPY_RESULT~pre_hos_days+season+age+IMMUNOSUPPRESSIVE+ANTIFUNGAL+INSTRUMENT_AIR+TUBE+NUTRITION_SUPPLY+TRAUMA+
                  NEUT+RBC+expectoration+hemopysis+cancer,family=binomial(link = "logit"),data,control = list(maxit=100))
summary(lgst_model)

lgst_model<-glm(data$XUEPY_RESULT~pre_hos_days+age+ANTIFUNGAL+TUBE+RBC +    INSTRUMENT_AIR+
                  NEUT+RBC,family=binomial(link = "logit"),data,control = list(maxit=100))
summary(lgst_model)

save(lgst_model,file="lgst_model.RData")

###画森林图
or_ci<-exp(confint(lgst_model) )   #输出OR的置信区间
or_value<-exp(coef(lgst_model)) #输出OR值
dat<-cbind(or_value,or_ci)
#dat<-read.csv('F:/侵袭性真菌病项目—罗燕萍主任/2-第二次提取数据/3-数据分析/2-分析结果/OR_CI.csv')

library(grid)
library(magrittr)
library(checkmate)
library(forestplot)
forestplot(rownames(dat),dat[,1],dat[,2],dat[,3],zero=1,graph.pos = 2,xlog = TRUE,lineheight=unit(1.2, "cm"))

##利用ROC选取最佳阈值
pre<-predict(lgst_model,type='response')
modelroc=roc(data$XUEPY_RESULT,pre)
plot(modelroc,print.auc=TRUE,auc.polygon=TRUE,grid=c(0.1,0.2),grid.col=c("green","red"),max.auc.polygon=TRUE,
     auc.polygon.col="skyblue",print.thres=TRUE)   #输出结果为0.477

##十折交叉
beta=1;Z=10;D=1
E=rep(0,Z);kap=rep(0,Z)
tpr=fpr=matrix(NA,nrow=Z,ncol=90);perfor=matrix(0,nrow=Z,ncol=5)
perfor_test=matrix(0,nrow=Z,ncol=5)
source("F:/mycode/Performance_iter.R")

for(i in 1:Z){
  m<-mm[[i]]   #另m为第i类的下标集 
  n1<-length(m)
  set.seed(1234)
  fit_model<-glm(XUEPY_RESULT~pre_hos_days+age+ANTIFUNGAL+INSTRUMENT_AIR+TUBE+
                   NEUT+RBC,family=binomial(link = "logit"),data[-m,],control = list(maxit=100))   #将m中下标集对应的行设置为测试集，利用训练集进行建模
  #####测试集
  pre<-predict(fit_model,data[-m,],type="response")
  pre_result<-ifelse(pre>0.477,1,0)
  z0=table(data[-m,1],pre_result)
  perfor_test[i,1:4]=Performance_iter(z0,beta)
  #计算AUC值
  pr<-prediction(pre,data[-m,1])
  perfor_test[i,5]=performance(pr,measure = 'auc')@y.values[[1]]
  
  #####验证集
  pre<-predict(fit_model,data[m,],type="response")
  pre_result<-ifelse(pre>0.477,1,0)
  z0=table(data[m,1],pre_result)
  perfor[i,1:4]=Performance_iter(z0,beta)
  kap[i]<-classAgreement(z0)$kappa
  E[i]<-sum(data[m,D]!=pre_result)/n1
  #计算AUC值
  pr<-prediction(pre,data[m,1])
  perfor[i,5]=performance(pr,measure = 'auc')@y.values[[1]]
  t<-length(performance(pr,'fpr','tpr')@x.values[[1]])
  tpr[i,1:t]=performance(pr,'fpr','tpr')@x.values[[1]]
  fpr[i,1:t]=performance(pr,'fpr','tpr')@y.values[[1]]
} 

##################绘制ROC曲线############
tpr_mean=fpr_mean=rep(NA,90)
tpr_sd=fpr_sd=nn=rep(NA,90)
for(j in 1:90){
  tpr_mean[j]=mean(na.omit(tpr[,j]))
  tpr_sd[j]=sd(na.omit(tpr[,j]))
  fpr_mean[j]=mean(na.omit(fpr[,j]))
  fpr_sd[j]=sd(na.omit(fpr[,j]))
  nn[j]=length(na.omit(tpr[,j]))
}

roc_mean2<-mean(perfor[,5])
tpr_mean2=tpr_mean;fpr_mean2=fpr_mean
#####添加均值的置信区间
tpr_upper=tpr_lower=rep(NA,length(na.omit(tpr_mean)))
for(i in 1:length(na.omit(tpr_mean))){
  tpr_upper[i]<-min(tpr_mean[i]+1.96*tpr_sd[i]/sqrt(nn[i]),1)
  tpr_lower[i]<-max(tpr_mean[i]-1.96*fpr_sd[i]/sqrt(nn[i]),0)
}

xx<-c(fpr_mean,rev(fpr_mean))
yy<-c(c(tpr_upper),rev(c(tpr_lower)))
plot(xx,yy,type='n',xlab='False positive rate',ylab="True positive rate",main='ROC curve of logisitic regression with 7 variables')
polygon(xx,yy,col='#DBDBDB',border ='#DBDBDB' )

#画一个对角线图
#画一个对角线图
x=y=seq(0,1,by=0.1)
lines(x,y,type='l',lty=2,lwd=2,col='grey')
for(i in 1:10){
  lines(fpr[i,],tpr[i,],type='l',col='grey')
}
lines(fpr_mean,tpr_mean,type='l',lty=1,lwd=2,col="red")
legend("center",cex=1.2,paste0(round(roc_mean2,3),' [',round(t.test(perfor[,5])$conf.int[1],3),' , ',round(t.test(perfor[,5])$conf.int[2],3),']'))

##评价指标
E;mean(E);sd(E)
kap;mean(kap);sd(E)
result<-paste("训练集结果如下：",
              "\nRecall=",round(colMeans(perfor_test)[1],3),"(",round(t.test(perfor_test[,1])$conf.int[1],3),",",round(t.test(perfor_test[,1])$conf.int[2],3),")",
              "\nSpecificity=",round(colMeans(perfor_test)[2],3),"(",round(t.test(perfor_test[,2])$conf.int[1],3),",",round(t.test(perfor_test[,2])$conf.int[2],3),")",
              "\nAccuracy=",round(colMeans(perfor_test)[3],3),"(",round(t.test(perfor_test[,3])$conf.int[1],3),",",round(t.test(perfor_test[,3])$conf.int[2],3),")",  
              "\nprecision=",round(colMeans(perfor_test)[4],3),"(",round(t.test(perfor_test[,4])$conf.int[1],3),",",round(t.test(perfor_test[,4])$conf.int[2],3),")",  
              "\nAUC=",round(colMeans(perfor_test)[5],3),"(",round(t.test(perfor_test[,5])$conf.int[1],3),",",round(t.test(perfor_test[,5])$conf.int[2],3),")","\n",sep="")   #定义Gmean为少数分类精确率和多数分类精确率的集合评价
cat(result)

result<-paste("验证集结果如下：",
              "\nRecall=",round(colMeans(perfor)[1],3),"(",round(t.test(perfor[,1])$conf.int[1],3),",",round(t.test(perfor[,1])$conf.int[2],3),")",
              "\nSpecificity=",round(colMeans(perfor)[2],3),"(",round(t.test(perfor[,2])$conf.int[1],3),",",round(t.test(perfor[,2])$conf.int[2],3),")",
              "\nAccuracy=",round(colMeans(perfor)[3],3),"(",round(t.test(perfor[,3])$conf.int[1],3),",",round(t.test(perfor[,3])$conf.int[2],3),")",  
              "\nprecision=",round(colMeans(perfor)[4],3),"(",round(t.test(perfor[,4])$conf.int[1],3),",",round(t.test(perfor[,4])$conf.int[2],3),")",  
              "\nAUC=",round(colMeans(perfor)[5],3),"(",round(t.test(perfor[,5])$conf.int[1],3),",",round(t.test(perfor[,5])$conf.int[2],3),")","\n",sep="")   #定义Gmean为少数分类精确率和多数分类精确率的集合评价
cat(result)



#################################################
#################xgboost#########################
#################################################

###构建xgboost模型
library(xgboost)
library(Matrix)

##################利用全数据集构建应用模型
set.seed(100)
trainset2<-as(as.matrix(data[,c(2:ncol(data))]),"dgCMatrix")
trainset3<-as.numeric(as.character(data[,1]))
trainset4<-list(data=trainset2,label=trainset3)
dtrain<-xgb.DMatrix(data=trainset4$data,label=trainset4$label)

##参数设置
param=list(
  objective="binary:logistic",
  booster='gbtree',
  max_depth=6,
  subsample=0.6,
  colsample_bytree=0.75,
  min_child_weight=1,
  alpha=1,
  gamma=0.1
  #early_stopping_rounds = 50
)
nround=10

set.seed(1234)
xgb_model <- xgb.train(params = param,data=dtrain,nrounds = nround,nthread=2)
save(xgb_model,file="xgb_model.RData")

#model<-xgb.dump(xgb_model,with_stats = T)   #显示计算过程，查看树结构
#library(DiagrammeR)
#xgb.plot.tree(model = xgb_model,feature_names = names)

names<-dimnames((data.matrix(data[,c(2:ncol(data))])))[[2]]
importance_matrix<-xgb.importance(names,model=xgb_model)
write.csv(importance_matrix,file='F:/侵袭性真菌病项目—罗燕萍主任/2-第二次提取数据/3-数据分析/2-分析结果/importance_shap.csv')

xgb.plot.importance(importance_matrix[,],top_n=10,left_margin = 20,rel_to_first = TRUE,cex.lab=1.5,cex=1)
###根据重要性图得到top20比较重要的变量。


####Visualizing the SHAP feature contribution to prediction dependencies on feature value.
col<-rgb(0,0,1,0.5)
xgb.plot.shap(trainset2,model=xgb_model,features = 'TUBE',col = col, pch = 16,plot_loess = FALSE) 
abline(h=0,col='red',lty=2,cex=4)
shap_value_tube<-xgb.plot.shap(trainset2,model=xgb_model,features = 'TUBE',plot=FALSE)
#shap_value_tube<-as.character(shap_value_tube)
#write(shap_value_tube,file='F:/侵袭性真菌病项目—罗燕萍主任/2-第二次提取数据/3-数据分析/2-分析结果/shap_value_tube.txt')
#tt<-read.table('F:/侵袭性真菌病项目—罗燕萍主任/2-第二次提取数据/3-数据分析/2-分析结果/shap_value_tube.txt')
tt<-as.matrix(cbind(shap_value_tube$data,shap_value_tube$shap_contrib))
write.csv(tt,file='F:/侵袭性真菌病项目—罗燕萍主任/2-第二次提取数据/3-数据分析/2-分析结果/shap_value_tube.csv')

xgb.plot.shap(trainset2,model=xgb_model,features = 'pre_hos_days',col = col, pch = 16,plot_loess = FALSE) 
abline(h=0,col='red',lty=2,cex=4)
abline(v=10,col='red',lty=2,cex=4)
text(6,0.05,"(10,0)")
shap_value_predays<-xgb.plot.shap(trainset2,model=xgb_model,features = 'pre_hos_days' ,plot=FALSE)
tt<-as.matrix(cbind(t(t(shap_value_predays$data)),shap_value_predays$shap_contrib))
write.csv(tt,file='F:/侵袭性真菌病项目—罗燕萍主任/2-第二次提取数据/3-数据分析/2-分析结果/shap_value_predays.csv')

xgb.plot.shap(trainset2,model=xgb_model,features = 'NEUT',col = col, pch = 16,plot_loess = FALSE) 
abline(h=0,col='red',lty=2,cex=4)
abline(v=0.71,col='red',lty=2,cex=4)
text(0.65,0.05,"(0.71,0)")
shap_value_NEUT<-xgb.plot.shap(trainset2,model=xgb_model,features = 'NEUT',plot=FALSE)
tt<-as.matrix(cbind(t(t(shap_value_NEUT$data)),shap_value_NEUT$shap_contrib))
write.csv(cbind(shap_value_NEUT$data@x,shap_value_NEUT$shap_contrib),file='F:/侵袭性真菌病项目—罗燕萍主任/2-第二次提取数据/3-数据分析/2-分析结果/shap_value_NEUT.csv')

xgb.plot.shap(trainset2,model=xgb_model,features = 'DPP',col = col, pch = 16,plot_loess = FALSE) 
abline(h=0,col='red',lty=2,cex=4)
abline(v=38,col='red',lty=2,cex=4)
text(30,0.05,"(38,0)")
shap_value_DPP<-xgb.plot.shap(trainset2,model=xgb_model,features = 'DPP',plot=FALSE)
tt<-as.matrix(cbind(t(t(shap_value_DPP$data)),shap_value_DPP$shap_contrib))
write.csv(cbind(shap_value_DPP$data@x,shap_value_DPP$shap_contrib),file='F:/侵袭性真菌病项目—罗燕萍主任/2-第二次提取数据/3-数据分析/2-分析结果/shap_value_DPP.csv')

xgb.plot.shap(trainset2,model=xgb_model,features = 'INSTRUMENT_AIR',col = col, pch = 16,plot_loess = FALSE) 
abline(h=0,col='red',lty=2,cex=4)
shap_value_INSTRUMENT<-xgb.plot.shap(trainset2,model=xgb_model,features = 'INSTRUMENT_AIR',plot=FALSE)
tt<-as.matrix(cbind(t(t(shap_value_INSTRUMENT$data)),shap_value_INSTRUMENT$shap_contrib))
write.csv(cbind(shap_value_INSTRUMENT$data@x,shap_value_INSTRUMENT$shap_contrib),file='F:/侵袭性真菌病项目—罗燕萍主任/2-第二次提取数据/3-数据分析/2-分析结果/shap_value_INSTRUMENT.csv')

xgb.plot.shap(trainset2,model=xgb_model,features = 'RBC',col = col, pch = 16,plot_loess = FALSE) 
abline(h=0,col='red',lty=2,cex=4)
abline(v=3.7,col='red',lty=2,cex=4)
text(3.4,-0.05,"(3.7,0)",cex=1.2)
shap_value_RBC<-xgb.plot.shap(trainset2,model=xgb_model,features = 'RBC',plot=FALSE)
tt<-as.matrix(cbind(t(t(shap_value_RBC$data)),shap_value_RBC$shap_contrib))
write.csv(cbind(shap_value_RBC$data@x,shap_value_RBC$shap_contrib),file='F:/侵袭性真菌病项目—罗燕萍主任/2-第二次提取数据/3-数据分析/2-分析结果/shap_value_RBC.csv')

xgb.plot.shap(trainset2,model=xgb_model,features = 'WBC',col = col, pch = 16,plot_loess = FALSE) 
abline(h=0,col='red',lty=2,cex=4)
abline(v=5,col='red',lty=2,cex=4)
#abline(v=10,col='red',lty=2,cex=4)
text(2,0.02,"(5,0)",cex=1)
shap_value_WBC<-xgb.plot.shap(trainset2,model=xgb_model,features = 'WBC',plot=FALSE)
tt<-as.matrix(cbind(t(t(shap_value_WBC$data)),shap_value_WBC$shap_contrib))
write.csv(cbind(shap_value_WBC$data@x,shap_value_WBC$shap_contrib),file='F:/侵袭性真菌病项目—罗燕萍主任/2-第二次提取数据/3-数据分析/2-分析结果/shap_value_WBC.csv')

xgb.plot.shap(trainset2,model=xgb_model,features = 'age',col = col, pch = 16,plot_loess = FALSE) 
abline(h=0,col='red',lty=2,cex=4)
abline(v=60,col='red',lty=2,cex=4)
#abline(v=10,col='red',lty=2,cex=4)
text(55,0.05,"(60,0)",cex=1)
shap_value_age<-xgb.plot.shap(trainset2,model=xgb_model,features = 'age',plot=FALSE)
tt<-as.matrix(cbind(t(t(shap_value_age$data)),shap_value_age$shap_contrib))
write.csv(cbind(shap_value_age$data@x,shap_value_age$shap_contrib),file='F:/侵袭性真菌病项目—罗燕萍主任/2-第二次提取数据/3-数据分析/2-分析结果/shap_value_age.csv')

xgb.plot.shap(trainset2,model=xgb_model,features = 'MPP',col = col, pch = 16,plot_loess = FALSE) 
abline(h=0,col='red',lty=2,cex=4)
abline(v=80,col='red',lty=2,cex=4)
abline(v=96,col='red',lty=2,cex=4)
text(85,0.02,"(80,0)",cex=1)
text(93,0,"(96,0)",cex=1)
shap_value_MPP<-xgb.plot.shap(trainset2,model=xgb_model,features = 'MPP',plot=FALSE)
tt<-as.matrix(cbind(t(t(shap_value_MPP$data)),shap_value_MPP$shap_contrib))
write.csv(cbind(shap_value_MPP$data@x,shap_value_MPP$shap_contrib),file='F:/侵袭性真菌病项目—罗燕萍主任/2-第二次提取数据/3-数据分析/2-分析结果/shap_value_MPP.csv')

xgb.plot.shap(trainset2,model=xgb_model,features = 'FG',col = col, pch = 16,plot_loess = FALSE) 
shap_value_FG<-xgb.plot.shap(trainset2,model=xgb_model,features = 'FG',plot=FALSE)
tt<-as.matrix(cbind(t(t(shap_value_FG$data)),shap_value_FG$shap_contrib))
write.csv(cbind(shap_value_FG$data@x,shap_value_FG$shap_contrib),file='F:/侵袭性真菌病项目—罗燕萍主任/2-第二次提取数据/3-数据分析/2-分析结果/shap_value_FG.csv')


##利用ROC选取最佳阈值
pre<-predict(xgb_model,trainset2)
modelroc=roc(trainset3,pre)
plot(modelroc,print.auc=TRUE,auc.polygon=TRUE,grid=c(0.1,0.2),grid.col=c("green","red"),max.auc.polygon=TRUE,
     auc.polygon.col="skyblue",print.thres=TRUE)   #输出结果为0.473

############十折交叉验证
beta=1
E=rep(0,Z);kap=rep(0,Z)
tpr=fpr=matrix(NA,nrow=Z,ncol=90);perfor=matrix(0,nrow=Z,ncol=5)
perfor_test=matrix(0,nrow=Z,ncol=5)
source("F:/mycode/Performance_iter.R")

for(i in 1:Z){
  m<-mm[[i]]   #另m为第i类的下标集
  n1<-length(m)
  ##训练集
  traindata2<-as(as.matrix(data[-m,c(2:39)]),"dgCMatrix")
  traindata3<-as.numeric(as.character(data[-m,1]))
  traindata4<-list(data=traindata2,label=traindata3)
  dtrain<-xgb.DMatrix(data=traindata4$data,label=traindata4$label)
  #验证集
  varifydata2<-as(as.matrix(data[m,c(2:39)]),"dgCMatrix")
  varifydata3<-as.numeric(as.character(data[m,1]))
  varifydata4<-list(data=varifydata2,label=varifydata3)
  dvarify<-xgb.DMatrix(data=varifydata4$data,label=varifydata4$label)
  #构建模型
  fit_model <- xgb.train(params = param,data=dtrain,nrounds = nround,nthread=2)
  #####训练集
  z<-predict(fit_model,traindata2)
  pre_result<-ifelse(z>0.473,1,0)
  confusion<-table(true=traindata3,pre=pre_result)
  perfor_test[i,1:4]=Performance_iter(confusion,beta)
  #计算AUC值
  pr<-prediction(z,traindata3)
  perfor_test[i,5]=performance(pr,measure = 'auc')@y.values[[1]]
  
  
  #####验证集
  z<-predict(fit_model,varifydata2)
  pre_result<-ifelse(z>0.473,1,0)
  confusion<-table(true=varifydata3,pre=pre_result)
  perfor[i,1:4]=Performance_iter(confusion,beta)
  kap[i]<-classAgreement(confusion)$kappa
  E[i]<-1-sum(diag(confusion))/n1
  #计算AUC值
  pr<-prediction(z,varifydata3)
  perfor[i,5]=performance(pr,measure = 'auc')@y.values[[1]]
  t<-length(performance(pr,'fpr','tpr')@x.values[[1]])
  tpr[i,1:t]=performance(pr,'fpr','tpr')@x.values[[1]]
  fpr[i,1:t]=performance(pr,'fpr','tpr')@y.values[[1]]
} 

##################绘制ROC曲线############
tpr_mean=fpr_mean=rep(NA,90)
tpr_sd=fpr_sd=nn=rep(NA,90)
for(j in 1:90){
  tpr_mean[j]=mean(na.omit(tpr[,j]))
  tpr_sd[j]=sd(na.omit(tpr[,j]))
  fpr_mean[j]=mean(na.omit(fpr[,j]))
  fpr_sd[j]=sd(na.omit(fpr[,j]))
  nn[j]=length(na.omit(tpr[,j]))
}

roc_mean3<-mean(perfor[,5])
tpr_mean3=tpr_mean;fpr_mean3=fpr_mean
#####添加均值的置信区间
tpr_upper=tpr_lower=rep(NA,length(na.omit(tpr_mean)))
for(i in 1:length(na.omit(tpr_mean))){
  tpr_upper[i]<-min(tpr_mean[i]+1.96*tpr_sd[i]/sqrt(nn[i]),1)
  tpr_lower[i]<-max(tpr_mean[i]-1.96*fpr_sd[i]/sqrt(nn[i]),0)
}

xx<-c(fpr_mean,rev(fpr_mean))
yy<-c(c(tpr_upper),rev(c(tpr_lower)))
plot(xx,yy,type='n',xlab='False positive rate',ylab="True positive rate",main='ROC curve of XGBoost model')
polygon(xx,yy,col='#DBDBDB',border ='#DBDBDB' )

#画一个对角线图
#画一个对角线图
x=y=seq(0,1,by=0.1)
lines(x,y,type='l',lty=2,lwd=2,col='grey')
for(i in 1:10){
  lines(fpr[i,],tpr[i,],type='l',col='grey')
}
lines(fpr_mean,tpr_mean,type='l',lty=1,lwd=2,col="red")
legend("center",cex=1.2,paste0(round(roc_mean3,3),' [',round(t.test(perfor[,5])$conf.int[1],3),' , ',round(t.test(perfor[,5])$conf.int[2],3),']'))


##评价指标
E;mean(E);sqrt(var(E));t.test(E)$conf.int
kap;mean(kap);sqrt(var(kap));t.test(kap)$conf.int
result<-paste("训练集结果如下：",
              "\nRecall=",round(colMeans(perfor_test)[1],3),"(",round(t.test(perfor_test[,1])$conf.int[1],3),",",round(t.test(perfor_test[,1])$conf.int[2],3),")",
              "\nSpecificity=",round(colMeans(perfor_test)[2],3),"(",round(t.test(perfor_test[,2])$conf.int[1],3),",",round(t.test(perfor_test[,2])$conf.int[2],3),")",
              "\nAccuracy=",round(colMeans(perfor_test)[3],3),"(",round(t.test(perfor_test[,3])$conf.int[1],3),",",round(t.test(perfor_test[,3])$conf.int[2],3),")",  
              "\nprecision=",round(colMeans(perfor_test)[4],3),"(",round(t.test(perfor_test[,4])$conf.int[1],3),",",round(t.test(perfor_test[,4])$conf.int[2],3),")",  
              "\nAUC=",round(colMeans(perfor_test)[5],3),"(",round(t.test(perfor_test[,5])$conf.int[1],3),",",round(t.test(perfor_test[,5])$conf.int[2],3),")","\n",sep="")   #定义Gmean为少数分类精确率和多数分类精确率的集合评价
cat(result)

result<-paste("验证集结果如下：",
              "\nRecall=",round(colMeans(perfor)[1],3),"(",round(t.test(perfor[,1])$conf.int[1],3),",",round(t.test(perfor[,1])$conf.int[2],3),")",
              "\nSpecificity=",round(colMeans(perfor)[2],3),"(",round(t.test(perfor[,2])$conf.int[1],3),",",round(t.test(perfor[,2])$conf.int[2],3),")",
              "\nAccuracy=",round(colMeans(perfor)[3],3),"(",round(t.test(perfor[,3])$conf.int[1],3),",",round(t.test(perfor[,3])$conf.int[2],3),")",  
              "\nprecision=",round(colMeans(perfor)[4],3),"(",round(t.test(perfor[,4])$conf.int[1],3),",",round(t.test(perfor[,4])$conf.int[2],3),")",  
              "\nAUC=",round(colMeans(perfor)[5],3),"(",round(t.test(perfor[,5])$conf.int[1],3),",",round(t.test(perfor[,5])$conf.int[2],3),")","\n",sep="")   #定义Gmean为少数分类精确率和多数分类精确率的集合评价
cat(result)


###########################################################################
#################根据提取的top10个变量进行建模##############################

indel<-c(1,2,5,10,11,18,19,20,23,25,26)
data_top10<-data[,indel]


set.seed(100)
trainset2<-as(as.matrix(data_top10[,c(2:ncol(data_top10))]),"dgCMatrix")
trainset3<-as.numeric(as.character(data_top10[,1]))
trainset4<-list(data=trainset2,label=trainset3)
dtrain<-xgb.DMatrix(data=trainset4$data,label=trainset4$label)

param<-list(objective="binary:logistic",
            booster='gbtree',
            max_depth=4,
            subsample=0.55,
            colsample_bytree=0.8,
            min_child_weight=4, 
            alpha=1,
            gamma=0.4
            #early_stopping_rounds = 5
)
nround=10

#set.seed(1234)
xgb_model_top10 <- xgb.train(params = param,data=dtrain,nrounds = nround,nthread=2,seed=1234)
save(xgb_model_top10,file="xgb_model_top10.RData")

##利用ROC选取最佳阈值
pre<-predict(xgb_model_top10,trainset2)
modelroc=roc(data_top10$XUEPY_RESULT,pre)
plot(modelroc,print.auc=TRUE,auc.polygon=TRUE,grid=c(0.1,0.2),grid.col=c("green","red"),max.auc.polygon=TRUE,
     auc.polygon.col="skyblue",print.thres=TRUE)   #输出结果为0.475

############十折交叉验证
beta=1
E=rep(0,Z);kap=rep(0,Z)
tpr=fpr=matrix(NA,nrow=Z,ncol=90);perfor=matrix(0,nrow=Z,ncol=5)
perfor_test=matrix(0,nrow=Z,ncol=5)
source("F:/mycode/Performance_iter.R")

for(i in 1:Z){
  m<-mm[[i]]   #另m为第i类的下标集
  n1<-length(m)
  ##训练集
  traindata2<-as(as.matrix(data_top10[-m,c(2:11)]),"dgCMatrix")
  traindata3<-as.numeric(as.character(data_top10[-m,1]))
  traindata4<-list(data=traindata2,label=traindata3)
  dtrain<-xgb.DMatrix(data=traindata4$data,label=traindata4$label)
  #验证集
  varifydata2<-as(as.matrix(data_top10[m,c(2:11)]),"dgCMatrix")
  varifydata3<-as.numeric(as.character(data_top10[m,1]))
  varifydata4<-list(data=varifydata2,label=varifydata3)
  dvarify<-xgb.DMatrix(data=varifydata4$data,label=varifydata4$label)
  #构建模型
  #set.seed(1234)
  fit_model <- xgb.train(params = param,data=dtrain,nrounds = nround,nthread=2)
  #####训练集
  z<-predict(fit_model,traindata2)
  pre_result<-ifelse(z>0.475,1,0)
  confusion<-table(true=traindata3,pre=pre_result)
  perfor_test[i,1:4]=Performance_iter(confusion,beta)
  #计算AUC值
  pr<-prediction(z,traindata3)
  perfor_test[i,5]=performance(pr,measure = 'auc')@y.values[[1]]
  
  
  ####验证集
  z<-predict(fit_model,varifydata2,type='response')
  pre_result<-ifelse(z>0.475,1,0)
  confusion<-table(true=varifydata3,pre=pre_result)
  perfor[i,1:4]=Performance_iter(confusion,beta)
  kap[i]<-classAgreement(confusion)$kappa
  E[i]<-1-sum(diag(confusion))/n1
  #计算AUC值
  pr<-prediction(z,varifydata3)
  perfor[i,5]=performance(pr,measure = 'auc')@y.values[[1]]
  t<-length(performance(pr,'fpr','tpr')@x.values[[1]])
  tpr[i,1:t]=performance(pr,'fpr','tpr')@x.values[[1]]
  fpr[i,1:t]=performance(pr,'fpr','tpr')@y.values[[1]]
} 

##################绘制ROC曲线############
tpr_mean=fpr_mean=rep(NA,90)
tpr_sd=fpr_sd=nn=rep(NA,90)
for(j in 1:90){
  tpr_mean[j]=mean(na.omit(tpr[,j]))
  tpr_sd[j]=sd(na.omit(tpr[,j]))
  fpr_mean[j]=mean(na.omit(fpr[,j]))
  fpr_sd[j]=sd(na.omit(fpr[,j]))
  nn[j]=length(na.omit(tpr[,j]))
}

roc_mean4<-mean(perfor[,5])
tpr_mean4=tpr_mean;fpr_mean4=fpr_mean
#####添加均值的置信区间
tpr_upper=tpr_lower=rep(NA,length(na.omit(tpr_mean)))
for(i in 1:length(na.omit(tpr_mean))){
  tpr_upper[i]<-min(tpr_mean[i]+1.96*tpr_sd[i]/sqrt(nn[i]),1)
  tpr_lower[i]<-max(tpr_mean[i]-1.96*fpr_sd[i]/sqrt(nn[i]),0)
}

xx<-c(fpr_mean,rev(fpr_mean))
yy<-c(c(tpr_upper),rev(c(tpr_lower)))
plot(xx,yy,type='n',xlab='False positive rate',ylab="True positive rate",main='ROC curve of XGBoost model with 10 variables')
polygon(xx,yy,col='#DBDBDB',border ='#DBDBDB' )

#画一个对角线图
x=y=seq(0,1,by=0.1)
lines(x,y,type='l',lty=2,lwd=2,col='grey')
for(i in 1:10){
  lines(fpr[i,],tpr[i,],type='l',col='grey')
}
lines(fpr_mean,tpr_mean,type='l',lty=1,lwd=2,col="red")
legend("center",cex=1.2,paste0(round(roc_mean4,3),' [',round(t.test(perfor[,5])$conf.int[1],3),' , ',round(t.test(perfor[,5])$conf.int[2],3),']'))


##评价指标
E;mean(E);sqrt(var(E));t.test(E)$conf.int
kap;mean(kap);sqrt(var(kap));t.test(kap)$conf.int
result<-paste("训练集结果如下：",
              "\nRecall=",round(colMeans(perfor_test)[1],3),"(",round(t.test(perfor_test[,1])$conf.int[1],3),",",round(t.test(perfor_test[,1])$conf.int[2],3),")",
              "\nSpecificity=",round(colMeans(perfor_test)[2],3),"(",round(t.test(perfor_test[,2])$conf.int[1],3),",",round(t.test(perfor_test[,2])$conf.int[2],3),")",
              "\nAccuracy=",round(colMeans(perfor_test)[3],3),"(",round(t.test(perfor_test[,3])$conf.int[1],3),",",round(t.test(perfor_test[,3])$conf.int[2],3),")",  
              "\nprecision=",round(colMeans(perfor_test)[4],3),"(",round(t.test(perfor_test[,4])$conf.int[1],3),",",round(t.test(perfor_test[,4])$conf.int[2],3),")",  
              "\nAUC=",round(colMeans(perfor_test)[5],3),"(",round(t.test(perfor_test[,5])$conf.int[1],3),",",round(t.test(perfor_test[,5])$conf.int[2],3),")","\n",sep="")   #定义Gmean为少数分类精确率和多数分类精确率的集合评价
cat(result)

result<-paste("验证集结果如下：",
              "\nRecall=",round(colMeans(perfor)[1],3),"(",round(t.test(perfor[,1])$conf.int[1],3),",",round(t.test(perfor[,1])$conf.int[2],3),")",
              "\nSpecificity=",round(colMeans(perfor)[2],3),"(",round(t.test(perfor[,2])$conf.int[1],3),",",round(t.test(perfor[,2])$conf.int[2],3),")",
              "\nAccuracy=",round(colMeans(perfor)[3],3),"(",round(t.test(perfor[,3])$conf.int[1],3),",",round(t.test(perfor[,3])$conf.int[2],3),")",  
              "\nprecision=",round(colMeans(perfor)[4],3),"(",round(t.test(perfor[,4])$conf.int[1],3),",",round(t.test(perfor[,4])$conf.int[2],3),")",  
              "\nAUC=",round(colMeans(perfor)[5],3),"(",round(t.test(perfor[,5])$conf.int[1],3),",",round(t.test(perfor[,5])$conf.int[2],3),")","\n",sep="")   #定义Gmean为少数分类精确率和多数分类精确率的集合评价
cat(result)



##########################################
######绘制四种方法的平均值结果############
x=y=seq(0,1,by=0.1)
plot(x,y,type='l',lty=2,lwd=2,col='grey',xlab='False positive rate',ylab="True positive rate",main='ROC curve')
#lines(fpr_mean1,tpr_mean1,type='l',lty=1,lwd=2,col=1)
lines(fpr_mean2,tpr_mean2,type='l',lty=1,lwd=2,col=2)
lines(fpr_mean3,tpr_mean3,type='l',lty=1,lwd=2,col=3)
lines(fpr_mean4,tpr_mean4,type='l',lty=1,lwd=2,col=4)
legend("bottomright",lty=rep(1,times=3),lwd=rep(2,times=3),col=c(1:3),
       c(paste0('logisitic(area=',round(roc_mean2,3),')'),
         paste0('xgboost(area=',round(roc_mean3,3),')'),
         paste0('xgboost+top10(area=',round(roc_mean4,3),')')))


