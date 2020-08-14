#dataset
WDBC

#data normalizing
normalise<-function(X){
  return((X-min(X))/(max(X)-min(X)))
}

for(i in 3:ncol(WDBC)){
  WDBC[,i]<-normalise(WDBC[,i])
}
print(WDBC)

#######################################FEATURE SELECTION##############################################

#store subset of features from all feature selection methods
feature_selection_results<-list()

############################################FS-RRC####################################################
#entropy, mutual information
library(DescTools)

#initialization
Slist<-c()
Slistnames<-c()
Scomp<-c()
Sbest<-c()

#symmetric uncertainty
symmetricUncertainty<-function(X,Y){
  return((2*MutInf(X,Y,base=2))/(Entropy(table(X),y=NULL,base=2)+Entropy(table(Y),y=NULL,base=2)))
}

#finding symmetric uncertainty for all the features
for(i in 3:ncol(WDBC)){
  su<-symmetricUncertainty(WDBC[,i],WDBC[,2])
  if(su>0){
    Slist[intToUtf8(i+62)]<-su
  }
}
print(Slist)

#sorting based on symmetric uncertainty
ordered_Slist<-order(as.numeric(Slist),decreasing=TRUE)
Slist<-Slist[ordered_Slist]
print(Slist)

SlistNew=Slist
print(SlistNew)

f<-Slist[1]
print(f)

#get next element in a list
getNextElement<-function(L,f){
  for(i in 1:length(L)){
    if(L[i]==f){
      if(i+1<=length(L))
        return(L[i+1])
      else
        return(NULL)
    }
  }
  return(NULL)
}

repeat {
  fNew<-getNextElement(Slist,f)
  repeat{
    if(symmetricUncertainty(WDBC[,(utf8ToInt(names(fNew))-64)],WDBC[,(utf8ToInt(names(f))-64)])>=
       symmetricUncertainty(WDBC[,(utf8ToInt(names(fNew))-64)],WDBC[,2]))
      Slist=Slist[Slist!=fNew]
    fNew<-getNextElement(Slist,fNew)
    if(is.null(fNew)==TRUE){
      break
    }
  }
  score<-c()
  for(fNew in names(SlistNew)){
    if(fNew!=f){
      score[fNew]<-(MutInf(WDBC[,(utf8ToInt(names(f))-64)],WDBC[,(utf8ToInt(fNew)-64)],base=2))
                            -(MutInf(WDBC[,(utf8ToInt(names(f))-64)],WDBC[,2],base=2))
    }
  }
  Scomp<-c(Scomp,names(score)[max(score)])
  f=getNextElement(Slist,f)
  if(is.null(f)==TRUE){
    break
  }
}
library(sets)
Slistnames<-names(Slist)
Sbest<-union(Slistnames,Scomp)
Sbest<-unique(Sbest)
SbestNew<-c()
for(i in Sbest){
  SbestNew<-c(SbestNew,utf8ToInt(i)-62)
}
print(SbestNew)
feature_selection_results[["FS-RRC"]]<-paste("V",SbestNew,sep = "")

#################################correlation (removing redundant features)##############################
# Subsetting the data and selecting only required variables
df <- WDBC[,3:ncol(WDBC)]
#df[,1]<-factor(gsub('B', '0', df[,1]))
#df[,1]<-factor(gsub('M', '1', df[,1]))
#df[,1]<-as.numeric(df[,1])

# Using corr function to generate correlation matrix
print(cor(df))

# Building correplot to visualize the correlartion matrix
library(corrplot)
corrplot(cor(df), method="number", is.corr=FALSE)

library(caret)
highlyCorrelated <- findCorrelation(cor(df), cutoff=0.5)
print(highlyCorrelated)
library(sets) 
lessCorrelated<-setdiff(c(1:30),highlyCorrelated)
feature_selection_results[["correlation"]]<-paste("V",lessCorrelated+2,sep = "")

############################Learning Vector Qunatization (rank by importance)###############################
library(mlbench)
library(caret)

# prepare training scheme
control <- trainControl(method="repeatedcv", number=10, repeats=3)
# train the model
model <- train(V2~., data=WDBC[,2:ncol(WDBC)], method="lvq", preProcess="scale", trControl=control)
# estimate variable importance
importance <- varImp(model, scale=FALSE)
# summarize importance
print(importance)
print(rownames(Sort(importance$importance,decreasing = TRUE)[1:7,]))
# plot importance
plot(importance)
feature_selection_results[["LVQ"]]<-rownames(Sort(importance$importance,decreasing = TRUE)[1:7,])

####################################Recursive Feature Elimination####################################
# load the library
library(mlbench)
library(caret)
# define the control using a random forest selection function
control <- rfeControl(functions=rfFuncs, method="cv", number=10)
# run the RFE algorithm
results <- rfe(WDBC[,3:ncol(WDBC)], WDBC[,2], sizes=c(3:ncol(WDBC)), rfeControl=control)
# summarize the results
print(results)
# list the chosen features
predictors(results)
# plot the results
plot(results, type=c("g", "o"))
feature_selection_results[["RFE"]]<-results$optVariables

##########################################Hypothesis Testing##########################################
df<-WDBC[,2:ncol(WDBC)]
df[,1]<-factor(gsub('B', '0', df[,1]))
df[,1]<-factor(gsub('M', '1', df[,1]))
df[,1]<-as.numeric(df[,1])

t.statistic<-c()
plotT<-c()
names(t.statistic)<-c()
for(i in 2:ncol(df)){
  print(t.test(df[,1], df[,i]))
  t.statistic<-c(t.statistic,t.test(df[,1], df[,i])$statistic)
  plotT<-c(plotT,t.test(df[,1], df[,i])$statistic)
  names(t.statistic)[i-1]<-i
}

print(t.statistic)
print(t.statistic[t.statistic>56])
plot((3:ncol(WDBC)),plotT,col=ifelse(plotT >56,'red','green'),xlab="features",ylab="t-statistic",pch=19)
feature_selection_results[["hypothesisTest"]]<-paste("V",as.numeric(names(t.statistic[t.statistic>56]))+1,sep="")

##############################################LASSO##################################################
library(glmnet)
`%ni%`<-Negate('%in%')

df<-WDBC[,2:ncol(WDBC)]
df[,1]<-factor(gsub('B', '0', df[,1]))
df[,1]<-factor(gsub('M', '1', df[,1]))
df[,1]<-as.numeric(df[,1])

x<-model.matrix(V2~.,data=df)
x=x[,-1]

glmnet1<-cv.glmnet(x=x,y=df[,2],type.measure='mse',nfolds=5,alpha=.5)

c<-coef(glmnet1,s='lambda.min',exact=TRUE)
inds<-which(c!=0)
variables<-row.names(c)[inds]
variables<-variables[variables %ni% '(Intercept)']
print(variables)
feature_selection_results[["LASSO"]]<-variables


feature_selection_results

#############################################CLASSIFIERS#############################################

#store performance of classifiers
accuracy<-list()

#################################################knn#################################################

library(class)
temp<-c()
for(i in 1:6){
  ##run knn function
  pr <- knn(WDBC[1:400,feature_selection_results[[i]]],WDBC[400:nrow(WDBC),feature_selection_results[[i]]],cl=WDBC[1:400,2],k=2)
  
  ##create confusion matrix
  tab <- table(pr,WDBC[400:nrow(WDBC),2])
  
  ##this function divides the correct predictions by total number of predictions that tell us how accurate teh model is.
  accuracyknn <- function(x){sum(diag(x)/(sum(rowSums(x)))) * 100}
  accuracyknn(tab)
  temp<-c(temp,accuracyknn(tab))
}

accuracy[["knn"]]<-temp

#################################################svm#################################################


# Fitting SVM to the Training set 
library(e1071) 
i=1
temp<-c()
for(i in 1:6){
  training_data<-cbind(WDBC[1:400,feature_selection_results[[i]]],WDBC[1:400,2:2])
  colnames(training_data)[ncol(training_data)]<-"V2"
  classifier = svm(formula =  V2~ ., 
                   data = training_data, 
                   type = 'C-classification', 
                   kernel = 'linear') 
  
  # Predicting the Test set results 
  y_pred = predict(classifier, newdata = WDBC[400:nrow(WDBC),feature_selection_results[[i]]])
  # Making the Confusion Matrix 
  cm = table(WDBC[400:nrow(WDBC),2:2], y_pred)
  temp<-c(temp,((cm[1,1]+cm[2,2])/(cm[1,1]+cm[2,2]+cm[1,2]+cm[2,1]))*100) 
}

accuracy[["svm"]]<-temp

##############################################Naive Bayes#################################################

#Loading the library
library(e1071)
i=1
temp<-c()
for(i in 1:6){
  #Fitting the Naive Bayes model
  training_data<-cbind(WDBC[1:400,feature_selection_results[[i]]],WDBC[1:400,2:2])
  colnames(training_data)[ncol(training_data)]<-"V2"
  Naive_Bayes_Model=naiveBayes(V2 ~., data=training_data)
  #Print the model summary
  print(Naive_Bayes_Model)
  
  #Prediction on the dataset
  NB_Predictions=predict(Naive_Bayes_Model,WDBC[400:nrow(WDBC),feature_selection_results[[i]]])
  #Confusion matrix to check accuracy
  cm = table(NB_Predictions,WDBC[400:nrow(WDBC),2:2])
  temp<-c(temp,((cm[1,1]+cm[2,2])/(cm[1,1]+cm[2,2]+cm[1,2]+cm[2,1]))*100) 
}
accuracy[["nb"]]<-temp


############################################comparing################################################

library("ggplot2")

dat <- data.frame(Accuracy = c(accuracy[[1]],accuracy[[2]],accuracy[[3]]),
                  Feature_Selection_Method= rep(c("FS-RRC","correlation","LVQ","RFE","hypothesisTest","LASSO"),3),
                  Classifier = c(rep("knn",6),rep("svm",6),rep("nb",6)))

ggplot(dat, aes(x=Feature_Selection_Method, y = Accuracy, group = Classifier,color = Classifier) )+
  geom_line() + labs(title = "FS-RRC vs Other FSMs")

