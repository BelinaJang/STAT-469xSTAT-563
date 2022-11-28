# Assignment #1-2: Logistic Regression, LDA, Classification tree and Elastic net

# The HIV Data

set.seed(0) #set the seed
library(MTPS)# data is stored in this package
library(tree)
data(HIV)
str(XX)
str(YY)

yBin <- as.matrix(YY)
cutoffs <- c(2,3,3,1.5,1.5) # cutoff value to be used to define drug resistance
for(ii in 1:5) yBin[,ii] <- (10^yBin[,ii] < cutoffs[ii])*1 
#smaller than cutoff -> 1:non-resistance - working
#bigger than cutoff -> 0:resistance - not working

y_df=as.data.frame(yBin)
names(y_df)[names(y_df) == '3TC'] <- 'Y3TC'
df=as.data.frame(XX)
#df = df + rnorm(228, 0, 0.01) # Add noise
drugs = list('ABC','Y3TC','AZT','D4T','DDI')

for (drug in drugs){
  cat("Drug: ", drug)
  # do this for each drug
  combined_data <- cbind(df, y_df[drug])
  
  # Cross Validation
  # Compare methods using the following 4 criteria: misclassification rate / precision/ recall / F1 score.
  
  nfold  = 5 #5-folds
  split = combined_data[drug] #maybe not necessary, can use ABC combined_data[ABC]
  
  #split data into resistance and non-resistance for stratified fold sampling
  #use stratified split for comparison of each drug. This ensures each fold has about the same proportion of resistance vs non-resistance.
  data.nres = combined_data[split==1,] #data where ABC is non-resistance
  data.res = combined_data[split==0,] #data where ABC is resistance
  
  # sample the index
  set.seed(0)
  idx.nres = sample(1:nfold, sum(split==1), replace=T) #non-resistance
  idx.res = sample(1:nfold, sum(split==0), replace=T) #resistance
  
  # One functions with different return values: precision, recall, misclassification rate, f1 score
  cv.accu = function(nfold, data.nres, data.res, idx.nres, idx.res, returntype) 
  {
    #copy code from above
    cM1= cM2 = cM3 = cM4 = matrix(0, 2, 2)
    
    for(ii in 1:nfold){
      cat("fold ", ii)
      data.train = rbind(data.nres[idx.nres!=ii,], data.res[idx.res!=ii,])
      data.test  = rbind(data.nres[idx.nres==ii,], data.res[idx.res==ii,])
      
      # logistic regression for drug
      # Put logistics first since it's causing an error so if there's an error then we discard this experiment
      model4 = glm(as.formula(paste(drug,"~.",sep="")), data=data.train, family="binomial")
      pred4 = try(predict(model4, newdata=data.test, type="response"))
      if(class(pred4)=="try-error") next
      cM4 = cM4 + table(pred4>0.5, data.test[,drug])
      
      # LDA for ABC
      model1 = lda(as.formula(paste(drug,"~.",sep="")), data=data.train)
      pred1 = predict(model1, newdata=data.test) # ASK this is not 0 1
      cM1 = cM1 + table(pred1$class, data.test[,drug])
      
      # DROP label for KNN
      #knn_train = subset(data.train, select= as.factor(paste("as.factor(-c(",drug,"))",sep="")))
      #knn_test = subset(data.test, select=as.factor(paste("as.factor(-c(",drug,"))",sep="")))
      
      # prediction for k=3
      #knn3 = knn(knn_train, knn_test, data.train[,drug], k=3) 
      #cMK3 = cMK3 + table(knn3,data.test[,drug])
      
      # Clasification Tree
      # fit the tree and plot it
      model2 = tree(as.formula(paste("as.factor(",drug,")~.",sep="")),data=data.train)
      pred2 = predict(model2, data.test, type="class") # Error: type "class" only for classification trees, fixed by using as.factor
      cM2 = cM2 + table(pred2, data.test[,drug])
      
      # Elastic net
      newY = model.matrix(as.formula(paste("~.-",drug,sep="")),data=data.train)
      model3 = cv.glmnet(newY, data.train[,drug], alpha = 1, family="binomial")
      newX = model.matrix(as.formula(paste("~.-",drug,sep="")),data=data.test)
      pred3 = predict(model3, s="lambda.1se", newx = newX, "response")
      cM3 = cM3 + table(pred3>0.5, data.test[,drug])
    }
    
    #Compare accuracy->prec,recall,f1,misclassification of the 3 models
    accu1 = sum(diag(cM1))/sum(cM1)
    accu2 = sum(diag(cM2))/sum(cM2)
    accu3 = sum(diag(cM3))/sum(cM3)
    accu4 = sum(diag(cM4))/sum(cM4)
    
    # Compare precision
    prec1 = cM1[1,1]/(cM1[1,1]+cM1[1,2]) #LDA
    prec2 = cM2[1,1]/(cM1[1,1]+cM2[1,2])
    prec3 = cM3[1,1]/(cM3[1,1]+cM3[1,2])
    prec4 = cM4[1,1]/(cM4[1,1]+cM4[1,2])
    
    # Compare recal
    recal1 = cM1[1,1]/(cM1[1,1]+cM1[2,1])
    recal2 = cM2[1,1]/(cM2[1,1]+cM2[2,1])
    recal3 = cM3[1,1]/(cM3[1,1]+cM3[2,1])
    recal4 = cM4[1,1]/(cM4[1,1]+cM4[2,1])
    
    # Compare f1
    f1_1 = (2*prec1*recal1)/(prec1+recal1) #same as (2*cM1[1,1])/(2*cM1[1,1]+cM1[2,1]+cM1[1,2])
    f1_2 = (2*prec2*recal2)/(prec2+recal2)
    f1_3 = (2*prec3*recal3)/(prec3+recal3)
    f1_4 = (2*prec4*recal4)/(prec4+recal4)
    
    # Compare misclassification
    misc1 = 1 - accu1
    misc2 = 1 - accu2
    misc3 = 1 - accu3
    misc4 = 1 - accu4
    
    if (returntype == 1){
      return(c(misc1, misc2, misc3, misc4))
    }
    if (returntype == 2){
      return(c(prec1, prec2, prec3, prec4))
    }
    if (returntype == 3){
      return(c(recal1, recal2, recal3, recal4))
    }
    if (returntype == 4){
      return(c(f1_1, f1_2, f1_3, f1_4))
    }
  } # work for accuracy, recall, f1, precision
  
  set.seed(0)
  nrep = 2
  
  idxmat.up = replicate(nrep, sample(1:nfold, sum(split==1), replace=T)) #ASK split is a 'double?' of ABC column (binary)
  idxmat.down = replicate(nrep, sample(1:nfold, sum(split==0), replace=T)) # what are we doing here?
  
  #Need separate matrix for misclassification, prec, recall, f1
  miscmat= matrix(NA, nrep, 4)
  precmat= matrix(NA, nrep, 4)
  recalmat= matrix(NA, nrep, 4)
  f1_mat= matrix(NA, nrep, 4)
  
  colnames(miscmat) = c("LDA", "Tree", "Elasticnet", "Logistic")
  colnames(precmat) = c("LDA", "Tree", "Elasticnet", "Logistic")
  colnames(recalmat) = c("LDA", "Tree", "Elasticnet", "Logistic")
  colnames(f1_mat) = c("LDA", "Tree", "Elasticnet", "Logistic")
  
  #for each repeat (each row ->each repeat?)
  for(jj in 1:nrep) miscmat[jj, ]= cv.accu(nfold, data.nres, data.res, idx.nres=idxmat.up[,jj], idx.res=idxmat.down[,jj], 1)
  for(jj in 1:nrep) precmat[jj, ]= cv.accu(nfold, data.nres, data.res, idx.nres=idxmat.up[,jj], idx.res=idxmat.down[,jj], 2)
  for(jj in 1:nrep) recalmat[jj, ]= cv.accu(nfold, data.nres, data.res, idx.nres=idxmat.up[,jj], idx.res=idxmat.down[,jj], 3)
  for(jj in 1:nrep) f1_mat[jj, ]= cv.accu(nfold, data.nres, data.res, idx.nres=idxmat.up[,jj], idx.res=idxmat.down[,jj], 4)
  # Please be patient, it takes very long
  
  # col(1,2,4) Wilcoxon test, you only need to pair-wise compare F1 scores of three models: LDA, Classification tree, Elastic net)
  ww12=wilcox.test(f1_mat[,1]-f1_mat[,2]) #between
  ww13=wilcox.test(f1_mat[,1]-f1_mat[,3]) 
  ww23=wilcox.test(f1_mat[,2]-f1_mat[,3]) 
  
  # plot misclassification rate / precision/ recall / F1 score
  boxplot(miscmat, ylab="Misclassification rate", main=paste0("for drug ",as.name(drug)))
  abline(h=seq(0, 1, 0.005), col = "gray", lty = "dotted")
  abline(v=seq(1, 3, 1.0), col = "gray", lty = "dotted")
  
  boxplot(precmat, ylab="Precision", main=paste0("for drug ",as.name(drug)))
  abline(h=seq(0, 1, 0.005), col = "gray", lty = "dotted")
  abline(v=seq(1, 3, 1.0), col = "gray", lty = "dotted")
  
  boxplot(recalmat, ylab="Recall", main=paste0("for drug ",as.name(drug)))
  abline(h=seq(0, 1, 0.005), col = "gray", lty = "dotted")
  abline(v=seq(1, 3, 1.0), col = "gray", lty = "dotted")
  
  boxplot(f1_mat, ylab="F1 score", main=paste0("for drug ",as.name(drug)))
  abline(h=seq(0, 1, 0.005), col = "gray", lty = "dotted")
  abline(v=seq(1, 3, 1.0), col = "gray", lty = "dotted")
  
  # Wilcox.test (1:LDA,2:Classification Tree,3:Elastic net)
  boxplot(f1_mat[,1]-f1_mat[,2], main=paste0("for drug ",as.name(drug)," P value = ", round(ww12$p.,12)),xlab="LDA-Tree", ylab="difference of f1-score")
  abline(h=seq(-1, 1, 0.0005), col = "gray", lty = "dotted")
  abline(h=0,col=2)
  
  boxplot(f1_mat[,1]-f1_mat[,3], main=paste0("for drug ",as.name(drug), "P value = ", round(ww13$p.,12)),xlab="LDA-Elastic", ylab="difference of f1-score")
  abline(h=seq(-1, 1, 0.0005), col = "gray", lty = "dotted")
  abline(h=0,col=2)
  
  boxplot(f1_mat[,2]-f1_mat[,3], main=paste0("for drug ",as.name(drug), "P value = ", round(ww23$p.,12)),xlab="Tree-Elastic", ylab="difference of f1-score")
  abline(h=seq(-1, 1, 0.005), col = "gray", lty = "dotted")
  abline(h=0,col=2)
  
  #############################################################################################################
}

