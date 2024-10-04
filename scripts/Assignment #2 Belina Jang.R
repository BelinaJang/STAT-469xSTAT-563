###########################################################################
# Assignment #2 STAT 469 
###########################################################################
# NOTE:
#  I used the sample answer of assignment 1(Doug Skye's code) posted on Brightspace
#  added two new models and made few changes to create the report
#  Below are mostly additions I made only (added parts marked as "-Belina added", minor changes such as changes in numbers of lines, abline and etc are NOT marked)
#  Please check Assignment1-2(usingloop).R to see the overall setup
      
# Clasification Tree - Belina added
      
# fit the tree and plot it
m = tree(as.factor(resist) ~ .,data=XXdf, subset=train)
p = predict(m, XXdf[test,], type="class") # Error: type "class" only for classification trees
cMs["Tree",] = cMs["Tree",] + table(p, yBin[test,c])
      
# Elastic net - Belina added
newY = model.matrix(~.-resist,data=XXdf[train,])
m = cv.glmnet(newY, yBin[train,c], alpha = 1, family="binomial")
newX = model.matrix(~.-resist,data=XXdf[test,])
p = predict(m, s="lambda.1se", newx = newX, "response")
cMs["Elasticnet",] = cMs["Elasticnet",] + table(p>0.5, yBin[test,c])
      
#m = cv.glmnet(XXdf[train,], yBin[train,c], alpha = 1, family="binomial")
#p = predict(m, s="lambda.1se", newx = XXdf[test,], "response")
#cMs["Elasticnet",] = cMs["Elasticnet",] + table(p>0.5, yBin[test,c])

###########################################################################

#check difference of f1 score between Elastic and LDA for AZT - Belina added
metrics.AZT = matrix(NA,49,2)
metrics.AZT.elastic = metrics[metrics$Mode=='Elasticnet' & metrics$Drug=='AZT',"F1"]
metrics.AZT.lda = metrics[metrics$Mode=='LDA' & metrics$Drug=='AZT',"F1"]
metrics.AZT = cbind(metrics.AZT.elastic,metrics.AZT.lda)
colnames(metrics.AZT) = c("Elastic","LDA")

ww12=wilcox.test(metrics.AZT[,1]-metrics.AZT[,2])
boxplot(metrics.AZT[,1]-metrics.AZT[,2], main=paste0("Figure 9: for AZT ", "  P value = ", round(ww12$p.,7)),xlab="Elasticnet-LDA", ylab="difference of f1-score")
abline(h=0,col=2)

