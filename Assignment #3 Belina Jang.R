###########################################################################
# Assignment #3 STAT 469 
###########################################################################
# NOTE:
#  I used the sample answer of assignment 1(Doug Skye's code) posted on Brightspace
#  added two new models and made few changes to create the report
#  Below are mostly additions I made only (added parts marked as "-Belina added", minor changes such as changes in numbers of lines, abline and etc are NOT marked)

cat("MTPS Experiment:")
for (e in 1:nexperiments) {
  cat(" ", e)
  # MTPS: no loop for drugs
  
  # confusion matrix accumulators
  cMsM = matrix(0, length(drugs), 4, dimnames=list(c(drugs),c("tn", "fp", "fn", "tp")))
  experiment.passes = TRUE
  
  # loop folds - Belina added
  for (f in 1:K) {
    test = folds.idx[e,]==f
    train =!test
    
    y.test.bin  <- yBin[test, ]
    y.train.bin <- yBin[-test, ]
    x.test.bin  <- XXdf[test, ]
    x.train.bin <- XXdf[-test, ]
    
    cat("MTPS: E=",e,"F=",f, "\n")
    
    # MTPS (residual stacking only) - Belina added
    fit.prs.std <- MTPS(xmat = x.train.bin, ymat=y.train.bin,
                        family = "binomial",
                        cv = FALSE, residual = TRUE,
                        method.step1 = rpart1,
                        method.step2 = lm1,
                        resid.type = "deviance", resid.std = TRUE) 
    pred.prs.std <- predict(fit.prs.std, x.test.bin)
    
    for (yy in 1 : ncol(yBin[test,])) {
      print(colnames(y.test.bin)[yy])
      print(table((pred.prs.std[,yy] > 0.5) * 1, y.test.bin[,yy]))
      cMsM[colnames(y.test.bin)[yy],] = cMsM[colnames(y.test.bin)[yy],] + table(pred.prs.std[,yy] > 0.5,y.test.bin[,yy])
    }
    
  }
 
# Elastic net - Belina added
set.seed(0)
newY = model.matrix(~.-resist,data=XXdf[train,])
m = cv.glmnet(newY, yBin[train,c], alpha = 1, family="binomial")
newX = model.matrix(~.-resist,data=XXdf[test,])
p = predict(m, s="lambda.1se", newx = newX, "response")
cMs["Elasticnet",] = cMs["Elasticnet",] + table(p>0.5, yBin[test,c])
      
# random forest (default setting) - Belina added
set.seed(0)
m = randomForest(as.factor(resist) ~ ., data=XXdf, subset = train)
p = predict(m,XXdf[test,])
cMs["RandomForest",] = cMs["RandomForest",] + table(p, yBin[test,c])

###########################################################################

#check difference of f1 score between Elastic and LDA - Belina added
metrics.f1 = matrix(NA,49,2)
metrics.f1.elastic = metrics[metrics$Mode=='Elasticnet',"F1"]
metrics.f1.random = metrics[metrics$Mode=='RandomForest',"F1"]
metrics.f1 = cbind(metrics.f1.elastic,metrics.f1.random)
colnames(metrics.f1) = c("Elastic","RandomForest")

ww12=wilcox.test(metrics.f1[,1]-metrics.f1[,2])
boxplot(metrics.f1[,1]-metrics.f1[,2], main=paste0("Figure 9: for F1-score", "  P value = ", round(ww12$p.,7)),xlab="Elasticnet-RandomForest", ylab="difference of f1-score")
abline(h=0,col=2)

