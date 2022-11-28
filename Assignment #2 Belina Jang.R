###########################################################################
# Assignment #2 STAT 469 
###########################################################################
# NOTE:
#  I used the sample answer of assignment 1(Doug Skye's code) posted on Brightspace
#  and added two new models and made few changes. (added parts marked as "-Belina added", minor changes such as changes in numbers of lines, abline and etc are NOT marked)

# Load the data
library(MTPS)
library(tree)
data(HIV)
str(XX)
str(YY)
setwd("~") # set working path, all files will be read from here and write to here unless otherwise specified

# Control parameters
K=5 # folds
nexperiments = 50

# Binary drug resistance matrix yBin
yBin <- as.matrix(YY)
cutoffs <- c(2,3,3,1.5,1.5) # cutoff means 1 is non-resistant
for(ii in 1:5) yBin[,ii] <- (10^yBin[,ii] < cutoffs[ii])*1

# Use Kmeans to cluster binary resistances
set.seed(0); kk=kmeans(yBin, centers=10, nstart=100)
ord=order(kk$cluster)

# Stratify using clusters and save all experiment fold indexes in "folds.idx"
set.seed(0); folds.idx=matrix(0, length(kk$cluster), nexperiments)
for (cl in 1:length(kk$size)) {
  fset = kk$cluster==cl
  folds.idx[fset,] = replicate(nexperiments, sample(1:K, sum(fset), replace=T))
}
folds.idx = t(folds.idx)

###########################################################################
# visualize folds and clusters

# I like "pdf" command and "dev.off" to be paired and right before and after plot
pdf(file="BoxPlotsDrugs.pdf", width=14, height=8)
layout(mat=matrix(1:2,2,1), heights=c(1,3))

par(mar=c(0,4.1,1.1,2.1))
image(1:nrow(yBin), 1:ncol(yBin),yBin[ord,], xlab="sample index", ylab="drug index")
abline(v=cumsum(head(kk$size,-1)),col="green",lwd=3)

par(mar=c(0.1,4.1,0.1,2.1))
image(1:ncol(folds.idx), 1:nrow(folds.idx),t(folds.idx)[ord,], xlab="sample index", ylab="fold index")
abline(v=cumsum(head(kk$size,-1)),col="green",lwd=3)
dev.off() #close device

###########################################################################
# Run all the experiments, collect metrics for named models
XXdf = as.data.frame(XX)
# You can use array to save many matrix, but when the dimension of matrix are small, 
# flatten the matrix into a vector is another good choice
metrics = data.frame(matrix(ncol = 7, nrow = 0))
colnames(metrics) = c("Model", "Drug", "experiment", "tn", "fp", "fn", "tp")
models = c("LDA","Tree","Elasticnet") #c("GLM","LDA",sprintf("KNN%02i",seq(2,10))) ##

# loop experiments
cat("Experiment:")
for (e in 1:nexperiments) {  # use "ee" to replace "e"
  cat(" ", e)
  # loop for drugs
  for (c in colnames(yBin)) {
    
    # confusion matrix accumulators
    cMs = matrix(0, length(models), 4, dimnames=list(c(models), c("tn", "fp", "fn", "tp")))
    XXdf$resist = yBin[,c]
    experiment.passes = TRUE
    
    # loop folds
    for (f in 1:K) {
      test = folds.idx[e,]==f
      train =!test
      
      #cat("E=",e,"C=",c,"F=",f, "\n")
      
      # Logistic Regression-deleted
      
      # LDA model
      m=lda(resist ~ ., data=XXdf, subset=train)
      p=predict(m, XXdf[test,])
      cMs["LDA",] = cMs["LDA",] + table(p$class,yBin[test,c])
      
      # KNN-deleted
      
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
    }
    
    # save good results
    if (experiment.passes)
      metrics = rbind(metrics, t(sapply(models, function(x) c(x, c, e, cMs[x,]), USE.NAMES = F)))
  }
}

###########################################################################
# Calculate MCR, Precision, Recall, and F1-score
colnames(metrics) = c("Model", "Drug", "experiment", "tn", "fp", "fn", "tp")
write.csv(metrics,"metrics.csv")
for (i in 3:7) {metrics[,i] = as.integer(metrics[,i])}

fillmetrics = function (metrics) {
  metrics$mcr = with(metrics, (fp+fn)/(tn+fp+fn+tp))
  metrics$acc = with(metrics, (tp+tn)/(tn+fp+fn+tp))
  metrics$precision = with(metrics, tp/(tp+fp))
  metrics$recall = with(metrics, tp/(tp+fn))
  metrics$F1 = with(metrics, 2*tp/(2*tp + fp + fn))
  return(metrics)
}

# Aggregate the metrics for all drugs
metrics.agg = aggregate(cbind(tn,fp,fn,tp) ~ Model + experiment, data=metrics, FUN=sum)
metrics = fillmetrics(metrics)
metrics.agg = fillmetrics(metrics.agg)

# count total valid experiments for each drug
drug.experiments = sapply(sort(colnames(yBin)), function(x) sum(metrics$Drug==x)/length(models), USE.NAMES = F)
all.experiments = sum(drug.experiments)

# Select KNN with highest F1-Median
#knn.f1.medians = sapply(sprintf("KNN%02i",seq(2,10)), function(x) median(metrics[metrics$Model==x,"F1"]))
#knn.f1.best = which.max(knn.f1.medians)

###########################################################################
# Plots "BoxPlotsDrugs.pdf" with details per drug per model
pal = hcl.colors(length(models), palette = "viridis", alpha = 0.6)
layout(mat=matrix(1:1,1,1))
par(las=2, mar=c(6,6,2,1), mgp=c(4.2,1,0))
#xlabel1 = sprintf("Model grouped by Drug: %s", paste0(sort(colnames(yBin)), collapse = " "))
cex.drug = 1.2
cex.exps = 1.0
cex.legd = 0.75
ncol.legd = 3

# MCR
boxplot(mcr ~ Model + Drug, data=metrics, horizontal = F, names=rep(sort(models),5),
        col = pal, rev = FALSE, fixup = TRUE, cex.lab=1.5,
        main = "Figure 1: Misclassification Rate of all Models grouped by Drug name",
        ylab = "Misclassification Rate", xlab = "")
abline(v=seq(0.5, 25, 4), col="green",lwd=2)
abline(h=seq(0, 0.25, 0.01), col = "gray", lty = "dotted")
abline(v=seq(1, 20, 1.0), col = "gray", lty = "dotted")
text(x=seq(2.5, 20, 4), y=0.16, sort(colnames(yBin)), cex=cex.drug)
text(x=seq(2.5, 20, 4), y=0.15, sprintf("%d experiments",drug.experiments), cex=cex.exps)
legend("topleft",sort(models),fill=pal ,title="Model", inset = c(0.04,0.01), cex=cex.legd, ncol=ncol.legd)

# Precision
boxplot(precision ~ Model + Drug, data=metrics, horizontal = F, names=rep(sort(models),5),
        col = pal, rev = FALSE, fixup = TRUE, cex.lab=1.5,
        main = "Figure 2: Precision of all Models grouped by Drug name",
        ylab = "Precision", xlab = "")
abline(v=seq(0.5, 25, 4), col="green",lwd=2)
abline(h=seq(0, 1, 0.01), col = "gray", lty = "dotted")
abline(v=seq(1, 20, 1.0), col = "gray", lty = "dotted")
text(x=seq(2.5, 20, 4), y=0.835, sort(colnames(yBin)), cex=cex.drug)
text(x=seq(2.5, 20, 4), y=0.825, sprintf("%d experiments",drug.experiments), cex=cex.exps)

# Recall
boxplot(recall ~ Model + Drug, data=metrics, horizontal = F, names=rep(sort(models),5),
        col = pal, rev = FALSE, fixup = TRUE, cex.lab=1.5,
        main = "Figure 3: Recall of all Models grouped by Drug name",
        ylab = "Recall", xlab = "")
abline(v=seq(0.5, 25, 4), col="green",lwd=2)
abline(h=seq(0, 1, 0.01), col = "gray", lty = "dotted")
abline(v=seq(1, 20, 1.0), col = "gray", lty = "dotted")
text(x=seq(2.5, 20, 4), y=0.81, sort(colnames(yBin)), cex=cex.drug)
text(x=seq(2.5, 20, 4), y=0.8, sprintf("%d experiments",drug.experiments), cex=cex.exps)
legend("bottomleft",sort(models),fill=pal,title="Model", inset = c(0.04,0.01), cex=cex.legd, ncol=ncol.legd)

# F1-score
boxplot(F1 ~ Model + Drug, data=metrics, horizontal = F, names=rep(sort(models),5),
        col = pal, rev = FALSE, fixup = TRUE, cex.lab=1.5,
        main = "Figure 4: F1-score of all Models grouped by Drug name",
        ylab = "F1-score", xlab = "")
abline(v=seq(0.5, 25, 4), col="green",lwd=2)
abline(h=seq(0, 1, 0.01), col = "gray", lty = "dotted")
abline(v=seq(1, 20, 1.0), col = "gray", lty = "dotted")
text(x=seq(2.5, 20, 4), y=0.88, sort(colnames(yBin)), cex=cex.drug)
text(x=seq(2.5, 20, 4), y=0.87, sprintf("%d experiments",drug.experiments), cex=cex.exps)


###########################################################################
# Plots "BoxPlotsAll.pdf" with details combined drugs per model

if (.Platform$GUI != "RStudio") {
  print("Setup PDF...")
  pdf(file="BoxPlotsAll.pdf", width=14, height=10)
}
layout(mat=matrix(1:2,1,2, byrow = T))
par(las=2, mar=c(6,6,2,1), mgp=c(4.2,1,0))
cex.legd = 1.0

# MCR
boxplot(mcr ~ Model, data=metrics.agg,horizontal = F, 
        col = pal, rev = FALSE, fixup = TRUE, cex.lab=1.5,
        main = sprintf("Figure 5: MCR of all Models all Drugs (%d exp)", all.experiments),
        ylab = "Misclassification Rate", xlab = "Model")
abline(h=seq(0, 0.2, 0.004), col = "gray", lty = "dotted")
abline(v=seq(1, 4, 1.0), col = "gray", lty = "dotted")

# Precision
boxplot(precision ~ Model, data=metrics.agg,horizontal = F, 
        col = pal, rev = FALSE, fixup = TRUE, cex.lab=1.5,
        main = sprintf("Figure 6: Precision of all Models all Drugs (%d exp)", all.experiments),
        ylab = "Precision", xlab = "Model")
abline(h=seq(0, 1, 0.004), col = "gray", lty = "dotted")
abline(v=seq(1, 4, 1.0), col = "gray", lty = "dotted")

# Recall
boxplot(recall ~ Model, data=metrics.agg,horizontal = F, 
        col = pal, rev = FALSE, fixup = TRUE, cex.lab=1.5,
        main = sprintf("Figure 7: Recall of all Models all Drugs (%d exp)", all.experiments),
        ylab = "Recall", xlab = "Model")
abline(h=seq(0, 1, 0.004), col = "gray", lty = "dotted")
abline(v=seq(1, 4, 1.0), col = "gray", lty = "dotted")

# F1-score
boxplot(F1 ~ Model, data=metrics.agg,horizontal = F, 
        col = pal, rev = FALSE, fixup = TRUE, cex.lab=1.5,
        main = sprintf("Figure 8: F1-Score of all Models all Drugs (%d exp)", all.experiments),
        ylab = "F1-score", xlab = "Model")
abline(h=seq(0, 1, 0.004), col = "gray", lty = "dotted")
abline(v=seq(1, 4, 1.0), col = "gray", lty = "dotted")

###########################################################################
# Wilcoxon Test: Table of P-Values pairwise comparison of F1-Scores
#
# Two types of comparisions:
# 1) Compare models to each other for each drug.
# 2) Compare drugs to each other in a given model.
#

models.wilcoxon = c("LDA","Tree","Elasticnet")
nm = length(models.wilcoxon)
nd = length(colnames(yBin))

# Wilcoxon stats by model
wwp.bymodel.f1 = list()
for (c in sort(colnames(yBin)))
  wwp.bymodel.f1[[c]] = matrix(NA, nm, nm, dimnames = list(sort(models.wilcoxon), sort(models.wilcoxon)))

for (m1 in sort(models.wilcoxon))
  for (m2 in sort(models.wilcoxon)) 
    if (m1 != m2 )
    {
      for (c in sort(colnames(yBin)))
        wwp.bymodel.f1[[c]][m1,m2] = wilcox.test(metrics[metrics$Model==m1 & metrics$Drug==c,"F1"],metrics[metrics$Model==m2 & metrics$Drug==c,"F1"], exact = F, paired = T)$p.value
    }

# Wilcoxon stats by drug
wwp.bydrug.f1 = list()
for (m1 in sort(models.wilcoxon))
  wwp.bydrug.f1[[m1]] = matrix(NA, nd, nd, dimnames = list(sort(colnames(yBin)), sort(colnames(yBin))))

for (d1 in sort(colnames(yBin)))
  for (d2 in sort(colnames(yBin))) 
    if (d1 != d2 )
    {
      for (m1 in sort(models.wilcoxon))
        wwp.bydrug.f1[[m1]][d1,d2] = wilcox.test(metrics[metrics$Model==m1 & metrics$Drug==d1,"F1"],metrics[metrics$Model==m1 & metrics$Drug==d2,"F1"], exact = F, paired = F)$p.value
    }

# make tables
library(knitr)

cat("Table 1 – Head of the metrics data frame\n")
print(head(metrics,4))

cat("Table 2 – Summary of accumulated Confusion Matrices\n")
print(aggregate(cbind(tn,fp,fn,tp) ~ Model, data=metrics, FUN=sum))

# plotting by model
tabnum = 3
for (ci in 1:length(colnames(yBin))) {
  c = sort(colnames(yBin))[ci]
  print(kable(wwp.bymodel.f1[[c]], caption = sprintf("%d: Wilcoxon Pairs of F1-Scores for Drug: %s (%d trials)", tabnum, c, drug.experiments[ci])))
  tabnum = tabnum + 1
}

# plotting by drug
for (m1 in sort(models.wilcoxon)) {
  print(kable(wwp.bydrug.f1[[m1]], caption = sprintf("%d: Wilcoxon Pairs of F1-Scores for Model: %s", tabnum, m1)))
  tabnum = tabnum + 1
}

#check difference of f1 score between Elastic and LDA for AZT - Belina added
metrics.AZT = matrix(NA,49,2)
metrics.AZT.elastic = metrics[metrics$Mode=='Elasticnet' & metrics$Drug=='AZT',"F1"]
metrics.AZT.lda = metrics[metrics$Mode=='LDA' & metrics$Drug=='AZT',"F1"]
metrics.AZT = cbind(metrics.AZT.elastic,metrics.AZT.lda)
colnames(metrics.AZT) = c("Elastic","LDA")

ww12=wilcox.test(metrics.AZT[,1]-metrics.AZT[,2])
boxplot(metrics.AZT[,1]-metrics.AZT[,2], main=paste0("Figure 9: for AZT ", "  P value = ", round(ww12$p.,7)),xlab="Elasticnet-LDA", ylab="difference of f1-score")
abline(h=0,col=2)

