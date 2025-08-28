library(future.apply) # install.packages("future.apply")
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(rlang)
# Get ps from H1 small mols. These are rchisq but then need to scrunch closer to y-axis
rm(list=ls())
N=1e4
propH1=0.05
nH1=round(propH1*N)
# tmp=rchisq(round(7*nH1/8), 1)
# tmp=tmp/max(tmp)
# tmp2=runif(nH1 - round(7*nH1/8))/500
# psH1=c(tmp, tmp2)
tmp=rexp(2*nH1)
psH1=tmp/(2* max(tmp))
psH1=psH1[order(psH1)]
# psH1=psH1[psH1<1]
psH1=psH1[1:nH1]^1.2
summary(psH1)
hist(psH1)
nH0=N-nH1
psH0=runif(nH0)
ps.df=data.frame(ps=c(psH1, psH0), 
  H1.ind=c(rep(1, nH1), rep(0, nH0)))
ps.df=ps.df[order(ps.df$ps),]
ps.df[1:20,]
print(summary(which(as.logical(ps.df$H1.ind))))
boxplot(ps.df$ps~ps.df$H1.ind)

# Get prec and recall 
precRecall.df=t(future_sapply(1:nrow(ps.df), function(x){
  tp=sum(ps.df$H1.ind[1:x])
  fp=sum(1-ps.df$H1.ind[1:x])
  fn=sum(ps.df$H1.ind[min((1+x), N):N])
  prec <- tp/(tp+fp)
  recall= tp/(tp+fn)
  return(c(prec, recall))
}))
dim(precRecall.df)
plot(precRecall.df[, 2], precRecall.df[, 1])
ps.df=data.frame(ps.df, precRecall.df)
colnames(ps.df)=c("pVal", "H1.ind", "precision", "recall")
dim(ps.df)
ps.df[1:5,]
g=ggplot(ps.df, aes(x = recall, y = precision)) +
  geom_path(linewidth = 1) +
  geom_point(size = 0.8, alpha = 0.4) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(
    title = "Screen Simulation Summary: Precision vs Recall",
    # subtitle = if (!is.null(auc_est)) sprintf("Trapezoidal PR-AUC ≈ %.3f", auc_est) else NULL,
    x = "Recall (TP / (TP + FN))",
    y = "Precision (TP / (TP + FP))"
  ) +
  theme_minimal(base_size = 12)
g
summary(ps.df)
