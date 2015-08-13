library(limma)
library(splines)

# 
# 
# ###
# ### Top DE genes in Control 
# ###
# ctr.top <- subset(topTable(fit.conti, coef= 2:(1+spline_df), number=nrow(exprs(eset))), adj.P.Val < pval_cutoff)
# # ctr.top <- topTable(fit.conti, coef= 2:(1+spline_df), number=nrow(exprs(eset.w)))
# # extract the control X Spline coefficient
# coef.X <- as.matrix(ctr.top[paste("X", 1:spline_df, sep="")] )
# # find out the fitted time dependence of the control 
# td.profile.X <- t( X %*% t(coef.X ))[, !duplicated(times)]
# # # figure out the 
# ctr.top$range <- apply(td.profile.X, 1, function(x) range(x)[2]-range(x)[1])
# ctr.top$Gene <- rownames(ctr.top)


diffExpress<-function(eset, spline_df=2, pval_cutoff=0.01, maxFC_cutoff=0.4){
  library(splines)
  times <- pData(eset)$time
  X <- ns(times, df=spline_df)
  # square root dose
  #TreatC <- sapply(pData(eset.w)$dose, function(x) sqrt(x))
  # log dose
  TreatC <- sapply(pData(eset)$dose, function(x) ifelse(x==0, 0, log10(x)+1))

  ###
  ### Top DE genes in treatment groups
  ###
  ## fit using Treat as continuous variable 
  design <- model.matrix(~ TreatC*X  )
  v <- voom(eset, design=design, plot=TRUE)
  fit.conti <- eBayes(lmFit(v, design), trend=T)
  #fit.conti <- eBayes(lmFit(log(exprs(eset)+0.5, 2), design), trend=T)
  ## figure out the column number of Spline interaction terms 
  TreatC_col <- grep("TreatC", colnames(fit.conti$t))
  #get the coeffecient F-stat
  treat.top <- subset(topTable(fit.conti, coef= TreatC_col, number=nrow(exprs(eset))), adj.P.Val < pval_cutoff)
  treat.top$Gene <- rownames(treat.top)
  
  # find out the fitted time dependence of the treated
  td.profile <- fit.conti$coefficients[treat.top$Gene, ] %*% t(design[!duplicated(pData(eset)$condi),])
  # figure out the max
  treat.top$maxFC <- apply(td.profile, 1, function(x) {y <- x-x[1]; y[which.max(abs(y))] })
  with(treat.top, plot(maxFC, log10(P.Value)))
  # cbind the Time response
  colnames(td.profile) <- pData(eset)$condi[!duplicated(pData(eset)$condi)]
  treat.top <- cbind(treat.top, td.profile)
  # remove genes with low expression change
  treat.top <- subset(treat.top, abs(maxFC)>maxFC_cutoff)

  # test concentration dependence at 1 day
  dose_vec <- pData(eset[, pData(eset)$time==24 &  pData(eset)$dose!=0])$dose
  exp_dose <- exprs(eset[, pData(eset)$time==24 &  pData(eset)$dose!=0])[treat.top$Gene, ]
  #dose_cor <- do.call(rbind, apply(exp_24, 1, function(x)   tidy(lm(log2(x+0.5) ~ dose_24))[2,c("statistic", "p.value")]))
  treat.top$dose.pval <- apply(exp_dose, 1, function(x)  dcov.test(log2(x + 0.5), dose_vec)$p.value )
  time_vec <- pData(eset[, pData(eset)$dose==3 ])$time
  exp_time <- exprs(eset[, pData(eset)$dose==3 ])[treat.top$Gene, ]
  treat.top$time.pval <- apply(exp_time, 1, function(x)  dcov.test(log2(x + 0.5), time_vec)$p.value )
  
  print(with(treat.top, table(dose.pval<0.1, time.pval<0.1)))
  #vec2df <-function(vec, keyname, value){
  #  df <- data.frame(V1=names(vec), V2=vec)
  #  colnames(df) <- c(keyname, value)
  #  return(df)
  #}
  #all.cluster <-  hclust(as.dist(1 - cor(t(as.matrix(eset[treat.top$Gene,])))), method="ward.D")
  #plot(all.cluster)
  #treat.top <- plyr::join(treat.top,  vec2df(sapply(cutree(all.cluster, h=100), function(x) paste("c", x, sep="")), "Gene", "cluster"), by="Gene")
  #print(table(treat.top$cluster))
  print(nrow(treat.top))
  return(treat.top)
}

