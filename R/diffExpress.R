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


diffExpress<-function(eset, spline_df=3, pval_cutoff=0.001, maxFC_cutoff=0.35){
  library(splines)
  times <- pData(eset)$time
  X <- ns(times, df=spline_df)
  # square root dose
  #TreatC <- sapply(pData(eset.w)$dose, function(x) sqrt(x))
  # log dose
  TreatC <- sapply(pData(eset)$dose, function(x) ifelse(x==0, 0, log10(x)+1))

  ## fit using Treat as continuous variable, note that -TreatC removes the concentration dependent shift
  design <- model.matrix(~ TreatC*X - TreatC)
  fit.conti <- eBayes(lmFit(exprs(eset), design))
  
  ###
  ### Top DE genes in treatment groups
  ###
  ## fit using Treat as continuous variable, note that -TreatC removes the concentration dependent shift
  design <- model.matrix(~ TreatC*X  - X )
  fit.conti <- eBayes(lmFit(exprs(eset), design))
  ## figure out the column number of Spline interaction terms 
  TreatCX_col <- grep("TreatC", colnames(fit.conti$t))
  #get the coeffecient F-stat
  treat.top <- subset(topTable(fit.conti, coef= TreatCX_col, number=nrow(exprs(eset))), adj.P.Val < pval_cutoff)
  treat.top$Gene <- rownames(treat.top)
  
  # extract the TreatX Spline coefficient
  coef.TreatCX <- treat.top[paste("TreatC.X", 1:spline_df, sep="")] 
  coef.TreatC <- treat.top$TreatC
  # find out the fitted time dependence of the treated
  td.profile <- ( t( X[!duplicated(times,),] %*% t(coef.TreatCX) ) + coef.TreatC ) * max(TreatC)
  # figure out the max
  treat.top$maxFC <- apply(td.profile, 1, function(x) x[which.max(abs(x))])
  hist(treat.top$maxFC)
  # cbind the Time response
  colnames(td.profile) <- paste("TC", 1:length(unique(times)), sep="")
  treat.top <- cbind(treat.top, td.profile)
  # remove genes with low expression change
  treat.top <- subset(treat.top, abs(maxFC)>maxFC_cutoff)

  all.cluster <-  hclust(as.dist(1 - cor(t(as.matrix(eset[treat.top$Gene,])))), method="ward.D")
  plot(all.cluster)
  #treat.top <- plyr::join(treat.top,  vec2df(sapply(cutree(all.cluster, k=6), function(x) paste("c", x, sep="")), "Gene", "cluster"), by="Gene")
  treat.top <- plyr::join(treat.top,  vec2df(sapply(cutree(all.cluster, h=100), function(x) paste("c", x, sep="")), "Gene", "cluster"), by="Gene")
  print(nrow(treat.top))
  print(table(treat.top$cluster))
  return(treat.top)
}

