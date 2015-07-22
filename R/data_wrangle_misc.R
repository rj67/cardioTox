#################################################################################################
# read the expression matrix
counts2DF <- function(counts_file){
  eset <- read.delim2(file=counts_file)
  
  # create df for experimental condition
  phenoData <- data.frame(ID = colnames(eset[2:ncol(eset)]), stringsAsFactors=F)
  phenoData$treatment <- sapply(phenoData$ID, function(x) strsplit(x, split=".", fixed=T)[[1]][1])
  phenoData$drug <- sapply( phenoData$treatment, function(x) gsub("[0-9]","", x))
  phenoData$dose <- sapply( phenoData$treatment, function(x) as.numeric(gsub("[A-Z]","", x, ignore.case=T)))
  phenoData$dose[is.na(phenoData$dose)] <- 1
  phenoData$time <- sapply(phenoData$ID, function(x) as.numeric(strsplit(strsplit(x, split=".", fixed=T)[[1]][2], split="_", fixed=T)[[1]][1]))
  phenoData$rep <- sapply(phenoData$ID, function(x) strsplit(strsplit(x, split=".", fixed=T)[[1]][2], split="_", fixed=T)[[1]][2])
  phenoData$condi <- with(phenoData, paste(treatment, time, sep="."))
  phenoData$treatment <- NULL
  rownames(phenoData) <- phenoData$ID
  
  # reshape eset into long data frame
  eset.l <- reshape2::melt(eset, id.vars="Gene", variable.name="ID", value.name="count")
  eset.l$logFC <- log(eset.l$count + 0.5, 2)
  eset.l %<>% plyr::join(., phenoData, by="ID")
  return(list(eset.l, phenoData))
}

# function to cast long df into matrix, then convert to ExpressionSet
df2Eset <- function(df, pData, value="count"){
  # recast the long count df into wide df
  eset.w <-  reshape2::dcast(df[c("Gene", "ID", value)], Gene~ID)
  rownames(eset.w) <- eset.w$Gene
  eset.w$Gene <- NULL
  eset <- ExpressionSet(assayData = as.matrix(eset.w) , phenoData = new("AnnotatedDataFrame", data=pData[colnames(eset.w),]))
  return(eset)
}

# # function quantile normalize expression Matrix
# qnEset <- function(eset){
#   library(preprocessCore)
#   eset.n <- normalize.quantiles(exprs(eset)) 
#   colnames(eset.n) <- colnames(exprs(eset))
#   rownames(eset.n) <- rownames(exprs(eset))
#   eset.n <- ExpressionSet(assayData = eset.n , phenoData = new("AnnotatedDataFrame", data=pData(eset)))
#   return(eset.n)
# }

# reshape normalized eset into long data frame
eset2DF <- function(eset, value="count"){
  eset.w <- as.data.frame(exprs(eset))
  eset.w$Gene <- rownames(eset.w)
  df <- reshape2::melt(eset.w, id.vars="Gene", variable.name="ID", value.name=value)
  if(value=="count"){
    df$logFC <- log(df$count + 0.5, 2)
  }
  df %<>% plyr::join(., pData(eset), by="ID")
  return(df)
}
# 
# # use the tools in EdgeR to normalize between array
# normalizeEdgeR <- function(eset){
#   library(edgeR)
#   # TMM scaling factor
#   # make sure eset is l
#   TMM <- calcNormFactors(DGEList(eset), method="TMM")$samples
#   exprs(eset) <- t(apply(exprs(eset), 1, function(x) x/(TMM$lib.size * TMM$norm.factors)*1e+6))
#   return(eset)
# }


# use the tools in EdgeR to normalize between array
normalizeDESeq <- function(eset){
  library(DESeq)
  # DESeq scaling factor
  norm.factors <- estimateSizeFactorsForMatrix(exprs(eset))
  exprs(eset) <- t(apply(exprs(eset), 1, function(x) x / norm.factors))
  return(eset)
}


# Combat remove batch effect
useComBat <- function(eset){
  library(sva)
  batch = factor(pData(eset)$rep==3)
  modcombat = model.matrix(~1, data=pData(eset))
  # eset is count data, take log2
  combat_edata = ComBat(dat=log(exprs(eset)+0.5, 2), batch=batch, mod=modcombat, par.prior=T, prior.plots=T)
  # convert back to count data
  eset <- ExpressionSet(assayData = 2^combat_edata-0.5 , phenoData = new("AnnotatedDataFrame", data=pData(eset)))
  return(eset)
}