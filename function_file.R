spliner <- function(vec, n.knots = 2) {
  ## Function to calculate bases for regression splines. Modified from code 
  ## provided in Crainiceanu, C., Ruppert, D. & Wand, M.P. Bayesian analysis for
  ## penalized spline regression using WinBUGS. J. Stat. Soft. 14, 1 24(2005).
  ## Parameter vec is a vector defining the raw data vector, n.knots defines the
  ## number of knots in the GAM.
  N         <- length(vec)
  x.time    <- c(vec)
  ZFE       <- cbind(rep(1,N), x.time)
  x.knots   <- quantile(unique(x.time), seq(0, 1, length = (n.knots+2))[-c(1,
                                                                           (n.knots+2))], na.rm = TRUE)
  Z_K       <- (abs(outer(x.time,x.knots,"-")))^3
  OMEGA.all <- (abs(outer(x.knots,x.knots,"-")))^3
  svd.OMEGA.all  <- svd(OMEGA.all)
  sqrt.OMEGA.all <- t(svd.OMEGA.all$v %*% (t(svd.OMEGA.all$u) *
                                             sqrt(svd.OMEGA.all$d)))
  Z.out     <- t(solve(sqrt.OMEGA.all, t(Z_K)))
  return(Z.out)
}

#MODIFIED

## modified version of https://github.com/cran/plotKML/blob/master/R/readGPX.R
track.readGPX.element <- function(file) {
  # element = "metadata", "wpt", "rte", "trk"
  
  ret <- xmlTreeParse(file, useInternalNodes = TRUE)
  # top structure: 
  top <- xmlRoot(ret)
  
  # check if there is any content:
  # tracks:
  ret <- NULL
  nu <- which(names(top) %in% element)
  for(c in seq_along(nu)){
    lst <- which(names(top[[nu[c]]]) %in% "trkseg")
    nm <- names(top[[nu[c]]][[lst[1]]][[1]])
    ret[[c]] <- list(NULL)
    for(i in seq_along(lst)) {
      trkpt <- top[[nu[c]]][[lst[i]]]
      ret[[c]][[i]] <- data.frame(NULL)
      ## get columns (http://www.topografix.com/GPX/1/1/#type_wptType)
      lon <- as.numeric(unlist(xmlSApply(trkpt, xmlGetAttr, "lon")))
      lat <- as.numeric(unlist(xmlSApply(trkpt, xmlGetAttr, "lat")))
      ret[[c]][[i]][1:length(lon),"lon"] <- lon
      ret[[c]][[i]][1:length(lat),"lat"] <- lat
      
      if(!length(nm)==0){
        for(j in 1:length(nm)){
          xm <- as.character(sapply(sapply(xmlChildren(trkpt), function(x) x[[nm[[j]]]]), xmlValue))
          ret[[c]][[i]][1:length(xm), nm[[j]]] <- xm 
        }
      } 
    }
    names(ret[[c]]) <- xmlValue(top[[nu[c]]][["name"]])
  }   
  
  return(ret)
}


#MODIFIED

# source("http://www.sthda.com/upload/rquery_cormat.r")

rquery.cormat.spearman <- function(x, type=c('lower', 'upper', 'full', 'flatten'),
                                   graph=TRUE, graphType=c("correlogram", "heatmap"),
                                   col=NULL, cor.method = "spearman", ...)
{
  library(corrplot)
  # Helper functions
  #+++++++++++++++++
  # Compute the matrix of correlation p-values
  cor.pmat <- function(x, cor.method) {
    mat <- as.matrix(x)
    n <- ncol(mat)
    p.mat<- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        tmp <- cor.test(mat[, i], mat[, j], method = "spearman", ...)
        p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      }
    }
    colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
    p.mat
  }
  # Get lower triangle of the matrix
  getLower.tri<-function(mat){
    upper<-mat
    upper[upper.tri(mat)]<-""
    mat<-as.data.frame(upper)
    mat
  }
  # Get upper triangle of the matrix
  getUpper.tri<-function(mat){
    lt<-mat
    lt[lower.tri(mat)]<-""
    mat<-as.data.frame(lt)
    mat
  }
  # Get flatten matrix
  flattenCorrMatrix <- function(cormat, pmat) {
    ut <- upper.tri(cormat)
    data.frame(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      cor  =(cormat)[ut],
      p = pmat[ut]
    )
  }
  # Define color
  if (is.null(col)) {
    col <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", 
                              "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
                              "#4393C3", "#2166AC", "#053061"))(200)
    col<-rev(col)
  }
  
  ?cor
  # Correlation matrix
  cormat<-signif(cor(x, use = "complete.obs", method = "spearman", ...),2)
  pmat<-signif(cor.pmat(x, ...),2)
  # Reorder correlation matrix
  ord<-corrMatOrder(cormat, order="hclust")
  cormat<-cormat[ord, ord]
  pmat<-pmat[ord, ord]
  # Replace correlation coeff by symbols
  sym<-symnum(cormat, abbr.colnames=FALSE)
  # Correlogram
  if(graph & graphType[1]=="correlogram"){
    corrplot(cormat, type=ifelse(type[1]=="flatten", "lower", type[1]),
             tl.col="black", tl.srt=45,col=col,...)
  }
  else if(graphType[1]=="heatmap")
    heatmap(cormat, col=col, symm=TRUE)
  # Get lower/upper triangle
  if(type[1]=="lower"){
    cormat<-getLower.tri(cormat)
    pmat<-getLower.tri(pmat)
  }
  else if(type[1]=="upper"){
    cormat<-getUpper.tri(cormat)
    pmat<-getUpper.tri(pmat)
    sym=t(sym)
  }
  else if(type[1]=="flatten"){
    cormat<-flattenCorrMatrix(cormat, pmat)
    pmat=NULL
    sym=NULL
  }
  list(r=cormat, p=pmat, sym=sym)
}
