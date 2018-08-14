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



invlogit <- function(x) exp(x)/(1 + exp(x))

sample.fixed.params <- function(model = model_test, nsample = 10000, param = "") {
  params <- names(model$marginals.fixed)
  which.param <- grep(param, params)
  
  if(sum(which.param)==0) stop("Check parameter names")
  dist <- matrix(ncol = length(which.param), nrow = nsample)
  for(i in 1:NROW(which.param)) {
    dist[,i] <- sample(model$marginals.fixed[[which.param[i]]][,1], nsample, TRUE,
                       model$marginals.fixed[[which.param[i]]][,2])
    
  } 
  dist <- data.frame(dist)
  names(dist) <- names(model$marginals.fixed)
  return(dist)
}


#### GENERAL BUILD OF AN EFFECT PLOT
effect.plot <- function(quantiles, lines, df, rug, with.lines, effect.of){
  
  gg <- 
    ggplot()+
    geom_ribbon(data = quantiles, aes(x = pred, ymin = lower*100, ymax = upper*100), fill = "gray")+
    geom_line(data = quantiles, aes(x = pred, y = mid*100), size = 1.5)+
    xlab(effect.of)+
    ylab("propability [%]")+
    theme_bw() 
  
  if(with.lines) gg <- gg + geom_line(data = lines, aes(x = pred, y = value * 100, group = variable), col = "red", alpha = .1)
  if(rug) gg <- gg + geom_point(data = df, aes(x = pred, y = zero), shape = "|", alpha = .1)
  
  if(!rug){
    gg <- grid.arrange(
      gg + theme(plot.margin=unit(c(-0.08,1,1,1), "cm")), 
      ggplot(df, aes(x = pred))+
        geom_histogram(bins = 100, col = "black", fill = "gray")+
        theme_void()+
        ggtitle(paste0("Effect of: ", effect.of))+
        theme(plot.margin=unit(c(.3,1,-0.08,1), "cm")),
      layout_matrix = matrix(c(2,rep(1,7),2,rep(1,7)), nc = 2))
  }
  
  
  return(gg)
}


#GENERATE CONDITIONALS
conditional.plot <- function(cond.vars = seqs[, names(seqs) %in% preds], effect.of = "VD", model = spatial_model, posterior.samples = 500, original.df, with.lines = F, rug = T){
  cond.vars[,-match(effect.of, names(cond.vars))] <- t(replicate(nrow(cond.vars),apply(cond.vars[,-match(effect.of, names(cond.vars))], 2, median), simplify = T))
  
  dist <- sample.fixed.params(model = model, nsample = posterior.samples)
  dist <- dist[,colnames(dist) %in% names(cond.vars)]
  dist <- dist[,match(names(dist), names(cond.vars))]
  
  m.cond.vars <- as.matrix(cond.vars)
  m.dist <- as.matrix(dist)
  
  out <- vector("list", length = nrow(m.dist))
  for(i in seq(nrow(m.dist))) out[[i]] <- as.vector(invlogit(m.cond.vars %*% m.dist[i,]))
  if(with.lines) gg.out.all <- reshape2::melt(data.frame(do.call(cbind, out), pred = original.seqs[,effect.of]), id.vars = "pred")
  gg.out.q <- data.frame(t(apply(do.call(cbind, out), 1, quantile, probs = c(0.025, 0.5, 0.975))), pred = original.seqs[,effect.of])
  names(gg.out.q)[1:3] <- c("lower", "mid", "upper")
  original.df <- data.frame(pred = original.df[,effect.of], zero = 0)
  
  gg <- effect.plot(quantiles = gg.out.q, lines = gg.out.all, df = original.df, rug = rug, with.lines = with.lines, effect.of = effect.of)
  
  return(gg)
}