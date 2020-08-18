##Load some important prerequisites for analysis. Script won't run without these!
{
  if (!requireNamespace("BiocManager", quietly = TRUE)) { install.packages("BiocManager") }
  BiocManager::install()

  if("tsne" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("tsne")}
  library("tsne")
    
  if("pheatmap" %in% rownames(installed.packages()) == FALSE) {install.packages("pheatmap")}
  library("pheatmap")
    
  if("MASS" %in% rownames(installed.packages()) == FALSE) {install.packages("MASS")}
  library("MASS")  
    
  if("cluster" %in% rownames(installed.packages()) == FALSE) {install.packages("cluster")}
  library("cluster")  
    
  if("mclust" %in% rownames(installed.packages()) == FALSE) {install.packages("mclust")}
  library("mclust")  
    
  if("flexmix" %in% rownames(installed.packages()) == FALSE) {install.packages("flexmix")}
  library("flexmix")  
  
  if("lattice" %in% rownames(installed.packages()) == FALSE) {install.packages('lattice')}
  library("lattice")
  
  if("fpc" %in% rownames(installed.packages()) == FALSE) {install.packages('fpc')}
  library("fpc")
    
  if("amap" %in% rownames(installed.packages()) == FALSE) {install.packages("amap")}
  library("amap")    

  if("RColorBrewer" %in% rownames(installed.packages()) == FALSE) {install.packages("RColorBrewer")}
  library("RColorBrewer")  
  
  if("locfit" %in% rownames(installed.packages()) == FALSE) {install.packages("locfit")}
  library("locfit")  

  if("vegan" %in% rownames(installed.packages()) == FALSE) {install.packages("vegan")}
  library("vegan")
  
  if("Rtsne" %in% rownames(installed.packages()) == FALSE) {install.packages("Rtsne")}
  library("Rtsne")
  
  if("scran" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("scran")}
  library("scran")
}

## load required packages.
require(tsne)
require(pheatmap)
require(MASS)
require(cluster)
require(mclust)
require(flexmix)
require(lattice)
require(fpc)
require(amap)
require(RColorBrewer)
require(locfit)
require(vegan)
require(Rtsne)
require(scran)


## class definition
SCseq <- setClass("SCseq", slots = c(expdata = "data.frame", ndata = "data.frame", fdata = "data.frame", distances = "matrix", tsne = "data.frame", cluster = "list", background = "list", out = "list", cpart = "vector", fcol = "vector", filterpar = "list", clusterpar = "list", outlierpar ="list" ))

setValidity("SCseq",
            function(object) {
              msg <- NULL
              if ( ! is.data.frame(object@expdata) ){
                msg <- c(msg, "input data must be data.frame")
              }else if ( nrow(object@expdata) < 2 ){
                msg <- c(msg, "input data must have more than one row")
              }else if ( ncol(object@expdata) < 2 ){
                msg <- c(msg, "input data must have more than one column")
              }else if (sum( apply( is.na(object@expdata),1,sum ) ) > 0 ){
                msg <- c(msg, "NAs are not allowed in input data")
              }else if (sum( apply( object@expdata,1,min ) ) < 0 ){
                msg <- c(msg, "negative values are not allowed in input data")
              }
              if (is.null(msg)) TRUE
              else msg
            }
            )

setMethod("initialize",
          signature = "SCseq",
          definition = function(.Object, expdata ){
            .Object@expdata <- expdata
            .Object@ndata <- expdata
            .Object@fdata <- expdata
            validObject(.Object)
            return(.Object)
          }
          )

setGeneric("filterdata", function(object, mintotal=3000, minexpr=5, minnumber=1, maxexpr=Inf, downsample=TRUE, sfn=FALSE, dsn=1, rseed=17000) standardGeneric("filterdata"))

setMethod("filterdata",
          signature = "SCseq",
          definition = function(object,mintotal,minexpr,minnumber,maxexpr,downsample,sfn,dsn,rseed) {
            if ( ! is.numeric(mintotal) ) stop( "mintotal has to be a positive number" ) else if ( mintotal <= 0 ) stop( "mintotal has to be a positive number" )
            if ( ! is.numeric(minexpr) ) stop( "minexpr has to be a non-negative number" ) else if ( minexpr < 0 ) stop( "minexpr has to be a non-negative number" )
            if ( ! is.numeric(minnumber) ) stop( "minnumber has to be a non-negative integer number" ) else if ( round(minnumber) != minnumber | minnumber < 0 ) stop( "minnumber has to be a non-negative integer number" )
            if ( ! ( is.numeric(downsample) | is.logical(downsample) ) ) stop( "downsample has to be logical (TRUE/FALSE)" )
            if ( ! ( is.numeric(sfn) | is.logical(sfn) ) ) stop( "sfn has to be logical (TRUE/FALSE)" )
            if ( ! is.numeric(dsn) ) stop( "dsn has to be a positive integer number" ) else if ( round(dsn) != dsn | dsn <= 0 ) stop( "dsn has to be a positive integer number" )
            object@filterpar <- list(mintotal=mintotal, minexpr=minexpr, minnumber=minnumber, maxexpr=maxexpr, downsample=downsample, dsn=dsn, sfn=sfn)
            object@ndata <- object@expdata[,apply(object@expdata,2,sum,na.rm=TRUE) >= mintotal]
            if ( downsample ){
              set.seed(rseed)
              object@ndata <- downsample(object@expdata,n=mintotal,dsn=dsn)
            }else if ( sfn ){
              d <- computeSumFactors(as.matrix( object@ndata ))
              object@ndata <- as.data.frame(t(t(object@ndata)/d)) + .1
            }else{
              x <- object@ndata
              ##object@ndata <- as.data.frame( t(t(x)/apply(x,2,sum))*median(apply(x,2,sum,na.rm=TRUE)) + .1 )
              object@ndata <- as.data.frame( t(t(x)/apply(x,2,sum))*min(apply(x,2,sum,na.rm=TRUE)) + .1 )
            }
            x <- object@ndata

            object@fdata <- x[apply(x>=minexpr,1,sum,na.rm=TRUE) >= minnumber,]
            x <- object@fdata
            object@fdata <- x[apply(x,1,max,na.rm=TRUE) < maxexpr,]

            x <- object@expdata[,names(object@ndata)]
            object@ndata <- as.data.frame( t(t(x)/apply(x,2,sum))*min(apply(x,2,sum,na.rm=TRUE)) + .1 )

            return(object)
          }
          )


downsample <- function(x,n,dsn){
  x <- round( x[,apply(x,2,sum,na.rm=TRUE) >= n], 0)
  nn <- min( apply(x,2,sum) )
  for ( j in 1:dsn ){
    z  <- data.frame(GENEID=rownames(x))
    rownames(z) <- rownames(x)
    initv <- rep(0,nrow(z))
    for ( i in 1:dim(x)[2] ){
      y <- aggregate(rep(1,nn),list(sample(rep(rownames(x),x[,i]),nn)),sum)
      na <- names(x)[i]
      names(y) <- c("GENEID",na)
      rownames(y) <- y$GENEID
      z[,na] <- initv
      k <- intersect(rownames(z),y$GENEID)
      z[k,na] <- y[k,na]
      z[is.na(z[,na]),na] <- 0
    }
    rownames(z) <- as.vector(z$GENEID)
    ds <- if ( j == 1 ) z[,-1] else ds + z[,-1]
  }
  ds <- ds/dsn + .1
  return(ds)
}

dist.gen <- function(x,method="euclidean", ...) if ( method %in% c("spearman","pearson","kendall") ) as.dist( 1 - cor(t(x),method=method,...) ) else if ( method %in% c("eupearson") ) sqrt(2*as.dist( 1 - cor(t(x),method="pearson",...) )) else if ( method %in% c("logpearson") ) as.dist( 1 - cor(log2(t(x)),method="pearson",...) ) else if ( method == "binary" )  dist( 1*(x > .1) ) else if ( method == "relentr") dist(x,FUN=relentr) else dist(x,method=method,...)


dist.gen.pairs <- function(x,y,...) dist.gen(t(cbind(x,y)),...)

binompval <- function(p,N,n){
  pval   <- pbinom(n,round(N,0),p,lower.tail=TRUE)
  pval[!is.na(pval) & pval > 0.5] <- 1-pval[!is.na(pval) & pval > 0.5]
  return(pval)
}

setGeneric("plotgap", function(object) standardGeneric("plotgap"))

setMethod("plotgap",
          signature = "SCseq",
          definition = function(object){
            if ( length(object@cluster$kpart) == 0 ) stop("run clustexp before plotgap")
            if ( sum(is.na(object@cluster$gap$Tab[,3])) > 0 ) stop("run clustexp with do.gap = TRUE first")
            plot(object@cluster$gap,ylim=c( min(object@cluster$gap$Tab[,3] - object@cluster$gap$Tab[,4]),  max(object@cluster$gap$Tab[,3] + object@cluster$gap$Tab[,4])))
          }
          )

setGeneric("plotjaccard", function(object) standardGeneric("plotjaccard"))

setMethod("plotjaccard",
          signature = "SCseq",
          definition = function(object){
            if ( length(object@cluster$kpart) == 0 ) stop("run clustexp before plotjaccard")
            if ( length(unique(object@cluster$kpart)) < 2 ) stop("only a single cluster: no Jaccard's similarity plot")
            barplot(object@cluster$jaccard,names.arg=1:length(object@cluster$jaccard),ylab="Jaccard's similarity")
          }
          )

setGeneric("plotoutlierprobs", function(object) standardGeneric("plotoutlierprobs"))

setMethod("plotoutlierprobs",
          signature = "SCseq",
          definition = function(object){
            if ( length(object@cpart) == 0 ) stop("run findoutliers before plotoutlierprobs")
            p <- object@cluster$kpart[ order(object@cluster$kpart,decreasing=FALSE)]
            x <- object@out$cprobs[names(p)]
            fcol <- object@fcol
            for ( i in 1:max(p) ){
              y <- -log10(x + 2.2e-16)
              y[p != i] <- 0
              if ( i == 1 ) b <- barplot(y,ylim=c(0,max(-log10(x + 2.2e-16))*1.1),col=fcol[i],border=fcol[i],names.arg=FALSE,ylab="-log10prob") else barplot(y,add=TRUE,col=fcol[i],border=fcol[i],names.arg=FALSE,axes=FALSE)
  }
            abline(-log10(object@outlierpar$probthr),0,col="black",lty=2)
            d <- b[2,1] - b[1,1]
            y <- 0
            for ( i in 1:max(p) ) y <- append(y,b[sum(p <=i),1] + d/2)
            axis(1,at=(y[1:(length(y)-1)] + y[-1])/2,lab=1:max(p))
            box()
          }
          )

setGeneric("plotbackground", function(object) standardGeneric("plotbackground"))

setMethod("plotbackground",
          signature = "SCseq",
          definition = function(object){
            if ( length(object@cpart) == 0 ) stop("run findoutliers before plotbackground")
            m <- apply(object@fdata,1,mean)
            v <- apply(object@fdata,1,var)
            fit <- locfit(v~lp(m,nn=.7),family="gamma",maxk=500)
            plot(log2(m),log2(v),pch=20,xlab="log2mean",ylab="log2var",col="grey")
            lines(log2(m[order(m)]),log2(object@background$lvar(m[order(m)],object)),col="red",lwd=2)
            lines(log2(m[order(m)]),log2(fitted(fit)[order(m)]),col="orange",lwd=2,lty=2)
            legend("topleft",legend=substitute(paste("y = ",a,"*x^2 + ",b,"*x + ",d,sep=""),list(a=round(coef(object@background$vfit)[3],2),b=round(coef(object@background$vfit)[2],2),d=round(coef(object@background$vfit)[1],2))),bty="n")
          }
          )

setGeneric("plotsensitivity", function(object) standardGeneric("plotsensitivity"))

setMethod("plotsensitivity",
          signature = "SCseq",
          definition = function(object){
            if ( length(object@cpart) == 0 ) stop("run findoutliers before plotsensitivity")
            plot(log10(object@out$thr), object@out$stest, type="l",xlab="log10 Probability cutoff", ylab="Number of outliers")
            abline(v=log10(object@outlierpar$probthr),col="red",lty=2)
          }
          )

setGeneric("diffgenes", function(object,cl1,cl2,mincount=5) standardGeneric("diffgenes"))

setMethod("diffgenes",
          signature = "SCseq",
          definition = function(object,cl1,cl2,mincount){
            part <- object@cpart
            cl1 <- c(cl1)
            cl2 <- c(cl2)
            if ( length(part) == 0 ) stop("run findoutliers before diffgenes")
            if ( ! is.numeric(mincount) ) stop("mincount has to be a non-negative number") else if (  mincount < 0 ) stop("mincount has to be a non-negative number")
            if ( length(intersect(cl1, part)) < length(unique(cl1)) ) stop( paste("cl1 has to be a subset of ",paste(sort(unique(part)),collapse=","),"\n",sep="") )
            if ( length(intersect(cl2, part)) < length(unique(cl2)) ) stop( paste("cl2 has to be a subset of ",paste(sort(unique(part)),collapse=","),"\n",sep="") )
            f <- apply(object@ndata[,part %in% c(cl1,cl2)],1,max) > mincount
            x <- object@ndata[f,part %in% cl1]
            y <- object@ndata[f,part %in% cl2]
            if ( sum(part %in% cl1) == 1 ) m1 <- x else m1 <- apply(x,1,mean)
            if ( sum(part %in% cl2) == 1 ) m2 <- y else m2 <- apply(y,1,mean)
            if ( sum(part %in% cl1) == 1 ) s1 <- sqrt(x) else s1 <- sqrt(apply(x,1,var))
            if ( sum(part %in% cl2) == 1 ) s2 <- sqrt(y) else s2 <- sqrt(apply(y,1,var))
            
            d <- ( m1 - m2 )/ apply( cbind( s1, s2 ),1,mean )
            names(d) <- rownames(object@ndata)[f]
            if ( sum(part %in% cl1) == 1 ){
              names(x) <- names(d)
              x <- x[order(d,decreasing=TRUE)]
            }else{
              x <- x[order(d,decreasing=TRUE),]
            }
            if ( sum(part %in% cl2) == 1 ){
              names(y) <- names(d)
              y <- y[order(d,decreasing=TRUE)]
            }else{
              y <- y[order(d,decreasing=TRUE),]
            }
            return(list(z=d[order(d,decreasing=TRUE)],cl1=x,cl2=y,cl1n=cl1,cl2n=cl2))
          }
          )

plotdiffgenes <- function(z,gene=g){
  if ( ! is.list(z) ) stop("first arguments needs to be output of function diffgenes")
  if ( length(z) < 5 ) stop("first arguments needs to be output of function diffgenes")
  if ( sum(names(z) == c("z","cl1","cl2","cl1n","cl2n")) < 5 ) stop("first arguments needs to be output of function diffgenes")
  if ( length(gene) > 1 ) stop("only single value allowed for argument gene")
  if ( !is.numeric(gene) & !(gene %in% names(z$z)) ) stop("argument gene needs to be within rownames of first argument or a positive integer number")
  if ( is.numeric(gene) ){ if ( gene < 0 | round(gene) != gene ) stop("argument gene needs to be within rownames of first argument or a positive integer number") }
  x <- if ( is.null(dim(z$cl1)) ) z$cl1[gene] else t(z$cl1[gene,])
  y <- if ( is.null(dim(z$cl2)) ) z$cl2[gene] else t(z$cl2[gene,])
  plot(1:length(c(x,y)),c(x,y),ylim=c(0,max(c(x,y))),xlab="",ylab="Expression",main=gene,cex=0,axes=FALSE)
  axis(2)
  box()
  u <- 1:length(x)
  rect(u - .5,0,u + .5,x,col="red")
  v <- c(min(u) - .5,max(u) + .5)
  axis(1,at=mean(v),lab=paste(z$cl1n,collapse=","))
  lines(v,rep(mean(x),length(v)))
  lines(v,rep(mean(x)-sqrt(var(x)),length(v)),lty=2)
  lines(v,rep(mean(x)+sqrt(var(x)),length(v)),lty=2)
  
  u <- ( length(x) + 1 ):length(c(x,y))
  v <- c(min(u) - .5,max(u) + .5)
  rect(u - .5,0,u + .5,y,col="blue")
  axis(1,at=mean(v),lab=paste(z$cl2n,collapse=","))
  lines(v,rep(mean(y),length(v)))
  lines(v,rep(mean(y)-sqrt(var(y)),length(v)),lty=2)
  lines(v,rep(mean(y)+sqrt(var(y)),length(v)),lty=2)
  abline(v=length(x) + .5)
}

setGeneric("plottsne", function(object,final=TRUE) standardGeneric("plottsne"))

setMethod("plottsne",
          signature = "SCseq",
          definition = function(object,final){
            if ( length(object@tsne) == 0 ) stop("run comptsne before plottsne")
            if ( final & length(object@cpart) == 0 ) stop("run findoutliers before plottsne")
            if ( !final & length(object@cluster$kpart) == 0 ) stop("run clustexp before plottsne")
            part <- if ( final ) object@cpart else object@cluster$kpart
            plot(object@tsne,xlab="Dim 1",ylab="Dim 2",pch=20,cex=1.5,col="lightgrey")
            for ( i in 1:max(part) ){
              if ( sum(part == i) > 0 ) text(object@tsne[part == i,1],object@tsne[part == i,2],i,col=object@fcol[i],cex=.75,font=4)
            }
          }
          )

setGeneric("plotlabelstsne", function(object,labels=NULL) standardGeneric("plotlabelstsne"))

setMethod("plotlabelstsne",
          signature = "SCseq",
          definition = function(object,labels){
            if ( is.null(labels ) ) labels <- names(object@ndata)
            if ( length(object@tsne) == 0 ) stop("run comptsne before plotlabelstsne")
            plot(object@tsne,xlab="Dim 1",ylab="Dim 2",pch=20,cex=1.5,col="lightgrey")
            text(object@tsne[,1],object@tsne[,2],labels,cex=.5)
          }
          )

setGeneric("plotsymbolstsne", function(object,types=NULL) standardGeneric("plotsymbolstsne"))

setMethod("plotsymbolstsne",
          signature = "SCseq",
          definition = function(object,types){
            if ( is.null(types) ) types <- names(object@fdata)
            if ( length(object@tsne) == 0 ) stop("run comptsne before plotsymbolstsne")
            if ( length(types) != ncol(object@fdata) ) stop("types argument has wrong length. Length has to equal to the column number of object@ndata")
            coloc <- rainbow(length(unique(types)))
            syms <- c()
            plot(object@tsne,xlab="Dim 1",ylab="Dim 2",pch=20,col="grey")
            for ( i in 1:length(unique(types)) ){
              f <- types == sort(unique(types))[i]
              syms <- append( syms, ( (i-1) %% 25 ) + 1 )
              points(object@tsne[f,1],object@tsne[f,2],col=coloc[i],pch=( (i-1) %% 25 ) + 1,cex=1)
            }
            legend("topleft", legend=sort(unique(types)), col=coloc, pch=syms)
          }
          )

setGeneric("plotexptsne", function(object,g,n="",logsc=FALSE) standardGeneric("plotexptsne"))

setMethod("plotexptsne",
          signature = "SCseq",
          definition = function(object,g,n="",logsc=FALSE){
            if ( length(object@tsne) == 0 ) stop("run comptsne before plottsne")
            if ( length(intersect(g,rownames(object@ndata))) < length(unique(g)) ) stop("second argument does not correspond to set of rownames slot ndata of SCseq object")
            if ( !is.numeric(logsc) & !is.logical(logsc) ) stop("argument logsc has to be logical (TRUE/FALSE)")
            if ( n == "" ) n <- g[1]
            l <- apply(object@ndata[g,] - .1,2,sum) + .1
            if (logsc) {
              f <- l == 0
              l <- log2(l)
              l[f] <- NA
            }
            mi <- min(l,na.rm=TRUE)
            ma <- max(l,na.rm=TRUE)
            #ColorRamp <- colorRampPalette(rev(brewer.pal(n = 7,name = "RdYlBu")))(100)
            ColorRamp <-  colorRampPalette((brewer.pal(n = 9,name = "YlOrRd")))(100)
            ColorRamp <-  colorRampPalette((brewer.pal(n = 9,name = "RdPu")))(100)
            ColorLevels <- seq(mi, ma, length=length(ColorRamp))
            v <- round((l - mi)/(ma - mi)*99 + 1,0)
            layout(matrix(data=c(1,3,2,4), nrow=2, ncol=2), widths=c(5,1,5,1), heights=c(5,1,1,1))
            par(mar = c(3,5,2.5,2))
            plot(object@tsne,xlab="Dim 1",ylab="Dim 2",main=n,pch=20,cex=0,col="grey")
            kk <- order(v,decreasing=F)
            for ( k in kk ){
              points(object@tsne[k,1],object@tsne[k,2],col=ColorRamp[v[k]],pch=20,cex=0.75)
            }
            par(mar = c(3,2.5,2.5,2))
            image(1, ColorLevels,
                  matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
                  col=ColorRamp,
                  xlab="",ylab="",
                  xaxt="n")
            layout(1)
          }
          )


plot.err.bars.y <- function(x, y, y.err, col="black", lwd=1, lty=1, h=0.1){
  arrows(x,y-y.err,x,y+y.err,code=0, col=col, lwd=lwd, lty=lty)
  arrows(x-h,y-y.err,x+h,y-y.err,code=0, col=col, lwd=lwd, lty=lty)
  arrows(x-h,y+y.err,x+h,y+y.err,code=0, col=col, lwd=lwd, lty=lty)
}

clusGapExt <-function (x, FUNcluster, K.max, B = 100, verbose = interactive(), method="euclidean",random=TRUE,
    ...) 
{
     stopifnot(is.function(FUNcluster), length(dim(x)) == 2, K.max >= 
        2, (n <- nrow(x)) >= 1, (p <- ncol(x)) >= 1)
    if (B != (B. <- as.integer(B)) || (B <- B.) <= 0) 
        stop("'B' has to be a positive integer")
    if (is.data.frame(x)) 
        x <- as.matrix(x)
    ii <- seq_len(n)
    W.k <- function(X, kk) {
        clus <- if (kk > 1) 
            FUNcluster(X, kk, ...)$cluster
        else rep.int(1L, nrow(X))
        0.5 * sum(vapply(split(ii, clus), function(I) {
            xs <- X[I, , drop = FALSE]
            sum(dist.gen(xs,method=method)/nrow(xs))
        }, 0))
    }
    logW <- E.logW <- SE.sim <- numeric(K.max)
    if (verbose) 
        cat("Clustering k = 1,2,..., K.max (= ", K.max, "): .. ", 
            sep = "")
    for (k in 1:K.max) logW[k] <- log(W.k(x, k))
    if (verbose) 
        cat("done\n")
    xs <- scale(x, center = TRUE, scale = FALSE)
    m.x <- rep(attr(xs, "scaled:center"), each = n)
    V.sx <- svd(xs)$v
    rng.x1 <- apply(xs %*% V.sx, 2, range)
    logWks <- matrix(0, B, K.max)
     if (random){
       if (verbose) 
         cat("Bootstrapping, b = 1,2,..., B (= ", B, ")  [one \".\" per sample]:\n", 
             sep = "")
       for (b in 1:B) {
         z1 <- apply(rng.x1, 2, function(M, nn) runif(nn, min = M[1], 
             max = M[2]), nn = n)
         z <- tcrossprod(z1, V.sx) + m.x
         ##z <- apply(x,2,function(m) runif(length(m),min=min(m),max=max(m)))
         ##z <- apply(x,2,function(m) sample(m))
         for (k in 1:K.max) {
           logWks[b, k] <- log(W.k(z, k))
         }
         if (verbose) 
           cat(".", if (b%%50 == 0) 
               paste(b, "\n"))
       }
       if (verbose && (B%%50 != 0)) 
         cat("", B, "\n")
       E.logW <- colMeans(logWks)
       SE.sim <- sqrt((1 + 1/B) * apply(logWks, 2, var))
     }else{
       E.logW <- rep(NA,K.max)
       SE.sim <- rep(NA,K.max)
     }
    structure(class = "clusGap", list(Tab = cbind(logW, E.logW, 
        gap = E.logW - logW, SE.sim), n = n, B = B, FUNcluster = FUNcluster))
}

clustfun <- function(x,clustnr=20,bootnr=50,metric="pearson",do.gap=FALSE,sat=TRUE,SE.method="Tibs2001SEmax",SE.factor=.25,B.gap=50,cln=0,rseed=17000,FUNcluster="kmedoids",distances=NULL,link="single")
{
  if ( clustnr < 2) stop("Choose clustnr > 1")
  if ( FUNcluster %in% c("kmedoids","clara") ){
    di  <- t(x)
    diM <- dist.gen(di,method=metric)
  }else{
    di <- dist.gen(t(x),method=metric)
  }
  if ( FUNcluster %in% c("clara") ){
    diCL <- as.data.frame( cmdscale(as.matrix(diM),k=nrow(di)-1) )
    rownames(diCL) <- names(x)
  }
  if ( nrow(di) - 1 < clustnr ) clustnr <- nrow(di) - 1
  if ( do.gap | sat | cln > 0 ){
    gpr <- NULL
    f <- if ( cln == 0 ) TRUE else FALSE
    if ( do.gap ){
      set.seed(rseed)
      if ( FUNcluster == "kmeans" )   gpr <- clusGapExt(as.matrix(di), FUN = kmeans, K.max = clustnr, B = B.gap, iter.max=100)
      if ( FUNcluster == "kmedoids" ) gpr <- clusGapExt(as.matrix(di), FUN = function(x,k) pam(dist.gen(x,method=metric),k), K.max = clustnr, B = B.gap, method=metric)
      if ( FUNcluster == "clara" ) gpr <- clusGapExt(diCL, FUN = function(x,k) clara(x,k), K.max = clustnr, B = B.gap, method="euclidean")
      if ( FUNcluster == "hclust" )   gpr <- clusGapExt(as.matrix(di), FUN = function(x,k){ y <- hclusterCBI(x,k,link=link,scaling=FALSE); y$cluster <- y$partition; y }, K.max = clustnr, B = B.gap) 
      if ( f ) cln <- maxSE(gpr$Tab[,3],gpr$Tab[,4],method=SE.method,SE.factor)
    }
    if ( sat ){
      if ( ! do.gap ){
        if ( FUNcluster == "kmeans" )   gpr <- clusGapExt(as.matrix(di), FUN = kmeans, K.max = clustnr, B = B.gap, iter.max=100, random=FALSE)
        if ( FUNcluster == "kmedoids" ) gpr <- clusGapExt(as.matrix(di), FUN = function(x,k) pam(dist.gen(x,method=metric),k), K.max = clustnr, B = B.gap, random=FALSE, method=metric)
        if ( FUNcluster == "clara" ) gpr <- clusGapExt(as.matrix(diCL), FUN = function(x,k) clara(x,k), K.max = clustnr, B = B.gap, random=FALSE, method="euclidean")
        if ( FUNcluster == "hclust" )   gpr <- clusGapExt(as.matrix(di), FUN = function(x,k){ y <- hclusterCBI(x,k,link=link,scaling=FALSE); y$cluster <- y$partition; y }, K.max = clustnr, B = B.gap, random=FALSE)
      }
      g <- gpr$Tab[,1]
      y <- g[-length(g)] - g[-1]
      mm <- numeric(length(y))
      nn <- numeric(length(y))
      for ( i in 1:length(y)){
        mm[i] <- mean(y[i:length(y)]) 
        nn[i] <- sqrt(var(y[i:length(y)]))
      }
      if ( f ) cln <- max(min(which( y - (mm + nn) < 0 )),1)
    }
    if ( cln <= 1 ) {
      clb <- list(result=list(partition=rep(1,dim(x)[2])),bootmean=1)
      names(clb$result$partition) <- names(x)
      return(list(x=x,clb=clb,gpr=gpr,di=if ( FUNcluster %in% c("kmedoids","clara") ) dist.gen(di,method=metric) else di))
    }
    if ( FUNcluster == "kmeans" ) clb <- clusterboot(di,B=bootnr,distances=FALSE,bootmethod="boot",clustermethod=kmeansCBI,krange=cln,scaling=FALSE,multipleboot=FALSE,bscompare=TRUE,seed=rseed)
    ##if ( FUNcluster == "kmedoids" ) clb <- clusterboot(diM,B=bootnr,bootmethod="boot",clustermethod=pamkCBI,k=cln,multipleboot=FALSE,bscompare=TRUE,seed=rseed)
    if ( FUNcluster == "kmedoids" ) clb <- clusterboot(diM,B=bootnr,bootmethod="boot",clustermethod=pamkdCBI,scaling=FALSE,diss=TRUE,k=cln,multipleboot=FALSE,bscompare=TRUE,seed=rseed)
    if ( FUNcluster == "clara" ) clb <- clusterboot(diCL,B=bootnr,bootmethod="boot",clustermethod=pamkdCBI,usepam=FALSE,scaling=FALSE,krange=cln,multipleboot=FALSE,bscompare=TRUE,seed=rseed)
    if ( FUNcluster == "hclust" ) clb <- clusterboot(di,B=bootnr,distances=FALSE,bootmethod="boot",clustermethod=hclusterCBI,k=cln,link=link,scaling=FALSE,multipleboot=FALSE,bscompare=TRUE,seed=rseed)
    return(list(x=x,clb=clb,gpr=gpr,di=if ( FUNcluster %in% c("kmedoids","clara") ) diM else di))
  }
}

setGeneric("clustexp", function(object,clustnr=20,bootnr=50,metric="pearson",do.gap=FALSE,sat=TRUE,SE.method="Tibs2001SEmax",SE.factor=.25,B.gap=50,cln=0,rseed=17000,FUNcluster="kmedoids",FSelect=FALSE) standardGeneric("clustexp"))

setMethod("clustexp",
          signature = "SCseq",
          definition = function(object,clustnr,bootnr,metric,do.gap,sat,SE.method,SE.factor,B.gap,cln,rseed,FUNcluster,FSelect) {
            if ( ! is.numeric(clustnr) ) stop("clustnr has to be a positive integer") else if ( round(clustnr) != clustnr | clustnr <= 0 ) stop("clustnr has to be a positive integer")
            if ( ! is.numeric(bootnr) ) stop("bootnr has to be a positive integer") else if ( round(bootnr) != bootnr | bootnr <= 0 ) stop("bootnr has to be a positive integer")
            if ( ! ( metric %in% c( "spearman","pearson","eupearson","logpearson","kendall","euclidean","maximum","manhattan","canberra","binary","minkowski","relentr") ) ) stop("metric has to be one of the following: spearman, pearson, eupearson, logpearson, kendall, euclidean, maximum, manhattan, canberra, binary, minkowski","relentr")
            if ( ! ( SE.method %in% c( "firstSEmax","Tibs2001SEmax","globalSEmax","firstmax","globalmax") ) ) stop("SE.method has to be one of the following: firstSEmax, Tibs2001SEmax, globalSEmax, firstmax, globalmax")
            if ( ! is.numeric(SE.factor) ) stop("SE.factor has to be a non-negative integer") else if  ( SE.factor < 0 )  stop("SE.factor has to be a non-negative integer")
            if ( ! ( is.numeric(do.gap) | is.logical(do.gap) ) ) stop( "do.gap has to be logical (TRUE/FALSE)" )
            if ( ! ( is.numeric(sat) | is.logical(sat) ) ) stop( "sat has to be logical (TRUE/FALSE)" )
            if ( ! is.numeric(B.gap) ) stop("B.gap has to be a positive integer") else if ( round(B.gap) != B.gap | B.gap <= 0 ) stop("B.gap has to be a positive integer")
            if ( ! is.numeric(cln) ) stop("cln has to be a non-negative integer") else if ( round(cln) != cln | cln < 0 ) stop("cln has to be a non-negative integer")          
            if ( ! is.numeric(rseed) ) stop("rseed has to be numeric")
            if ( !do.gap & !sat & cln == 0 ) stop("cln has to be a positive integer or either do.gap or sat has to be TRUE")
            if ( ! ( FUNcluster %in% c("kmeans","hclust","kmedoids","clara") ) ) stop("FUNcluster has to be one of the following: kmeans, hclust,kmedoids","clara")
            if ( ! ( is.numeric(FSelect) | is.logical(FSelect) ) ) stop( "FSelect has to be logical (TRUE/FALSE)" )
            object@clusterpar <- list(clustnr=clustnr,bootnr=bootnr,metric=metric,do.gap=do.gap,sat=sat,SE.method=SE.method,SE.factor=SE.factor,B.gap=B.gap,cln=cln,rseed=rseed,FUNcluster=FUNcluster,FSelect=FSelect)

            if ( FSelect ){
              x <- object@fdata
              m <- apply(x,1,mean)
              v <- apply(x,1,var )

              ml <- log2(m)
              vl <- log2(v)
              f <- ml > -Inf & vl > -Inf
              ml <- ml[f]
              vl <- vl[f]
              mm <- -8
              repeat{
                fit <- lm(vl ~ ml + I(ml^2)) 
                if( coef(fit)[3] >= 0 | mm >= -1){
                  break
                }
                mm <- mm + .5
                f <- ml > mm
                ml <- ml[f]
                vl <- vl[f]
              }
              lvar <- function(x) 2**(coef(fit)[1] + log2(x)*coef(fit)[2] + coef(fit)[3] * log2(x)**2)
              vln <- log2(v)  - log2(sapply(m,FUN=lvar))
              n <- names(vln)[vln>0]
              ##n <- names(head(vln[order(vln,decreasing=T)],round(length(vln)/5,0)))
            }else{
              n <- rownames(object@fdata)
            }
            
            y <- clustfun(object@fdata[n,],clustnr,bootnr,metric,do.gap,sat,SE.method,SE.factor,B.gap,cln,rseed,FUNcluster)
            object@cluster   <- list(kpart=y$clb$result$partition, jaccard=y$clb$bootmean, gap=y$gpr, clb=y$clb, features=n)
            object@distances <- as.matrix( y$di )
            set.seed(111111)
            object@fcol <- sample(rainbow(max(y$clb$result$partition)))
            return(object)
          }
          )

setGeneric("findoutliers", function(object,outminc=5,outlg=2,probthr=1e-3,thr=2**-(1:40),outdistquant=.95) standardGeneric("findoutliers"))

setMethod("findoutliers",
          signature = "SCseq",
          definition = function(object,outminc,outlg,probthr,thr,outdistquant) {
            if ( length(object@cluster$kpart) == 0 ) stop("run clustexp before findoutliers")
            if ( ! is.numeric(outminc) ) stop("outminc has to be a non-negative integer") else if ( round(outminc) != outminc | outminc < 0 ) stop("outminc has to be a non-negative integer")
            if ( ! is.numeric(outlg) ) stop("outlg has to be a non-negative integer") else if ( round(outlg) != outlg | outlg < 0 ) stop("outlg has to be a non-negative integer")
            if ( ! is.numeric(probthr) ) stop("probthr has to be a number between 0 and 1") else if (  probthr < 0 | probthr > 1 ) stop("probthr has to be a number between 0 and 1")
            if ( ! is.numeric(thr) ) stop("thr hast to be a vector of numbers between 0 and 1") else if ( min(thr) < 0 | max(thr) > 1 ) stop("thr hast to be a vector of numbers between 0 and 1")
            if ( ! is.numeric(outdistquant) ) stop("outdistquant has to be a number between 0 and 1") else if (  outdistquant < 0 | outdistquant > 1 ) stop("outdistquant has to be a number between 0 and 1")
            
            object@outlierpar <- list( outminc=outminc,outlg=outlg,probthr=probthr,thr=thr,outdistquant=outdistquant )
            ### calibrate background model
            m <- log2(apply(object@fdata,1,mean))
            v <- log2(apply(object@fdata,1,var))
            f <- m > -Inf & v > -Inf
            m <- m[f]
            v <- v[f]
            mm <- -8
            repeat{
              fit <- lm(v ~ m + I(m^2)) 
              ##if( coef(fit)[3] >= 0 | mm >= 3){
              if( coef(fit)[3] >= 0 | mm >= -1){
                break
              }
              mm <- mm + .5
              f <- m > mm
              m <- m[f]
              v <- v[f]
            }
            object@background <- list()
            object@background$vfit <- fit
            object@background$lvar <- function(x,object) 2**(coef(object@background$vfit)[1] + log2(x)*coef(object@background$vfit)[2] + coef(object@background$vfit)[3] * log2(x)**2)
            object@background$lsize <- function(x,object) x**2/(max(x + 1e-6,object@background$lvar(x,object)) - x)

            ### identify outliers
            out   <- c()
            stest <- rep(0,length(thr))
            cprobs <- c()
            outgene <- list()
            for ( n in 1:max(object@cluster$kpart) ){     
              if ( sum(object@cluster$kpart == n) == 1 ){ cprobs <- append(cprobs,.5); names(cprobs)[length(cprobs)] <- names(object@cluster$kpart)[object@cluster$kpart == n]; next }
              x <- object@fdata[,object@cluster$kpart == n]
              x <- x[apply(x,1,max) > outminc,]
              z <- t( apply(x,1,function(x){ apply( cbind( pnbinom(round(x,0),mu=mean(x),size=object@background$lsize(mean(x),object)) , 1 - pnbinom(round(x,0),mu=mean(x),size=object@background$lsize(mean(x),object)) ),1, min) } ) )
              cp <- apply(z,2,function(x){ y <- p.adjust(x,method="BH"); y <- y[order(y,decreasing=FALSE)]; return(y[outlg]);})
              f <- cp < probthr
              cprobs <- append(cprobs,cp)
              if ( sum(f) > 0 ) out <- append(out,names(x)[f])
              for ( j in 1:length(thr) )  stest[j] <-  stest[j] + sum( cp < thr[j] )
              fg <- apply(z,1,min) < probthr
              outgene[[n]] <- if ( sum(fg) > 0 ) z[fg,] else 0
            }
            object@out <-list(out=out,stest=stest,thr=thr,cprobs=cprobs,outgene=outgene)

            ### cluster outliers
            clp2p.cl <- c()
            cols <- names(object@fdata)
            cpart <- object@cluster$kpart
            di   <- as.data.frame(object@distances)
            for ( i in 1:max(cpart) ) {
              tcol <- cols[cpart == i]
              if ( sum(!(tcol %in% out)) > 1 ) clp2p.cl <- append(clp2p.cl,as.vector(t(di[tcol[!(tcol %in% out)],tcol[!(tcol %in% out)]])))
            }
            clp2p.cl <- clp2p.cl[clp2p.cl>0]
              
            cadd  <- list()
            if ( length(out) > 0 ){
              if (length(out) == 1){
                cadd <- list(out)
              }else{
                n <- out
                m <- as.data.frame(di[out,out])
                
                for ( i in 1:length(out) ){
                  if ( length(n) > 1 ){
                    o   <- order(apply(cbind(m,1:dim(m)[1]),1,function(x)  min(x[1:(length(x)-1)][-x[length(x)]])),decreasing=FALSE)
                    m <- m[o,o]
                    n <- n[o]          
                    f <- m[,1] < quantile(clp2p.cl,outdistquant) | m[,1] == min(clp2p.cl)
                    ind <- 1  
                    if ( sum(f) > 1 ) for ( j in 2:sum(f) ) if ( apply(m[f,f][j,c(ind,j)] > quantile(clp2p.cl,outdistquant) ,1,sum) == 0 ) ind <- append(ind,j)
                    cadd[[i]] <- n[f][ind]
                    g <- ! n %in% n[f][ind]
                    n <- n[g]
                    m <- m[g,g]
                    if ( sum(g) == 0 ) break
          
                  }else if (length(n) == 1){
                    cadd[[i]] <- n
                    break
                  }
                }
              }
    
              for ( i in 1:length(cadd) ){
                cpart[cols %in% cadd[[i]]] <- max(cpart) + 1
              }
            }

            ### determine final clusters
            for ( i in min(cpart):max(cpart) ){
              if ( sum(cpart == i) == 0 ) next
              f <- cols[cpart == i]
              d <- object@fdata[object@cluster$features,]
              if ( length(f) == 1 ){
                cent <- d[,f]
                md <- f
              }else{
                if ( object@clusterpar$FUNcluster %in% c("kmedoids","clara") ){
                  ##x <- apply(as.matrix(dist.gen(t(d[,f]),method=object@clusterpar$metric)),2,mean)
                  x <- apply(object@distances[f,f],2,mean)
                  md <- f[which(x == min(x))[1]]
                  cent <- d[,md]
                }else{
                  cent <- apply(d[,f],1,mean)
                }
              }
              if ( i == min(cpart) ) dcent <- data.frame(cent) else dcent <- cbind(dcent,cent)
              if ( object@clusterpar$FUNcluster %in% c("kmedoids","clara") ){
                if ( i == min(cpart) ) tmp <- data.frame(object@distances[,md]) else tmp <- cbind(tmp,object@distances[,md])
              }else{
                if ( i == min(cpart) ) tmp <- data.frame(apply(d,2,dist.gen.pairs,y=cent,method=object@clusterpar$metric)) else tmp <- cbind(tmp,apply(d,2,dist.gen.pairs,y=cent,method=object@clusterpar$metric))
              }
            }
            cpart <- apply(tmp,1,function(x) order(x,decreasing=FALSE)[1])
            
            for  ( i in max(cpart):1){if (sum(cpart==i)==0) cpart[cpart>i] <- cpart[cpart>i] - 1 }

            object@cpart <- cpart

            set.seed(111111)
            object@fcol <- sample(rainbow(max(cpart)))
            return(object)
          }
        )


setGeneric("comptsne", function(object,rseed=15555,sammonmap=FALSE,initial_cmd=TRUE,fast=FALSE,perplexity=30,...) standardGeneric("comptsne"))

setMethod("comptsne",
          signature = "SCseq",
          definition = function(object,rseed,sammonmap,initial_cmd,fast,perplexity,...){
            if ( length(object@cluster$kpart) == 0 ) stop("run clustexp before comptsne")
            set.seed(rseed)
            di <- if ( object@clusterpar$FUNcluster %in% c("kmedoids","clara")) as.dist(object@distances) else dist.gen(as.matrix(object@distances))
            if ( sammonmap ){
              object@tsne <- as.data.frame(sammon(di,k=2)$points)
            }else{
              if ( fast ){
                ts <- if ( initial_cmd ) Rtsne(di,dims=2,initial_config=cmdscale(di,k=2),perplexity=perplexity,...)$Y else Rtsne(di,dims=2,perplexity=perplexity,...)$Y
              }else{
                ts <- if ( initial_cmd ) tsne(di,k=2,initial_config=cmdscale(di,k=2),perplexity=perplexity,...) else tsne(di,k=2,perplexity=perplexity,...)
              }
              object@tsne <- as.data.frame(ts)
            }
            return(object)
          }
          )


setGeneric("clustdiffgenes", function(object,pvalue=.01) standardGeneric("clustdiffgenes"))

setMethod("clustdiffgenes",
          signature = "SCseq",
          definition = function(object,pvalue){
            if ( length(object@cpart) == 0 ) stop("run findoutliers before clustdiffgenes")
            if ( ! is.numeric(pvalue) ) stop("pvalue has to be a number between 0 and 1") else if (  pvalue < 0 | pvalue > 1 ) stop("pvalue has to be a number between 0 and 1")
            cdiff <- list()
            x     <- object@ndata
            y     <- object@expdata[,names(object@ndata)]
            part  <- object@cpart
            for ( i in 1:max(part) ){
              if ( sum(part == i) == 0 ) next
              m <-  if ( sum(part != i) > 1 ) apply(x[,part != i],1,mean) else x[,part != i]
              n <-  if ( sum(part == i) > 1 ) apply(x[,part == i],1,mean) else x[,part == i]
              no <- if ( sum(part == i) > 1 ) median(apply(y[,part == i],2,sum))/median(apply(x[,part == i],2,sum)) else sum(y[,part == i])/sum(x[,part == i])
              m <- m*no
              n <- n*no
              pv <- binompval(m/sum(m),sum(n),n)
              d <- data.frame(mean.ncl=m,mean.cl=n,fc=n/m,pv=pv)[order(pv,decreasing=FALSE),]
              cdiff[[paste("cl",i,sep=".")]] <- d[d$pv < pvalue,]
            }
            return(cdiff)
          }
          )

setGeneric("plotsaturation", function(object,disp=FALSE) standardGeneric("plotsaturation"))

setMethod("plotsaturation",
          signature = "SCseq",
          definition = function(object,disp){
            if ( length(object@cluster$gap) == 0 ) stop("run clustexp before plotsaturation")
            g <- object@cluster$gap$Tab[,1]
            y <- g[-length(g)] - g[-1]
            mm <- numeric(length(y))
            nn <- numeric(length(y))
            for ( i in 1:length(y)){
              mm[i] <- mean(y[i:length(y)]) 
              nn[i] <- sqrt(var(y[i:length(y)]))
            }
            cln <- max(min(which( y - (mm + nn) < 0 )),1)
            x <- 1:length(y)
            if (disp){
              x <- 1:length(g)
              plot(x,g,pch=20,col="grey",xlab="k",ylab="log within cluster dispersion")
              f <- x == cln
              points(x[f],g[f],col="blue")
            }else{
              plot(x,y,pch=20,col="grey",xlab="k",ylab="Change in log within cluster dispersion")
              points(x,mm,col="red",pch=20)
              plot.err.bars.y(x,mm,nn,col="red")
              points(x,y,col="grey",pch=20)
              f <- x == cln
              points(x[f],y[f],col="blue")
            }
          }
          )

setGeneric("plotsilhouette", function(object) standardGeneric("plotsilhouette"))

setMethod("plotsilhouette",
          signature = "SCseq",
          definition = function(object){
            if ( length(object@cluster$kpart) == 0 ) stop("run clustexp before plotsilhouette")
            if ( length(unique(object@cluster$kpart)) < 2 ) stop("only a single cluster: no silhouette plot")
            kpart <- object@cluster$kpart
            distances  <- if ( object@clusterpar$FUNcluster == "kmedoids" ) as.dist(object@distances) else dist.gen(object@distances)
            si <- silhouette(kpart,distances)
            plot(si)
          }
          )

compmedoids <- function(x,part,metric="pearson"){
  m <- c()
  for ( i in sort(unique(part)) ){
    f <- names(x)[part == i]
    if ( length(f) == 1 ){
      m <- append(m,f)
    }else{
      y <- apply(as.matrix(dist.gen(t(x[,f]),method=metric)),2,mean)
      m <- append(m,f[which(y == min(y))[1]])
    }
  }
  m
}

setGeneric("clustheatmap", function(object,final=FALSE,hmethod="single") standardGeneric("clustheatmap"))

setMethod("clustheatmap",
          signature = "SCseq",
          definition = function(object,final,hmethod){
            if ( final & length(object@cpart) == 0 ) stop("run findoutliers before clustheatmap")
            if ( !final & length(object@cluster$kpart) == 0 ) stop("run clustexp before clustheatmap")
            x <- object@fdata  
            part <- if ( final ) object@cpart else object@cluster$kpart
            na <- c()
            j <- 0
            for ( i in 1:max(part) ){
              if ( sum(part == i) == 0 ) next
              j <- j + 1
              na <- append(na,i)
              d <- x[,part == i]
              if ( sum(part == i) == 1 ) cent <- d else cent <- apply(d,1,mean)
              if ( j == 1 ) tmp <- data.frame(cent) else tmp <- cbind(tmp,cent)
            }
            names(tmp) <- paste("cl",na,sep=".")
            ld <- if ( object@clusterpar$FUNcluster == "kmedoids" ) dist.gen(t(tmp),method=object@clusterpar$metric) else dist.gen(as.matrix(dist.gen(t(tmp),method=object@clusterpar$metric)))
            if ( max(part) > 1 )  cclmo <- hclust(ld,method=hmethod)$order else cclmo <- 1
            q <- part
            for ( i in 1:max(part) ){
              q[part == na[cclmo[i]]] <- i
            }
            part <- q
            di <-  if ( object@clusterpar$FUNcluster == "kmedoids" ) object@distances else as.data.frame( as.matrix( dist.gen(t(object@distances)) ) )
            pto <- part[order(part,decreasing=FALSE)]
            ptn <- c()
            for ( i in 1:max(pto) ){ pt <- names(pto)[pto == i]; z <- if ( length(pt) == 1 ) pt else pt[hclust(as.dist(t(di[pt,pt])),method=hmethod)$order]; ptn <- append(ptn,z) }
            col <- object@fcol
            mi  <- min(di,na.rm=TRUE)
            ma  <- max(di,na.rm=TRUE)
            layout(matrix(data=c(1,3,2,4), nrow=2, ncol=2), widths=c(5,1,5,1), heights=c(5,1,1,1))
            ColorRamp   <- colorRampPalette(brewer.pal(n = 7,name = "RdYlBu"))(100)
            ColorLevels <- seq(mi, ma, length=length(ColorRamp))
            if ( mi == ma ){
              ColorLevels <- seq(0.99*mi, 1.01*ma, length=length(ColorRamp))
            }
            par(mar = c(3,5,2.5,2))
            image(as.matrix(di[ptn,ptn]),col=ColorRamp,axes=FALSE)
            abline(0,1)
            box()
            
            tmp <- c()
            for ( u in 1:max(part) ){
              ol <- (0:(length(part) - 1)/(length(part) - 1))[ptn %in% names(x)[part == u]]
              points(rep(0,length(ol)),ol,col=col[cclmo[u]],pch=15,cex=.75)
              points(ol,rep(0,length(ol)),col=col[cclmo[u]],pch=15,cex=.75)
              tmp <- append(tmp,mean(ol))
            }
            axis(1,at=tmp,lab=cclmo)
            axis(2,at=tmp,lab=cclmo)
            par(mar = c(3,2.5,2.5,2))
            image(1, ColorLevels,
                  matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
                  col=ColorRamp,
                  xlab="",ylab="",
                  xaxt="n")
            layout(1)
            return(cclmo)
          }
          )






## class definition
Ltree <- setClass("Ltree", slots = c(sc = "SCseq", ldata = "list", entropy = "vector", trproj = "list", par = "list", prback = "data.frame", prbacka = "data.frame", ltcoord = "matrix", prtree = "list", sigcell = "vector", cdata = "list"  ))

setValidity("Ltree",
            function(object) {
              msg <- NULL
              if ( class(object@sc)[1] != "SCseq" ){
                msg <- c(msg, "input data must be of class SCseq")
              }
              if (is.null(msg)) TRUE
              else msg
            }
            )

setMethod("initialize",
          signature = "Ltree",
          definition = function(.Object, sc ){
            .Object@sc <- sc
            validObject(.Object)
            return(.Object)
          }
          )

setGeneric("compentropy", function(object) standardGeneric("compentropy"))

setMethod("compentropy",
          signature = "Ltree",
          definition = function(object){
            probs   <- t(t(object@sc@ndata)/apply(object@sc@ndata,2,sum))
            object@entropy <- -apply(probs*log(probs)/log(nrow(object@sc@ndata)),2,sum)
            return(object)
          }            
          )


compproj <- function(pdiloc,lploc,cnloc,mloc,d=NULL){
  pd    <- data.frame(pdiloc)
  k     <- paste("X",sort(rep(1:nrow(pdiloc),length(mloc))),sep="")
  pd$k  <- paste("X",1:nrow(pdiloc),sep="")
  pd    <- merge(data.frame(k=k),pd,by="k")
 
  if ( is.null(d) ){
    cnv   <- t(matrix(rep(t(cnloc),nrow(pdiloc)),nrow=ncol(pdiloc)))
    pdcl  <- paste("X",lploc[as.numeric(sub("X","",pd$k))],sep="")
    rownames(cnloc) <- paste("X",mloc,sep="")
    pdcn  <- cnloc[pdcl,]
    v     <- cnv - pdcn
  }else{
    v    <- d$v
    pdcn <- d$pdcn
  }
  w <- pd[,names(pd) != "k"] - pdcn
  
  h <- apply(cbind(v,w),1,function(x){
    x1 <- x[1:(length(x)/2)];
    x2 <- x[(length(x)/2 + 1):length(x)];
    x1s <- sqrt(sum(x1**2)); x2s <- sqrt(sum(x2**2)); y <- sum(x1*x2)/x1s/x2s; return( if (x1s == 0 | x2s == 0 ) NA else y ) }) 
  
  rma <- as.data.frame(matrix(h,ncol=nrow(pdiloc)))
  names(rma) <- unique(pd$k)
  pdclu  <- lploc[as.numeric(sub("X","",names(rma)))]
  pdclp  <- apply(t(rma),1,function(x) if (sum(!is.na(x)) == 0 ) NA else mloc[which(abs(x) == max(abs(x),na.rm=TRUE))[1]])
  pdclh  <- apply(t(rma),1,function(x) if (sum(!is.na(x)) == 0 ) NA else x[which(abs(x) == max(abs(x),na.rm=TRUE))[1]])
  pdcln  <-  names(lploc)[as.numeric(sub("X","",names(rma)))]
  names(rma) <- pdcln
  rownames(rma) <- paste("X",mloc,sep="")
  res    <- data.frame(o=pdclu,l=pdclp,h=pdclh)
  rownames(res) <- pdcln
  return(list(res=res[names(lploc),],rma=as.data.frame(t(rma[,names(lploc)])),d=list(v=v,pdcn=pdcn)))
}

pdishuffle <- function(pdi,lp,cn,m,all=FALSE){
  if ( all ){
    d <- as.data.frame(pdi)
    y <- t(apply(pdi,1,function(x) runif(length(x),min=min(x),max=max(x))))
    names(y)    <- names(d)
    rownames(y) <- rownames(d)
    return(y)
  }else{
    fl <- TRUE
    for ( i in unique(lp)){
      if ( sum(lp == i) > 1 ){
        y <-data.frame( t(apply(as.data.frame(pdi[,lp == i]),1,sample)) )
      }else{
        y <- as.data.frame(pdi[,lp == i])
      }
      names(y) <- names(lp)[lp == i]
      rownames(y) <- names(lp)
      z <- if (fl) y else cbind(z,y)
      fl <- FALSE
    }
    z <- t(z[,names(lp)])
    return(z)
  }
}

setGeneric("projcells", function(object,cthr=0,nmode=FALSE) standardGeneric("projcells"))

setMethod("projcells",
          signature = "Ltree",
          definition = function(object,cthr,nmode){
            if ( ! is.numeric(cthr) ) stop( "cthr has to be a non-negative number" ) else if ( cthr < 0 ) stop( "cthr has to be a non-negative number" )
            if ( ! length(object@sc@cpart == 0) ) stop( "please run findoutliers on the SCseq input object before initializing Ltree" )
            if ( !is.numeric(nmode) & !is.logical(nmode) ) stop("argument nmode has to be logical (TRUE/FALSE)")
         
            object@par$cthr  <- cthr
            object@par$nmode <- nmode
            
            lp <- object@sc@cpart
            ld <- object@sc@distances
            n  <- aggregate(rep(1,length(lp)),list(lp),sum)
            n  <- as.vector(n[order(n[,1],decreasing=FALSE),-1])
            m  <- (1:length(n))[n>cthr]
            f  <- lp %in% m
            lp <- lp[f]
            ld <- ld[f,f]

            pdil <- object@sc@tsne[f,]
            cnl  <- aggregate(pdil,by=list(lp),median)
            cnl  <- cnl[order(cnl[,1],decreasing=FALSE),-1]

            pdi <- suppressWarnings( cmdscale(as.matrix(ld),k=ncol(ld)-1) )
            cn <- as.data.frame(pdi[compmedoids(object@sc@fdata[,names(lp)],lp),])
            rownames(cn) <- 1:nrow(cn)

            x <- compproj(pdi,lp,cn,m)
            res <- x$res
           
            if ( nmode ){
              rma <- x$rma
              z <- paste("X",t(as.vector(apply(cbind(lp,ld),1,function(x){ f <- lp != x[1]; lp[f][which(x[-1][f] == min(x[-1][f]))[1]] }))),sep="")
              k <- apply(cbind(z,rma),1,function(x) (x[-1])[names(rma) == x[1]])
              rn <- res
              rn$l <- as.numeric(sub("X","",z))
              rn$h <- as.numeric(k)
              res <- rn
              x$res <- res
            }

            object@ldata  <- list(lp=lp,ld=ld,m=m,pdi=pdi,pdil=pdil,cn=cn,cnl=cnl)
            object@trproj <- x
            return(object)
          }
          )





setGeneric("projback", function(object,pdishuf=2000,nmode=FALSE,fast=FALSE,rseed=17000) standardGeneric("projback"))

setMethod("projback",
          signature = "Ltree",
          definition = function(object,pdishuf,nmode,fast,rseed){
            if ( !is.numeric(nmode) & !is.logical(nmode) ) stop("argument nmode has to be logical (TRUE/FALSE)")
            if ( !is.numeric(fast) & !is.logical(fast) ) stop("argument fast has to be logical (TRUE/FALSE)")
            if ( ! is.numeric(pdishuf) ) stop( "pdishuf has to be a non-negative integer number" ) else if ( round(pdishuf) != pdishuf | pdishuf < 0 ) stop( "pdishuf has to be a non-negative integer number" )
            if ( length(object@trproj) == 0 ) stop("run projcells before projback")
            object@par$pdishuf  <- pdishuf
            object@par$rseed    <- rseed
            object@par$fast     <- fast
            if ( ! nmode & ! fast ){
              set.seed(rseed)
              for ( i in 1:pdishuf ){
                cat("pdishuffle:",i,"\n")
                x <- compproj(pdishuffle(object@ldata$pdi,object@ldata$lp,object@ldata$cn,object@ldata$m,all=TRUE),object@ldata$lp,object@ldata$cn,object@ldata$m,d=object@trproj$d)$res
                y <- if ( i == 1 ) t(x) else cbind(y,t(x))
              }    
              ##important
              object@prback <- as.data.frame(t(y))
              
              x <- object@prback
              x$n <- as.vector(t(matrix(rep(1:(nrow(x)/nrow(object@ldata$pdi)),nrow(object@ldata$pdi)),ncol=nrow(object@ldata$pdi))))
              object@prbacka <- aggregate(data.frame(count=rep(1,nrow(x))),by=list(n=x$n,o=x$o,l=x$l),sum)
            }
            return( object )
          }
          )





setGeneric("lineagetree", function(object,pthr=0.01,nmode=FALSE,fast=FALSE) standardGeneric("lineagetree"))

setMethod("lineagetree",
          signature = "Ltree",
          definition = function(object,pthr,nmode,fast){
            if ( !is.numeric(nmode) & !is.logical(nmode) ) stop("argument nmode has to be logical (TRUE/FALSE)")
            if ( !is.numeric(fast) & !is.logical(fast) ) stop("argument fast has to be logical (TRUE/FALSE)")
            if ( length(object@trproj) == 0 ) stop("run projcells before lineagetree")
            if ( max(dim(object@prback)) == 0 & ! nmode & ! fast  ) stop("run projback before lineagetree")
            if ( ! is.numeric(pthr) ) stop( "pthr has to be a non-negative number" ) else if ( pthr < 0 ) stop( "pthr has to be a non-negative number" )
            
            object@par$pthr <- pthr
            cnl    <- object@ldata$cnl
            pdil   <- object@ldata$pdil
            cn    <- object@ldata$cn
            pdi   <- object@ldata$pdi
            m      <- object@ldata$m
            lp     <- object@ldata$lp
            res    <- object@trproj$res
            rma    <- object@trproj$rma
            prback <- object@prback
            
            cm <- as.matrix(dist(cnl))*0
            linl <- list()
            linn <- list()
            for ( i in 1:length(m) ){
              for ( j in i:length(m) ){
                linl[[paste(m[i],m[j],sep=".")]] <- c()
                linn[[paste(m[i],m[j],sep=".")]] <- c()
              }
            }
            sigcell <- c()
            for ( i in 1:nrow(res) ){
              a <- which( m == res$o[i])
              if ( sum( lp == m[a] ) == 1 ){
                k <- t(cnl)[,a]
                k <- NA
                sigcell <- append(sigcell, FALSE)
              }else{
                b <- which(m == res$l[i] )
                h <- res$h[i]
                if ( nmode | fast ){
                  sigcell <- append(sigcell, FALSE)
                }else{
                  f <- prback$o == m[a] & prback$l == m[b]
                  if ( sum(f) > 0 ){
                    ql <- quantile(prback[f,"h"],pthr,na.rm=TRUE)
                    qh <- quantile(prback[f,"h"],1 - pthr,na.rm=TRUE)
                  }else{
                    ql <- 0
                    qh <- 0
                  }
                  sigcell <- if (is.na(h) ) append(sigcell, NA) else if ( h > qh |  h < min(0,ql) ) append(sigcell, TRUE) else append(sigcell, FALSE)
                }
                if ( !is.na(res$h[i]) ){
                  w <- t(pdil)[,i] - t(cnl)[,a]
                  v <- t(cnl)[,b] - t(cnl)[,a]
                  
                  wo <- t(pdi)[,i] - t(cn)[,a]
                  vo <-  t(cn)[,b] - t(cn)[,a]
                  df <-( h*sqrt(sum(wo*wo)) )/sqrt(sum(vo*vo))*v
                  k <- df + t(cnl)[,a]
                  cm[a,b] <- cm[a,b] + 1
                  so <- m[sort(c(a,b))]
                  dfl <-  sign(h)*sqrt( sum( df*df ) )/sqrt(sum(v*v))
                  if ( a > b ) dfl <-  1 - dfl
                  linn[[paste(so[1],so[2],sep=".")]] <- append( linn[[paste(so[1],so[2],sep=".")]], rownames(pdi)[i] )
                  linl[[paste(so[1],so[2],sep=".")]] <- append( linl[[paste(so[1],so[2],sep=".")]], dfl ) 
                }else{
                  k <- t(cnl)[,a]
                  for ( j in unique(lp[lp != m[a]]) ){
                    b <- which(j == m)
                    so <- m[sort(c(a,b))]
                    dfl <- 0
                    if ( a > b ) dfl <-  1 - dfl
                    linn[[paste(so[1],so[2],sep=".")]] <- append( linn[[paste(so[1],so[2],sep=".")]], rownames(pdi)[i] )
                    linl[[paste(so[1],so[2],sep=".")]] <- append( linl[[paste(so[1],so[2],sep=".")]], dfl ) 
                  }
                }
              }
              lt <- if ( i == 1 ) data.frame(k) else cbind(lt,k)
            }
            lt <- t(lt)
            cm <- as.data.frame(cm)
            names(cm) <- paste("cl",m,sep=".")
            rownames(cm) <- paste("cl",m,sep=".")
            lt <- as.data.frame(lt)
            rownames(lt) <- rownames(res)
            object@ltcoord <- as.matrix(lt)
            object@prtree  <- list(n=linn,l=linl)
            object@cdata$counts <- cm
            names(sigcell) <- rownames(res)
            object@sigcell <- sigcell

            return( object )
          }
          )





setGeneric("comppvalue", function(object,pethr=0.01,nmode=FALSE,fast=FALSE) standardGeneric("comppvalue"))

setMethod("comppvalue",
          signature = "Ltree",
          definition = function(object,pethr,nmode,fast){
            if ( !is.numeric(nmode) & !is.logical(nmode) ) stop("argument nmode has to be logical (TRUE/FALSE)")
            if ( !is.numeric(fast) & !is.logical(fast) ) stop("argument fast has to be logical (TRUE/FALSE)")
            if ( length(object@prtree) == 0 ) stop("run lineagetree before comppvalue")
            if ( ! is.numeric(pethr) ) stop( "pethr has to be a non-negative number" ) else if ( pethr < 0 ) stop( "pethr has to be a non-negative number" )
            object@par$pethr <- pethr
            cm <- object@cdata$counts
            cmpv   <- cm*NA
            cmpvd  <- cm*NA
            cmbr   <- cm*NA
            cmpvn  <- cm*NA
            cmpvnd <- cm*NA
            cmfr   <- cm/apply(cm,1,sum)
            if ( nmode ){
              N <- apply(cm,1,sum) + 1
              N0 <- sum(N) - N
              n0 <- t(matrix(rep(N,length(N)),ncol=length(N)))
              p <- n0/N0
              n <- cm
              k <- cbind(N,p,n)
              cmpv   <- apply(k,1,function(x){N <- x[1]; p <- x[2:( ncol(cm) + 1 )];  n <- x[( ncol(cm) + 2 ):( 2*ncol(cm) + 1)]; apply(cbind(n,p),1,function(x,N) binom.test(x[1],N,min(1,x[2]),alternative="g")$p.value,N=N)})
              cmpvd   <- apply(k,1,function(x){N <- x[1]; p <- x[2:( ncol(cm) + 1 )];  n <- x[( ncol(cm) + 2 ):( 2*ncol(cm) + 1)]; apply(cbind(n,p),1,function(x,N) binom.test(x[1],N,min(1,x[2]),alternative="l")$p.value,N=N)})
              cmpvn  <- cmpv
              cmpvnd <- cmpvd
              cmbr   <- as.data.frame(n0/N0*N)
              names(cmbr)    <- names(cm)
              rownames(cmbr) <- rownames(cm)
            }else if ( fast ){
              set.seed(object@par$rseed)
              p <- sort(unique(object@ldata$lp))
              for ( i in 1:length(p) ){
                mm <- 1
                dloc <- object@trproj$d
                dloc$v <- as.data.frame(matrix(rep(as.matrix(dloc$v),mm),ncol=ncol(dloc$v)))
                dloc$pdcn <- as.data.frame(matrix(rep(as.matrix(dloc$pdcn),mm),ncol=ncol(dloc$pdcn)))
                lploc <- rep(p[i] + object@ldata$lp*0,mm)
                names(lploc) <- 1:length(lploc)
                x <- compproj(pdishuffle(matrix(rep(t(object@ldata$pdi),mm),ncol=ncol(object@ldata$pdi)),object@ldata$lp,object@ldata$cn,object@ldata$m,all=TRUE),lploc,object@ldata$cn,object@ldata$m,d=dloc)$res

                #x <- compproj(pdishuffle(object@ldata$pdi,object@ldata$lp,object@ldata$cn,object@ldata$m,all=TRUE),p[i] + object@ldata$lp*0,object@ldata$cn,object@ldata$m,d=object@trproj$d)$res
                z <-merge(data.frame(p=p),aggregate(rep(1,nrow(x)),by=list(x$l),sum),by.x="p",by.y="Group.1",all.x=T)
                z$x[is.na(z$x)] <- 0
                pr <- z$x/sum(z$x)
                if ( i == 1 ) d <- data.frame(z) else d[,i] <- pr
              }
              N <- apply(cm,1,sum) + 1
              k <- cbind(N,d,cm)
              cmpv   <- apply(k,1,function(x){N <- x[1]; p <- x[2:( ncol(cm) + 1 )];  n <- x[( ncol(cm) + 2 ):( 2*ncol(cm) + 1)]; apply(cbind(n,p),1,function(x,N) binom.test(x[1],N,min(1,x[2]),alternative="g")$p.value,N=N)})
              cmpvd   <- apply(k,1,function(x){N <- x[1]; p <- x[2:( ncol(cm) + 1 )];  n <- x[( ncol(cm) + 2 ):( 2*ncol(cm) + 1)]; apply(cbind(n,p),1,function(x,N) binom.test(x[1],N,min(1,x[2]),alternative="l")$p.value,N=N)})
              ##cmpv   <- apply(k,1,function(x){N <- x[1]; p <- x[2:( ncol(cm) + 1 )];  n <- x[( ncol(cm) + 2 ):( 2*ncol(cm) + 1)]; apply(cbind(n,p),1,function(x,N) fisher.test(matrix(c(1e6,round(1e6*x[2],0),N,x[1]),ncol=2),alternative="greater")$p.value,N=N)})
              ##cmpvd   <- apply(k,1,function(x){N <- x[1]; p <- x[2:( ncol(cm) + 1 )];  n <- x[( ncol(cm) + 2 ):( 2*ncol(cm) + 1)]; apply(cbind(n,p),1,function(x,N) fisher.test(matrix(c(1e6,round(1e6*x[2],0),N,x[1]),ncol=2),alternative="less")$p.value,N=N)})
 
              cmpvn  <- cmpv
              cmpvnd <- cmpvd
              cmbr   <- as.data.frame(d*N)
              names(cmbr)    <- names(cm)
              rownames(cmbr) <- rownames(cm)
            }else{
              for ( i in 1:nrow(cm) ){
                for ( j in 1:ncol(cm) ){
                  f <- object@prbacka$o == object@ldata$m[i] & object@prbacka$l == object@ldata$m[j]
                  x <- object@prbacka$count[f]
                  if ( sum(f) < object@par$pdishuf ) x <- append(x,rep(0, object@par$pdishuf - sum(f)))
                  cmbr[i,j]   <- if ( sum(f) > 0 ) mean(x) else 0
                  cmpv[i,j]   <- if ( quantile(x,1 - pethr) < cm[i,j] ) 0 else 0.5
                  cmpvd[i,j]  <- if ( quantile(x,pethr) > cm[i,j] ) 0 else 0.5
                  cmpvn[i,j]  <- sum( x >= cm[i,j])/length(x)
                  cmpvnd[i,j] <- sum( x <= cm[i,j])/length(x)
                }
              }
            }

            diag(cmpv)   <- .5
            diag(cmpvd)  <- .5
            diag(cmpvn)  <- NA
            diag(cmpvnd) <- NA

            object@cdata$counts.br <- cmbr
            object@cdata$pv.e <- cmpv
            object@cdata$pv.d <- cmpvd
            object@cdata$pvn.e <- cmpvn
            object@cdata$pvn.d <- cmpvnd

            m    <- object@ldata$m
            linl <- object@prtree$l
            ls   <- as.data.frame(matrix(rep(NA,length(m)**2),ncol=length(m)))
            names(ls) <- rownames(ls) <- paste("cl",m,sep=".")
            for ( i in 1:( length(m) - 1 )){
              for ( j in (i + 1):length(m) ){
                na <- paste(m[i],m[j],sep=".")
                if ( na %in% names(linl) &  min(cmpv[i,j],cmpv[j,i],na.rm=TRUE) < pethr ){
                  y <- sort(linl[[na]])
                  nn <- ( 1 - max(y[-1] - y[-length(y)]) )
                }else{
                  nn <- 0
                }
                ls[i,j] <- nn
              }
            }
            object@cdata$linkscore <- ls

            return(object)
          }
          )

setGeneric("plotlinkpv", function(object) standardGeneric("plotlinkpv"))

setMethod("plotlinkpv",
          signature = "Ltree",
          definition = function(object){
            if ( length(object@cdata) <= 0 ) stop("run comppvalue before plotlinkpv")
            pheatmap(-log2(object@cdata$pvn.e + 1/object@par$pdishuf/2))
          }
          )

setGeneric("plotlinkscore", function(object) standardGeneric("plotlinkscore"))

setMethod("plotlinkscore",
          signature = "Ltree",
          definition = function(object){
            if ( length(object@cdata) <= 0 ) stop("run comppvalue before plotlinkscore")
            pheatmap(object@cdata$linkscore,cluster_rows=FALSE,cluster_cols=FALSE)
          }
          )

setGeneric("plotmapprojections", function(object) standardGeneric("plotmapprojections"))

setMethod("plotmapprojections",
          signature = "Ltree",
          definition = function(object){
            if ( length(object@cdata) <= 0 ) stop("run comppvalue before plotmapprojections")
         
            cent <- object@sc@fdata[object@sc@cluster$features,compmedoids(object@sc@fdata[object@sc@cluster$features,],object@sc@cpart,metric=object@sc@clusterpar$metric)]
            ##dc <- as.data.frame(1 - cor(cent))
            dc <- as.data.frame(as.matrix(dist.gen(t(cent),method=object@sc@clusterpar$metric)))
            names(dc) <- sort(unique(object@sc@cpart))
            rownames(dc) <- sort(unique(object@sc@cpart))
            trl <- spantree(dc[object@ldata$m,object@ldata$m])

            u <- object@ltcoord[,1]
            v <- object@ltcoord[,2]
            cnl <- object@ldata$cnl
            plot(u,v,cex=1.5,col="grey",pch=20,xlab="Dim 1",ylab="Dim 2")
            for ( i in unique(object@ldata$lp) ){ f <- object@ldata$lp == i; text(u[f],v[f],i,cex=.75,font=4,col=object@sc@fcol[i]) }
            points(cnl[,1],cnl[,2])
            text(cnl[,1],cnl[,2],object@ldata$m,cex=2)
            for ( i in 1:length(trl$kid) ){
              lines(c(cnl[i+1,1],cnl[trl$kid[i],1]),c(cnl[i+1,2],cnl[trl$kid[i],2]),col="black")
            }
          }
          )


setGeneric("plotmap", function(object) standardGeneric("plotmap"))

setMethod("plotmap",
          signature = "Ltree",
          definition = function(object){
            if ( length(object@cdata) <= 0 ) stop("run comppvalue before plotmap")
         
            cent <- object@sc@fdata[,compmedoids(object@sc@fdata,object@sc@cpart)]
            dc <- as.data.frame(1 - cor(cent))
            names(dc) <- sort(unique(object@sc@cpart))
            rownames(dc) <- sort(unique(object@sc@cpart))
            trl <- spantree(dc[object@ldata$m,object@ldata$m])

  
            u <- object@ldata$pdil[,1]
            v <- object@ldata$pdil[,2]
            cnl <- object@ldata$cnl
            plot(u,v,cex=1.5,col="grey",pch=20,xlab="Dim 1",ylab="Dim 2")
            for ( i in unique(object@ldata$lp) ){ f <- object@ldata$lp == i; text(u[f],v[f],i,cex=.75,font=4,col=object@sc@fcol[i]) }
            points(cnl[,1],cnl[,2])
            text(cnl[,1],cnl[,2],object@ldata$m,cex=2)
            for ( i in 1:length(trl$kid) ){
              lines(c(cnl[i+1,1],cnl[trl$kid[i],1]),c(cnl[i+1,2],cnl[trl$kid[i],2]),col="black")
            }
          }
          )


setGeneric("plottree", function(object,showCells=TRUE,nmode=FALSE,scthr=0) standardGeneric("plottree"))

setMethod("plottree",
          signature = "Ltree",
          definition = function(object,showCells,nmode,scthr){
            if ( length(object@cdata) <= 0 ) stop("run comppvalue before plotmap")
            if ( !is.numeric(nmode) & !is.logical(nmode) ) stop("argument nmode has to be logical (TRUE/FALSE)")
            if ( !is.numeric(showCells) & !is.logical(showCells) ) stop("argument showCells has to be logical (TRUE/FALSE)")
            if ( ! is.numeric(scthr) ) stop( "scthr has to be a non-negative number" ) else if ( scthr < 0 | scthr > 1 ) stop( "scthr has to be a number between 0 and 1" )
           
         
            ramp <- colorRamp(c("darkgreen", "yellow", "brown"))
            mcol <- rgb( ramp(seq(0, 1, length = 101)), max = 255)
            co <- object@cdata
            fc <- (co$counts/( co$counts.br + .5))*(co$pv.e < object@par$pethr)
            fc <- fc*(fc > t(fc)) + t(fc)*(t(fc) >= fc)
            fc <- log2(fc + (fc == 0))

            k <- -log10(sort(unique(as.vector(t(co$pvn.e))[as.vector(t(co$pv.e))<object@par$pethr])) + 1/object@par$pdishuf)
            if (length(k) == 1) k <- c(k - k/100,k)
            mlpv <- -log10(co$pvn.e + 1/object@par$pdishuf)
            diag(mlpv) <- min(mlpv,na.rm=TRUE)
            dcc <- t(apply(round(100*(mlpv - min(k))/(max(k) - min(k)),0) + 1,1,function(x){y <- c(); for ( n in x ) y <- append(y,if ( n < 1 ) NA else mcol[n]); y }))


            cx <- c()
            cy <- c()
            va <- c()
            m <- object@ldata$m
            for ( i in 1:( length(m) - 1 ) ){
              for ( j in ( i + 1 ):length(m) ){
                if ( min(co$pv.e[i,j],co$pv.e[j,i],na.rm=TRUE) < object@par$pethr ){
                  if ( mlpv[i,j] > mlpv[j,i] ){
                    va <- append(va,dcc[i,j])
                  }else{
                    va <- append(va,dcc[j,i])
                  }
                  cx <- append(cx,i)
                  cy <- append(cy,j)
                }
              }
            }



            cnl <- object@ldata$cnl
            u <- object@ltcoord[,1]
            v <- object@ltcoord[,2]
            layout( cbind(c(1, 1), c(2, 3)),widths=c(5,1,1),height=c(5,5,1))
            par(mar = c(12,5,1,1))

            if ( showCells ){
              plot(u,v,cex=1.5,col="grey",pch=20,xlab="Dim 1",ylab="Dim 2")
              if ( !nmode ) points(u[object@sigcell],v[object@sigcell],col="black")
            }else{
              plot(u,v,cex=0,xlab="Dim 1",ylab="Dim 2")
            }
    
            if ( length(va) > 0 ){
              f <- order(va,decreasing=TRUE)
              for ( i in 1:length(va) ){
                if ( object@cdata$linkscore[cx[i],cy[i]] > scthr ){
                  if ( showCells ){
                    lines(cnl[c(cx[i],cy[i]),1],cnl[c(cx[i],cy[i]),2],col=va[i],lwd=2)
                  }else{
                    ##nn <- min(10,fc[cx[i],cy[i]])
                    lines(cnl[c(cx[i],cy[i]),1],cnl[c(cx[i],cy[i]),2],col=va[i],lwd=5*object@cdata$linkscore[cx[i],cy[i]])
                  }
                }
              }
            }



            en <- aggregate(object@entropy,list(object@sc@cpart),median)
            en <- en[en$Group.1 %in% m,]
    
            mi <- min(en[,2],na.rm=TRUE)
            ma <- max(en[,2],na.rm=TRUE)
            w <- round((en[,2] - mi)/(ma - mi)*99 + 1,0)
            ramp <- colorRamp(c("red","orange", "pink","purple", "blue"))
            ColorRamp <- rgb( ramp(seq(0, 1, length = 101)), max = 255)
            ColorLevels <- seq(mi, ma, length=length(ColorRamp))
            if ( mi == ma ){
              ColorLevels <- seq(0.99*mi, 1.01*ma, length=length(ColorRamp))
            }
            for ( i in m ){
              f <- en[,1] == m
              points(cnl[f,1],cnl[f,2],cex=5,col=ColorRamp[w[f]],pch=20)
            }
            text(cnl[,1],cnl[,2],m,cex=1.25,font=4,col="white")
            par(mar = c(5, 4, 1, 2))
            image(1, ColorLevels,
                  matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
                  col=ColorRamp,
                  xlab="",ylab="",
                  xaxt="n")
            coll <- seq(min(k), max(k), length=length(mcol))
            image(1, coll,
                  matrix(data=coll, ncol=length(mcol),nrow=1),
                  col=mcol,
                  xlab="",ylab="",
                  xaxt="n")
            layout(1)
          }
          )


setGeneric("plotdistanceratio", function(object) standardGeneric("plotdistanceratio"))

setMethod("plotdistanceratio",
          signature = "Ltree",
          definition = function(object){
            if ( length(object@ldata) <= 0 ) stop("run projcells before plotdistanceratio")
            l <- as.matrix(dist(object@ldata$pdi))
            z <- (l/object@ldata$ld)
            hist(log2(z),breaks=100,xlab=" log2 emb. distance/distance",main="")
          }
          )


setGeneric("getproj", function(object,i) standardGeneric("getproj"))

setMethod("getproj",
          signature = "Ltree",
          definition = function(object,i){
            if ( length(object@ldata) <= 0 ) stop("run projcells before plotdistanceratio")
            if ( ! i %in% object@ldata$m )  stop(paste("argument i has to be one of",paste(object@ldata$m,collapse=",")))
            x <- object@trproj$rma[names(object@ldata$lp)[object@ldata$lp == i],]
            x <- x[,names(x) != paste("X",i,sep="")]
            f <- !is.na(x[,1])
            x <- x[f,]
            if ( nrow(x) > 1 ){
              y <- x
              y <- as.data.frame(t(apply(y,1,function(x) (x - mean(x))/sqrt(var(x)))))
            }
            names(x) = sub("X","cl.",names(x))
            names(y) = sub("X","cl.",names(y))
            return(list(pr=x,prz=y))
          }
          )


setGeneric("projenrichment", function(object) standardGeneric("projenrichment"))

setMethod("projenrichment",
          signature = "Ltree",
          definition = function(object){
            if ( length(object@cdata) <= 0 ) stop("run comppvalue before plotmap")
            
            ze <- ( object@cdata$pv.e < object@par$pethr | object@cdata$pv.d < object@par$pethr) * (object@cdata$counts + .1)/( object@cdata$counts.br + .1 )
            pheatmap(log2(ze + ( ze == 0 ) ),cluster_rows=FALSE,cluster_cols=FALSE)
          }
          )





setGeneric("compscore", function(object,nn=1,scthr=0) standardGeneric("compscore"))

setMethod("compscore",
          signature = "Ltree",
          definition = function(object,nn,scthr){
            if ( length(object@cdata) <= 0 ) stop("run comppvalue before compscore")
            if ( ! is.numeric(nn) ) stop( "nn has to be a non-negative integer number" ) else if ( round(nn) != nn | nn < 0 ) stop( "nn has to be a non-negative integer number" )
            if ( ! is.numeric(scthr) ) stop( "scthr has to be a non-negative number" ) else if ( scthr < 0 | scthr > 1 ) stop( "scthr has to be a number between 0 and 1" )
            scm <- as.matrix(object@cdata$linkscore)
            diag(scm) <- rep(0,ncol(scm))
            for ( i in 1:ncol(scm)) scm[i:nrow(scm),i] <- scm[i,i:ncol(scm)]

            x <- object@cdata$counts*(object@cdata$pv.e < object@par$pethr)*(scm > scthr)>0
            y <- x | t(x)
            
            if ( max(y) > 0 ){
              z <- apply(y,1,sum)
              nl <- list()
              n <- list()
              for ( i in 1:nn ){
                if ( i == 1 ){
                  n[[i]] <- as.list(apply(y,1,function(x) grep(TRUE,x)))
                  nl <- data.frame( apply(y,1,sum) )
                }
                if ( i > 1 ){
                  v <- rep(0,nrow(nl))
                  n[[i]] <- list()
                  for ( j in 1:length(n[[i-1]]) ){
                    cl <- n[[i-1]][[j]]
                    if ( length(cl) == 0 ){
                      n[[i]][[paste("d",j,sep="")]] <- NA
                      v[j] <- 0
                    }else{
                      k  <- if ( length(cl) > 1 ) apply(y[cl,],2,sum) > 0 else if ( length(cl) == 1 ) y[cl,]
                      n[[i]][[paste("d",j,sep="")]] <- sort(unique(c(cl,grep(TRUE,k))))
                      v[j] <- length(n[[i]][[paste("d",j,sep="")]])
                    }
                  }
                  names(n[[i]]) <- names(z)
                  nl <- cbind(nl,v)
          
                }
              }
              x <- nl[,i]
              names(x) <- rownames(nl)
            }else{
              x <- rep(0,length(object@ldata$m))
              names(x) <- paste("cl",object@ldata$m,sep=".")
            }
            
            v <- aggregate(object@entropy,list(object@sc@cpart),median)
            v <- v[v$Group.1 %in% object@ldata$m,]
            w <- as.vector(v[,-1])
            names(w) <- paste("cl.",v$Group.1,sep="")
            w <- w - min(w)
            
            return(list(links=x,entropy=w,StemIDscore=x*w))
          }
          )




setGeneric("plotscore", function(object,nn=1,scthr=0) standardGeneric("plotscore"))

setMethod("plotscore",
          signature = "Ltree",
          definition = function(object,nn,scthr){
            if ( length(object@cdata) <= 0 ) stop("run comppvalue before plotscore")
            if ( ! is.numeric(scthr) ) stop( "scthr has to be a non-negative number" ) else if ( scthr < 0 | scthr > 1 ) stop( "scthr has to be a number between 0 and 1" )
            x <- compscore(object,nn,scthr)
            layout(1:3)
            barplot(x$links,names.arg=sub("cl\\.","",object@ldata$m),xlab="Cluster",ylab="Number of links",cex.names=1)
            barplot(x$entropy,names.arg=sub("cl\\.","",object@ldata$m),xlab="Cluster",ylab="Delta-Entropy",cex.names=1)
            barplot(x$StemIDscore,names.arg=sub("cl\\.","",object@ldata$m),xlab="Cluster",ylab="Number of links * Delta-Entropy",cex.names=1)
            layout(1)
          }
          )


setGeneric("branchcells", function(object,br) standardGeneric("branchcells"))

setMethod("branchcells",
          signature = "Ltree",
          definition = function(object,br){
            if ( length(object@ldata) <= 0 ) stop("run projcells before branchcells")
            msg <- paste("br needs to be list of length two containing two branches, where each has to be one of", paste(names(object@prtree$n),collapse=","))
            if ( !is.list(br) ) stop(msg) else if ( length(br) != 2 ) stop(msg) else if ( ! br[[1]] %in% names(object@prtree$n) | ! br[[2]] %in% names(object@prtree$n) ) stop(msg)

             
            n <- list()
            scl <- object@sc
            k <- c()
            cl <- intersect( as.numeric(strsplit(br[[1]],"\\.")[[1]]), as.numeric(strsplit(br[[2]],"\\.")[[1]]))
            if ( length(cl) == 0 ) stop("the two branches in br need to have one cluster in common.")
                      
            for ( i in 1:length(br) ){
              f <- object@sc@cpart[ object@prtree$n[[br[[i]]]] ] %in% cl
              if ( sum(f) > 0 ){
                n[[i]] <- names(object@sc@cpart[ object@prtree$n[[br[[i]]]] ])[f]
                k <- append(k, max( scl@cpart ) + 1)
                scl@cpart[n[[i]]] <- max( scl@cpart ) + 1
              }else{
                stop(paste("no cells on branch",br[[i]],"fall into clusters",cl))
              }
            }
            set.seed(111111)
            scl@fcol <- sample(rainbow(max(scl@cpart)))
            z <- diffgenes(scl,k[1],k[2])
            return( list(n=n,scl=scl,k=k,diffgenes=z) )
          }
          )



pamk <- function (data, krange = 2:10, criterion = "asw", usepam = TRUE, 
    scaling = FALSE, alpha = 0.001, diss = inherits(data, "dist"), 
    critout = FALSE, ns = 10, seed = NULL, ...) 
{
    ddata <- as.matrix(data)
    if (!identical(scaling, FALSE)) 
        sdata <- scale(ddata, scale = scaling)
    else sdata <- ddata
    cluster1 <- 1 %in% krange
    critval <- numeric(max(krange))
    pams <- list()
    for (k in krange) {
        if (usepam) 
            pams[[k]] <- pam(as.dist(sdata), k, diss = TRUE)
        else pams[[k]] <- clara(as.dist(sdata), k, diss = TRUE)
        if (k != 1) 
            critval[k] <- switch(criterion, asw = pams[[k]]$silinfo$avg.width, 
                multiasw = distcritmulti(sdata, pams[[k]]$clustering, 
                  seed = seed, ns = ns)$crit.overall, ch = ifelse(diss, 
                  cluster.stats(sdata, pams[[k]]$clustering)$ch, 
                  calinhara(sdata, pams[[k]]$clustering)))
        if (critout) 
            cat(k, " clusters ", critval[k], "\n")
    }
    k.best <- if ( length(krange) == 1 ) krange else (1:max(krange))[which.max(critval)]
    if (cluster1) {
        if (diss) 
            cluster1 <- FALSE
        else {
            cxx <- dudahart2(sdata, pams[[2]]$clustering, alpha = alpha)
            critval[1] <- cxx$p.value
            cluster1 <- cxx$cluster1
        }
    }
    if (cluster1) 
        k.best <- 1
    out <- list(pamobject = pams[[k.best]], nc = k.best, crit = critval)
    out
}

pamkdCBI <- function (data, krange = 2:10, k = NULL, criterion = "asw", usepam = TRUE, 
    scaling = TRUE, diss = inherits(data, "dist"), ...) 
{
    if (!is.null(k)) 
        krange <- k
    c1 <- pamk(as.dist(data), krange = krange, criterion = criterion, 
        usepam = usepam, scaling = scaling, diss = diss, ...)
    partition <- c1$pamobject$clustering
    cl <- list()
    nc <- c1$nc

    for (i in 1:nc) cl[[i]] <- partition == i
    out <- list(result = c1, nc = nc, clusterlist = cl, partition = partition, 
        clustermethod = "pam/estimated k", criterion = criterion)
    out
}
