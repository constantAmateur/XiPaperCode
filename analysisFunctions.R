#' Functions to help with downstream analysis


#############
# Libraries #
#############

library(Seurat)
library(miloR)
library(SingleCellExperiment)
library(Matrix)
library(gsl)
library(distr)
library(bettermc)
import('Code/scQC.R')

#############
# Functions #
#############


#' Plot fractions
#'
#' Plots X-Inactivation fractions per cluster relative to the model fit value.  Essentially a fix for the stupidity of FeaturePlot.
#'
#' @param srat The seurat object with StateX and clusters
#' @param tau The fit vaule of tau.  If NULL, raw fractions plotted.
#' @param useConfInt If TRUE, any value with confidence interval overlapping tau set to 0. 
plotFractions = function(srat,tau=NULL,useConfInt=TRUE){
  dd = cbind(srat@meta.data,srat@reductions$umap@cell.embeddings)
  #Get the per-cluster fraction
  tmp = table(dd$seurat_clusters,dd$stateX)
  class(tmp) ='matrix'
  tmp = data.frame(tmp)
  tmp$tot = rowSums(tmp)
  tmp$frac = tmp$Maternal/tmp$tot
  tmp$fracLow = sapply(seq(nrow(tmp)),function(e) binom.test(tmp$Maternal[e],tmp$tot[e],0.5)$conf.int[1])
  tmp$fracHigh = sapply(seq(nrow(tmp)),function(e) binom.test(tmp$Maternal[e],tmp$tot[e],0.5)$conf.int[2])
  tmp$cluster = rownames(tmp)
  if(is.null(tau) || !useConfInt){
    tmp$fracDiff = tmp$frac - ifelse(is.null(tau),0,tau)
  }else{
    tmp$fracDiff = ifelse(tmp$frac>tau,
                          ifelse(tmp$fracLow>tau,
                                 tmp$fracLow-tau,
                                 0),
                          ifelse(tmp$fracHigh<tau,
                                 tmp$fracHigh-tau,
                                 0)
                          )
  }
  #Shove back into dd
  dd$fracDiff = tmp$fracDiff[match(dd$seurat_clusters,tmp$cluster)]
  srat$fracDiff = tmp$fracDiff[match(srat@meta.data$seurat_clusters,tmp$cluster)]
  #Now plot it
  if(is.null(tau)){
    gg = ggplot(dd,aes(UMAP_1,UMAP_2))+geom_point(aes(colour=fracDiff),size=0.0)+scale_color_gradient(low='red',high='blue')
  }else{
    gg = ggplot(dd,aes(UMAP_1,UMAP_2))+geom_point(aes(colour=fracDiff),size=0.0)+scale_color_gradient2(midpoint=0,low='red',mid='grey',high='blue')
  }
  plot(gg)
  srat@misc$plotFractions = gg
  srat@misc$clustStats = tmp[order(tmp$fracDiff),]
  return(srat)
}


#' Run Milo
#'
#' Pretend that the stateX labels are the two conditions and look for differences in neighbourhoods.  This doesn't work by default, mainly as we have only one sample and so a negative-binomial model doesn't much work.
#'
#' @param srat The annotated Seurat object.
#' @param tau Maternal fraction to test against.
#' @param resultsPassthrough Pass through these metadata columns to the results objects from Seurat meta.data.
#' @param prop What fraction of cells to randomly sample?
#' @param k How many neighbours to use in graph.
#' @param d How many PCs to use.
#' @param logFC If FALSE, uses linear difference instead of logFC (label remains logFC to not break milo).
#' @param useGraphFDR Use the graph FDR.
runMilo = function(srat,tau=0.5,resultsPassthrough='seurat_clusters',prop=0.2,k=30,d=50,logFC=FALSE,useGraphFDR=TRUE){
  milo = Milo(as.SingleCellExperiment(srat[,!is.na(srat@meta.data$stateX)]))
  milo = buildGraph(milo,k=k,d=d)
  milo = makeNhoods(milo,prop=prop,k=k,d=d,refined=TRUE)
  tmp = quantile(colSums(nhoods(milo)))
  message(sprintf('Neighbourhoods have median size %d (range %d to %d)',round(tmp[3]),tmp[1],tmp[5]))
  if(tmp[3]<20)
    stop("Neighbourhoods to small, increase prop or k or both")
  #Now generate counts in neighbourhoods
  milo = countCells(milo, meta.data = data.frame(milo@colData), sample="stateX")
  milo = calcNhoodDistance(milo, d=d)
  #Do our own test because of reasons
  dd = data.frame(milo@nhoodCounts)
  dd$tot = rowSums(dd)
  dd$frac = dd$Maternal/dd$tot
  dd$pVal = sapply(seq(nrow(dd)),function(e) binom.test(dd$Maternal[e],dd$tot[e],tau)$p.value)
  #This needs to mimic the results format exactly
  res = data.frame(logFC = if(logFC){log(dd$frac)-log(tau)}else{dd$frac-tau},
                   logCPM = 14, #Dummy values
                   `F` = 0.2, #Dummy values
                   PValue = dd$pVal,
                   FDR = p.adjust(dd$pVal,method='BH'),
                   Nhood = as.numeric(rownames(dd)))
  if(!useGraphFDR){
    res$SpatialFDR = res$FDR
  }else{
    tmp = graphSpatialFDR(x.nhoods = nhoods(milo),
                          graph = graph(milo),
                          k=milo@.k,
                          pvalues = res$PValue[order(res$Nhood)],
                          indices = nhoodIndex(milo),
                          distances = nhoodDistances(milo),
                          reduced.dimensions = reducedDim(milo,'PCA'))
    res$SpatialFDR[order(res$Nhood)] = tmp
  }
  #Make thing for visualising
  milo = buildNhoodGraph(milo)
  #Add in metadata
  for(tgt in resultsPassthrough)
    res = annotateNhoods(milo, res, coldata_col = tgt)
  #Now return the two objects for you to with what you will
  return(list(milo=milo,res=res))
}

#' Standardise X data
#'
#' Given Seurat object with annotation and X-Inactivation predictions, flatten it into a data.frame with various standard properties to be used downstream.
#'
#' @param srat Seurat object.  Must have 'states', 'donor', and an annotation meta.data.
#' @param annotCol Name of annotation column.
standardiseX = function(srat,annotCol='annot'){
  #Get a distribution across age that accounts for sample variability
  dd = split(srat@meta.data,srat@meta.data[,annotCol])
  dd = lapply(dd,function(e) aggregate(cbind(matCount=states,patCount=1-states)~donor,data=e,FUN=sum,na.rm=TRUE))
  dd = cbind(do.call(rbind,dd),cellType = rep(names(dd),sapply(dd,nrow)))
  rownames(dd) = NULL
  dd$tot = dd$matCount+dd$patCount
  dd$postA = 1 + dd$matCount
  dd$postB = 1 + dd$patCount
  dd$matFrac = dd$matCount/dd$tot
  #Get median fracs by donor
  baseRates = split(dd,dd$donor)
  baseRates = lapply(baseRates,function(e) unlist(lapply(seq(nSamps),function(ee) median(rbeta(nrow(e),e$postA,e$postB)))))
  dd$matFracRel = dd$matFrac-sapply(baseRates,median)[dd$donor]
  return(dd)
}

#' Sample posterior fractions
#'
#' Sample from posterior for all observations and constuct a nrow(dd) x nSamps matrix giving various derived quantities: the sampled fractions, the baseRate for each, and the relative fraction.
#'
#' The base rate calculation mode is set using \code{baseToUse}.  If set to median, or mean, the base rate is taken to be the median or mean of the cell types, under the assumption that most populations won't have a clonal expansion.  If 'fixed', then a constant value (such as 0.5) is assumed.  By default all cell types are used for this calculation, but this can be restricted using the \code{basePops} parameter.  If 'tau', use the global value.
#'
#' @param dd Clonality data, in flat format, as produced by \code{\link{standardiseX}}.
#' @param baseToUse How to calculate the base rate to assume for populations of cells with no clonal expansion.  Options are 'median','mean', and 'fixed'.  See details.
#' @param basePops This specifiecs the populations (cell types) that should be used when calculating the base rate.  If NULL, all are used.
#' @param meanTrim If \code{baseToUse} is 'mean', this specifies how much the data will be trimmed before the mean calculated.  Ignored otherwise.
#' @param fixedVal If \code{baseToUse} is 'fixed', this specifiecs the fixed base rate to assume.  Ignored otherwise.
#' @param nSamps How many samples to draw.
#' @return A list containing the three sampled matricies.
sampPostFracs = function(dd,baseToUse=c('median','mean','fixed','tau'),basePops=NULL,meanTrim=0,fixedVal=0.5,nSamps=1000){
  baseToUse = match.arg(baseToUse)
  cc = unique(dd$cellType)
  nn = unique(dd$donor)
  if(is.null(basePops))
    basePops = cc
  #Simulate across everything
  tmp = matrix(rbeta(nSamps*nrow(dd),dd$postA,dd$postB),ncol=nSamps)
  #For each simulation, work out the base rate 
  w = which(dd$cellType %in% basePops)
  if(baseToUse=='fixed'){
    baseRate=matrix(fixedVal,nrow=length(nn),ncol=nSamps)
    rownames(baseRate) = nn
  }else if(baseToUse=='median'){
    baseRate = do.call(rbind,lapply(split(seq(nrow(tmp)[w]),dd$donor[w]),function(e) apply(tmp[e,,drop=FALSE],2,median,na.rm=TRUE)))
  }else if(baseToUse=='mean'){
    baseRate = do.call(rbind,lapply(split(seq(nrow(tmp)[w]),dd$donor[w]),function(e) apply(tmp[e,,drop=FALSE],2,mean,trim=meanTrim,na.rm=TRUE)))
  }else if(baseToUse=='tau'){
    baseRate = do.call(rbind,lapply(split(seq(nrow(tmp)[w]),dd$donor[w]),function(e) colSums(dd[e,'tot']*tmp[e,,drop=FALSE])/sum(dd[e,'tot'])))
  }
  #Expand to match
  baseRate = baseRate[match(dd$donor,rownames(baseRate)),]
  return(list(sampFrac = tmp,baseRate=baseRate,relFrac = tmp-baseRate))
}

#' Calculate clonal fraction
#'
#' Calculates clonal fraction from X-Inactivation data.
#'
#' The formula for the clonal population is ifelse(x>r,(x-r)/(1-r),(r-x)/r), where x is the observed maternal fraction and r is the assumed base rate that would be observed if there were no clonal expansion.  The function \code{\link{sampPostFracs}} specifies how r is calculated.
#'
#' @param dd Clonality data, in flat format, as produced by \code{\link{standardiseX}}.
#' @param alpha Confidence parameter.  The default 0.05 maps to 95% confidence interval.
#' @param ar Aspect ratio for plot.
#' @param doPlot Should we make a plot?
#' @param returnFull Return all the raw sampled data.
#' @param ... Passed to \code{\link{sampPostFracs}}
calcClonal = function(dd,alpha=0.05,ar=4/3,doPlot=TRUE,returnFull=FALSE,...){
  samps = sampPostFracs(dd,...)
  sampFrac = samps$sampFrac
  baseRate = samps$baseRate
  cc = sort(unique(dd$cellType))
  nn = sort(unique(dd$donor))
  #Convert to clonality estimates
  cloneFrac = ifelse(sampFrac>=baseRate,(sampFrac-baseRate)/(1-baseRate),(baseRate-sampFrac)/baseRate)
  #Now plot one at a time
  if(doPlot){
    #Medians to sort by
    oo = sort(sapply(split(apply(cloneFrac,1,median),dd$cellType),median))
    cc = names(oo)
    x = gridDimensions(length(cc),ar)
    y = x[2]
    x = x[1]
    labIdxs = gridLayout(length(cc),x,y)
    for(i in seq(y)){
      for(j in seq(x)){
        tmp = setGridMargins(x,y,i,j,labIdxs)
        xaxt = tmp$xaxt
        yaxt = tmp$yaxt
        ii = tmp$ii
        if(ii>length(cc))
          next
        cType = cc[ii]
        w = which(dd$cellType==cType)
        tt = t(cloneFrac[w,])
        ttt = matrix(NA,ncol=length(nn),nrow=nSamps)
        colnames(ttt) = nn
        ttt[,match(colnames(tt),nn)] = tt
        boxplot(ttt,main=cType,xaxt=xaxt,yaxt=yaxt,ylim=c(0,1),las=2)
        par(xpd=FALSE)
        ttt = apply(ttt,1,median,na.rm=TRUE)
        if(!all(is.na(ttt)))
          abline(h=quantile(ttt,c(alpha/2,0.5,1-(alpha/2))),lty=c(2,1,2),col='red')
      }
    }
  }
  dd$clonalFrac = apply(cloneFrac,1,median)
  dd$clonalFracLow = apply(cloneFrac,1,quantile,alpha/2)
  dd$clonalFracHigh = apply(cloneFrac,1,quantile,1-(alpha/2))
  if(returnFull)
    return(list(dd=dd,cloneFrac=cloneFrac))
  return(dd)
}

#' Sample from exact sampling distribution for pearson correlation
#'
#' Given a pearson correlation coefficient \code{r} estimated from \code{n} observations, randomly sample \code{nSamps} from the sampling distribution.  
#'
#' The standard \code{\link{cor.test}} simply uses Fisher's transformation to approximate the sampling distribution.  But this breaks down for small n or strong correlation, exactly where the uncertainty matters.  This function implements the exact sampling distribution pdf.  As such, it will obviously be slower, but more accurate.
#'
#' @param r Observed pearson correlation.
#' @param n Number of paired observations from which r was derived.
#' @param nSamps Number of random samples to draw.
#' @return \code{nSamps} from the sampling distribution of the correlation coefficient.
sampleExactConfRho = function(r,n,nSamps=1000){
  #Given r - measured correlation, and n - number of observations
  vu = n-1
  #See here for definition of sampling distribution https://en.wikipedia.org/wiki/Pearson_correlation_coefficient#Using_the_exact_distribution
  d = function(x){
    logp = log(vu) + log(vu-1) + lgamma(vu-1)-lgamma(vu+.5)-.5*log(2*pi)+(vu-1)/2*log(1-r*r) +(vu-2)/2*log(1-x*x) + (1-2*vu)/2*log(1-x*r) + log(hyperg_2F1(3/2,-1/2,vu+0.5,(1+r*x)/2))
    return(exp(logp))
  }
  dist =AbscontDistribution(d=d,low1=-1,up1=1,low=-1,up=1)
  return(r(dist)(nSamps))
}

#' Calculate population correlations
#'
#' Calculates the posterior distribution for the correlation of all populations.  By default this is done using fractions relative to a sample specific base rate, but this can be changed by altering the parameters passed to \code{\link{sampPostFracs}}.
#'
#' If \code{sampleRho} is set (the default), then the output will take into account the uncertainty due to the number of individuals used to determine the correlation.  That is, the function will draw random samples from the sampling distribution of the correlation coefficient.  This uses the exact formulation of the sampling distribution (see \code{\link{sampleExactConfRho}}) and as such is much slower.  But should not be disabled unless the number of individuals is very large (>100).
#'
#' The function \code{sumFun} determines how the random samples (over the uncertainty in the correlation coefficient) are summarised into a value used to colour each cell.  This is also used to determine the ordering of rows/columns.
#'
#' @param dd Clonality data, in flat format, as produced by \code{\link{standardiseX}}.
#' @param alpha Confidence parameter.  The default 0.05 maps to 95% confidence interval.
#' @param ar Aspect ratio for plot.
#' @param doPlot Should we make a plot?
#' @param colFun Function to use to generate colours for histogram.  Should take a vector of correlations and produce a vector of colours.
#' @param sumFun Summary function that operates on the random samples of the correlation (over the uncertainty in the correlation) and produces a value used to colour things.
#' @param sampleRho Should we include variability due to sampling of the correlation coefficient?
#' @param verbose Be verbose?
#' @param mcParams Overwrite default parameters controlling the call to \code{\link{mclapply}}
#' @param ... Passed to \code{\link{sampPostFracs}}
#' @return A list containing all the samples from the posterior of the pairwise correlations and a table summarising the correlations.
calcCor = function(dd,doPlot=TRUE,colFun=circlize::colorRamp2(c(-1,0,1),c('blue','white','red')),sumFun=median,sampleRho=TRUE,verbose=1,mcParams=list(),...){
  samps = sampPostFracs(dd,...)
  cc = unique(dd$cellType)
  nn = unique(dd$donor)
  idx = cbind(match(dd$cellType,cc),match(dd$donor,nn))
  corVals = list()
  for(i in seq(nSamps)){
    tt = matrix(NA,nrow=length(cc),ncol=length(nn))
    rownames(tt) = cc
    colnames(tt) = nn
    tt[idx]= samps$relFrac[,i]
    corVals[[i]] = cor(t(tt),use='pairwise.complete.obs')
  }
  #Get number of datapoints for mutual pairs
  nPairs = matrix(NA,nrow=length(cc),ncol=length(cc))
  rownames(nPairs)=colnames(nPairs)=cc
  for(a in cc){
    for(b in cc){
      nPairs[match(a,cc),match(b,cc)] = sum(table(dd[dd$cellType %in% c(a,b),'donor'])==2)
    }
  }
  #Work out which indicies we need to run on.  As symmetric and diagonals trivial, only need to operate on the upper triangle
  w = which(nPairs>=4,arr.ind=TRUE)
  w = w[w[,2]>w[,1],]
  corVals = lapply(corVals,function(e) {tmp=e;tmp[seq_along(tmp)]=NA;tmp[w]=e[w];tmp})
  if(sampleRho){
    #The slow bit, sample over sampling distribution...
    defParams = list(X = corVals,
                    FUN = function(e){
                      for(i in seq(nrow(w))){
                        e[w[i,,drop=FALSE]] = sampleExactConfRho(e[w[i,,drop=FALSE]],nPairs[w[i,,drop=FALSE]],1)
                      }
                      return(e)},
                    mc.allow.error=TRUE,
                    mc.allow.fatal=TRUE,
                    mc.silent=verbose<1,
                    mc.retry.silent = verbose<1,
                    mc.dump.frames='no',
                    mc.shm.ipc=FALSE,
                    mc.share.copy=FALSE,
                    mc.share.vectors=FALSE,
                    mc.progress= interactive() && verbose,
                    mc.preschedule=TRUE,
                    mc.retry=10,
                    mc.cores=nParallel)
    for(nom in names(mcParams))
      defParams[[nom]] = mcParams[[nom]]
    corVals = do.call(mclapply,defParams)
  }
  #Make back into symmetric
  corVals = lapply(corVals,function(e) {e[lower.tri(e)]=t(e)[lower.tri(e)];e})
  ###corVals = list()
  ###for(i in seq(nSamps)){
  ###  if(verbose)
  ###    message(sprintf('Random sample %d of %d',i,nSamps))
  ###  tt = matrix(NA,nrow=length(cc),ncol=length(nn))
  ###  rownames(tt) = cc
  ###  colnames(tt) = nn
  ###  tt[idx]= samps$relFrac[,i]
  ###  corVals[[i]] = cor(t(tt),use='pairwise.complete.obs')
  ###  #Slow as anything.  Fold in variability due to sampling
  ###  for(ii in seq(ncol(corVals[[i]]))){
  ###    for(jj in seq(nrow(corVals[[i]]))){
  ###      #Don't touch the nonsense self correlation
  ###      if(ii==jj)
  ###        next
  ###      n = sum(table(dd[dd$cellType %in% c(rownames(corVals[[i]])[jj],colnames(corVals[[i]])[ii]),]$donor)==2)
  ###      r = corVals[[i]][ii,jj]
  ###      #If we have too few for a confidence interval, we have too few to say anything
  ###      if(n<4){
  ###        corVals[[i]][ii,jj]=NA
  ###      }else{
  ###        if(sampleRho)
  ###          corVals[[i]][ii,jj] = sampleExactConfRho(r,n,1)
  ###      }
  ###    }
  ###  }
  ###}
  #Drop any that are uninformative
  tmp = matrix(NA,nrow=length(cc),ncol=length(cc),dimnames=list(cc,cc))
  for(i in seq_along(cc)){
    for(j in seq_along(cc)){
      tmp[i,j] = sumFun(sapply(corVals,function(e) e[i,j]))
    }
  }
  good = which(colSums(!is.na(tmp))>1 & rowSums(!is.na(tmp))>1)
  if(length(good)==0)
    stop("Catestrophic failure, all NAs")
  cc = cc[good]
  tmp = tmp[good,good]
  corVals = lapply(corVals,function(e) e[good,good])
  #p(arctanh(rho)==x) = dnorm(mean=arctanh(r),sd=1/sqrt(n-3))
  #p(rho==x) = dnorm(mean
  #Flatten into table and plot if we're plotting
  aa = list()
  hc=NULL
  oo=seq_along(cc)
  if(doPlot){
    x=y=length(cc)
    labIdxs=gridLayout(x*y,x,y)
    hc = hclust(dist(tmp))
    oo = hc$order
  }
  corVals = lapply(corVals,function(e) e[oo,oo])
  cc = cc[oo]
  for(i in seq_along(cc)){
    iType=cc[i]
    for(j in seq_along(cc)){
      jType=cc[j]
      tt = unlist(lapply(corVals,function(e) e[i,j]))
      aa[[paste0(iType,'_',jType)]] = tt
      #Plot if we're going to
      if(doPlot){
        tmp = setGridMargins(x,y,i,j,labIdxs,bottomMar=1,leftMar=1)
        xaxt = tmp$xaxt
        yaxt = tmp$yaxt
        flag = TRUE
        if(all(is.na(tt)) || i==j){
          plot.new()
          flag=FALSE
        }else{
          #Plot twice, overwriting
          for(ii in seq(2)){
            hist(tt,
                 main='',
                 xaxt='n',
                 yaxt='n',
                 xlim=c(-1,1),
                 xlab='',
                 ylab='',
                 breaks=50,
                 border=NA,
                 col='black'
            )
            if(ii==1){
              rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = colFun(sumFun(tt)))
              par(new=TRUE)
            }
          }
        }
        par(xpd=NA)
        if(yaxt=='s')
          title(ylab=iType,line=0)
        if(xaxt=='s')
          title(xlab=jType,line=0)
        par(xpd=FALSE)
        if(flag)
          abline(v=c(-1,0,1),lty=c(2,1,2),col='black')
      }
      par(bg='white')
    }
  }
  #Clean up
  aa = aa[!duplicated(sapply(lapply(strsplit(names(aa),'_'),sort),paste,collapse='_'))]
  aa = aa[order(sapply(aa,median))]
  bb = lapply(aa,quantile,na.rm=TRUE)
  bb = do.call(rbind,bb)
  bb = bb[!duplicated(sapply(lapply(strsplit(rownames(bb),'_'),sort),paste,collapse='_')),]
  bb = bb[order(-bb[,3]),]
  return(list(allCorSamps=aa,corTable=bb,hc=hc,corVals=corVals))
}

