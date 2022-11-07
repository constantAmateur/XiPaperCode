library(ComplexHeatmap)

#' Make standard plots for a given tissue
#' 
#' Makes standard set of plots from a Seurat object with annotation and state columns.
#'
#' The logic for the confidence interval of the sampling distribution is as follows.  It can be shown that X^2 = ((n-1)S^2)/sig^2 is chi-squared distributed with n-1 degrees of freedom, where S^2 is the variance estimated from the data, n is the number of data points, and sig is the true standard deviation.  From this it follows that,
#' p(X^2==x) = dchisq(x,n-1)
#' p(((n-1)S^2)/sig^2 = x) = dchisq(x,n-1)
#' p(sig^2 = ((n-1)S^2)/x) = dchisq(x,n-1)
#' If we now set z=((n-1)S^2)/x
#' p(sig^2 = z) = dchisq(((n-1)S^2)/z,n-1)
#' We can then use this formula to sample from the distribution of the true variance.  Of course, this whole thing works on the assumption that n are normally distributed, which is probably only true-ish, but for large enough n the central limit theorem will save us anyway.
#'
#' @param srat The seurat object.
#' @param plotBase The base filename to use to save plots.  Plot names will be formed by appending _plotType.pdf.
#' @param dropSamps Drop these samples for plotting.
#' @param dropCells Drop these cell types for plotting.
#' @param annotCol Column containing annotation.
#' @param meanTrim When using trimmed mean, how much to trim.
#' @param sampleVariance For rough population size estimates, should we account for the variability due to the limited number of individuals by sampling from the sampling distribution of the variance?
makeStandardPlots = function(srat,plotBase,dropSamps=c(),dropCells=c(),annotCol='annot',meanTrim=0.1,sampleVariance=TRUE){
  dd = standardiseX(srat,annotCol=annotCol)
  dd = dd[!dd$cellType %in% dropCells,]
  dd = dd[!dd$donor %in% dropSamps,]
  cc = unique(dd$cellType)
  nn = unique(dd$donor)
  ###############################
  # Cell type specific clonality
  #pdf(paste0(plotBase,'_clonality_median.pdf'))
  #  dd = calcClonal(dd,baseToUse='median')
  #dev.off()
  #pdf(paste0(plotBase,'_clonality_tau.pdf'))
  #  dd = calcClonal(dd,baseToUse='tau')
  #dev.off()
  #pdf(paste0(plotBase,'_clonality_mean.pdf'))
  #  dd = calcClonal(dd,baseToUse='mean',meanTrim=meanTrim)
  #dev.off()
  pdf(paste0(plotBase,'_clonality_fixed.pdf'))
    dd = calcClonal(dd,baseToUse='fixed')
  dev.off()
  ##############
  # Correlation
  pdf(paste0(plotBase,'_correlation_fixed.pdf'),width=14,height=14)
    fixedCor = calcCor(dd,baseToUse='fixed')
  dev.off()
  #pdf(paste0(plotBase,'_correlation_relative_median.pdf'),width=14,height=14)
  #  relativeCor = calcCor(dd,baseToUse='median')
  #dev.off()
  #pdf(paste0(plotBase,'_correlation_relative_tau.pdf'),width=14,height=14)
  #  relativeCor = calcCor(dd,baseToUse='tau')
  #dev.off()
  #pdf(paste0(plotBase,'_correlation_relative_mean.pdf'),width=14,height=14)
  #  relativeCor = calcCor(dd,baseToUse='mean',meanTrim=meanTrim)
  #dev.off()
  ###########
  # Pop size
  pdf(paste0(plotBase,'_pop_size_rough.pdf'))
    postSamps = sampPostFracs(dd,baseToUse='fixed')
    x=list()
    for(i in seq_along(cc)){
      cType=cc[i]
      w = which(dd$cellType==cType)
      x[[cType]] = apply(postSamps$sampFrac[w,],2,var)
      if(sampleVariance)
        x[[cType]] = (length(w)-1)*x[[cType]]/rchisq(length(x[[cType]]),length(w)-1)
      x[[cType]] = 1/(4*x[[cType]])
    }
    par(mar=c(5,8,3,1))
    x = x[order(sapply(x,median))]
    #Add in global
    tmp = rbeta(nSamps*length(srat@misc$fitStats$taus),
                srat@misc$fitStats$taus*srat@misc$fitStats$nCells+1,
                (1-srat@misc$fitStats$taus)*srat@misc$fitStats$nCells+1)
    tmp = matrix(tmp,nrow=length(srat@misc$fitStats$taus))
    x[['Global']] = apply(tmp,2,var)
    if(sampleVariance)
      x[['Global']] = (nrow(tmp)-1)*x[['Global']]/rchisq(length(x[['Global']]),nrow(tmp)-1)
    x[['Global']] = 1/(4*x[['Global']])
    boxplot(x,
            ylim=c(0,32),
            horizontal=TRUE,
            main='Inferred pop size',
            las=1
            )
    #Add in points that show the naiive estimates of pop-size with no variability accounted for
    a=1/(4*sapply(split(dd$matCount/dd$tot,dd$cellType),var))
    a['Global'] = 1/(4*var(srat@misc$fitStats$taus))
    a = a[names(x)]
    points(a,seq_along(a),pch=19,col='red')
    abline(v=a['Global'],col='red')
  dev.off()
  #######
  # Misc
  #Raw fractions
  pdf(paste0(plotBase,'_raw_fractions.pdf'),width=14,height=14)
    postSamps = sampPostFracs(dd,baseToUse='fixed')
    #Want to split by sample/celltype combination.  Which is each row...
    par(mfrow=gridDimensions(length(cc)))
    par(mar=c(4,5,1,0.5))
    x = par('mfrow')[1]
    y = par('mfrow')[2]
    nBottom = ((length(cc)-1) %% y)+1
    for(i in seq_along(cc)){
      cType=cc[i]
      w = which(dd$cellType==cType)
      tmp = matrix(NA,ncol=length(nn),nrow=nSamps)
      colnames(tmp) = nn
      xx = t(postSamps$sampFrac[w,])
      tmp[,match(dd$donor[w],nn)]=xx
      xaxt = ifelse(ceiling(i/y)==x || (ceiling(i/y)==x-1 && (((i-1) %% y)+1) > nBottom),'s','n')
      boxplot(tmp,
              las=2,
              xaxt=xaxt,
              ylim=c(0,1),
              main=cType
              )
      abline(h=0.5,col='red')
    }
  dev.off()
  #Pairwise population size
  out = matrix(NA,nrow=length(cc),ncol=length(cc),dimnames=list(cc,cc))
  for(i in seq_along(cc)){
    for(j in seq_along(cc)){
      w = which(srat@meta.data[,annotCol] %in% cc[c(i,j)])
      w=unlist(lapply(split(w,srat@meta.data[w,annotCol]),sample,min(table(srat@meta.data[w,annotCol]))))
      out[i,j]= 1/(4*var(sapply(split(srat@meta.data$states[w],srat@meta.data$donor[w]),mean,na.rm=TRUE),na.rm=TRUE))
    }
  }
  #a = ((out-outer(diag(out),rep(1,nrow(out)))))
  #Heatmap(sign(a)*log(abs(a)+1),cluster_rows=FALSE,cluster_columns=FALSE)
  pdf(paste0(plotBase,'_pairwise_pop_size.pdf'))
    hm = Heatmap(out,
                 name='PopSize'
                 )
    draw(hm)
  dev.off()
  ####Get a maximally mixed estimate
  ###w = seq(nrow(srat@meta.data))
  ###w=unlist(lapply(split(w,srat@meta.data$annotFull[w]),sample,min(table(srat@meta.data$annotFull[w]))))
  ###x = 1/(4*var(sapply(split(srat@meta.data$states[w],srat@meta.data$donor[w]),mean,na.rm=TRUE),na.rm=TRUE))  
  return(dd)
}

