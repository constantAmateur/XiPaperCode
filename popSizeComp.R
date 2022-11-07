#' Use the fit models to infer the effective population bottleneck size (presumably no cells contributing from epiblast)
#'
#' This calls and pulls this information from multiple tissues

#############
# Libraries #
#############

library(ggplot2)
library(reshape2)
library(inactiveXX)
library(miloR)

##########
# Params #
##########

tgtFiles = c(Lung='Results/HCA_lung.RDS',
             Vagina = 'Results/HCA_vag.RDS',
             Oral = 'Results/oralMucosa.RDS',
             PBMC = 'Results/HCA_AIDA.RDS',
             Kidney = 'Results/kidney.RDS',
             Placenta = 'Results/placenta_F.RDS')
annotCols = c(Lung='annotCl',
              Kidney = 'annotFull')
custFilt = list(PBMC = function(e) {ageCut=35
                tmp = e[,e@meta.data$donorAge < ageCut]
                w = which(tmp@misc$fitStats$ages<ageCut)
                tmp@misc$fitStats = lapply(tmp@misc$fitStats,function(e) e[w])
                return(tmp)},
                Kidney = function(e) e$mKid
                )
sampleVariance=TRUE
plotDir='Results/plots'
#Alias together similar cell types
aliases = list(Tcells=c('TCells','Tcell','Tcells'),
               VascEndo = c('Endothelium','VascularEndo'),
               Erythroid = c('Late erythroid','Erythroid'),
               Plasma=c('Plasma','Plasma cells'),
               NK = c('NK','NK cells'),
               MNP=c('Macrophages','Classical monocytes','Macrophage','Macro'))
aliasFlat = setNames(rep(names(aliases),lengths(aliases)),unlist(aliases))
#Cell types in placenta that are definitely placental in origin
placentaPops = c('EVT','SCT','VCT','Hofbauer')


#################
# Call pop size #
#################

out=list()
dds = list()
for(nom in names(tgtFiles)){
  srat = readRDS(tgtFiles[nom])
  if(!is.null(custFilt[[nom]]))
    srat = custFilt[[nom]](srat)
  #Get annotation column
  annot=annotCols[nom]
  if(is.na(annot))
    annot='annot'
  #Standardise
  dds[[nom]] = standardiseX(srat,annotCol=annot)
  dd = dds[[nom]]
  cc = unique(dd$cellType)
  nn = unique(dd$donor)
  #Get pop size
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
   #Add in points that show the naiive estimates of pop-size with no variability accounted for
  a=1/(4*sapply(split(dd$matCount/dd$tot,dd$cellType),var))
  a['Global'] = 1/(4*var(srat@misc$fitStats$taus))
  a = a[names(x)]
  #Save them
  out[[nom]] = list(pointEstimates = a,uncertaintySamples = x)
}

###############
# Plot global #
###############

#Tissue specific global.  Calculate epiblast and hypoblast separately.
x = lapply(out[!names(out) %in% c('Placenta')],function(e) e$uncertaintySamples$Global)
#Add in an average of the definitive hypoblast tissues in the placenta
srat = readRDS(tgtFiles['Placenta'])
tmp = srat@meta.data[which(srat@meta.data$annot %in% placentaPops),]
cnts = aggregate(cbind(nMat=states,nPat=1-states) ~ donor,FUN=sum,data=tmp,na.rm=TRUE)
tmp = rbeta(nSamps*nrow(cnts),
            1 + cnts$nMat,
            1 + cnts$nPat)
tmp = matrix(tmp,nrow=nrow(cnts))
x[['Placenta']] = apply(tmp,2,var)
if(sampleVariance)
  x[['Placenta']] = (nrow(cnts)-1)*x[['Placenta']]/rchisq(length(x[['Placenta']]),nrow(cnts)-1)
x[['Placenta']] = 1/(4*x[['Placenta']])
#Now plot
xx = seq_along(x)
pdf(file.path(plotDir,'globalPopSize.pdf'),width=4,height=4)
plot(NA,
     type='n',
     xlab='Tissue',
     ylab='Population size',
     xlim=c(0.5,length(x)+0.5),
     ylim=c(0,32),
     xaxt='n',
     bty='n',
     )
abline(h=c(mean(sapply(x,median)),median(x[['Placenta']])),col=c('darkred','darkgreen'))
points(xx,sapply(x,median),pch=16,cex=2)
axis(1,xx,labels=names(x))
# Add error bars
arrows(x0=xx, y0=sapply(x,quantile,0.05), x1=xx, y1=sapply(x,quantile,0.95), code=3, angle=90, length=0.1)
dev.off()

######################
# Cell type pop size #
######################

#Get shared cell type list 
cTypes = lapply(out,function(e) names(e$uncertaintySamples))
noms = unique(unlist(cTypes))
cTypes = data.frame(row.names=noms,do.call(cbind,lapply(cTypes,function(e) noms%in%e)))
shared = setdiff(rownames(cTypes[rowSums(cTypes)>1,]),'Global')
x = lapply(out,function(e) e$uncertaintySamples[setdiff(names(e$uncertaintySamples),'Global')])
#Filter out the non-hypoblasts from placenta
x$Placenta = x$Placenta[placentaPops]
#Cell type annotation: leukocytes, endothelial, epithelial, mesenchymal, nervous
cellGroups = list(Lung = list(Epithelial = c("Goblet (nasal)","SMG mucous","Goblet (bronchial)","AT1",'AT2',"Ionocyte","Multiciliated (non-nasal)","Basal"),
                              Endothelial = c("LymphEndo","VascularEndo" ),
                              Leukocytes = c("Erythroid","Monocytes","Plasma","Classical monocytes","Tcells","B cells","Mast cells","NK cells","DC2","Alveolar macrophages"),
                              Mesenchymal =c("Smooth muscle","Fibroblasts")),
                  Vagina = list(Leukocytes = c("Erythroid","Plasma","BCells","Mast","Macro","TCells"),
                                Endothelial = c("LymphEndo","Endothelium"),
                                Epithelial =c("Myoepithelium","Epithelium"),
                                Mesenchymal = c("Fibroblasts","SmoothMuscle"),
                                Neuronal =c("Neurons")),
                  Oral = list(Endothelial = c("LymphEndo","VascularEndo"),
                              Leukocytes = c("Erythroid","Plasma","Bcell","Mast","NK","Macrophages","Tcells"),
                              Epithelial = c("Epithelium"),
                              Mesenchymal = c("Perivascular","Fibroblasts" )),
                  PBMC = list(Leukocytes = c("Late erythroid", "Megakaryocytes/platelets", "DC1", "HSC/MPP", 
"Memory CD4+ cytotoxic T cells", "Plasma cells", "pDC", "NK cells", 
"DC2", "Cytotoxic T cells", "Non-classical monocytes", "Naive B cells", 
"Regulatory T cells", "gdT", "Classical monocytes", "MAIT cells", 
"Tcm/Naive cytotoxic T cells", "Memory B cells", "Tem/Effector helper T cells", 
"Tcm/Naive helper T cells", "CD16+ NK cells")),
                  Kidney = list(Leukocytes = c("Erythroid","Bcell","NK","Tcell","Macrophage"),
                                Epithelial = c("Podocyte","DistalTubules","LoH","IntercalatedCell","PT",'Pelvis'),
                                Mesenchymal = c("Myofibroblast"),
                                Endothelial = c("Endothelium")),
                  Placenta = list(Mesenchymal = c( "EVT","SCT","VCT"),
                                  Leukocytes = c("Hofbauer"))
                  )
cellGroups = lapply(cellGroups, function(e) setNames(rep(names(e),lengths(e)),unlist(e)))
#Sort based on this information
#x = setNames(lapply(names(x),function(e) x[[e]][order(cellGroups[[e]][names(x[[e]])],sapply(x[[e]],median))]),names(x))
pdf(file.path(plotDir,'cellSpecificPopSize.pdf'),width=7,height=4)
gap=1
cols = list(globalQuants = paste0('#000000',as.hexmode(round(255*0.2))),
            globalMedian = paste0('#000000',as.hexmode(round(255*0.8))))
#Create a background showing the global inferred size for each tissue typ
par(mar=c(7,5,1,1))
plot(NA,
     xlim=c(0.5,sum(lengths(x))+gap*(length(x)-1)+0.5),
     ylim=c(0,32),
     xaxt='n',
     xlab='',
     ylab='Population size',
     bty='n',
     type='n')
left=0.5
i=1
ats=c()
labs=c()
for(tgt in names(x)){
  quants = quantile(out[[tgt]]$uncertaintySamples$Global,c(0.05,0.5,0.95))
  right = left+length(x[[tgt]])
  rect(xleft=left,
       xright = right,
       ybottom=quants[1],
       ytop=quants[3],
       border=NA,
       col=cols$globalQuants)
  lines(c(left,right),rep(quants[2],2),col=cols$globalMedian,lwd=1)
  xx = seq(left+0.5,right,1)
  ats = c(ats,xx)
  labs = c(labs,names(x[[tgt]]))
  points(xx,sapply(x[[tgt]],median),pch=16,cex=1,col='black')#match(cellGroups[[tgt]][names(x[[tgt]])],unique(unlist(cellGroups))))
  # Add error bars
  arrows(x0=xx, y0=sapply(x[[tgt]],quantile,0.05), x1=xx, y1=sapply(x[[tgt]],quantile,0.95), code=3, angle=90, length=0.05)
  left = right+gap
}
axis(1,at=ats,labels=labs,las=2,cex.axis=0.4,gap.axis=-1)
dev.off()
popSizes=x

#########################
# Pop size vs clonality #
#########################

pdf(file.path(plotDir,'popSizeVsCloneSize.pdf'),width=5,height=4)
plot(NA,
     xlim=c(1,32),
     ylim=c(0,1),
     type='n',
     xlab='PopSize',
     ylab='Clonal fraction')
dClone = list()
for(tissue in names(dds)){
  if(tissue=='Placenta')
    next
  dd = dds[[tissue]]
  samps = sampPostFracs(dd,baseToUse='tau')
  #Geth the pop Size estimate
  s = split(seq(nrow(dd)),dd$cellType)
  popSize = sapply(s,function(e) 1/(4*apply(samps$sampFrac[e,],2,var)))
  cloneFrac = ifelse(samps$sampFrac>=samps$baseRate,
                     samps$relFrac/(1-samps$baseRate),
                     -samps$relFrac/samps$baseRate)
  cloneSize = sapply(s,function(e) apply(cloneFrac[e,],2,median))
  for(i in seq_along(s)){
    xx = median(popSize[,i])
    yy = median(cloneSize[,i])
    #points(xx,yy,pch=19,cex=.5)
    arrows(x0=xx,
           x1=xx,
           y0=quantile(cloneSize[,i],0.05),
           y1=quantile(cloneSize[,i],0.95),
           lwd=0.1,
           col='#00000095',
           code=3,
           angle=90,
           length=0.1)
    arrows(x0=quantile(popSize[,i],0.05),
           x1=quantile(popSize[,i],0.95),
           lwd=0.1,
           col='#00000095',
           y0=yy,
           y1=yy,
           code=3,
           angle=90,
           length=0.1)
  }
  dClone[[length(dClone)+1]] = data.frame(cellType = names(s),
                                          tissue=tissue,
                                          popSize=apply(popSize,2,median),
                                          cloneSize=apply(cloneSize,2,median))
}
dClone = do.call(rbind,dClone)
dClone$cellTypeFixed = ifelse(dClone$cellType %in% unlist(aliases),aliasFlat[dClone$cellType],dClone$cellType)
points(dClone$popSize,dClone$cloneSize,pch=19,cex=0.5)
w = (dClone$popSize<6)
low =lm(cloneSize ~ popSize,data=dClone[w,]) 
abline(low)
high = lm(cloneSize~popSize,data=dClone[!w,])
abline(high)
dev.off()


#####################
# Common cell types #
#####################

#Do the same thing as above, but by cell type
toGet = rownames(cTypes[rowSums(cTypes)>1,])
toGet = c(toGet,unlist(aliases),names(aliases))
x = lapply(out,function(e) {tmp = e$uncertaintySamples[intersect(toGet,names(e$uncertaintySamples))];names(tmp) = ifelse(names(tmp) %in% unlist(aliases),aliasFlat[names(tmp)],names(tmp));tmp})
#Flatten and reorder
tmp = list()
for(i in x)
  tmp = c(tmp,i)
tissues = rep(names(x),lengths(x))
#Work out what the median pop size is by type and sort by it
typeSize = sort(sapply(split(sapply(tmp,median),names(tmp)),median))
o = order(factor(names(tmp),levels=names(typeSize)))
tmp = tmp[o]
tissues = tissues[o]
#Plot them
gap=5
#Create a background showing the global inferred size for each tissue typ
pdf(file.path(plotDir,'popSizeSharedCells.pdf'),width=7,height=4)
par(mar=c(7,5,1,1))
plot(NA,
     xlim=c(0.5,length(tmp)+gap*(length(unique(names(tmp)))-1)+0.5),
     ylim=c(0,32),
     xaxt='n',
     xlab='',
     ylab='Population size',
     bty='n',
     type='n')
left=0.5
ats=c()
labs=c()
for(tgt in unique(names(tmp))){
  right = left + sum(names(tmp)==tgt)
  xx = seq(left+0.5,right,1)
  w = which(names(tmp)==tgt)
  ats = c(ats,median(xx))
  labs = c(labs,tgt)
  points(xx,sapply(tmp[w],median),pch=16,cex=1)#,col=match(tissues[w],unique(tissues)))
  # Add error bars
  arrows(x0=xx, y0=sapply(tmp[w],quantile,0.05), x1=xx, y1=sapply(tmp[w],quantile,0.95), code=3, angle=90, length=0.05)
  left = right+gap
}
axis(1,at=ats,labels=labs,las=2)
dev.off()
###############
# Raw material
cloneFrac = lapply(dds,calcClonal,doPlot=FALSE,alpha=0.1,returnFull=TRUE)
#Make global versions of both bits
dd = do.call(rbind,lapply(cloneFrac,function(e) e$dd))
cc = do.call(rbind,lapply(cloneFrac,function(e) e$cloneFrac))
#Fix cell type
dd$cellTypeFixed = ifelse(dd$cellType %in% unlist(aliases),aliasFlat[dd$cellType],dd$cellType)
dd$tissue = gsub('\\.[0-9]+$','',rownames(dd))
#Drop placenta, those not well formulated, and those that are not in at least 4 tissues
nTissues = sapply(split(dd$tissue,dd$cellTypeFixed),function(e) length(unique(e)))
w = which(!dd$tissue %in% c('Placenta') & dd$cellTypeFixed %in% toGet & dd$cellTypeFixed %in% names(nTissues)[nTissues>=4])
dd = dd[w,]
cc = cc[w,]
##############################
# Per cell type distributions
#Get median, and quantiles to plot for each 
cellGap = 0.2
tissueGap= 0.05
inkSize=20
x = dd
#Plot by celltype, then by tissue, with gaps
shortLabs=c(Kidney='K',Oral='O',Placenta='P',PBMC='B',Lung='L',Vagina='V')
tmp = lapply(split(x,x$cellTypeFixed),function(e) e[order(e$tissue,e$clonalFrac),])
#Order by average across cell type
tmp = tmp[names(sort(-sapply(tmp,function(e) median(e$clonalFrac))))]
pdf(file.path(plotDir,'cellLevelClonality.pdf'),width=7,height=4)
plot(NA,
     xlim=c(0,length(tmp)),
     ylim=c(0,1),
     bty='n',
     xaxt='n',
     ylab='Clonality',
     xlab='',
     type='n')
left=0
ats=c()
labs=c()
#Loop through cell types 
left = cellGap/2
for(cell in names(tmp)){
  #And over tissue
  tt = tmp[[cell]]
  right = left+1-cellGap
  med = sapply(split(tt$clonalFrac,tt$tissue),median)
  #Work out the spacing between points and position them
  xpos = (right-left - (length(med)-1)*tissueGap)/(nrow(tt))
  xpos = left + (seq(nrow(tt))-1)*xpos + xpos*0.5 + (as.numeric(factor(tt$tissue,levels=unique(tt$tissue)))-1)*tissueGap
  #xpos = (right-left)/(nrow(tt)+1)
  #xpos = seq(left+xpos,right-xpos,xpos)
  lineCol = toupper(as.hexmode(round(255*min(1,inkSize/nrow(tt)))))
  lineCol= paste0('#000000',ifelse(nchar(lineCol)==1,'0',''),lineCol)
  for(i in seq_along(xpos))
    lines(c(xpos[i],xpos[i]),c(tt$clonalFracLow[i],tt$clonalFracHigh[i]),lwd=0.1,col=lineCol)
  points(xpos,tt$clonalFrac,cex=0.1,pch=19)
  #lines(c(left,right),rep(median(tt$clonalFrac),2),lwd=1)
  #Also make tissue type specific medians
  for(i in seq_along(med))
    lines(range(xpos[tt$tissue==names(med)[i]]),rep(med[i],2),lwd=1)
    #lines(c(left,right),rep(med[i],2),lwd=0.2,col=match(names(med)[i],names(shortLabs)))
  #Construct axis ticks for each
  t2=sapply(split(xpos,tt$tissue),median)
  axis(1,at=t2,labels=shortLabs[names(t2)],lwd=0,lwd.ticks=1,line=0,gap.axis=-1)
  ats=c(ats,0.5*(left+right))
  labs = c(labs,cell)
  left=right+cellGap
}
axis(1,at=ats,labels=labs,lwd=0,lwd.ticks=0,line=1,gap.axis=-1)
dev.off()
########################
# Per-tissue histograms
nBins=20
tissueGap=1
cellgap=4
tissueSize=8
maxFrac=0.4
shortLabs=c(Kidney='K',Oral='O',Placenta='P',PBMC='B',Lung='L',Vagina='V')
tmp = lapply(split(dd,dd$cellTypeFixed),function(e) lapply(split(e,e$tissue),function(ee) ee[order(ee$clonalFrac),]))
tmp = tmp[names(sort(-sapply(split(dd$clonalFrac,dd$cellTypeFixed),median)))]
plot(NA,
     xlim=c(0,sum(lengths(tmp)*tissueSize+(lengths(tmp)-1)*tissueGap)+(length(tmp)-1)*cellGap),
     ylim=c(0,1),
     bty='n',
     xaxt='n',
     ylab='Clonality',
     xlab='',
     type='n')
#Loop over fixed cell types
left=0
ats=c()
labs=c()
for(cType in names(tmp)){
  wCell = dd$cellTypeFixed==cType
  #And tissue type
  tissues = unique(dd$tissue[wCell])
  cellLeft = left
  tissAts=c()
  tissLabs=c()
  for(tissue in tissues){
    #if(cType=="Fibroblasts" && tissue=='Oral')
    #  break
    right = left+tissueSize
    w = wCell & dd$tissue==tissue
    #Calculate the histogram for each
    x = hist(cc[w,],breaks=seq(0,1,length.out=nBins+1),plot=FALSE)
    tissAts=c(tissAts,(right+left)*0.5)
    tissLabs=c(tissLabs,tissue)
    #Make a line that we're drawing from
    lines(c(left,left),c(0,1))
    #Do the plot
    rect(xleft=rep(left,length(x$counts)),
         #xright = left + (x$counts/sum(w)/ncol(cc))*tissueSize/maxFrac,
         xright = left + (x$counts/max(x$counts))*tissueSize,
         ybottom=x$breaks[-length(x$breaks)],
         ytop=x$breaks[-1],
         col='grey')
    lines(c(left,right),rep(median(cc[w,]),2),lwd=1,col='blue')
    lines(rep(0.5*(right+left),2),quantile(cc[w,],c(0.25,0.75)),lwd=1,col='blue')
    left = right+tissueGap
    #readline(prompt="Press [enter] to continue")
  }
  lines(c(cellLeft,left),rep(median(cc[wCell,]),2),lwd=2,col='red')
  lines(rep(0.5*(cellLeft+left),2),quantile(cc[wCell,],c(0.25,0.75)),lwd=2,col='red')
  axis(1,at=tissAts,labels=shortLabs[tissLabs],lwd=0,lwd.ticks=1,cex.axis=0.4,gap.axis=-1)
  ats=c(ats,(cellLeft+left)*0.5)
  labs = c(labs,cType)
  left = left+cellGap
}
axis(1,at=ats,labels=labs,lwd=0,lwd.ticks=1,line=3,gap.axis=-1)
################################
# Per-sample with uncertainties
tissueSize=8
tissueGap=1
cellGap = 6
inkSize=8
x = dd
#Plot by celltype, then by tissue, with gaps
shortLabs=c(Kidney='K',Oral='O',Placenta='P',PBMC='B',Lung='L',Vagina='V')
x = x[order(x$cellTypeFixed,x$tissue),]
tmp = lapply(split(x,x$cellTypeFixed),function(e) lapply(split(e,e$tissue),function(ee) ee[order(ee$clonalFrac),]))
#Order by average across cell type
tmp = tmp[names(sort(-sapply(split(x$clonalFrac,x$cellTypeFixed),median)))]
plot(NA,
     xlim=c(0,sum(lengths(tmp)*tissueSize+(lengths(tmp)-1)*tissueGap)+(length(tmp)-1)*cellGap),
     ylim=c(0,1),
     bty='n',
     xaxt='n',
     ylab='Clonality',
     xlab='',
     type='n')
left=0
ats=c()
labs=c()
#Loop through cell types 
for(cell in names(tmp)){
  #And over tissue
  cellLeft = left
  tissAts=c()
  tissLabs=c()
  for(tissue in names(tmp[[cell]])){
    tt = tmp[[cell]][[tissue]]
    right = left+tissueSize
    tissAts=c(tissAts,(right+left)*0.5)
    tissLabs=c(tissLabs,tissue)
    xpos = (right-left)/(nrow(tt)+1)
    xpos = seq(left+xpos,right-xpos,xpos)
    for(i in seq_along(xpos))
      lines(c(xpos[i],xpos[i]),c(tt$clonalFracLow[i],tt$clonalFracHigh[i]),lwd=0.1,col=paste0('#000000',as.hexmode(round(255*min(1,inkSize/nrow(tt))))))
    points(xpos,tt$clonalFrac,cex=0.1,pch=19)
    lines(c(left,right),rep(median(tt$clonalFrac),2),lwd=1)
    left=right+tissueGap
  }
  lines(c(cellLeft,left),rep(median(unlist(lapply(tmp[[cell]],function(e) (e$clonalFrac)))),2),lwd=2)
  axis(1,at=tissAts,labels=shortLabs[tissLabs],lwd=0,lwd.ticks=1,cex.axis=0.4,gap.axis=-1)
  ats=c(ats,(cellLeft+left)*0.5)
  labs = c(labs,cell)
  left = left+cellGap
}
axis(1,at=ats,labels=labs,lwd=0,lwd.ticks=1,line=3,gap.axis=-1)


######################
# Correlation kidney #
######################

#Get correlation coefficients
rhos = list()
nom='Kidney'
for(tgt in c('fixed','mean','median')){
  message(nom)
  rhos[[tgt]] = calcCor(dds[[nom]],baseToUse=tgt)
}
#Fancy plot, define parameters
nBins=50
noms = c(Tubules='Podocyte',
         Tubules='PT',
         Tubules='LoH',
         Tubules='DistalTubules',
         Tubules='IntercalatedCell',
         Stroma='Myofibroblast',
         Stroma='Endothelium',
         Leukocytes='Bcell',
         Leukocytes='Tcell',
         Leukocytes='Macrophage',
         Leukocytes='NK',
         Leukocytes='Erythroid')
xBig = rhos$mean$corVals
#Make a matrix that we use for clustering, etc
xAvg = matrix(NA,nrow=nrow(xBig[[1]]),ncol=ncol(xBig[[1]]),dimnames=dimnames(xBig[[1]]))
for(i in seq(nrow(xAvg))){
  for(j in seq(ncol(xAvg))){
    xAvg[i,j] = mean(unlist(lapply(xBig,function(e) e[i,j])))
  }
}
cell_fun = function(j,i,x,y,width,height,fill){
  #White out the area
  #grid.rect(x = x, y = y, width = width, height = height,gp = gpar(col = NA, fill = 'white'))
  #Grab the raw data
  dat = unlist(lapply(xBig,function(e) e[i,j]))
  #Make histogram
  brks = seq(-1,1,length.out=nBins)
  t1 = hist(dat,breaks=brks,plot=FALSE)
  heights = t1$counts/max(t1$counts)
  #Draw histogram
  for(ii in seq_along(t1$counts)){
    grid.rect(x=x - width/2 + width*((ii-1)/nBins),
              y=y-height/2,
              just=c(0,0),
              width= width/nBins,
              height=heights[ii]*height,
              gp = gpar(col=NA,fill='black'))
  }
  #Vertical line
  grid.lines(x=c(x,x),y=c(y-height/2,y+height/2),gp=gpar(lty=2))
}
hm = Heatmap(xAvg,
             col = circlize::colorRamp2(c(-1,0,1),c('blue','white','red')),
             name='Correlation',
             row_order = noms,
             column_order = noms,
             column_split = names(noms[match(colnames(xAvg),noms)]),
             row_split = names(noms[match(colnames(xAvg),noms)]),
             cell_fun=cell_fun)
pdf(file.path(plotDir,'kidneyCellCorrelation.pdf'),width=7,height=7)
draw(hm)
dev.off()

######################
# Correlation global #
######################

#Get correlation coefficients
rhos = list()
for(nom in names(dds)){
  message(nom)
  rhos[[nom]] = calcCor(dds[[nom]],baseToUse='fixed',fixedVal=0)
}
#Flatten.  Will double count everything, so watch out for that...
x = melt(lapply(rhos,function(e) e$corVals))
#Get rid of duplicate entries (keep only upper-tri)
toKeep=c()
for(tgt in unique(x$L1)){
  cc = unique(as.character(x$Var1[x$L1==tgt]))
  toKeep = c(toKeep,which(match(x$Var1,cc)>match(x$Var2,cc) & x$L1==tgt)) 
}
x = x[toKeep,]
x = x[!is.na(x$value),]
#Fix up cell type names
x$Var1 = as.character(x$Var1)
x$Var2 = as.character(x$Var2)
x$Var1 = ifelse(x$Var1 %in% names(aliasFlat),aliasFlat[x$Var1],x$Var1)
x$Var2 = ifelse(x$Var2 %in% names(aliasFlat),aliasFlat[x$Var2],x$Var2)
#Fix up placenta
x = x[x$L1!='Placenta' | (x$Var1 %in% placentaPops  & x$Var2 %in% placentaPops),]
#Allocation of cell types to germ layers.
germAllocation = list(Kidney=list(mesoderm=c("Podocyte", "LoH", "DistalTubules", "PT", "VascEndo", "IntercalatedCell", "MNP", "NK", "Bcell", "Tcells", "Myofibroblast", "Erythroid")),
                      Lung=list(endoderm=c("Basal", "Multiciliated (non-nasal)", "AT1",  "Goblet (nasal)", "Ionocyte", "Goblet (bronchial)", "SMG mucous", "AT2"),
                                mesoderm=c("Smooth muscle", "Tcells", "Fibroblasts", "VascEndo", "MNP", "Alveolar macrophages", "LymphEndo", "DC2", "Monocytes", "B cells","Plasma", "NK", "Mast cells" )),
                      Oral=list(mesoderm=c("Plasma", "MNP", "NK", "Tcells", "LymphEndo", "VascEndo", "Fibroblasts", "Mast", "Bcell", "Erythroid",'Epithelium','Perivascular')),
                      PBMC=list(mesoderm=c("Tem/Effector helper T cells", "MNP", "Tcm/Naive cytotoxic T cells", "Tcm/Naive helper T cells", "CD16+ NK cells", "Memory B cells", "DC2", "Naive B cells", "Cytotoxic T cells", "Non-classical monocytes", "NK", "Megakaryocytes/platelets", "pDC", "Regulatory T cells", "gdT", "MAIT cells", "Erythroid", "HSC/MPP", "Memory CD4+ cytotoxic T cells", "Plasma", "DC1")),
                      Vagina = list(mesoderm = c("SmoothMuscle", "MNP", "Tcells", "Mast", "VascEndo", "Plasma", "LymphEndo", "Myoepithelium", "Epithelium", "BCells", "Erythroid", "Fibroblasts"),
                                    ectoderm = "Neurons"))
x$germLayer1 = NA
x$germLayer2 = NA
for(tgt in unique(x$L1)){
  if(tgt=='Placenta')
    next
  w = which(x$L1==tgt)
  tmp = germAllocation[[tgt]]
  tmp = setNames(rep(names(tmp),lengths(tmp)),unlist(tmp))
  x$germLayer1[w] = tmp[x$Var1[w]]
  x$germLayer2[w] = tmp[x$Var2[w]]
}
x$compType=ifelse(x$germLayer1==x$germLayer2,'Within','Across')
x$mark = paste0(x$Var1,'_',x$Var2)
x$fullMark = paste0(x$L1,':',x$mark)
#Excluded cell types
excludeTypes = c('Erythroid','Plasma')
pdf(file.path(plotDir,'allCorrelations.pdf'),width=5,height=4)
#Now do some plotting
xx = x[!(x$Var1 %in% excludeTypes | x$Var2 %in% excludeTypes) & x$L1!='Placenta',]
tissueGap=0.2
plot(NA,
     xlim=c(0,length(unique(xx$L1))),
     ylim=c(-1,1),
     ylab='Correlation coefficient',
     xlab='Tissue',
     xaxt='n',
     type='n',
     bty='n')
left= tissueGap/2
#Order as we'd like
t1 = sapply(split(xx$value,xx$fullMark),median)
xx = xx[order(xx$L1,xx$compType,t1[xx$fullMark]),]
tgts = unique(xx$L1)
for(tgt in tgts){
  right = left+(1-tissueGap)
  w = which(xx$L1==tgt)
  t1 = xx[w,]
  t2 = as.data.frame(t(sapply(split(t1$value,t1$mark),quantile,c(0.05,0.5,0.95))))
  t2$compType = t1$compType[match(rownames(t2),t1$mark)]
  #Preserve ordering
  t2 = t2[unique(t1$mark),]
  xpos = seq(left,right,length.out=nrow(t2))
  for(i in seq(nrow(t2)))
    lines(c(xpos[i],xpos[i]),c(t2[i,1],t2[i,3]),
          lwd=0.1,
          col=ifelse(is.na(t2$compType[i]),
                     'black',
                     ifelse(t2$compType[i]=='Within',
                            'darkgrey',
                            'darkgrey'
                            )
                     )
          )
  points(xpos,t2[,2],pch=19,cex=0.1,col=ifelse(t2[,2]>0,'black','lightgrey'))
  #The grand median(s)
  t3 = unique(t2$compType)
  if(length(t3)>1){
    for(t4 in t3){
      lines(range(xpos[t2$compType==t4]),rep(median(t1$value[t1$compType==t4]),2),col='black',lwd=2)
    }
  }else{
    lines(c(left,right),rep(median(t1$value),2),col='black',lwd=2)
  }
  left = right+tissueGap
}
abline(h=0,lty=2,lwd=0.5)
axis(1,at=seq_along(tgts)-0.5,labels=tgts)
dev.off()


################
# Milo example #
################

#Define annotation meta-data for fetal samples
annotGrp = list("MYELOID"= c("DC2", "DC1", "NEUTROPHIL", "PROMONOCYTE",  "PRE_DC2", "MAST_CELL", "MOP", "EOSINOPHIL_BASOPHIL", "DC_PROGENITOR", "PDC", "MYELOCYTE", "MONOCYTE_I_CXCR4", "CYCLING_PDC", "MIGRATORY_DC",  "AS_DC",   "CYCLING_DC", "MONOCYTE_II_CCR2", "MONOCYTE_III_IL1B", "MACROPHAGE_PROLIFERATING", "MACROPHAGE_MHCII_HIGH", "MACROPHAGE_LYVE1_HIGH", "MACROPHAGE_PERI", "LMPP_MLP",  "LANGERHANS_CELLS", "MACROPHAGE_ERY", "MACROPHAGE_IRON_RECYCLING", "MACROPHAGE_KUPFFER_LIKE", "MACROPHAGE_TREM2"), 
                "B CELLS"= c("LARGE_PRE_B", "PRE_PRO_B", "PRO_B", "CYCLING_B", "SMALL_PRE_B", "MATURE_B", "B1", "IMMATURE_B", "LATE_PRO_B", "PLASMA_B"), 
                "ERYTHROID CELLS"= c("EARLY_ERY", "EARLY_MK", "PROMYELOCYTE", "MID_ERY", "LATE_ERY", "LATE_MK","YS_ERY", "CYCLING_YS_ERY"),
                "ILC"= c("ILC3", "ILC2", "CYCLING_ILC"),
                "OTHER"= c("DOUBLET_IMMUNE_FIBROBLAST", "LOW_Q_INCONSISTENT", "DOUBLET_LYMPHOID_MACROPHAGE", "LOW_QUALITY", "HIGH_MITO", "DOUBLETS_FIBRO_ERY", "DOUBLET_ENDOTHELIUM_ERYTHROCYTE", "DOUBLET_ERY_B", "LOW_QUALITY_MACROPHAGE", "LOW_QUALITY_MID_ERY_(HIGH_RIBO)", "PLACENTAL_CONTAMINANTS", "DOUBLET","DOUBLET_VSMC_ERYTHROCYTE"),
                "NK/T CELLS"= c("CYCLING_T", "CD4+T", "CD8+T", "TREG", "NK", "CYCLING_NK", "DP(P)_T", "CD8AA", "ABT(ENTRY)", "DP(Q)_T", "DN(P)_T", "DN(early)_T", "DN(Q)_T", "TYPE_1_INNATE_T", "TYPE_3_INNATE_T"),
                "ENDOTHELIUM"=c("ENDOTHELIUM_I","ENDOTHELIUM_II","ENDOTHELIUM_III", "ENDOTHELIUM_IV","ENDOTHELIUM_V"), 
                "FIBROBLAST"=c("CYCLING_FIBROBLAST_I", "FIBROBLAST_IV", "CYCLING_FIBROBLAST_II",  "FIBROBLAST_I","FIBROBLAST_II","FIBROBLAST_III", "FIBROBLAST_V", "FIBROBLAST_VI", "FIBROBLAST_VII", "FIBROBLAST_VIII","FIBROBLAST_IX", "FIBROBLAST_X", "FIBROBLAST_XI","FIBROBLAST_XII","FIBROBLAST_XIII", "FIBROBLAST_XIV", "FIBROBLAST_XV", "FIBROBLAST_XVI","FIBROBLAST_XVII"),
                "MUSCLE" = c("MYOFIBROBLAST","MYOFIBROBLAST_I","SMOOTH_MUSCLE","MUSCLE_SATELLITE","SKELETAL_MUSCLE"),
                "MELANOCYTE" = c("MELANOCYTE",  "KERATINOCYTE"),
                "NEURONAL" = c("GLIAL","NEURON"),
                "MESOTHELIUM" = c("MESOTHELIUM"),
                "HEPATOCYTE" = c("HEPATOCYTE-LIKE","HEPATOCYTE_I", "HEPATOCYTE_II"),
                "EPITHELIUM" = c("CYCLING_EPITHELIUM","EPITHELIUM_I", "EPITHELIUM_II"),
                "CAJAL" = c("INTERSTITIAL_CELLS_OF_CAJAL"),
                "BONE" = c("OSTEOBLAST",  "CHONDROCYTE","OSTEOCLAST"),
                "NEPHRON" = c("DEVELOPING_NEPHRON_I", "DEVELOPING_NEPHRON_II"),
                "ENTEROENDOCRINE" = c("ENTEROENDOCRINE_I",  "ENTEROENDOCRINE_II"),
                "PERICYTE" = c("VSMC_PERICYTE",    "VSMC_PERICYTE_I",     "VSMC_PERICYTE_II",  "VSMC_PERICYTE_III"),
                "MESENCHYME"= c("MESENCHYMAL_LYMPHOID_TISSUE_ORGANISER"),
                "YS_STROMA" = c("YS_STROMA"),
                "PROGENITORS"= c("MEMP", "GMP", "HSC_MPP", "MEP", "CMP", "CYCLING_MEMP", "CYCLING_MPP")
)
mDat = read.table('Data/PAN.A01.v01.20210429.sample_metadata.csv',sep=',',header=TRUE)
cellDat = read.table('Data/PAN.A01.v01.entire_data_normalised_log.20210429.full_obs.annotated.clean.csv',sep=',',header=TRUE)
cellDat$cellID = paste0(cellDat$file,'_',gsub('.*-','',cellDat$X))
cellDat$annotHigh = setNames(rep(names(annotGrp),lengths(annotGrp)),unlist(annotGrp))[cellDat$anno_lvl_2_final_clean]
#Load the relevant sample
big = readRDS('Results/fitEM/F45_fit.RDS')
#Get annotation
srat = big$srat
m = match(srat@meta.data$orig.ident,mDat$file)
srat@meta.data$organ = mDat$organ[m]
srat@meta.data$cellSort = mDat$Sort_id[m]
srat@meta.data$techMethod = mDat$method[m]
m = match(colnames(srat),cellDat$cellID)
srat@meta.data$annot = cellDat$anno_lvl_2_final_clean[m]
srat@meta.data$annotHigh = cellDat$annotHigh[m]
#Make the call of maternal/paternal more leniant
srat@meta.data$stateX = ifelse(is.na(srat@meta.data$stateProbs) | abs(srat@meta.data$stateProbs-0.5)<0.4,
                               NA,
                               ifelse(srat@meta.data$stateProbs>0.5,
                                      'Maternal',
                                      'Paternal'
                                      )
                               )
#Average within clusters
srat = FindClusters(srat,resolution=10)
x = aggregate(stateProbs ~ seurat_clusters,data=srat@meta.data,FUN=mean,na.rm=TRUE)
srat@meta.data$stateByClust = x$stateProbs[match(srat@meta.data$seurat_clusters,x$seurat_clusters)]
#Implied rho
tau = big$fit$tau
srat@meta.data$rho = ifelse(srat@meta.data$stateByClust>tau,
                            (tau - srat@meta.data$stateByClust)/(tau-1),
                            -1*(tau - srat@meta.data$stateByClust)/(tau)#Make estimates negative to show that the observed fraction was less than the population average.
                            )
#Save for table construction
saveRDS(srat,'Results/fitEM/F45_fit_withAnnot.RDS')
#Now plot the averaged out 
dd = cbind(srat@meta.data,srat@reductions$umap@cell.embeddings)
#Drop ones without annotation
dd = dd[!is.na(dd$annotHigh),]
#First, make a plot of the high level organ type
for(i in seq(2)){
  if(i==1){
    pdf(file.path(plotDir,'F45_organ.pdf'),width=4,height=4)
  }else{
    png(file.path(plotDir,'F45_organ.png'),width=4,height=4,units='in',res=300)
  }
  cmap = c(LI='#855C75',KI='#D9AF6B',BM='#AF6458','#736F4C',SP='#526A83','#625377',TH='#68855C','#9C9C5E','#A06177',SK='#8C785D','#467378','#7C7C7C')
  plot(dd$UMAP_1,dd$UMAP_2,
       xlab='UMAP1',
       ylab='UMAP2',
       xaxt='n',
       yaxt='n',
       cex=0.1,
       pch=19,
       col = cmap[dd$organ])
  if(i==1){
    par(xpd=TRUE)
    legend(x='topleft',
           legend = unique(dd$organ),
           col = cmap[unique(dd$organ)],
           pt.cex=3,
           inset=c(-.15,0),
           pch=19)
  }
  dev.off()
}
#Plot the implied clone size
pdf(file.path(plotDir,'F45_rho.pdf'),width=4,height=4)
divScheme = c('#D6E0FF','#111111','#C2EFB4')
divScheme = c('darkgreen','lightgrey','darkblue')
rhoPts = c(-.3,0,.3)
cmap = circlize::colorRamp2(rhoPts,divScheme)
plot(dd$UMAP_1,dd$UMAP_2,
     xlab='UMAP1',
     ylab='UMAP2',
     cex=0.1,
     pch=19,
     col=cmap(dd$rho))
t1 = seq(rhoPts[3],rhoPts[1],length.out=5)
legend(x='topleft',
       legend=t1,
       col=cmap(t1),
       pch=19,
       pt.cex=3)
dev.off()
#Or just plot the raw difference
pdf(file.path(plotDir,'F45_tau.pdf'),width=4,height=4)
tauPts = c(-.2,0,.2)
cmap = circlize::colorRamp2(tauPts,divScheme)
plot(dd$UMAP_1,dd$UMAP_2,
     xlab='UMAP1',
     ylab='UMAP2',
     cex=0.1,
     pch=19,
     col=cmap(dd$stateByClust-tau))
t1 = seq(tauPts[1],tauPts[3],length.out=5)
legend(x='topleft',
       legend=round(t1+tau,2),
       col=cmap(t1),
       pch=19,
       pt.cex=3)
dev.off()
#Now run milo
sratPass = srat[,!is.na(srat@meta.data$annot)]
tau = mean(sratPass@meta.data$stateX=='Maternal',na.rm=TRUE)
tauCell = split(sratPass@meta.data,sratPass@meta.data$annot)
tauCell = sapply(tauCell,function(e) mean(e$stateX=='Maternal',na.rm=TRUE))
tauCell = median(tauCell,na.rm=TRUE)
milo = runMilo(sratPass,tau=tauCell,resultsPassthrough = c('seurat_clusters','organ','cellSort','techMethod','annot','annotHigh'))
milo$res$rho = ifelse(milo$res$logFC>0,
                      -milo$res$logFC/(tauCell-1),
                      milo$res$logFC/tauCell)
milo$res$rhoAbs = abs(milo$res$rho)
#Hack the milo plot
for(i in seq(2)){
  if(i==1){
    pdf(file.path(plotDir,'F45graphEnrich.pdf'),width=5,height=3)
  }else{
    png(file.path(plotDir,'F45graphEnrich.png'),width=5,height=3,units='in',res=300)
  }
  pl = plotNhoodGraphDA(milo$milo, milo$res, alpha=0.1,res_column='logFC',node_stroke=0.05)
  #Change colour scale
  pl = pl + scale_fill_gradient2(low=divScheme[1],mid=divScheme[2],high=divScheme[3])
  t1 =  ggplot_build(pl)
  t1$data[[1]]$edge_width = t1$data[[1]]$edge_width/10
  plot(ggplot_gtable(t1))
  dev.off()
}
#Swarm by organ
pdf(file.path(plotDir,'F45swarmOrgan.pdf'),width=4,height=4)
pl = plotDAbeeswarm(milo$res, group.by = "organ")
pl = pl + scale_color_gradient2(low=divScheme[1],mid=divScheme[2],high=divScheme[3])
t1 =  ggplot_build(pl)
t1$data[[1]]$size = 1.0
t1$data[[1]]$alpha = 1.0
t1$data[[1]]$stroke = 0
t1$data[[1]]$size = ifelse(t1$data[[1]]$colour=='grey50',1,3)
plot(ggplot_gtable(t1))
dev.off()
#Swarm by celltype
#Work out which cell types have at least a few significant regions
x = split(milo$res$FDR,milo$res$annot)
x = sapply(x,function(e) mean(e<0.05))
x = x[x>0.1]
cTypes = unique(milo$res$annot)
biggest = "DOUBLET_LYMPHOID_MACROPHAGE"
nCols=3
maxTypes = ceiling(length(cTypes)/nCols)
cTypeGroups = split(cTypes,floor((seq_along(cTypes)-1)/maxTypes))
for(j in seq_along(cTypeGroups)){
  cTypeGroup = cTypeGroups[[j]]
  subMilo = milo
  subMilo$res = subMilo$res[subMilo$res$annot %in% cTypeGroup,]
  if(!biggest %in% cTypeGroup){
    message(sprintf("Replacing %s with %s in group %i",cTypeGroup[1],biggest,j))
    subMilo$res$annot[subMilo$res$annot==cTypeGroup[1]]=biggest
  }
  pl = plotDAbeeswarm(subMilo$res, group.by = "annot")
  pl = pl + scale_color_gradient2(low=divScheme[1],mid=divScheme[2],high=divScheme[3])
  t1 =  ggplot_build(pl)
  t1$data[[1]]$size = 1.0
  t1$data[[1]]$alpha = 1.0
  t1$data[[1]]$stroke = 0
  t1$data[[1]]$size = ifelse(t1$data[[1]]$colour=='grey50',1,3)
  pdf(file.path(plotDir,paste0('F45swarmCellType_',j,'.pdf')),width=7,height=7)
  plot(ggplot_gtable(t1))
  dev.off()
}




######################
# Some illustrations #
######################

#Make an animated gif of population size
Ns = seq(2,32)
Ns = c(Ns,rev(Ns))
tfiles=c()
for(N in Ns){
  tfiles = c(tfiles,tempfile())
  png(tfiles[length(tfiles)],width=5,height=4,units='in',res=300)
  hist(rbeta(10000,(N-1)/2,(N-1)/2),
       xlim=c(0,1),
       n=100,
       main=sprintf('N=%d',N),
       xlab='Maternal %',
       ylab='',
       xaxt='n',
       yaxt='n',
       bty='n',
       col='black',
       border=NA)
  axis(1,at=seq(0,1,.2),labels=seq(0,100,20),xpd=NA)
  dev.off()
}
system(sprintf('convert -delay 10 -loop 0 %s %s',paste(tfiles,collapse=' '),file.path(plotDir,'popSize.gif')))
unlink(tfiles)
#What to the distribution as you change effective population size?
Ns = seq(4,32,length.out=10)
par(mfrow=c(1,length(Ns)))
par(mar=c(0,0,0,0))
for(N in Ns){
  hist(rbeta(10000,(N-1)/2,(N-1)/2),
       n=100,
       main='',
       xlab='Maternal %',
       ylab='',
       xaxt='n',
       yaxt='n',
       bty='n',
       col='black',
       border=NA)
}

#Make an animated gif of population size
N = 12
rhos = seq(0,1,.025)
rhos = c(rhos,rev(rhos))
tfiles=c()
for(rho in rhos){
  tfiles = c(tfiles,tempfile())
  png(tfiles[length(tfiles)],width=5,height=4,units='in',res=300)
  r = rbeta(10000,(N-1)/2,(N-1)/2)
  hist(rho*(runif(10000)>r)+(1-rho)*r,
       xlim=c(0,1),
       n=100,
       main=sprintf('rho=%.02f',rho),
       xlab='Maternal %',
       ylab='',
       xaxt='n',
       yaxt='n',
       bty='n',
       col='black',
       border=NA)
  axis(1,at=seq(0,1,.2),labels=seq(0,100,20),xpd=NA)
  dev.off()
}
system(sprintf('convert -delay 10 -loop 0 %s %s',paste(tfiles,collapse=' '),file.path(plotDir,'rhoShift.gif')))
unlink(tfiles)


#The same for clone at fixed N
N=12
rhos = seq(0,1,.1)
par(mfrow=c(1,length(rhos)))
par(mar=c(0,0,0,0))
for(rho in rhos){
  r = rbeta(10000,(N-1)/2,(N-1)/2)
  hist(rho*(runif(10000)>r)+(1-rho)*r,
       n=100,
       main='',
       xlab='',
       ylab='',
       xaxt='n',
       yaxt='n',
       bty='n',
       col='black',
       border=NA)
}
#Alternative in one plot
xs = seq(0,1,.001)
ys = lapply(rhos,function(rho) dbeta(xs/(1-rho),(N-1)/2,(N-1)/2)*(1-(xs/(1-rho))) + dbeta((xs-rho)/(1-rho),(N-1)/2,(N-1)/2)*((xs-rho)/(1-rho)))
plot(NA,
     type='n',
     xlim=c(0,1),
     ylim=c(0,max(unlist(ys),na.rm=TRUE)),
     xlab='',
     ylab='',
     main=11)
for(i in seq_along(ys))
  lines(xs,ys[[i]],lwd=1+i*1)


