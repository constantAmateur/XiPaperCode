#' Plots validating the accuracy of the method

#############
# Libraries #
#############

library(inactiveXX)

##########
# Params #
##########

srcDir = 'Results/fitEM'
goldData = c('Wilms2','Wilms3','GOSH025')
#These are all of them
cancerData = c('GOSH016',"GOSH023",  "GOSH026", "GOSH028", "GOSH029", "GOSH032", "GOSH033", "IMP672",  "RCC2", "VHL_RCC")
#But only a handful are not shit.  Filter by those that 
cancerData = c("GOSH028", "GOSH029", "GOSH032", "IMP672", "VHL_RCC", "RCC2")
#Threshold for calling something as definitively in one state or the other
pCut=0.95
#What sample to do the downsampling experiment on
downsampleTgt = 'GOSH025'
#How many downsampling reps to do
nRepsDownsample = 5
#Downsampling target cell numbers.  Values above those given ignored.
downCells = c(100,200,500,1000,2000)
#What ratios to test
testRatios = c(0.01,0.02,0.05,0.1,0.2,0.5)
refGenomeRNA = 'Data/refGenomeRNA.fa'
outDir = 'Results/fitEM'
plotDir = 'Results/plots'
stateCols = c('#902D41','#133C55')
#Jitter.  As fraction of x,y range
jitter=c(.01,.01)

#############
# Functions #
#############

#' Plot solution distribution
#'
#' Given a fit object, with multiple random starts, plot the distribution of population skew in X-inactivation (tau) versus log-likelihood (Q).  Add some jitter to this so that overlapping points don't get too obscured.
#'
#' @param f The fit object.
#' @param jitter Fraction of plot range to jitter by in x,y direction.
#' @return Nothing, makes a plot.
plotSolutionDistribution = function(f,jitter=c(.01,.01)){
  taus = sapply(f$allFits,function(e) e$tau)
  Qs = sapply(f$allFits,function(e) e$Q)
  d = data.frame(Q=Qs,tau=taus)
  plot(x=taus+rnorm(length(taus),0,sd=abs(diff(range(taus)))*jitter[1]),
       y=Qs+rnorm(length(Qs),0,sd=abs(diff(range(Qs)))*jitter[2]),
       xlim=c(0,1),
       xlab='Maternal %',
       ylab='Log likelihood',
       main='Solution distribution',
       pch=19,
       cex=0.1)
  #gg = ggplot(d,aes(tau,Q)) +
  #  geom_hex(bins=10) +
  #  geom_point(size=1.0,color='lightgrey',alpha=0.1) + 
  #  theme_bw() +
  #  xlim(0,1) + 
  #  xlab('Fractional maternal') +
  #  ylab('Log likelihood') +
  #  ggtitle('Solution distribution')
  #return(gg)
}





###########################
# Load gold standard fits #
###########################

outs = list()
for(src in goldData){
  message('\n\n\n')
  message(sprintf("Stats for case %s",src))
  #Load the fit object
  dat = readRDS(file.path(srcDir,paste0(src,'_fit.RDS')))
  coal = dat$fitList$FromRNA
  gold = dat$fitList$FromDNA_phased
  ##################
  # Cell comparison
  comp = merge(coal$cellSummary,gold$cellSummary,
               all=TRUE,
               by='cell',
               suffixes=c('Coal','Gold'))
  #Keep only cells with data in both
  compFull = comp[!is.na(comp$matCountGold) & !is.na(comp$matCountCoal),] 
  compFullHighQual = compFull[(compFull$highConfCallGold & compFull$highConfCallCoal),]
  cells = list(all=comp,full=compFull,qual=compFullHighQual)
  message('\n')
  message('Cells')
  message('All')
  print(table(gold=compFull$statesGold>0,coal=compFull$statesCoal>0))
  message(sprintf('Bad frac %f',sum((compFull$statesGold>0) != (compFull$statesCoal>0))/nrow(compFull)))
  message('High quality')
  print(table(gold=compFullHighQual$statesGold>0,coal=compFullHighQual$statesCoal>0))
  message(sprintf('Bad frac %f',sum((compFullHighQual$statesGold>0) != (compFullHighQual$statesCoal>0))/nrow(compFull)))
  #################
  # SNP comparison
  gold$snpSummary$genotype = gold$genotype[gold$snpSummary$SNP]
  coal$snpSummary$genotype = coal$genotype[coal$snpSummary$SNP]
  gold$snpSummary$SNP = toupper(gold$snpSummary$SNP)
  comp = merge(coal$snpSummary,gold$snpSummary,
               all=TRUE,
               by='SNP',
               suffixes=c('Coal','Gold'))
  #Keep only cells with data in both
  compFull = comp[!is.na(comp$matCountGold) & !is.na(comp$matCountCoal),] 
  compFullHighQual = compFull[(compFull$highConfCallGold & compFull$highConfCallCoal),]
  snps = list(all=comp,full=compFull,qual=compFullHighQual)
  message('\n')
  message('SNPs')
  message('All')
  print(table(gold=compFull$genotypeGold,coal=compFull$genotypeCoal))
  message(sprintf('Bad frac %f',sum(compFull$genotypeCoal != compFull$genotypeGold)/nrow(compFull)))
  message('High quality')
  print(table(gold=compFullHighQual$genotypeGold,coal=compFullHighQual$genotypeCoal))
  message(sprintf('Bad frac %f',sum(compFullHighQual$genotypeCoal != compFullHighQual$genotypeGold)/nrow(compFull)))
  outs[[src]]=list(cell=cells,snp=snps,srat=dat$srat)
}

######################
# Example tau-Q plot #
######################

pdf(file.path(plotDir,'solnDistrn.pdf'),width=3,height=4)
plotSolutionDistribution(coal)
dev.off()


#################################
# Summarise predictive accuracy #
#################################

pdf(file.path(plotDir,'validationAccuracy.pdf'),width=6,height=4)
par(mfcol=c(2,length(goldData)))
cols=list(correct='#4D4D4D',
          error='#AEAEAF',
          uncalled='#E6E6E6')
for(src in goldData){
  #Margins for top row
  par(mar=c(1,5,3,0))
  x = outs[[src]]$cell$full
  #Keep only the ones we're really sure of
  x = x[x$highConfCallGold,]
  #Order by the quality of the coal call
  x = x[order(-abs(x$stateLogitsCoal)),]
  #Construct the ROC curve values with end-points
  d = data.frame(TP = cumsum((x$stateProbsCoal>0.5)==(x$stateProbsGold>0.5)),
                 FP = cumsum((x$stateProbsCoal>0.5)!=(x$stateProbsGold>0.5)),
                 Uncalled = seq(nrow(x),1,-1),
                 stateLogit = x$stateLogitsCoal)
  barplot(t(d[,1:3]),
          main=src,
          las=2,
          col=unlist(cols[c('correct','error','uncalled')]),
          border=NA,
          space=0,
          ylab='# Cells')
  w = max(which(abs(d$stateLogit)>=abs(log(1-pCut)-log(pCut))))
  abline(v=w,col='red')
  #Work out the uncalled fraction and error rate at this cut-off
  errRate = c(fracUncalled = d$Uncalled[w]/nrow(d),
              fracErrors = d$FP[w]/(d$FP[w]+d$TP[w]),
              fracCorrect = d$TP[w]/(d$FP[w]+d$TP[w]))
  #And for the full allocation
  errRateFull = c(fracUncalled=0,
                  fracErrors = d$FP[nrow(d)]/nrow(d),
                  fracCorrect = d$TP[nrow(d)]/nrow(d))
  text(x=1,
       y=nrow(d)*1,
       adj=c(0,1),
       xpd=NA,
       sprintf('State prob>%g: %d %% uncalled, %d %% correct',pCut,round(100*errRate['fracUncalled']),round(100*errRate['fracCorrect'])))
  text(x=1,
       y=nrow(d)*.8,
       adj=c(0,1),
       xpd=NA,
       sprintf('All allocated: %d %% uncalled, %d %% correct',round(100*errRateFull['fracUncalled']),round(100*errRateFull['fracCorrect'])))
  #And for the bottom row
  par(mar=c(5,5,0,0))
  #And now the same, but for SNPs
  x = outs[[src]]$snp$full
  #Keep only the ones we're really sure of
  x = x[x$highConfCallGold,]
  #Order by the quality of the coal call
  x = x[order(x$totCoal),]
  #Split by 2-5, 5+
  tmp = split(x,cut(x$totCoal,breaks = c(0,5,10,20,Inf)))
  d = data.frame(coverage = names(tmp),
                 nSNPs = sapply(tmp,nrow),
                 nCorrect = sapply(tmp,function(e) sum(e$genotypeCoal==e$genotypeGold)),
                 nError = sapply(tmp,function(e) sum(e$genotypeCoal!=e$genotypeGold)))
  barplot(t(d[,3:4]/d$nSNPs),
          col = unlist(cols[c('correct','error')]),
          border=NA,
          las=1,
          ylab='Fraction of  SNPs',
          xlab='Total SNP coverage')
}
dev.off()

##############
# Downsample #
##############

#Gosh025 downsample stats
dat = readRDS(file.path(srcDir,paste0(downsampleTgt,'_fit.RDS')))
gold = dat$fitList$FromDNA_phased
gold$snpSummary$genotype = gold$genotype[gold$snpSummary$SNP]
gold$snpSummary$SNP = toupper(gold$snpSummary$SNP)
cellNames = gold$cellSummary$cell[(gold$cellSummary$highConfCall)]
downStats = data.frame(nCells = rep(downCells,each=nRepsDownsample),
                       nCorrectCells = NA,
                       nFailedCells = NA,
                       nUncalledCells = NA,
                       nCorrectSNPs = NA,
                       nFailedSNPs = NA,
                       repNumber = rep(seq(nRepsDownsample),length(downCells)))
ii=1
dsamps = list()
for(nCells in downCells){
  message(sprintf('-----------------------------\n- Downsampling to %d cells -\n-----------------------------',nCells))
  dsamps[[nCells]] = list()
  for(i in seq(nRepsDownsample)){
    message(sprintf('------------------------------\n- Rep %d of %d',i,nRepsDownsample))
    #If we're asking for more cells than exist, do one final report and finish
    finalRun = nCells>=length(cellNames)
    if(finalRun){
      fit = dat$fitList$FromRNA
    }else{
      #Get random set of cells to use
      cellsToUse = sample(cellNames,nCells)
      #Call SNPs
      message('Calling SNPs')
      snps = hetSNPsFromRNA(dat$bams10X,refGenomeRNA,
                            cellsToUse = cellsToUse,
                            errRate=0.01,
                            alpha=0.1,
                            outputs = file.path(outDir,sprintf('DOWNSAMPLE_%s_scRNA_1kSNPs_XCnts_from_%s_rep_%d',downsampleTgt,names(dat$bams10X),i)),
                            nParallel=nParallel,
                            skipIfExists=FALSE,
                            nMaxSim=1000,
                            useVartrix=FALSE)
      #Get X counts
      XCnts = filterCountsX(snps)
      #Do the fit
      message('Fitting model')
      fit = inferInactiveX(XCnts,nParallel=nParallel,nStarts=1000,verbose=1,tauDiffWarnOnly=TRUE)
      dsamps[[nCells]][[i]] = fit
      finalRun=FALSE
    }
    #Match with gold standard stuff
    coal = fit
    cells = merge(coal$cellSummary,gold$cellSummary,
                  all=TRUE,
                  by='cell',
                  suffixes=c('Coal','Gold'))
    #Keep only those that are complete in both
    cells = cells[!is.na(cells$matCountGold) & !is.na(cells$matCountCoal),]
    #And they should all (by construction) be high quality in gold, but check this
    if(!finalRun && !all(cells$highConfCallGold))
      stop('Sanity check fail')
    downStats$nUncalledCells[ii] =  sum(abs(cells$stateLogitsCoal)<abs(log(1-pCut)-log(pCut)))
    downStats$nCorrectCells[ii] =  sum(abs(cells$stateLogitsCoal)>=abs(log(1-pCut)-log(pCut)) & (cells$stateLogitsCoal>0)==(cells$stateLogitsGold>0))
    downStats$nFailedCells[ii] =  sum(abs(cells$stateLogitsCoal)>=abs(log(1-pCut)-log(pCut)) & (cells$stateLogitsCoal>0)!=(cells$stateLogitsGold>0))
    #Same for snps
    coal$snpSummary$genotype = coal$genotype[coal$snpSummary$SNP]
    snps = merge(coal$snpSummary,gold$snpSummary,
                 all=TRUE,
                 by='SNP',
                 suffixes=c('Coal','Gold'))
    #Keep only cells with data in both
    snps = snps[!is.na(snps$matCountGold) & !is.na(snps$matCountCoal),] 
    #And those that we are confident we know accurately.
    snps = snps[snps$highConfCallGold,]
    #Get pCut accuracy of cells for 
    downStats$nCorrectSNPs[ii] = sum(snps$genotypeCoal==snps$genotypeGold)
    downStats$nFailedSNPs[ii] = sum(snps$genotypeCoal!=snps$genotypeGold)
    print(downStats[ii,,drop=FALSE])
    if(finalRun){
      downStats$nCells[ii]=length(cellNames)
      downStats = downStats[seq(ii),]
      break
    }
    ii=ii+1
  }
}
#Where nFailed>nCorrect, can (and should) swap definitions of genotype to get better performance.
w = downStats$nCorrectSNPs < downStats$nFailedSNPs
tmp = downStats
downStats$nCorrectCells[w] = tmp$nFailedCells[w]
downStats$nFailedCells[w] = tmp$nCorrectCells[w]
downStats$nCorrectSNPs[w] = tmp$nFailedSNPs[w]
downStats$nFailedSNPs[w] = tmp$nCorrectSNPs[w]
#Drop the final row
d = downStats[-nrow(downStats),]
d$cellErrRate = d$nFailedCells/(d$nFailedCells+d$nCorrectCells)
#Make two plots to show that the cell error rate goes down as the number of SNPs goes up
pdf(file.path(plotDir,'validationDownsample.pdf'),width=4,height=5)
layout(matrix(c(1,2),nrow=2),heights=c(1,1))
par(mar=c(0.5,5,2,1))
#First plot the error rate
nCells = unique(d$nCells)
plot(NA,
   xlim=c(0.5,length(nCells)+0.5),
   ylim=c(0.001,max(d$cellErrRate)),
   log='y',
   ylab='error rate',
   las=1,
   xaxt='n',
   type='n')
points(x = rep(seq_along(nCells),each=nRepsDownsample)+rnorm(nrow(d),sd=0.1),
     y = d$cellErrRate,
     pch=19,
     cex=0.5)
for(i in seq_along(nCells)){
m = median(d$cellErrRate[d$nCells==nCells[i]])
lines(i + c(-0.2,.2),c(m,m),lwd=2)
}
#axis(1,at=seq_along(nCells),labels=nCells)
#Then below do the number of detected SNPs
par(mar=c(3,5,0,1))
plot(NA,
   xlim=c(0.5,length(nCells)+0.5),
   ylim=c(1,max(d$nCorrectSNPs)),
   ylab='# genotyped SNPs',
   las=1,
   xaxt='n',
   type='n')
points(x = rep(seq_along(nCells),each=nRepsDownsample)+rnorm(nrow(d),sd=0.1),
     y = d$nCorrectSNPs,
     pch=19,
     cex=0.5)
for(i in seq_along(nCells)){
m = median(d$nCorrectSNPs[d$nCells==nCells[i]])
lines(i + c(-0.2,.2),c(m,m),lwd=2)
}
axis(1,at=seq_along(nCells),labels=nCells)
title(xlab='# Cells',line=2)
dev.off()

##############
# Ratio test #
##############

#GOSH025 extreme ratio tests
dat = readRDS(file.path(srcDir,paste0(downsampleTgt,'_fit.RDS')))
gold = dat$fitList$FromDNA_phased
gold$snpSummary$genotype = gold$genotype[gold$snpSummary$SNP]
gold$snpSummary$SNP = toupper(gold$snpSummary$SNP)
cellNames = gold$cellSummary$cell[(gold$cellSummary$highConfCall)]
#Split by state
cellNames = split(cellNames,gold$cellSummary[cellNames,'stateLogits']>0)
cellNames = cellNames[order(lengths(cellNames))]
nStates = sort(as.numeric(table(gold$cellSummary$stateLogits[gold$cellSummary$highConfCall]>0)))
#Target number of each cell type for all ratios
nTgtCells = list()
for(testRatio in testRatios){
  #Decrease the lower one
  if(testRatio < nStates[1]/sum(nStates)){
    nTgtCells[[length(nTgtCells)+1]] = c(floor(nStates[2]*testRatio/(1-testRatio)),nStates[2])
  }else{
    nTgtCells[[length(nTgtCells)+1]] = c(nStates[1],ceiling(nStates[1]*(1-testRatio)/testRatio))
  }
}
ii=1
rsamps = list()
for(i in seq_along(testRatios)){
  testRatio = testRatios[i]
  message(sprintf('-----------------------------\n- Testing ratio  %g -\n-----------------------------',testRatio))
  rsamps[[sprintf('%.02f',testRatio)]] = list()
  #Which index to downsample
  sampTgt = which(nTgtCells[[i]]<nStates)
  for(j in seq(nRepsDownsample)){
    message(sprintf('------------------------------\n- Rep %d of %d',j,nRepsDownsample))
    #Downsample the required one then merge
    cellsToUse = cellNames
    cellsToUse[[sampTgt]] = sample(cellsToUse[[sampTgt]],nTgtCells[[i]][sampTgt])
    cellsToUse = unlist(cellsToUse,use.names=FALSE)
    #Call SNPs
    message('Calling SNPs')
    snps = hetSNPsFromRNA(dat$bams10X,refGenomeRNA,
                          cellsToUse = cellsToUse,
                          errRate=0.01,
                          alpha=0.1,
                          outputs = file.path(outDir,sprintf('RATIOTEST_%s_scRNA_1kSNPs_XCnts_from_%s_rep_%d',sprintf('%.02f',testRatio),names(dat$bams10X),j)),
                          nParallel=nParallel,
                          skipIfExists=FALSE, #As we haven't set a seed, the random samples will be random and we need to redo this every time...
                          nMaxSim=1000,
                          useVartrix=FALSE)
    #Get X counts
    XCnts = filterCountsX(snps)
    #Do the fit
    message('Fitting model')
    fit = inferInactiveX(XCnts,nParallel=nParallel,nStarts=1000,verbose=1,tauDiffWarnOnly=TRUE)
    taus = sapply(fit$allFits,function(e) e$tau)
    Qs = sapply(fit$allFits,function(e) e$Q)
    message(sprintf('Best fit tau=%g',fit$tau))
    plot(taus,Qs,pch=19)
    rsamps[[i]][[j]] = fit
  }
}
ratioStats = data.frame(testRatio = rep(testRatios,each=nRepsDownsample),
                        repNo = rep(seq(nRepsDownsample),length(testRatios)),
                        nCellsSmall = rep(sapply(nTgtCells,`[`,1),each=nRepsDownsample),
                        nCellsLarge = rep(sapply(nTgtCells,`[`,2),each=nRepsDownsample),
                        tau = NA,
                        tauDiff = NA)
ii=1
for(i in seq_along(rsamps)){
  #par(mfcol=c(2,3))
  message(sprintf('Ratio %s',names(rsamps)[i]))
  for(j in seq_along(rsamps[[i]])){
    taus = sapply(rsamps[[i]][[j]]$allFits,function(e) e$tau)
    Qs = sapply(rsamps[[i]][[j]]$allFits,function(e) e$Q)
    o = order(Qs)
    taus = taus[o]
    Qs = Qs[o]
    tau = tail(taus,n=1)
    ratioStats$tau[ii]=tau
    #WLOG make tau <0.5
    if(tau>0.5)
      tau = 1-tau
    isExtreme = taus>=tau | taus<=(1-tau)
    w = tail(seq_along(taus),n=100)
    tauDiff = pmin(1-taus,taus)-tau
    ratioStats$tauDiff[ii] = mean(tauDiff[w])/tau
    #plot(tauDiff[w],Qs[w],pch=19,
    #     main=sprintf('ratio %s,tau=%.02f, test %d',names(rsamps)[i],tau,j))
    #Average of the difference between top 100 fits and best fit
    message(sprintf('Average tau diff is %g',mean(tauDiff[w])/tau))
    #plot(taus,Qs,
    #     pch=19,
    #     main=sprintf('ratio %s,tau=%.02f, test %d',names(rsamps)[i],tau,j),
    #     col=ifelse(isExtreme,'red','black')
    #     )
    #print(table(tail(isExtreme,n=300)))
    ii=ii+1
  }
}
#Make a plot summarising things
pdf(file.path(plotDir,'validationRatioTest.pdf'),width=4,height=4)
par(mar=c(5,5,1,1))
with(ratioStats,
     plot(testRatio,pmin(tau,1-tau),
          log='xy',
          col = 'black',#circlize::colorRamp2(c(0,.2),c('black','red'))(abs(tauDiff)>0.05),
          xlim=c(0.01,0.5),
          ylim=c(0.01,0.5),
          pch=ifelse(tauDiff>0.05,3,19),
          cex=0.6,
          xlab='True ratio',
          ylab='Inferred ratio'
          )
     )
abline(0,1,col='black',lty=2)
dev.off()
#Show an example of a failure to converge on a value of tau across samples
pdf(file.path(plotDir,'validationRatioTestBadExample.pdf'),width=3,height=4)
plotSolutionDistribution(rsamps$`0.01`[[1]])
dev.off()
#And one well converged
pdf(file.path(plotDir,'validationRatioTestGoodExample.pdf'),width=3,height=4)
plotSolutionDistribution(rsamps$`0.10`[[1]])
dev.off()

#############################
# Cancer samples validation #
#############################

#Samples to use
#RCC2 - Good number of clear cancer cells (~300)
#VHL_RCC - Only ~50 tumour cells, but clearly definable.
#IMP672 - Wilms. The best one.  Thousands of tumour cells in two blobs.
#GOSH032 - Ewings. Not many cells here, but the data is clear.  Also a good example of near clonal endothelium, which may or may not be desirable
#GOSH029 - MRT.  Dominated by the tumour cells (thousands of them), but extremely clean signal.
#GOSH028 - RMS, so can only use once Nathan has published.  Would be usable but nothing really amazing to add that the others don't.  I guess the main benefit would be that there are populations of cells of unclear origin around the periphery of the tumour blob.  This shows that while the main tumour blob is extremely homogenous, the peripheral populations must have at least some non-tumour cells.
#Rough annotation 
annots=list()
annots[['GOSH028']] = list(NK = c(19, 16), T = 16, MNP = 17, RBC = 21:22, Unknown = 15, 
    Tum1 = c(0, 10, 8, 6, 7, 3, 5, 4, 14, 23, 1, 2, 18, 11), TumUnknown = c(20, 
    9), TumUnknown2 = 13, TumUnknown3 = 12)
annots[['GOSH029']] = list(T = c(27, 19, 20), MNP = c(7, 8, 18, 25), RBC = 22, Tum1 = c(0, 
2, 5, 10, 15, 13, 3, 9, 1, 4, 6, 11, 16, 29), U17 = 17, U26 = 26, 
    U28 = 28, U14 = 14, U12 = 12, U21 = 21, NK = 23)
annots[['GOSH032']] = list(T=3,MNP = c(6,8), Endo = 7,RBC=c(0,5), Tum1 = c( 2, 1, 4, 9,10)) 
annots[['IMP672']] = list(T = 14, MNP = 11, U12 = 12, U15 = 15, U14 = 14, Tum1 = c(2, 
4, 5), Tum2 = c(8, 7, 9, 13, 6, 3, 10, 0, 1))
annots[['RCC2']] = list(T = c(23, 9, 13, 2, 1, 26), NK = c(17, 15, 4), MNP = c(7, 
10, 11), B = 27, Mast = 30, RBC = c(14, 18), Tum1 = 12, TumAdj = c(8, 
24, 20), Endo = c(6, 22), Fibro = 19, Pelvis = c(16, 25), Plasma = 21, 
    Tub1 = 31, Tub2 = 32, Unknown = 29, Tub = c(0, 3, 5, 28), 
    `Stem?` = 25)
annots[['VHL_RCC']] = list(NK = c(11, 12), T = c(3, 7, 6), MNP = 8, Endo = c(18, 10, 
13), Mast = 21, Pelvis = 20, Tum1 = 15, B = 22, Fibro = 17, Unknown1 = 19, 
    Prolif = 16, Tub = c(5, 9, 2, 0, 4, 1, 14))
outs = list()
for(src in cancerData){
  message('\n\n\n')
  message(sprintf("Stats for case %s",src))
  #Load the fit object
  dat = readRDS(file.path(srcDir,paste0(src,'_fit.RDS')))
  srat = dat$srat
  if(!is.null(annots[[src]])){
    srat@meta.data$annot = setNames(rep(names(annots[[src]]),lengths(annots[[src]])),unlist(annots[[src]]))[as.character(srat@active.ident)]
  }else{
    srat@meta.data$annot = as.character(srat@active.ident)
  }
  #Save it out
  saveRDS(srat,file.path(srcDir,paste0(src,'_fit_withAnnot.RDS')))
  cells = dat$fitList$FromRNA$cellSummary
  x = cells[cells$highConfCall,]
  x$annot = srat@meta.data[x$cell,'annot']
  tmp = split(factor(ifelse(x$stateLogits>=0,'M','P'),levels=c('M','P')),srat@meta.data[x$cell,'annot'])
  tmp = lapply(tmp,function(e) table(e)/length(e))
  tmp = do.call(cbind,tmp)
  barplot(tmp)
  table(x$annot)
  outs[[src]] = list(cells=cells,mDat=srat@meta.data,x=x)
}
noms = unique(unlist(lapply(outs,function(e) e$x$annot)))
noms = noms[!is.na(noms)]
#keepNoms = c('MNP','NK','T','Endo','Tum1','Tum2','B','Fibro','Tub','Mast')
keepNoms = c('Tum1','Tum2','Endo','Fibro','Tub','Mast','MNP','NK','T','B')
#NOTE: The expansion factor for last row needs to be tweaked to ensure a consistent height
pdf(file.path(plotDir,'validatingCancer.pdf'),width=3,height=4)
layout(matrix(seq_along(outs),ncol=1),
       heights=ifelse(seq_along(outs)==length(outs),2,1))
plotAll=FALSE
aoeu=list()
for(i in seq_along(outs)){
  if(i==length(outs)){
    par(mar=c(3,5,.2,1))
  }else{
    par(mar=c(.2,5,.2,1))
  }
  #No cut-offs, just sum probabilities
  if(plotAll){
    x = outs[[i]]$cells
    x$annot = outs[[i]]$mDat[x$cell,'annot']
    x$mProb = (1+exp(-x$stateLogits))**-1
    x$pProb = (1+exp(x$stateLogits))**-1
    x = aggregate(cbind(mProb,pProb) ~ annot,FUN=sum,data=x)
    x = x[match(keepNoms,x$annot),]
    x = t(as.matrix(x[,2:3]))
    x = t(t(x)/colSums(x))
    colnames(x) = keepNoms
    tmp = x 
  }else{
    #Just look at high-confidence cells
    x = outs[[i]]$x
    tmp = split(factor(ifelse(x$stateLogits>=0,'M','P'),levels=c('M','P')),factor(x$annot,levels=keepNoms))
    tmp = lapply(tmp,function(e) table(e)/length(e))
    tmp = do.call(cbind,tmp)
    aoeu[[i]]=tmp
  }
  barplot(tmp,
          col=stateCols,
          xaxt=ifelse(i==length(outs),'s','n'),
          yaxt='n',
          xlab=ifelse(i==length(outs),'cellType',''),
          ylab=paste0(names(outs)[i],'\nMat Frac')
          )
  axis(2,at=c(0,1),labels=c(0,1),las=1)
}
dev.off()

