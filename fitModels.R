#' Fit the EM model to each dataset to determine X genotype and which is inactivated.
#deterministicAnnealing.R

#############
# Libraries #
#############

library(inactiveXX)
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(ggExtra)
import('Code/scQC.R',as='scQC')

##########
# Params #
##########

outDir = 'Results/fitEM'
refGenomeDNA = 'Data/refGenomeDNA.fa'
refGenomeRNA = 'Data/refGenomeRNA.fa'
liftChain = 'Data/hg19ToHg38_noChr.over.chain'
X1kFile = 'Data/X1k.RDS'
gtf = 'Data/GRCh38.gtf'
regionsToUse=c('Exonic','Intronic')
#Using vartrix may be slightly more accurate due to remapping, but is much slower particularly for donors with many samples
useVartrix = FALSE

#############
# Functions #
#############


##############
# Fit models #
##############

###############
#### Prep SNPs
####Definition from https://academic.oup.com/gbe/article/5/10/1863/522116, which is hg19
strata = c(0,2.78,5.04,8.43,30.62,55.78,75.53,99.98,130.82,145.73,155.72)*1e6
#These are the boundaries lifted over and the PAR definitions from a hg38 source
strata = c(0,2.781479,5.12,8.46,30.60,55.75,76.31,100.73,131.69,146.65,155.701383)*1e6
strata = GRanges('X',IRanges(strata,width=c(diff(strata),1e9)))
names(strata) = paste0('strata',seq_along(strata))
names(strata)[1] = 'PAR1'
names(strata)[length(strata)]='PAR2'
counterLab=ifelse(useVartrix,'using_Vartrix','using_alleleCounter')
##########
# Process
for(nom in unique(unlist(areUsed))){
  #Is this one that needs to be split by demultiplexing?
  dat = scDat[scDat$donor==nom,]
  #dat  = dataSources[[nom]]
  if(any(dat$ID %in% demux$sampleID)){
    demuxSamples=sort(unique(demux$donorID[demux$sampleID %in% dat$ID]))
    toDemux=TRUE
  }else{
    demuxSamples=NA
    toDemux=FALSE
  }
  for(demuxSample in demuxSamples){
    #There are no samples where we care about the maternal component, so don't waste time on them...
    if(!is.na(demuxSample) && demuxSample=='M')
      next
    dnom = nom
    if(toDemux)
      dnom = paste0(dnom,'___',demuxSample)
    message(sprintf("Fitting model for sample %s",dnom))
    #It's now too much to save everything in one big file.  So save the model fit one at a time
    outFit = file.path(outDir,sprintf('%s_fit.RDS',dnom))
    doIt = !file.exists(outFit)
    #Don't process this one if it's already been fit
    if(doIt){
      out = list()
      #If the fit is not converged, set this
      out$warnFlag=FALSE
      ##########################
      # Basic Seurat processing
      #Construct a Seurat object that integrates this information
      message("Constructing Seurat object")
      cellsToUse = getCellBarcodes(setNames(dat$stagedPathCnts,dat$ID))
      if(toDemux)
        cellsToUse = cellsToUse[cellsToUse %in% gsub('-[0-9]+$','',demux$cellID[demux$donorID==demuxSample])]
      srat = setNames(dat$stagedPathCnts,dat$ID)
      srat = Read10X(srat)
      srat = srat[,cellsToUse]
      #Do clustering
      srat = scQC$quickCluster(srat)
      ############
      # Call SNPs
      snpList = list()
      bams10X = setNames(dat$stagedPathBAM,dat$ID)
      out$bams10X = bams10X
      #Always get it from scRNA
      message("Getting hSNPs from scRNA")
      for(ii in seq(2)){
        snpList[['FromRNA']] = tryCatch(
                                        {
                                          hetSNPsFromRNA(bams10X,refGenomeRNA,
                                                cellsToUse = cellsToUse,
                                                errRate=0.01,
                                                alpha=0.1,
                                                outputs = file.path(outDir,sprintf('%s_scRNA_1kSNPs_Xcnts_from_%s_%s.tsv',dnom,names(bams10X),counterLab)),
                                                nParallel=nParallel,
                                                skipIfExists=ii==1,
                                                nMaxSim=1000,
                                                useVartrix=useVartrix)
                                        },
                                        error=function(cond) {
                                          message("Fetch from file failed, recalculating from BAM file.") 
                                          return(NULL)
                                        })
        #Only proceed to the next if counts are the old type, with -N on the end, or some other failure.
        if(!(is.null(snpList[['FromRNA']]) || any(grepl('-[0-9]$',snpList$FromRNA@metadata$cnts$cellID))))
          break
      }
      #And also from DNA if available
      dna = datMap[datMap$donor==nom & grepl('DNA',datMap$dataType),]
      #If multiplexing, filter out unless appropriate
      if(!is.na(demuxSample) && demuxSample!='F')
        dna = dna[NULL,]
      if(any(dna$dataType=='normDNA')){
        x = dna[dna$dataType=='normDNA',]
        #Get from BAM
        message("Getting hSNPs from DNA")
        hSNPs = findHetSNPs(x$mappedDataPath,refGenomeDNA,
                            skipIfExists=TRUE,
                            outVCF = file.path(outDir,sprintf('%s_normBAM_hetSNPs.vcf',dnom)),
                            nParallel=nParallel)
        #Add in strata
        tmp = hSNPs
        o = findOverlaps(tmp,strata)
        tmp$strata = NA
        tmp$strata[queryHits(o)] = names(strata)[subjectHits(o)]
        tmp$strata = factor(tmp$strata,levels=names(strata))
        #Convert genome
        tmp = changeGenomeVersion(tmp,liftChain)
        tmp = annotateSNPs(tmp,gtf=gtf)
        tmp = tmp[as.character(seqnames(tmp))=='X']
        #o = findOverlaps(tmp,X1k@metadata$strata)
        #tmp$strata[queryHits(o)] = names(X1k@metadata$strata)[subjectHits(o)]
        #tmp$strata = factor(tmp$strata,levels=names(X1k@metadata$strata))
        snpList[['FromDNA']] = tmp
        #Phase if we can
        if(any(dna$dataType %in% c('mumDNA','dadDNA'))){
          message("Phasing hSNPs from DNA")
          #Fix up path to mum/dad DNA
          if('mumDNA' %in% dna$dataType){
            mum=dna$mappedDataPath[dna$dataType=='mumDNA']
          }else{
            mum=NULL
          }
          if('dadDNA' %in% dna$dataType){
            dad=dna$mappedDataPath[dna$dataType=='dadDNA']
          }else{
            dad=NULL
          }
          hSNPs = phaseSNPsFromParents(hSNPs,refGenomeDNA,mum,dad,
                                       outBase=file.path(outDir,paste0(dnom,'_phasing')),
                                       nParallel=nParallel)
          hSNPs = hSNPs[!is.na(hSNPs$altIsMum)]
          #Add in strata
          tmp = hSNPs
          o = findOverlaps(tmp,strata)
          tmp$strata = NA
          tmp$strata[queryHits(o)] = names(strata)[subjectHits(o)]
          tmp$strata = factor(tmp$strata,levels=names(strata))
          #Change genome
          tmp = changeGenomeVersion(hSNPs,liftChain)
          tmp = annotateSNPs(tmp,gtf=gtf)
          tmp = tmp[as.character(seqnames(tmp))=='X']
          #o = findOverlaps(tmp,X1k@metadata$strata)
          #tmp$strata[queryHits(o)] = names(X1k@metadata$strata)[subjectHits(o)]
          #tmp$strata = factor(tmp$strata,levels=names(X1k@metadata$strata))
          snpList[['FromDNA_phased']] = tmp
        }
      }
      out$snpList = snpList
      #############
      # Fit models
      #Do one fit per source of hSNPs.  That is scRNA (always), germline DNA (sometimes), parental DNA (rarely)
      out$XCntsList = list()
      out$fitList = list()
      srat@misc$gentypeList = list()
      for(snpType in names(snpList)){
        #Process using each source of SNPs in turn
        message(sprintf('Fitting model using het SNPs %s',snpType))
        #Check if we can just load the counts from the object
        if(!is.null(snpList[[snpType]]@metadata$cnts)){
          #If we're reusing the calculated counts from when we called het SNPs, we need to annotate
          XCnts = filterCountsX(snpList[[snpType]])
        }else{
          #Check for old version and repeat if needed
          for(ii in seq(2)){
            XCnts = tryCatch({
              #If we have het SNPs from DNA, we need to calculate the counts in the scRNA
              if(useVartrix){
                XCnts = vartrixCnts(snpList[[snpType]],bams10X,cellsToUse,refGenomeRNA,
                                           outputs=file.path(outDir,sprintf('%s_scRNA_X_cnts_from_%s_using_snps_%s_%s.tsv',dnom,names(bams10X),snpType,counterLab)),
                                           nParallel=nParallel)
              }else{
                XCnts = getAllelicExpression(snpList[[snpType]],refGenomeRNA,bams10X,
                                           outputs=file.path(outDir,sprintf('%s_scRNA_X_cnts_from_%s_using_snps_%s_%s.tsv',dnom,names(bams10X),snpType,counterLab)),
                                           skipIfExists=ii==1,
                                           nParallel=nParallel)
              }
              XCnts
            },
            error = function(cond){
              message("Fetching from file failed, going back to BAM file")
              return(NULL)
            }
            )
            #Stop if we've got the new one.  Otherwise have to go back to BAM
            if(!is.null(XCnts) && any(!grepl('-[0-9]+$',XCnts$cellID)))
              break
          }
          #Redundant in Vartrix mode, but needed otherwise
          XCnts = XCnts[XCnts$cellID %in% cellsToUse]
          XCnts = filterCountsX(XCnts)
        }
        out$XCntsList[[snpType]] = XCnts
        #Fit!  Put the failures in the sinBin and proceed.
        fit = inferInactiveX(XCnts,nParallel=nParallel,nStarts=1000,tauDiffWarnOnly=TRUE)
        if(fit$warnFlag){
          #This needs to still progress and save, but we need to retain some obvious marker it hasn't worked.
          out$warnFlag=TRUE
        }
        out$fitList[[snpType]] = fit
        out$fit = fit
        #Store
        #Add in all useful fit metadata do this twice, in an overwritey way that will ultimately store the highest quality results and in a per-approach way.
        #Haven't bothered working out which way around mat/pat should be here
        tmp = fit$cellSummary
        tmp = tmp[,c('matCount','patCount','tot','offCount','matFrac','badFrac','stateProbs','stateXi')]
        for(cNom in colnames(tmp))
          srat@meta.data[,cNom] = tmp[colnames(srat),cNom]
        srat@meta.data$nSNPs = as.numeric(table(fit$dd$cell)[colnames(srat)])
        srat@misc$genotype = fit$genotype
        srat@misc$tau = fit$tau
        #Now store in the more aware of origin way
        for(cNom in colnames(tmp))
          srat@meta.data[,paste0(cNom,'_',snpType)] = tmp[colnames(srat),cNom]
        srat@meta.data[,paste0('nSNPs_',snpType)] = as.numeric(table(fit$dd$cell)[colnames(srat)])
        srat@misc[[paste0('genotype_',snpType)]] = fit$genotype
        srat@misc[[paste0('tau_',snpType)]] = fit$tau
        ###############
        # Plot summary
        pdf(file.path(outDir,sprintf('%s_fit_%s.pdf',dnom,snpType)))
        {
          #New version
          #Solution overview.  ggMarginal draws blank page, can't be bothered fixing
          if(snpType!='FromDNA_phased'){
            plotSolutions(fit)
          }
          #Bad SNPs
          x = fit$snpSummary
          x = x[!is.na(x$offCount),]
          x$strata = X1k$strata[match(gsub('_.*','',x$SNP),as.character(X1k))]
          #Get the lowest badFrac consistent with 95% confidence interval
          x$badFracLow = sapply(seq(nrow(x)),function(e) if(x$tot[e]==0){NA}else{binom.test(x$offCount[e],x$totOff[e],0.5)$conf.int[1]})
          boxplot(split(x$badFracLow,x$strata),
                  main='Genomic distribution of bad SNPs',
                  ylab='badFracLow',
                  xlab='',
                  ylim=c(0,0.25)
          )
          #Plot the model (high hcoverage)
          plotModelHeatmap(fit,minCov=20)
          #All of them
          plotModelHeatmap(fit,minCov=0)
          #Marginal on cells
          plotModelHeatmap(fit,summariseBy='cell')
          #And SNPs
          plotModelHeatmap(fit,summariseBy='SNP')
        }
        dev.off()
      }
      out$srat = srat
      out$mDat = dat
      #########
      # Finish
      saveRDS(out,outFit)
    }
  }
}
