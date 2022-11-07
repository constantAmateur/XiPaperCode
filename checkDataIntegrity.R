#' Check integrity of data
#'
#' Do some basic checks that the data is consistent and as we expect.  Will ultimately replace most of the pre-processing in fitModel.
#'
#' External variables:
#' @param dataSources/scDat Defined by running pipeline.

#############
# Libraries #
#############

library(alleleIntegrator)
library(Seurat)

##########
# Params #
##########

outBase = 'Results'
liftChain = 'Data/hg19ToHg38_noChr.over.chain'
xistCutOff=0.1
rps4y1CutOff=0.05

#############
# Functions #
#############


##################
# Check all data #
##################

for(nom in names(dataSources)){
  message('--------------------------------------------------')
  message(sprintf("Checking processed data for sample %s",nom))
  dat = dataSources[[nom]]
  dat$checksPassed=TRUE
  mDat = scDat[scDat$ID %in% dat$IDs,,drop=FALSE]
  if(nrow(mDat)==0){
    dat$checksPassed=FALSE
    warning(sprintf("Could not find scDat entry for sample %s",nom))
    dataSources[[nom]]=dat
    next
  }
  bams10X = setNames(mDat$stagedPathBAM,mDat$ID)
  mtx10X = setNames(file.path(mDat$stagedPathCnts,'barcodes.tsv.gz'),mDat$ID)
  #Do we have bams as we expect?
  #bams10X = file.path(dat$lustrePathByID,'possorted_genome_bam.bam')
  #w = file.exists(bams10X)
  #alt = file.path(dat$lustrePathByID,'outs','possorted_genome_bam.bam')
  #bams10X[!w] = alt[!w]
  #names(bams10X) = names(dat$lustrePathByID)
  if(!all(file.exists(bams10X))){
    dat$checksPassed=FALSE
    warning(sprintf("BAM files missing in sample %s",nom))
    dataSources[[nom]]=dat
    next
  }
  if(length(bams10X)==0){
    dat$checksPassed=FALSE
    warning(sprintf("No BAM files found for sample %s",nom))
    next
  }
  dat$bams10X = bams10X
  #And filtered counts?
  #mtx10X = file.path(dirname(bams10X),'filtered_gene_bc_matrices/GRCh38/barcodes.tsv')
  #w = file.exists(mtx10X)
  #alt = file.path(dirname(bams10X),'filtered_feature_bc_matrix/barcodes.tsv.gz')
  #mtx10X[!w] = alt[!w]
  #names(mtx10X) = names(bams10X)
  dat$mtx10X = mtx10X
  if(!all(file.exists(mtx10X))){
    dat$checksPassed=FALSE
    warning(sprintf("Counts missing in sample %s",nom))
    dataSources[[nom]]=dat
    next
  }
  #Load the barcodes that we will use
  cellsToUse=list()
  for(i in seq_along(mtx10X)){
    cellsToUse[[i]] = paste0(names(mtx10X)[i],'_',read.table(mtx10X[i],sep='\t',header=FALSE)[,1])
  }
  cellsToUse = unlist(cellsToUse,use.names=FALSE)
  dat$cellsToUse = cellsToUse
  #Check counts can be loaded
  srat = Read10X(setNames(dirname(mtx10X),names(mtx10X)))
  ############
  # Check sex
  gns=c('XIST','RPS4Y1')
  sexGns = srat[gns,]
  sexGns = do.call(cbind,lapply(split(colnames(srat),gsub('_[ACGT]+(-[0-9])?$','',colnames(srat))),function(e) rowMeans(sexGns[,e,drop=FALSE]>0)))
  dat$sexGeneExpression = sexGns
  #We should just base it on the absence of male markers
  #dat$isFemale = sexGns['XIST',]>xistCutOff & sexGns['RPS4Y1',]<rps4y1CutOff
  dat$isFemale = sexGns['RPS4Y1',]<rps4y1CutOff
  if(!all(dat$isFemale)){
    warning(sprintf("Not all samples in %s look female",nom))
    print(dat$sexGeneExpression)
    dat$checksPassed=FALSE
  }
  ###########################
  # Check barcode sharedness
  bcodes = colnames(srat)
  srcs = gsub('_[ACGT]+(-[0-9])?$','',bcodes)
  bcodes = gsub('.*_([ACGT]+(-[0-9])?)$','\\1',bcodes)
  #Report a broad summary
  olap = split(bcodes,srcs)
  message(sprintf("Across %d mappings we found an average of %g cells (5th and 95th quantiles %g and %g, range %g to %g)",length(bams10X),mean(lengths(olap)),quantile(lengths(olap),0.05),quantile(lengths(olap),0.95),min(lengths(olap)),max(lengths(olap))))
  #This magic is from here https://stackoverflow.com/questions/24614391/intersect-all-possible-combinations-of-list-elements
  olap = crossprod(table(stack(olap)))
  #Convert to fractional
  cnts = diag(olap)
  frac = t(olap/cnts)
  #Make a flat version of the frac
  dd = frac
  dd[lower.tri(dd)]=NA
  diag(dd)=NA
  dd = na.omit(data.frame(as.table(dd)))
  colnames(dd)=c('A','B','fracOlap')
  dd$cntsA = cnts[dd$A]
  dd$cntsB = cnts[dd$B]
  dd$cntsAB = dd$fracOlap*dd$cntsB
  #Now need a criteria for what is too much
  nPool=737000 #This should really be 3million for newer chemistries, but if I'm going to just use one criteria better to be conservative and use the smaller pool
  dd$pVal = phyper(dd$cntsAB-1,dd$cntsA,nPool-dd$cntsA,dd$cntsB,lower.tail = FALSE)
  dd$qVal = p.adjust(dd$pVal,method='BH')
  dd = dd[order(dd$pVal),]
  dat$barcodeComparison = dd
  if(any(dd$qVal<0.05)){
    warning("Significant barcode overlap detected.  Likely barcode swapping or sample mixup.")
    print(dd)
    dat$checksPassed=FALSE
  }
  #############
  # Genotyping
  #This be slooooooow, but only have to do it once and can be compared across everything
  pdf(file.path(outBase,'genotypes',sprintf('genotypeCheck_%s.pdf',nom)))
  genoCheck = matchBAMs(BAMs = bams10X,
                        refGenomes = file.path(refGenome,'fasta','genome.fa'), #refGenome should be predefined from prepData.R
                        outputs = file.path(outBase,'genotypes',sprintf('genotypeCheck_%s___%s.tsv',nom,gsub('-','_',names(bams10X)))),
                        liftOvers = liftChain,
                        is10X=TRUE,
                        nMaxSim=6,
                        nChunks=8,
                        nParallel=nParallel)
  dev.off()
  #Convert genotype consistency to useful format
  ibs = genoCheck$ibs$ibs
  ibs[lower.tri(ibs)]=NA
  diag(ibs)=NA
  ibs = na.omit(data.frame(as.table(ibs)))
  colnames(ibs) = c('A','B','IBS')
  dat$genoCheck = ibs
  dat$genoCheckFull = genoCheck
  if(dat$checksPassed){
    message("All checks passed.")
  }else{
    message("Some data checks failed for this sample.")
  }
  #Make exceptions on a case by case basis
  passAnyway = c(GOSH028="Needs remapping combining the two sequencing runs in one.  But have used data in this state elsehere without issue",
                 GOSH029="Needs remapping combining the two sequencing runs in one.  But have used data in this state elsehere without issue",
                 IMP672= "Needs remapping combining the two sequencing runs in one.  But have used data in this state elsehere without issue",
                 GOSH016="Technically fails, but trivial numbers shared.",
                 InfALL1="Multiplexed.  Trivial numbers shared.",
                 CongAML="Trivial numbers shared despite reaching significance.",
                 F17="Male fetus, keeping just for the one decidua channel.",
                 F20="Something is really fucked with this, as 100k cells called in one channel.  But fails due to overlapping barcodes which don't seem a big concern.",
                 F25="Male fetus, keeping for decidua channel.",
                 F27="Male fetus, keeping for decidua channel.",
                 F36="Male fetus, keeping for decidua channel.",
                 F39="Male fetus, keeping for decidua channel.",
                 FR6="Signficiant barcode overlap, but numbers are trivial.",
                 F73="Signifgicant barcode overlap, but numbers are trivial.",
                 F72="Signifgicant barcode overlap, but numbers are trivial.")
  if(nom %in% names(passAnyway)){
    message(sprintf("Passing sample anyway because: %s",passAnyway[nom]))
    dat$checksPassed=TRUE
  }
  dataSources[[nom]] = dat
}
