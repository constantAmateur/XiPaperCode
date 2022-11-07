#' Demultiplex mapped data
#'
#' Takes prepared data and demultiplexes and organises channels that have more than one donor per channel.  Ultimately, this should produce an extra object in dataSources, 'donorByCell', indicating which cells belong to which donor.
#'
#' External variables:
#' @param dataSources/scDat Defined by running prepData.R

#############
# Libraries #
#############

library(alleleIntegrator)
library(Seurat)

##########
# Params #
##########

outBase = 'Results/souporcell'

#############
# Functions #
#############

#' De multiplex mixed channel
#'
#' If a 10X channel contains multiple samples, identify which barcodes belong to which samples, without prior genotype information.
#' 
#' @param cellrangerBAM Path to the cellranger BAM to use.
#' @param bcodeFile Path to file containing barcodes.
#' @param k How many samples to demultiplex.
#' @param finDir Directory to store output files.
#' @param varFile VCF containing the common variants to use.
#' @param nRestarts How many restarts to do.
#' @param skipIfExists Skip if output files found.
#' @param nParallel How many threads to use.
#' @return Something...
deMultiplex = function(cellrangerBAM,bcodeFile,k,finDir,varFile=normalizePath('~/src/common_variants_grch38.vcf'),fasta=file.path(refGenome,'fasta','genome.fa'),nRestarts=500,skipIfExists=TRUE,nParallel=1){
  tmpBase=tempdir()
  outDir = tempfile()
  dir.create(outDir)
  #outDir = file.path(tmpBase,'outs')
  #finDir = file.path(normalizePath(dirname(cellrangerBAM)),'souporcellOut')
  outFiles = c('clusters.tsv','clusters.err','ambient_rna.txt','cluster_genotypes.vcf')
  if(!file.exists(finDir) || !skipIfExists){
    #Decide if we have chr or not...
    hasChr = seqinfo(Rsamtools::BamFile(cellrangerBAM))
    hasChr = any(grepl('^chr',seqnames(hasChr)))
    vcfHasChr = readLines(varFile,n=1000)
    vcfHasChr = length(grep('^chr',gsub('\t.*','',grep('^#',vcfHasChr,invert = TRUE,value=TRUE))))>0
    genomeHasChr = system(sprintf("grep '>' %s",fasta),intern=TRUE)
    genomeHasChr = any(grepl('^>chr',genomeHasChr))
    #Pull image
    system(sprintf('cd %s && singularity pull shub://wheaton5/souporcell',tmpBase))
    #Move things there, singularity struggles otherwise
    file.copy(cellrangerBAM,tmpBase)
    file.copy(paste0(cellrangerBAM,'.bai'),tmpBase)
    cellrangerBAM = basename(cellrangerBAM)
    #bcodeFile = file.path(dirname(cellrangerBAM),'filtered_feature_bc_matrix','barcodes.tsv.gz')
    #if(!file.exists(bcodeFile))
    #  bcodeFile = file.path(dirname(cellrangerBAM),'filtered_gene_bc_matrices','GRCh38','barcodes.tsv')
    file.copy(bcodeFile,tmpBase)
    bcodeFile = basename(bcodeFile)
    #Can we just copy the genome or do we need to modify it?
    if(hasChr==genomeHasChr){
      file.copy(fasta,tmpBase)
    }else{
      #Do we need to add it?
      if(hasChr){
        #Yes
        system(sprintf("sed 's/^>/>chr/g' %s > %s",fasta,file.path(tmpBase,'genome.fa')))
      }else{
        #Nope, remove it
        system(sprintf("sed 's/^>chr/>/g' %s > %s",fasta,file.path(tmpBase,'genome.fa')))
      }
    }
    #Won't copy this as it's fast to generate and another thing to modify
    #file.copy(file.path(refGenome,'fasta','genome.fa.fai'),tmpBase)
    #Can we just copy the VCF or do we need to do something
    if(hasChr==vcfHasChr){
      file.copy(varFile,file.path(tmpBase,'commonVariants.vcf'))
    }else{
      #Do we need to add it
      if(hasChr){
        #Yes
        system(sprintf("awk '{if($0 !~ /^#/) print \"chr\"$0; else print $0}' %s > %s",varFile,file.path(tmpBase,'commonVariants.vcf')))
      }else{
        #Nope, remove it
        system(sprintf("awk '{gsub(/^chr/,\"\"); print}' %s > %s",varFile,file.path(tmpBase,'commonVariants.vcf')))
      }
    }
    #Run pipeline
    cmd = sprintf('cd %s && singularity exec souporcell_latest.sif souporcell_pipeline.py --bam %s --barcodes %s -f %s -t %d -o %s -k %d -p 2 --restarts %d --common_variants %s --skip_remap true',
                  tmpBase,
                  cellrangerBAM,
                  bcodeFile,
                  'genome.fa',
                  nParallel,
                  outDir,
                  k,
                  nRestarts,
                  'commonVariants.vcf')
    system(cmd)
    #Check it worked.
    tgtFiles =file.path(outDir,outFiles)
    if(!all(file.exists(tgtFiles)))
      stop("Demultiplexing failed!")
    #Move files
    dir.create(finDir,recursive=TRUE)
    for(tgtFile in tgtFiles)
      system(sprintf('mv %s %s',tgtFile,finDir))
    #Delete tmp file
    system(sprintf('rm %s/*',tmpBase))
  }
  #Load and return output
  tgtFiles = file.path(finDir,outFiles)
  clusts = read.table(tgtFiles[1],sep='\t',header=TRUE)
  amb = readChar(tgtFiles[3],1e6)
  gtPath = tgtFiles[4]
  rho = as.numeric(gsub('\\%','',gsub('ambient RNA estimated as ','',amb)))/100
  return(list(clusts=clusts,rho=rho,gtPath=gtPath))
}



#' Match demultiplexed genotype
#'
#' If demultiplexing has been performed on multiples channels containing the same samples, it is necessary to work out which sample matches with which.  This uses the souporcell genotype vcf to do this matching very roughly.
#'
#' @param vcfs A collection of vcfs generated by souporcell with at least some samples overlapping.  Should be uniquely named.
#' @param minForMatch If the genotype agrees on fewer than this many SNPs, don't match samples even if they're the closest of the options.
#' @param minSNPs Need coverage of at least this many SNPs before we consider a match valid.
matchGenotypes = function(vcfs,minForMatch=0.5,minSNPs=500){
  noms = names(vcfs)
  if(is.null(noms) || any(duplicated(noms)))
    stop("Invalid vcfs parameter, should be named vector.")
  #Read them all in and get sample names
  dats = list()
  for(nom in noms){
    dats[[nom]] = read.table(vcfs[nom],sep='\t',header=FALSE)
    #Load the header
    hdr = readLines(vcfs[nom],n=100)
    hdr = hdr[max(grep('^#',hdr))]
    hdr = strsplit(gsub('^#','',hdr),'\t')[[1]]
    colnames(dats[[nom]]) = hdr
    #Drop the common columns, to be left with samples
    samps = hdr[!hdr %in% c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")]
    #Extract genotype string for each
    for(samp in samps){
      dats[[nom]][,paste0('GT_',samp)] = sapply(strsplit(dats[[nom]][,samp],':'),`[`,1)
    }
    #Convert to GRanges
    tmp = GRanges(dats[[nom]]$CHROM,IRanges(dats[[nom]]$POS,width=1))
    mcols(tmp) = dats[[nom]][,paste0('GT_',samps),drop=FALSE]
    dats[[nom]] = tmp
  }
  #Now one file at a time, one sample at a time
  samps = lapply(dats,function(e) gsub('^GT_','',colnames(mcols(e))))
  sampsFlat = paste0(rep(names(samps),lengths(samps)),':',unlist(samps,use.names=FALSE))
  out = matrix(NA,ncol=length(sampsFlat),nrow=length(sampsFlat))
  colnames(out) = rownames(out) = sampsFlat
  dd=list()
  for(i in seq_along(sampsFlat)){
    for(j in seq_along(sampsFlat)){
      a = sampsFlat[i]
      b = sampsFlat[j]
      aFile = gsub(':.*','',a)
      bFile = gsub(':.*','',b)
      aSamp = gsub('.*:','',a)
      bSamp = gsub('.*:','',b)
      #Match the two vcfs
      m = match(dats[[aFile]],dats[[bFile]])
      aGT = mcols(dats[[aFile]])[!is.na(m),paste0('GT_',aSamp)]
      bGT = mcols(dats[[bFile]])[m[!is.na(m)],paste0('GT_',bSamp)]
      #Now which of these have a genotype in both?
      w = which(aGT!='./.' & bGT!='./.')
      #Store summary stats
      dd[[length(dd)+1]] = data.frame(IDA=a,
                                      IDB=b,
                                      fileA=aFile,
                                      fileB=bFile,
                                      sampleA=aSamp,
                                      sampleB=bSamp,
                                      nA = length(dats[[aFile]]),
                                      nB = length(dats[[bFile]]),
                                      nCmn = length(aGT),
                                      nWithDatA = sum(aGT!='./.'),
                                      nWithDatB = sum(bGT!='./.'),
                                      nWithDat = length(w),
                                      nMatch = sum(aGT[w]==bGT[w]),
                                      matchFrac = sum(aGT[w]==bGT[w])/length(w))
      out[i,j] = dd[[length(dd)]]$matchFrac
    }
  }
  hm=out
  dd = do.call(rbind,dd)
  #Now for each sample, find the likely match in each file
  out = list()
  for(sampFlat in sampsFlat){
    for(nom in noms){
      tmp = dd[dd$IDA == sampFlat & dd$fileB==nom ,]
      #Work out if any match
      tmp$isMatch = FALSE
      if(any(tmp$nWithDat>minSNPs)){
        tt = max(tmp$matchFrac[tmp$nWithDat>minSNPs])
        if(!is.na(tt) && tt>minForMatch)
          tmp$isMatch[match(tt,tmp$matchFrac)]=TRUE
      }
      out[[length(out)+1]] = tmp
    }
  }
  out = do.call(rbind,out)
  #Get sample groupings across files
  grps = out[out$isMatch,]
  grps = unique(lapply(split(grps$IDB,grps$IDA),sort))
  return(list(matchStats=out,sampleGroupings=grps,matchMatrix=hm))
}

###################################
# Demultiplex maternal/fetal data #
###################################

#Load pre-done demultilplexing, thank you Roser!
demux = read.table('Data/demultiplexTable.tsv',sep='\t',header=TRUE)
demux$sampleID = gsub('_[ACGT]+(-[0-9]+)?$','',demux$cellID)
#Use pan-fetal immune cells where there are gaps
tmp = read.table('Data/panFetalMaternalBarcodes.csv',sep= ',',header=TRUE)
tmp = gsub('-','_',tmp$x)
tmp = gsub('_1$','-1',tmp)
tmp = paste0(tmp,ifelse(grepl('-1$',tmp),'','-1'))
#Fix up the gaps in F37 from pan-fetal.  These are the FCAImm samples and deliberately don't include any entries for FCA747067 as it's low quality
dat = scDat[scDat$donor=='F37',]
dat = dat[grep('FCAImm',dat$ID),]
bcodes = Read10X(setNames(dat$stagedPathCnts,dat$ID))
bcodes = paste0(colnames(bcodes),'-1')
demux = rbind(demux,data.frame(cellID = bcodes,
                               donorID = ifelse(bcodes %in% tmp,'M','F'),
                               dataset = 'Placenta',
                               sampleID = gsub('_[ACGT]*-1$','',bcodes)))
#This concludes the sanger portion of our presentation.  On to HCA

