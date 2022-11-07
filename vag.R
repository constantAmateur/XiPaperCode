#' Comparitive anylsis of lung samples

#############
# Libraries #
#############


##########
# Params #
##########

#Sources of data.  NA means all channel
srcs = scDat[grep('HCA_VagWall',scDat$donor),]
#Keep only the females that have a fit object
srcs = srcs[file.exists(sprintf('Results/fitEM/%s_fit.RDS',srcs$donor)),]


#############
# Load fits #
#############

tgtFile = 'Results/HCA_vag.RDS'
if(file.exists(tgtFile)){
  srat = readRDS(tgtFile)
}else{
  fits = list()
  for(nom in unique(srcs$donor)){
    if(nom %in% names(fits))
      next
    message(sprintf("Loading fit for %s",nom))
    tgt = sprintf('Results/fitEM/%s_fit.RDS',nom)
    if(file.exists(tgt))
      fits[[nom]] = readRDS(sprintf('Results/fitEM/%s_fit.RDS',nom))$fitList$FromRNA
  }
  taus = sapply(fits,function(e) e$tau)
  nSNPs = sapply(fits,function(e) length(e$genotype))
  nCells = sapply(fits,function(e) length(e$states))
  ###################
  # Construct seurat
  srat = setNames(srcs$stagedPathCnts,srcs$ID)
  srat = Read10X(srat)
  srat = scQC$quickCluster(srat)
  srat = FindClusters(srat,resolution=2)
  #Store the basics
  srat@misc$fitStats = list(taus=taus,nSNPs=nSNPs,nCells=nCells)
  #Annotate with broad markers from the paper (https://www.nature.com/articles/s41467-020-20358-y)
  annot = list(Epithelium=c(33,39,49,57,35),
               Fibroblasts=c(42,31,48,3,51,17,25,28,16,26,4,23,0,58,41,21,45,38,14,2,54),
               SmoothMuscle=c(11,5,24,53,37,6,9),
               Myoepithelium=c(13),
               Endothelium=c(10,15,18),
               LymphEndo=c(29),
               Macro=c(19,20,46,32,22),
               TCells=c(1,36,8,7,55,12),
               BCells=c(34,56),
               Plasma=c(40,44),
               Erythroid=c(47),
               Mast=c(27),
               Neurons=c(43)
               )
  #The un-annotated cluster is pure shit, so it should be ignored
  srat@meta.data$annot = setNames(rep(names(annot),lengths(annot)),unlist(annot))[as.character(srat@active.ident)]
  #Get the state mappings
  states = unlist(lapply(fits,function(e) e$states))
  names(states) = gsub('.*\\.','',names(states))
  names(states) = gsub('-1$','',names(states))
  srat@meta.data$orig.ident = gsub('_[ACGT]*(-1)?$','',colnames(srat))
  srat@meta.data$donor = scDat$donor[match(srat@meta.data$orig.ident,scDat$ID)]
  srat@meta.data$states = states[gsub('-1$','',colnames(srat))]
  srat@meta.data$stateX = ifelse(is.na(srat@meta.data$states) | abs(srat@meta.data$states-0.5)<0.4,NA,ifelse(srat@meta.data$states>0.5,'Maternal','Paternal'))
  #Save it out
  saveRDS(srat,tgtFile)
}

############
# Analysis #
############

t1 = nParallel
nParallel = min(12,t1)
dd = makeStandardPlots(srat,'Results/plots/vag')
nParallel = t1






