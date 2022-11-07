
##########
# Params #
##########

srcs = list('P13'=NA,
            'Hrv98'=NA,
            'Hrv99'=NA,
            'Hrv100'=NA,
            'F37'=c('FCA7474068','FCA7474069'),
            'F40'=c('FCA7511884','FCA7511885','FCA7511886'))


#############
# Load fits #
#############

tgtFile = 'Results/placenta_F.RDS'
if(file.exists(tgtFile)){
  srat = readRDS(tgtFile)
}else{
  fits = list()
  for(nom in names(srcs)){
    message(sprintf("Loading fit for %s",nom))
    fits[[nom]] = readRDS(sprintf('Results/fitEM/%s___F_fit.RDS',nom))$fitList$FromRNA
  }
  taus = sapply(fits,function(e) e$tau)
  nSNPs = sapply(fits,function(e) length(e$genotype))
  nCells = sapply(fits,function(e) length(e$states))
  #Construct seurat object
  srat = lapply(names(srcs),function(e) if(length(srcs[e])==1 && is.na(srcs[e])){scDat[scDat$donor==e,,drop=FALSE]}else{scDat[match(srcs[[e]],scDat$ID),]})
  srat = do.call(rbind,srat)
  srat = Read10X(setNames(srat$stagedPathCnts,srat$ID))
  #Keep only those that are fetal
  srat = srat[,colnames(srat) %in% gsub('-[0-9]+$','',demux$cellID)[demux$donorID=='F']]
  srat = scQC$quickCluster(srat)
  srat = FindClusters(srat,resolution=2)
  srat@meta.data$orig.ident = gsub('_[ACGT]+$','',colnames(srat))
  srat@meta.data$donor = scDat$donor[match(srat@meta.data$orig.ident,scDat$ID)]
  annot = list(Hofbauer=c(19,28,45,14,21,48),
               Erythroid=c(33),
               Crap=c(52,49),
               Macrophages=c(13,16,32),
               Bcell=51,
               Tcell=47,
               NK=c(43,30),
               Endothelium=c(53,26,36),
               Epithelium=c(46,23),
               Stroma=c(38,15),
               Fibroblasts=c(40,24,34,6,37,9,41),
               Perivascular=c(12,29),
               VCT=c(27,35,50,10,17,39,3,4,25,2,18,11),
               SCT=c(31,1),
               EVT=c(8,44,20,22,5,0,7,42))
  srat@meta.data$annot = setNames(rep(names(annot),lengths(annot)),unlist(annot))[as.character(srat@active.ident)]
  srat@meta.data$annot[srat@meta.data$annot=='Crap']=NA
  states = setNames(unlist(lapply(fits,function(e) e$cellSummary$stateXi)),unlist(lapply(fits,function(e) e$cellSummary$cell)))
  srat@meta.data$stateX = ifelse(states[colnames(srat)]=='Undetermined',NA,states[colnames(srat)])
  states = unlist(lapply(fits,function(e) e$states))
  names(states) = gsub('.*\\.','',names(states))
  srat@meta.data$states = states[colnames(srat)]
  #srat@meta.data$stateX = ifelse(is.na(srat@meta.data$states) | abs(srat@meta.data$states-0.5)<0.4,NA,ifelse(srat@meta.data$states>0.5,'Maternal','Paternal'))
  srat@misc$fitStats = list(taus=taus,nSNPs=nSNPs,nCells=nCells)
  srat@misc$fits = fits
  saveRDS(srat,tgtFile)
}

t1=nParallel
nParallel = min(12,t1)
dd = makeStandardPlots(srat,'Results/plots/placenta')
nParallel = t1


