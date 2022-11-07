#' Tissue specific processing of oral mucosa.


##########
# Params #
##########

srcs = scDat[which(scDat$dataSrc=='HCA_OralMucosa'),]

#############
# Load Fits #
#############

tgtFile = 'Results/oralMucosa.RDS'
if(file.exists(tgtFile)){
  srat = readRDS(tgtFile)
}else{
  fits = list()
  for(nom in unique(srcs$donor)){
    message(sprintf("Loading fit for %s",nom))
    fits[[nom]] = readRDS(sprintf('Results/fitEM/%s_fit.RDS',nom))$fitList$FromRNA
  }
  taus = sapply(fits,function(e) e$tau)
  nSNPs = sapply(fits,function(e) length(e$genotype))
  nCells = sapply(fits,function(e) length(e$states))
  #Construct seurat object
  srat = Read10X(setNames(srcs$stagedPathCnts,srcs$ID))
  colnames(srat) = ifelse(grepl('-[0-9]+$',colnames(srat)),colnames(srat),paste0(colnames(srat),'-1'))
  #Do any QCing of the cells
  srat = scQC$quickCluster(srat)
  srat = FindClusters(srat,resolution=2)
  srat@meta.data$orig.ident = gsub('_[ACGT]+-[0-9]+$','',colnames(srat))
  srat@meta.data$donor = scDat$donor[match(srat@meta.data$orig.ident,scDat$ID)]
  #Annotate the thing again...
  annot = list(Tcells=c(14,9,5,13,42),
               NK = c(26,43),
               Erythroid = c(8,27,55,51,36,52,56,45,58,32,11,17,31,39),
               Macrophages = c(33,20,57,37,44),
               Mast=c(61,19,62),
               Bcell=38,
               Plasma=c(21,18),
               VascularEndo=c(4,0,3,40,15,12,10),
               LymphEndo=c(24,59),
               Perivascular=c(50,6,28,30),
               Crap=c(54,49),
               Fibroblasts=c(34,41,1,2,29,25,7,16),
               Epithelium=c(47,35,23,22)
               )
  srat@meta.data$annot = setNames(rep(names(annot),lengths(annot)),unlist(annot))[as.character(srat@active.ident)]
  srat@meta.data$annot[srat@meta.data$annot=='Crap']=NA
  states = unlist(lapply(fits,function(e) e$states))
  names(states) = paste0(gsub('.*\\.','',names(states)),'-1')
  srat@meta.data$states = states[colnames(srat)]
  srat@meta.data$stateX = ifelse(is.na(srat@meta.data$states) | abs(srat@meta.data$states-0.5)<0.4,NA,ifelse(srat@meta.data$states>0.5,'Maternal','Paternal'))
  srat@misc$fitStats = list(taus=taus,nSNPs=nSNPs,nCells=nCells)
  saveRDS(srat,tgtFile)
}

############
# Analysis #
############

t1=nParallel
nParallel = min(12,t1)
dd = makeStandardPlots(srat,'Results/plots/oralMucosa')
nParallel = t1



