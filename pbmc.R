#' Do stuff with the PBMC data
#'
#' Starting with AIDA

#############
# Libraries #
#############



##########
# Params #
##########

ageBreaks=c(10,30,40,60,80)
ageCut=35
srcs = scDat[scDat$dataSrc=='HCA_AIDA',]
srcs = srcs[!is.na(srcs$stagedPathCnts),]

##############
# Processing #
##############

############
# Load fits
tgtFile = 'Results/HCA_AIDA.RDS'
if(file.exists(tgtFile)){
  srat = readRDS(tgtFile)
}else{
  fits = list()
  for(nom in unique(srcs$donor)){
    message(sprintf("Loading fit for %s",nom))
    tgt = sprintf('Results/fitEM/%s_fit.RDS',nom)
    if(file.exists(tgt))
      fits[[nom]] = readRDS(sprintf('Results/fitEM/%s_fit.RDS',nom))$fitList$FromRNA
  }
  taus = sapply(fits,function(e) e$tau)
  nSNPs = sapply(fits,function(e) length(e$genotype))
  nCells = sapply(fits,function(e) length(e$states))
  ages = setNames(as.numeric(gsub(' year','',dataSources$HCA_AIDA$donor_organism.organism_age[match(names(fits),dataSources$HCA_AIDA$humanDonorName)])),names(fits))
  oldAnnot = readRDS('Results/HCA_AIDA_old.RDS')@meta.data
  ###################
  # Construct seurat
  srat = setNames(srcs$stagedPathCnts,srcs$ID)
  srat = Read10X(srat)
  colnames(srat) = ifelse(grepl('-[1-9]+$',colnames(srat)),colnames(srat),paste0(colnames(srat),'-1'))
  srat = scQC$quickCluster(srat)
  srat = FindClusters(srat,resolution=2)
  #Store the basics
  srat@misc$fitStats = list(taus=taus,ages=ages,nSNPs=nSNPs,nCells=nCells)
  #Bring over old annotation
  srat@meta.data$donor = scDat$donor[match(gsub('_[ACGT]+-1$','',colnames(srat)),scDat$ID)]
  srat@meta.data$orig.ident = gsub('_[ACGT]+-1$','',colnames(srat))
  srat@meta.data$barcode = gsub('(.*)_([ACGT]+-1)$','\\2',colnames(srat))
  m = match(paste0(srat@meta.data$donor,'_',srat@meta.data$barcode),rownames(oldAnnot))
  srat@meta.data$annotCell = oldAnnot$annot[m]
  #Do majority voting by cluster. Dumb summary with no checks
  srat@meta.data$annot = sapply(split(srat@meta.data$annotCell,srat@meta.data$seurat_clusters),function(e) names(-sort(-round(table(e)/length(e)*100))[1]))[srat@meta.data$seurat_clusters]
  ##Annotate with celltypist
  #annot = runCelltypist(srat@assays$RNA@counts,model='Immune_All_Low.pkl')
  ##Stick cell level calls into srat 
  #srat@meta.data$annotCell = annot$labMat$predicted_labels[match(colnames(srat),annot$labMat$X)]
  ##Do majority voting by cluster. Dumb summary with no checks
  #srat@meta.data$annot = sapply(split(srat@meta.data$annotCell,srat@meta.data$seurat_clusters),function(e) names(-sort(-round(table(e)/length(e)*100))[1]))[srat@meta.data$seurat_clusters]
  ##Fix up the annotation using Laura's input
  #srat@meta.data$annot[srat@meta.data$RNA_snn_res.1 %in% c(18)]='gdT'
  #srat@meta.data$annot[srat@meta.data$RNA_snn_res.1 %in% c(21)]='NKT'
  #Get the state mappings
  states = unlist(lapply(names(fits),function(e) setNames(fits[[e]]$states,paste0(e,gsub('.*_','_',names(fits[[e]]$states))))))
  m = match(paste0(srat@meta.data$donor,'_',srat@meta.data$barcode),names(states))
  srat@meta.data$states = states[m]
  srat@meta.data$stateX = ifelse(is.na(srat@meta.data$states) | abs(srat@meta.data$states-0.5)<0.4,NA,ifelse(srat@meta.data$states>0.5,'Maternal','Paternal'))
  srat@meta.data$donorAge = ages[srat@meta.data$donor]
  #Save it out
  saveRDS(srat,tgtFile)
}

############
# Analysis #
############

t1 = nParallel
nParallel=min(8,t1)
#Drop those that have a reasonable chance of containing clonal haem
tmp = srat[,srat@meta.data$donorAge<ageCut]
w = which(tmp@misc$fitStats$ages<ageCut)
tmp@misc$fitStats = lapply(tmp@misc$fitStats,function(e) e[w])
dd = makeStandardPlots(tmp,'Results/plots/pbmcs')
nParallel=t1

