#' Comparitive anylsis of lung samples

#############
# Libraries #
#############


##########
# Params #
##########

#Sources of data.  NA means all channel
srcs = scDat[grep('HCA_[CP]F',scDat$donor),]
#Keep only the females that have a fit object
srcs = srcs[file.exists(sprintf('Results/fitEM/%s_fit.RDS',srcs$donor)),]

#############
# Load fits #
#############

tgtFile = 'Results/HCA_lung.RDS'
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
  #Annotate with celltypist
  annot = runCelltypist(srat@assays$RNA@counts,model='Human_Lung_Atlas.pkl')
  #Stick cell level calls into srat 
  srat@meta.data$annotCell = annot$labMat$predicted_labels[match(colnames(srat),annot$labMat$X)]
  #Do majority voting by cluster. Purge some of the shit clusters
  #Get the counts for each cluster
  clMap = setNames(rep(NA,length(levels(srat@active.ident))),levels(srat@active.ident))
  x = table(srat@meta.data$annotCell,srat@meta.data$seurat_clusters)
  class(x)='matrix'
  x = t(t(x)/colSums(x))
  maxCl = setNames(rownames(x)[apply(x,2,which.max)],colnames(x))
  w = which(apply(x,2,max)>0.75)
  clMap[w] = maxCl[w]
  w = which(is.na(clMap) & lengths(apply(x,2,function(e) e[e>0.1]))<=2)
  clMap[w] = maxCl[w]
  #Do this and manually make a choiec for the rest
  #w = which(is.na(clMap))
  #apply(x[,w],2,function(e) sort(e[e>0]*100))
  clMap[as.character(42)]='Fibroblasts'
  clMap[as.character(c(19,38))]='Macrophages'
  clMap[as.character(37)]='Plasma'
  clMap[as.character(c(70,67))]='Erythroid'
  clMap[as.character(12)]=maxCl[as.character(12)]
  clMap[as.character(17)]=maxCl[as.character(17)]
  clMap[as.character(25)]=maxCl[as.character(25)]
  #These ones are a bit iffy, but they'll all get merged below anyway
  clMap[as.character(21)]='Basal'
  clMap[as.character(52)]=maxCl[as.character(52)]
  clMap[as.character(49)]=maxCl[as.character(49)]
  clMap[as.character(78)]=NA
  clMap[as.character(75)]=NA
  srat@meta.data$annotCl = clMap[as.character(srat@active.ident)]
  #srat@meta.data$annot = sapply(split(srat@meta.data$annotCell,srat@meta.data$seurat_clusters),function(e) names(-sort(-round(table(e)/length(e)*100))[1]))[srat@meta.data$seurat_clusters]
  #Smush some of the similar things
  annot = list(Tcells = c('CD4 T cells','CD8 T cells','T cells proliferating'),
               VascularEndo = c("EC aerocyte capillary","EC arterial","EC general capillary","EC venous systemic"),
               LymphEndo = 'Lymphatic EC mature',
               Monocytes = c("Monocyte-derived Mφ",'Macrophages'),
               Basal =c('Suprabasal','SMG serous (bronchial)','Transitional Club-AT2','Basal resting'),
               Plasma = c('Plasma','Plasma cells'))
  annot = setNames(rep(names(annot),lengths(annot)),unlist(annot))
  srat@meta.data$annotCl = ifelse(srat@meta.data$annotCl %in% names(annot),annot[srat@meta.data$annotCl],srat@meta.data$annotCl)
  #Get the state mappings
  states = unlist(lapply(fits,function(e) e$states))
  names(states) = gsub('.*\\.','',names(states))
  srat@meta.data$donor = scDat$donor[match(srat@meta.data$orig.ident,scDat$ID)]
  srat@meta.data$states = states[colnames(srat)]
  srat@meta.data$stateX = ifelse(is.na(srat@meta.data$states) | abs(srat@meta.data$states-0.5)<0.4,NA,ifelse(srat@meta.data$states>0.5,'Maternal','Paternal'))
  #Save it out
  saveRDS(srat,tgtFile)
}

############
# Analysis #
############

#nParallel=12
t1 = nParallel
nParallel = min(12,nParallel)
dd = makeStandardPlots(srat,'Results/plots/lung',annotCol='annotCl')
nParallel = t1

