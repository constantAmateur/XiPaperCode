#' Comparitive anylsis of kidney samples

#############
# Libraries #
#############


##########
# Params #
##########

#Sources of data.  NA means all channel
srcs = list('HCA_KidCortex005'=NA,
            "HCA_kidney002" = NA, 
            "HCA_kidney007" = NA, 
            "HCA_kidney008" = NA, 
            "HCA_kidney009" = NA, 
            "HCA_kidney010" = NA, 
            "HCA_kidney012" = NA, 
            "HCA_kidney013" = NA,
            "HCA_kidney014" = NA, 
            "HCA_kidney015" = NA,
            "HCA_kidney017" = NA,
            "HCA_kidney018" = NA,
            'RCC2'=c('4602STDY6976426','4602STDY6976427','4602STDY6976428'),
            'VHL_RCC'=c('4602STDY6949178','4602STDY6949179','4602STDY6949180','4602STDY6949181'),
            'Wilms2'=c('4602STDY7018926','4602STDY7018923','4602STDY7018924','4602STDY7018925'),
            'Wilms3'=c('4602STDY7090428','4602STDY7090429'),
            'KidTransplant' = c('4602STDY6949184','4602STDY6949185','4602STDY6949187','4602STDY6949188','4602STDY6949186','4602STDY6949189'),
            'F35'=c('FCAImmP7462242','FCAImmP7462243'),
            'F41'=c('FCAImmP7555859'),
            'F45'=c('FCAImmP7579214','FCAImmP7579215','FCAImmP7579225'))



#############
# Load fits #
#############

tgtFile = 'Results/kidney.RDS'
if(file.exists(tgtFile)){
  srats = readRDS(tgtFile)
  fKid = srats$fKid
  mKid = srats$mKid
}else{
  fits = list()
  for(nom in names(srcs)){
    message(sprintf("Loading fit for %s",nom))
    fits[[nom]] = readRDS(sprintf('Results/fitEM/%s_fit.RDS',nom))$fitList$FromRNA
  }
  taus = sapply(fits,function(e) e$tau)
  nSNPs = sapply(fits,function(e) length(e$genotype))
  nCells = sapply(fits,function(e) length(e$states))
  #########################
  # Construct fetal object
  srcsF = srcs[grep('^F[0-9]+$',names(srcs))]
  srat = scDat[match(unlist(srcsF),scDat$ID),]
  srat = Read10X(setNames(srat$stagedPathCnts,srat$ID))
  colnames(srat) = paste0(colnames(srat),'-1')
  mtFrac = colSums(srat[grep("^MT-",rownames(srat)),])/colSums(srat)
  srat = srat[,mtFrac<0.6]
  srat = scQC$quickCluster(srat)
  srat = FindClusters(srat,resolution=2)
  srat@meta.data$donor = gsub('_[^_]+$','',gsub('_[ACGT]+(-[0-9])?$','',colnames(srat)))
  srat@meta.data$donor = scDat$donor[(match(srat@meta.data$donor,scDat$ID))]
  states = unlist(lapply(fits,function(e) e$states))
  names(states) = gsub('.*\\.','',names(states))
  srat@meta.data$states = states[colnames(srat)]
  srat@meta.data$stateX = ifelse(is.na(srat@meta.data$states) | abs(srat@meta.data$states-0.5)<0.4,NA,ifelse(srat@meta.data$states>0.5,'Maternal','Paternal'))
  annot = runCelltypist(srat@assays$RNA@counts,model='Data/Fetal_Full_v3.pkl')
  annotMY = runCelltypist(srat@assays$RNA@counts,model='Data/fKid.pkl')
  srat@meta.data$annotCell = annot$labMat$predicted_labels[match(colnames(srat),annot$labMat$X)]
  #Do majority voting by cluster. Dumb summary with no checks
  srat@meta.data$annot = sapply(split(srat@meta.data$annotCell,srat@meta.data$seurat_clusters),function(e) names(-sort(-round(table(e)/length(e)*100))[1]))[srat@meta.data$seurat_clusters]
  srat@meta.data$annotMYCell = annotMY$labMat$predicted_labels[match(colnames(srat),annotMY$labMat$X)]
  srat@meta.data$annotMY = sapply(split(srat@meta.data$annotMYCell,srat@meta.data$seurat_clusters),function(e) names(-sort(-round(table(e)/length(e)*100))[1]))[srat@meta.data$seurat_clusters]
  x = do.call(rbind,lapply(unique(srat@meta.data$seurat_clusters),function(e) srat@meta.data[srat@meta.data$seurat_clusters==e,c('annot','annotMY','seurat_clusters')][1,]))
  x$seurat_clusters = as.character(x$seurat_clusters)
  rownames(x) = NULL
  #This is the cellGen Mapping annotation
  finalAnn = list(Tcell = c(24,46),
                  NK = c(11),
                  Bcell=c(30),
                  Mast=40,
                  MK=41,
                  Erythroid=45,
                  Neutrophil=38,
                  Macrophage=c(4,3,1),
                  cDC1 = 43,
                  cDC2 = c(22,36,31,15),
                  MSC = c(14),
                  ICa = c(10,7,16,23,8,25,6,33),
                  ICb = c(35,2),
                  CapMes = c(0,12,9,13,18,39,17,29),
                  Pod=c(26),
                  ErPrT=c(21,34,5),
                  SSBm.d=c(28),
                  DTLH=42,
                  CnT=c(27,20),
                  UBCD=44,
                  Endothelium=19)
  srat@meta.data$annotFull = setNames(rep(names(finalAnn),lengths(finalAnn)),as.character(unlist(finalAnn)))[as.character(srat@meta.data$seurat_clusters)]
  m = unique(srat@meta.data$donor)
  srat@misc$fitStats = list(taus=taus[m],nSNPs=nSNPs[m],nCells=nCells[m])
  fKid = srat
  ##########################
  # Construct mature object
  srcsM = srcs[grep('^F[0-9]+$',names(srcs),invert=TRUE)]
  srat = lapply(names(srcsM),function(e) if(length(srcsM[e])==1 && is.na(srcsM[e])){scDat[scDat$donor==e,,drop=FALSE]}else{scDat[match(srcsM[[e]],scDat$ID),]})
  srat = do.call(rbind,srat)
  srat = Read10X(setNames(srat$stagedPathCnts,srat$ID))
  colnames(srat) = ifelse(grepl('-1$',colnames(srat)),colnames(srat),paste0(colnames(srat),'-1'))
  mtFrac = colSums(srat[grep("^MT-",rownames(srat)),])/colSums(srat)
  srat = srat[,mtFrac<0.6]
  srat = scQC$quickCluster(srat)
  srat = FindClusters(srat,resolution=2)
  srat@meta.data$donor = gsub('_[^_]+$','',gsub('_[ACGT]+(-[0-9])?$','',colnames(srat)))
  srat@meta.data$donor = scDat$donor[(match(srat@meta.data$donor,scDat$ID))]
  states = unlist(lapply(fits,function(e) e$states))
  names(states) = gsub('.*\\.','',names(states))
  srat@meta.data$states = states[colnames(srat)]
  srat@meta.data$stateX = ifelse(is.na(srat@meta.data$states) | abs(srat@meta.data$states-0.5)<0.4,NA,ifelse(srat@meta.data$states>0.5,'Maternal','Paternal'))
  annot= runCelltypist(srat@assays$RNA@counts,model='Data/Mature_Full_v3.pkl')
  annotMY = runCelltypist(srat@assays$RNA@counts,model='Data/fKid.pkl')
  annotFetalKid = runCelltypist(srat@assays$RNA@counts,model='Data/Fetal_Full_v3.pkl')
  srat@meta.data$annotCell = annot$labMat$predicted_labels[match(colnames(srat),annot$labMat$X)]
  #Do majority voting by cluster. Dumb summary with no checks
  srat@meta.data$annot = sapply(split(srat@meta.data$annotCell,srat@meta.data$seurat_clusters),function(e) names(-sort(-round(table(e)/length(e)*100))[1]))[srat@meta.data$seurat_clusters]
  srat@meta.data$annotMYCell = annotMY$labMat$predicted_labels[match(colnames(srat),annotMY$labMat$X)]
  srat@meta.data$annotMY = sapply(split(srat@meta.data$annotMYCell,srat@meta.data$seurat_clusters),function(e) names(-sort(-round(table(e)/length(e)*100))[1]))[srat@meta.data$seurat_clusters]
  srat@meta.data$annotFetKidCell = annotFetalKid$labMat$predicted_labels[match(colnames(srat),annotFetalKid$labMat$X)]
  srat@meta.data$annotFetKid = sapply(split(srat@meta.data$annotFetKidCell,srat@meta.data$seurat_clusters),function(e) names(-sort(-round(table(e)/length(e)*100))[1]))[srat@meta.data$seurat_clusters]
  x = do.call(rbind,lapply(unique(srat@meta.data$seurat_clusters),function(e) srat@meta.data[srat@meta.data$seurat_clusters==e,c('annot','annotMY','seurat_clusters')][1,]))
  x$seurat_clusters = as.character(x$seurat_clusters)
  rownames(x) = NULL
  #Annotation of cellGen Mapped version
  finalAnn = list(Bcell = 46,
                  Tcell = c(26,25),
                  NK=c(31),
                  Macrophage=c(30,32),
                  Pelvis=c(44),
                  Podocyte=35,
                  Endothelium=c(13,50,24,39,37),
                  Myofibroblast=c(21),
                  IntercalatedCell=c(27,43,51,18,40),
                  DistalTubules=c(29,9,16,19,23,12),
                  LoH=c(33,14,15,20,42),
                  Erythroid=c(36,34),
                  PT=c(17,22,2,1,0,52,8,5,6,11,4,7,3,10,47,53)
                  )
  srat@meta.data$annotFull = setNames(rep(names(finalAnn),lengths(finalAnn)),as.character(unlist(finalAnn)))[as.character(srat@meta.data$seurat_clusters)]
  m = unique(srat@meta.data$donor)
  srat@misc$fitStats = list(taus=taus[m],nSNPs=nSNPs[m],nCells=nCells[m])
  mKid = srat
  saveRDS(list(mKid=mKid,fKid=fKid),tgtFile)
}

############
# Analysis #
############

t1 = nParallel
nParallel = min(t1,12)
dd = makeStandardPlots(mKid,'Results/plots/kidneyAdult',annotCol='annotFull')
#Not enough samples to bother with this
#dd = makeStandardPlots(fKid,'Results/plots/kidneyFetal',annotCol='annotFull')
nParallel=t1

