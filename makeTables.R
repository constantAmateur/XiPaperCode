#' Generate supplementary tables for paper

######################
# Sample level table #
######################
sampTab = datMap[datMap$donor %in% unlist(areUsedAnalysis),c('ID','donor','dataType','dataSrc')]
colnames(sampTab)[1] = 'sampleID'
#Drop spatial
sampTab = sampTab[sampTab$dataType!='spatial',]
#Drop RMS
sampTab = sampTab[sampTab$donor!='GOSH028',]
#Extra study level details
sampTab$dataSrc[sampTab$donor %in% c("RCC2","VHL_RCC","Wilms2","Wilms3","KidTransplant")] = 'SANGER_kidney'
sampTab$dataSrc[sampTab$donor %in% c('GOSH028')] = 'SANGER_RMS'
sampTab$dataSrc[sampTab$donor %in% c('GOSH029','IMP672')] = 'SANGER_kidney2'
sampTab$dataSrc[sampTab$donor %in% c('GOSH025')] ='SANGER_neuroblastoma'
sampTab$dataSrc[sampTab$donor %in% c('GOSH032')] = 'SANGER_alleleIntegrator'
sampTab$dataSrc[sampTab$donor %in% c("F37","F40","F45","F41","F35")] = 'SANGER_dev'
sampTab$dataSrc[sampTab$donor %in% c("Hrv98","Hrv99","Hrv100","P13")] = 'SANGER_placenta'
#Add in extra details about each study.
studyDetails = read.table('Data/studyDetails.tsv',sep='\t',header=TRUE)
sampTab = merge(sampTab,studyDetails,by='dataSrc')
#Add a column saying what grouping this sample is used for
tt = data.frame(donor = unlist(areUsedAnalysis),usedFor=rep(names(areUsedAnalysis),lengths(areUsedAnalysis)))
tt = sapply(split(tt$usedFor,tt$donor),function(e) paste(unique(e),collapse=','))
sampTab$usedFor = tt[sampTab$donor]
#Now need fine grained disease, etc information
sampTab$disease = NA
sampTab$age = NA
sangerMani = read.delim('Data/sampleManifestSanger.tsv',sep='\t',header=TRUE)
for(dSrc in unique(sampTab$dataSrc)){
  w = which(sampTab$dataSrc==dSrc)
  if(grepl('^HCA_',dSrc)){
    m = match(sampTab$sampleID[w],dataSources[[dSrc]]$sequencing_input.biomaterial_core.biomaterial_id)
    sampTab$disease[w] = dataSources[[dSrc]]$donor_organism.diseases[m]
    sampTab$age[w] = dataSources[[dSrc]]$donor_organism.organism_age[m]
  }else{
    m = match(sampTab$sampleID[w],sangerMani$sangerSampleID1)
    if(any(is.na(m))){
      m[is.na(m)] = match(sampTab$sampleID[w[is.na(m)]],sangerMani$sangerSampleID2)
    }
    sampTab$disease[w] = sangerMani$Diagnosis[m]
  }
  if(dSrc=='SANGER_dev'){
    t1 = read.table('Data/PAN.A01.v01.20210429.sample_metadata.csv',sep=',',header=TRUE)
    t1$organ = c(SP='Spleen',TH='Thymus',SK='Skin',BM='BoneMarrow',LI='Liver',KI='Kidney',YS='YolkSac',GU='Gut',MLN='MLN')[t1$organ]
    m = match(sampTab$sampleID[w],t1$X)
    sampTab$tissue[w] = ifelse(is.na(m),'Placenta',t1$organ[m])
    sampTab$age[w] = ifelse(is.na(m),NA,paste0(t1$age[m],' PCW'))
  }
}
#Final fix up
sampTab$disease[sampTab$usedFor %in% c('bigIndividual,kidney','pbmc','placenta','kidney')]='Normal'
sampTab$disease[grepl('Mature Kidney',sampTab$disease)]='Normal'
sampTab$disease[sampTab$disease %in% c('Mature Ureter','Wilms Kidney Cortex DNA','Wilms Paternal DNA','Wilms Maternal DNA','Kidney Cortex DNA','Ewings Sarcoma normal DNA','Wilms without treatment, Normal DNA','MRT Normal DNA','Alveolar RMS Normal DNA')]='Normal'
sampTab$disease[is.na(sampTab$disease) & sampTab$usedFor=='goldStandard,kidney']='Normal'
sampTab$disease[is.na(sampTab$disease) & sampTab$usedFor=='goldStandard']='Normal'
sampTab$disease[sampTab$disease=='Wilms without treatment']='Wilms'
sampTab$disease[sampTab$disease=='Neuroblastoma SC']='Neuroblastoma'
sampTab$disease[sampTab$disease=='Normal']='normal'
sampTab$disease[grep('ccRCC',sampTab$disease)]='ccRCC'
write.table(sampTab,'Results/sampleTable.tsv',sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE)


####################
# Cell level table #
####################
#The main lot of samples
tgtFiles = c(Lung='Results/HCA_lung.RDS',
             Vagina = 'Results/HCA_vag.RDS',
             Oral = 'Results/oralMucosa.RDS',
             PBMC = 'Results/HCA_AIDA.RDS',
             Kidney = 'Results/kidney.RDS',
             Placenta = 'Results/placenta_F.RDS')
cellTab = list()
for(nom in names(tgtFiles)){
  message(nom)
  srat = readRDS(tgtFiles[nom])
  if(nom=='Kidney')
    srat = srat$mKid
  t1 = srat@meta.data
  #Keep only those that we're using.  The sample table should determine this
  t1 = t1[(t1$orig.ident %in% sampTab$sampleID) & (t1$donor %in% sampTab$donor),]
  annotTgt = 'annot'
  if(nom=='Kidney')
    annotTgt='annotFull'
  if(nom=='Lung')
    annotTgt='annotCl'
  cellTab[[nom]] = data.frame(row.names = rownames(t1),
                              cellID = rownames(t1),
                              sampleID = t1$orig.ident,
                              donor = t1$donor,
                              donorGroup = nom,
                              annot = t1[,annotTgt],
                              stateProbXi = t1$states,
                              stateXi = t1$stateX,
                              goldStateXi = NA)
}
#The big boy
srat = readRDS('Results/fitEM/F45_fit_withAnnot.RDS')
t1 = srat@meta.data
t1$donor = sampTab$donor[match(t1$orig.ident,sampTab$sampleID)] 
t1 = t1[(t1$orig.ident %in% sampTab$sampleID) & (t1$donor %in% sampTab$donor),]
nom = 'highCellNumberIndividual'
cellTab[[nom]] = data.frame(row.names = rownames(t1),
                              cellID = rownames(t1),
                              sampleID = t1$orig.ident,
                              donor = t1$donor,
                              donorGroup = nom,
                              annot = t1$annot,
                              stateProbXi = t1$stateProbs,
                              stateXi = t1$stateXi,
                              goldStateXi = NA)
#Gold standard samples
goldData = c('Wilms2','Wilms3','GOSH025')
for(nom in goldData){
  message(nom)
  dat = readRDS(file.path('Results/fitEM',paste0(nom,'_fit.RDS')))
  srat = dat$srat
  t1 = srat@meta.data
  t1$donor = sampTab$donor[match(t1$orig.ident,sampTab$sampleID)] 
  t1 = t1[(t1$orig.ident %in% sampTab$sampleID) & (t1$donor %in% sampTab$donor),]
  cellTab[[nom]] = data.frame(row.names = rownames(t1),
                              cellID = rownames(t1),
                              sampleID = t1$orig.ident,
                              donor = t1$donor,
                              donorGroup = 'GoldStandard',
                              annot = NA,
                              stateProbXi = t1$stateProbs_FromRNA,
                              stateXi = t1$stateXi_FromRNA,
                              goldStateXi = t1$stateXi_FromDNA_phased)
}
#Cancer samples
cancerData = c("GOSH029", "GOSH032", "IMP672", "VHL_RCC", "RCC2")
for(nom in cancerData){
  message(nom)
  srat = readRDS(file.path('Results/fitEM',paste0(nom,'_fit_withAnnot.RDS')))
  t1 = srat@meta.data
  t1$orig.ident = gsub('_[ACGT]+$','',rownames(t1))
  t1$donor = sampTab$donor[match(t1$orig.ident,sampTab$sampleID)] 
  t1 = t1[(t1$orig.ident %in% sampTab$sampleID) & (t1$donor %in% sampTab$donor),]
  cellTab[[nom]] = data.frame(row.names = rownames(t1),
                              cellID = rownames(t1),
                              sampleID = t1$orig.ident,
                              donor = t1$donor,
                              donorGroup = 'CancerValidation',
                              annot = t1$annot,
                              stateProbXi = t1$stateProbs_FromRNA,
                              stateXi = t1$stateXi_FromRNA,
                              goldStateXi = NA)
}
cellTab = do.call(rbind,cellTab)
#Drop the non-annotated ones (except in gold standard)
cellTab = cellTab[!(is.na(cellTab$annot) & cellTab$donorGroup=='GoldStandard'),]
write.table(cellTab,'Results/cellTable.tsv',sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE)

