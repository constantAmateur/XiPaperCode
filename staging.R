#' Stages files needed from irods
#'
#' Everything else can be run whereever, but some data depends on being able to 'iget' from the farm.  So anything for which that is true needs to be run in here.  As such, it should preceed all other scripts.
#'
#' NOTE: This really needs to be reworked to function when not downloading everything and just generating the IDs of things to be mapped.  Currently it relies on files being downloaded from irods in order to list them...
#'
#' Extrenal dependencies (i.e., variables that must be defined elsewhere before running):
#' @param forceReftch Should we refetch things that exist.
#' @return Adds dataSources to variables.


#############
# Libraries #
#############

##########
# Params #
##########

dstBase = 'Data/mappedData'
liveServer = file.exists('/usr/bin/iget')

#############
# Functions #
#############

#' Fetch data from ID
#'
#' Uses the getCellRangerFromIrods.sh code to pull 10X data from irods to somewhere on lustre.
#'  
#' @param ids The IDs of the samples that need fetching.
#' @param dst The folder to save them too on lustre.  Created if doesn't exist.
#' @param srcs Folders to check for mapped cellranger data on irods.  Can be a vector of possible locations.  If NULL, checks just /archive/HCA/10X.
#' @param getAll If TRUE, gets everything get everything (including the BAM file), while FALSE gets everything but the BAMs to save space.
#' @param skipIfExists Should we skip fetching if files exist.  If \code{getAll=TRUE} checks for BAM file, otherwise just for a suitable folder.
#' @return Nothing, but hopefully the data will get downloaded.
fetchDataFromIrods = function(ids,dst,srcs=NULL,getAll=TRUE,skipIfExists=TRUE){
  #Make destination
  dir.create(dst,showWarnings=FALSE,recursive=TRUE)
  #Are we getting the cram files?
  ###getCrams = !is.null(srcs) && srcs=='UNMAPPED'
  #Check on existing files
  if(skipIfExists){
    tmp = list.dirs(dst,recursive=FALSE)
    tmp = tmp[tmp!=dst]
    #Now check each ID
    tgts = lapply(ids,function(e) grep(e,tmp,value=TRUE))
    #Check each ID for either just the directory, BAM adds too many complications...
    isFound = lengths(tgts)>0
    ###if(getAll && !getCrams){
    ###  #Check the second format to skip those that have been remapped 
    ###  isFound = sapply(tgts,function(e) any(file.exists(file.path(e,'possorted_genome_bam.bam'))))
    ###  isFound = isFound | sapply(tgts,function(e) any(file.exists(file.path(e,'outs','possorted_genome_bam.bam'))))
    ###}else{
    ###  #If we're getting crams, skip if just the folder exists...
    ###  isFound = lengths(tgts)>0
    ###}
    if(any(lengths(tgts)>1))
      warning(sprintf("Multiple matching folder found for IDs: %s",paste(tgts[lengths(tgts)>1],collapse=', ')))
    #Drop those that are already found
    if(all(isFound))
      return()
    ids = ids[!isFound]
  }
  #If we get here, we haven stuff to do
  if(!liveServer){
    stop('Files yet to be staged, but running on server incapable of staging (i.e., no iget).  Rerun on staging server first.')
  }
  #Now search for hits
  if(is.null(srcs))
    srcs = c('/archive/HCA/10X/')
  #For each, check which IDs are present 
  found = rep(FALSE,length(ids))
  toRun = list()
  for(src in srcs){
    ret = system2('ils',src,stdout=TRUE)
    tmp = sapply(ids,function(e) any(grepl(e,ret)))
    w = which(tmp & !found)
    found = found | tmp
    if(length(w)>0)
      toRun[[src]] = ids[w]
    if(all(found))
      break
  }
  if(!all(found)){
    message(sprintf("Attempting to get crams as could not find mapped data for IDs: %s",paste(ids[!found],collapse=', ')))
    #We're getting crams, do the imeta pull
    for(id in ids[!found]){
      tDir = file.path(dst,id)
      dir.create(tDir)
      system(sprintf('cd %s && %s %s',normalizePath(tDir),normalizePath('Code/irodsCramFetch.sh'),id))
    }
  }
  if(any(found)){
    message(sprintf("Downloading mapped data found for IDs: %s",paste(ids[found],collapse=', ')))
    #Now get the things that we want
    if(getAll){
      allFlag='-a'
    }else{
      allFlag=''
    }
    tmp = getwd()
    setwd(dst)
    for(i in seq_along(toRun)){
      system(sprintf('getCellrangerFromIrods.sh %s %s %s',allFlag,names(toRun[i]),paste(toRun[[i]],collapse = ' ')))
    }
    setwd(tmp)
    #Fix up any that have starsolo
    #Get the paths for each ID
    x = list.files(dst,full.names=TRUE)
    for(id in ids){
      #Get the directories matching this
      tgts = grep(id,x,value=TRUE)
      #Check each in turn
      for(tgt in tgts){
        dd = list.files(tgt)
        if('starsolo' %in% dd)
          unlink(file.path(tgt,'starsolo'),recursive=TRUE)
        #Be lazy.  Use UNIX
        if('outs' %in% dd){
          tmp = file.path(tgt,'outs')
          system(sprintf('mv %s/* %s && rmdir %s',tmp,tgt,tmp))
        }
      }
    }
  }
}



########################
# Previously used data #
########################

#Big definition of where to get what
#srcs refers to a path on irods to get the 10X bams.
#DNA BAMs can either be a lustre filepath (not preferred), an irods ID (OK I guess), or a CASM project and ID (in format project/id).
#F15 is useless,  there are like 2 read on X across all samples.
#fAdr_techComp and the two baby adrenals are good data.  But I can't get hold of the irods permissions or the cellranger BAMs.
tt = list(
                   PA513=list(IDs=c('5028STDY7152014','5028STDY7152015','5028STDY7152016','5028STDY7152017'),
                              normBAM = '1798/PD37517b',
                              srcs=NULL
                              ),
                   RCC2 = list(IDs=c('4602STDY6976426','4602STDY6976427','4602STDY6976428','4602STDY6976422','4602STDY6976423','4602STDY6976424'),
                               normBAM = '1711/PD37228b',
                               srcs=NULL,
                               notes='Last ID in manifest is garbarge, hence not included. Also renal hilum.'
                               ),
                   VHL_RCC = list(IDs=c('4602STDY6949178','4602STDY6949179','4602STDY6949180','4602STDY6949181','4602STDY6949182','4602STDY6949183'),
                                  normBAM = '1711/PD36793b',
                                  srcs=NULL
                                  ),
                   Wilms2 = list(IDs = c('4602STDY7018926','4602STDY7018923','4602STDY7018924','4602STDY7018925','4602STDY7018920','4602STDY7018921','4602STDY7018922'),
                                 normBAM = '1711/PD37276c',
                                 dadBAM = '1711/PD37286b',
                                 mumBAM = '1711/PD37287b',
                                 srcs=NULL,
                                 notes='Skipping A&V and renal pelvis as they suck.'
                                 ),
                   Wilms3 = list(IDs = c('4602STDY7090428','4602STDY7090429','4602STDY7090425','4602STDY7090426','4602STDY7090427','4602STDY7090431'),
                                 normBAM = '1711/PD37272c',
                                 mumBAM = 'Data/stagedDNA/PD37288b.v1.sample.dupmarked.bam',
                                 #mumBAM = '1711/PD37288b',
                                 dadBAM = 'Data/stagedDNA/PD37289b.v1.sample.dupmarked.bam',
                                 #dadBAM = '1711/PD37289b',
                                 srcs=NULL,
                                 Notes = 'Dropping last ID as no coverage.'
                                 ),
                   GOSH028 = list(IDs = c('CG_SB_NB8113359','CG_SB_NB8113360','CG_SB_NB8113361','CG_SB_NB8368298','CG_SB_NB8368299'),
                                  normBAM = '1607/PD46696b',
                                  srcs=c('/seq/30855/cellranger/','/seq/illumina/runs/31/31667/cellranger')
                                  ),
                   GOSH026 = list(IDs = c('4602STDY8004990','4602STDY8004991','4602STDY8004992'),
                                  normBAM = '1607/PD46692b',
                                  srcs= '/seq/illumina/cellranger/'
                                  ),
                   GOSH029 = list(IDs = c('CG_SB_NB8350518','CG_SB_NB8350519','CG_SB_NB8350520','CG_SB_NB8368300','CG_SB_NB8368301'),
                                  normBAM = '1985/PD47704b',
                                  srcs = c('/seq/31387/cellranger','/seq/illumina/runs/31/31667/cellranger')
                                  ),
                   IMP672 = list(IDs = c('CG_SB_NB8350521','CG_SB_NB8350522','CG_SB_NB8350523','CG_SB_NB8368876','CG_SB_NB8368877'),
                                 normBAM = '1985/PD48777b',
                                 srcs= c('/seq/31387/cellranger/','/seq/illumina/runs/31/31667/cellranger/')
                                 ),
                   #GOSH016 = list(IDs = c('4602STDY7685338','4602STDY7685339','4602STDY7685335','4602STDY7685336','4602STDY7685337'),
                   GOSH016 = list(IDs = c('4602STDY7685338','4602STDY7685339','4602STDY7685340','4602STDY7685341','4602STDY7685342'),
                                  normBAM = '1985/PD42187b',
                                  srcs = '/seq/illumina/cellranger/'
                                  ),
                   GOSH023 = list(IDs = c('4602STDY7843576','4602STDY7843577','4602STDY7843578'),
                                  normBAM = '1607/PD42752d',
                                  srcs = '/seq/28881/cellranger/'
                                  ),
                   GOSH025 = list(IDs = c('4602STDY8004894','4602STDY8004902','4602STDY8004910'),
                                  normBAM = 'Data/stagedDNA/PD46693b.v1.sample.dupmarked.bam',
                                  #normBAM = '1607/PD46693b',
                                  mumBAM = 'Data/stagedDNA/PD46695b.v1.sample.dupmarked.bam',
                                  #mumBAM = '1607/PD46695b',
                                  dadBAM = 'Data/stagedDNA/PD46694b.v1.sample.dupmarked.bam',
                                  #dadBAM = '1607/PD46694b',
                                  srcs= '/seq/30063/cellranger/'
                                  ),
                   GOSH032 = list(IDs = c('CG_SB_NB8715409','CG_SB_NB8715410','CG_SB_NB8715411'),
                                  normBAM = '1607/PD47706b',
                                  srcs = '/seq/illumina/runs/33/33315/cellranger/'
                                  ),
                   GOSH033 = list(IDs = c('CG_SB_NB8715412','CG_SB_NB8715413','CG_SB_NB8715414'),
                                  normBAM = '1607/PD47708b',
                                  srcs = '/seq/illumina/runs/33/33315/cellranger'
                                  ),
                   InfALL1 = list(IDs = c('4602STDY7920960','4602STDY7920961','4602STDY7920962','4602STDY7920963','4602STDY7920964'),
                                  srcs = '/seq/29622/cellranger/',
                                  notes = "This is multiplexed with a boy and girl together."
                                  ),
                   CongAML = list(IDs = c('4602STDY7920965','4602STDY7920966'),
                                  srcs = '/seq/29622/cellranger/'
                                  ),
                   ClassicalHodgkins = list(IDs =c('4710STDY6879626','4710STDY6879627','4710STDY6879628','4710STDY6879629'),
                                            srcs = '/seq/illumina/cellranger'
                                            ),
                   fAdr_late = list(IDs = c('WSSS_F_Adr8768489','WSSS_F_Adr8768490'),
                                    srcs = '/seq/illumina/runs/33/33701/cellranger/'
                                    ),
                   IMP582 = list(IDs = c('4602STDY7733090','4602STDY7733091','4602STDY7733092'),
                                 srcs = '/seq/27838/cellranger/'
                                 ),
                   KidTransplant = list(IDs = c('4602STDY6949190','4602STDY6949184','4602STDY6949185','4602STDY6949187','4602STDY6949188','4602STDY6949186','4602STDY6949189'),
                                         srcs=NULL
                                         ),
                   F19 = list(IDs = c("FCAImmP7241240", "FCAImmP7241241","FCAImmP7241242", "FCAImmP7241243","FCA7167219", "FCA7167223", "FCA7167224", "FCA7167231"),
                              normBAM = 'Data/stagedDNA/PD38222b.v1.sample.dupmarked.bam',
                              #normBAM = "1836/PD38222b",
                              mumBAM = 'Data/stagedDNA/PD38223b.v1.sample.dupmarked.bam',
                              #mumBAM = '1836/PD38223b',
                              srcs = NULL),
                   F20 = list(IDs = c("FCA7167230", "FCA7167221", "FCA7167222","FCA7167226", "FCA7167225", "FCA7167232"),
                              srcs = NULL,
                              notes = "Fetus is female"),
                   F26 = list(IDs = c("FCA7196222", "FCA7196223","FCA7196228"),#From the genotyping, have no idea what FCA7196221 is
                              srcs = NULL,
                              notes = "Fetus is female.  Not in Roser's manifest."),
                   F37 = list(IDs = c("FCA7474067", "FCA7474068", "FCA7474069"),
                              srcs = NULL,
                              notes = "Fetus is female"),
                   F40 = list(IDs = c("FCA7511881", "FCA7511882", "FCA7511883","FCA7511884", "FCA7511885", "FCA7511886"),
                              srcs = NULL,
                              notes = "Fetus is female"),
                   Hrv43 = list(IDs = c('Pla_HDBR10142767','Pla_HDBR10142768','Pla_HDBR10701667','Pla_HDBR10701666'),#Multiomes 'Pla_HDBR10142863','Pla_HDBR10142864','Pla_HDBR10142865'
                                #visiumIDs = c('Pla_HDBR9518710'),
                                mapType=rep(c('scRNA','snRNA'),c(3,1)),
                                srcs = NULL,
                                notes = "Roser's data."),
                   Hrv46 = list(IDs = c('Pla_HDBR10142769','Pla_HDBR10142770','Pla_HDBR10701668'),
                                #visiumIDs=c('Pla_HDBR9518711'),
                                srcs = NULL,
                                notes = "Roser's data."),
                   H2 = list(IDs = c('Pla_HDBR8624430','Pla_HDBR8624431'),
                             srcs = NULL,
                             notes = "Roser's data."),
                   H7and9 = list(IDs = c('Pla_HDBR8768477','Pla_HDBR8715512','Pla_HDBR8715514'),
                                 srcs = NULL,
                                 notes = "Rosers's data.  Two samples mutliplexed together."),
                   Hrv98 = list(IDs = 'Pla_HDBR10917730',
                                srcs = NULL,
                                notes = "Roser's data."),
                   Hrv99 = list(IDs = 'Pla_HDBR10917731',
                                srcs = NULL,
                                notes = "Roser's data."),
                   Hrv100 = list(IDs = 'Pla_HDBR10917733',
                                 srcs = NULL,
                                 notes = "Roser's data."),
                   P17 = list(IDs=c('6044STDY8711544','WSSS_PLA8810748','WSSS_PLA8810749','Pla_Camb10691974'),
                              #visiumIDs=c('WS_PLA_S9101768','WS_PLA_S9101771'),
                              srcs = NULL,
                              mapType='snRNA',
                              notes = "Roser's data."),
                   P13 = list(IDs = c('WSSS_PLA8764121','WSSS_PLA8764122','WSSS_PLA8810750','WSSS_PLA8810751','Pla_Camb10691970','Pla_Camb10691971','WS_PLA_S9101764','WS_PLA_S9101765','WS_PLA_S9101766','WS_PLA_S9101767','Pla_Camb9779195','Pla_Camb9518737'),#These ones are the multiome, need custom code ,'Pla_Camb10714919','Pla_Camb10714920'),
                              srcs = c('/seq/illumina/runs/34/34914/spaceranger/','/seq/illumina/runs/34/34882/spaceranger/','/seq/illumina/runs/36/36173/spaceranger/','/seq/illumina/runs/36/36645/spaceranger/'),
                              mapType=rep(c('snRNA','spatial'),c(6,6)),
                              notes = "Roser's data."),
                   P14 = list(IDs = c('WSSS_PLA8811068','WSSS_PLA8811069','WSSS_PLA8811070','Pla_Camb10691972'),#Multiome 'Pla_Camb10714918'
                              #visiumIDs=c('WS_PLA_S9101769','WS_PLA_S9101770'),
                              srcs = NULL,
                              mapType = 'snRNA',
                              notes = "Roser's data."),
                   P34 = list(IDs='Pla_Camb10691975',
                              srcs = NULL,
                              mapType = 'snRNA',
                              notes = "Roser's data.  Typo on sheet for this one?"),
                   WholeEmbryo1 = list(IDs=c("WS_wEMB10202336", "WS_wEMB10202337", "WS_wEMB10202338", "WS_wEMB10202339", "WS_wEMB10202340", "WS_wEMB10202341", "WS_wEMB10202342", "WS_wEMB10202343", "WS_wEMB10202344", "WS_wEMB10202345", "WS_wEMB10202346", "WS_wEMB10202347", "WS_wEMB10202348", "WS_wEMB10202349", "WS_wEMB10202350", "WS_wEMB10202351", "WS_wEMB10202352", "WS_wEMB10202353", "WS_wEMB10202354", "WS_wEMB10202355", "WS_wEMB10202356", "WS_wEMB10202357", "WS_wEMB10202358", "WS_wEMB10202359", "WS_wEMB10202360", "WS_wEMB10202361", "WS_wEMB10202362", "WS_wEMB10202363", "WS_wEMB10202364", "WS_wEMB10202365", "WS_wEMB10202366", "WS_wEMB10202367", "WS_wEMB10202368", "WS_wEMB10202369", "WS_wEMB10202370", "WS_wEMB10202371", "WS_wEMB10202372", "WS_wEMB10202373", "WS_wEMB10202374", "WS_wEMB10202375", "WS_wEMB10202376", "WS_wEMB10202377", "WS_wEMB10202378", "WS_wEMB10202379", "WS_wEMB10202380", "WS_wEMB11031919", "WS_wEMB11031920", "WS_wEMB11031921", "WS_wEMB11031922", "WS_wEMB11031923", "WS_wEMB11031924", "WS_wEMB11031925", "WS_wEMB11031926", "WS_wEMB11031927", "WS_wEMB11031928", "WS_wEMB11031929", "WS_wEMB11031930", "WS_wEMB11031931", "WS_wEMB11031932", "WS_wEMB11031933", "WS_wEMB11031934", "WS_wEMB11031935", "WS_wEMB11031936", "WS_wEMB11031937", "WS_wEMB11031938", "WS_wEMB11031939", "WS_wEMB11031940", "WS_wEMB11031941", "WS_wEMB11031942", "WS_wEMB11031943", "WS_wEMB11031944", "WS_wEMB11031945", "WS_wEMB11031946", "WS_wEMB11031947", "WS_wEMB11031948", "WS_wEMB11031949", "WS_wEMB11031950", "WS_wEMB11031951", "WS_wEMB11031952", "WS_wEMB11031953", "WS_wEMB11031954", "WS_wEMB11031955", "WS_wEMB11031956", "WS_wEMB11031957", "WS_wEMB11031958", "WS_wEMB11031959", "WS_wEMB11031960", "WS_wEMB11031961", "WS_wEMB11031962", "WS_wEMB11031963", "WS_wEMB11031964", "WS_wEMB11031965", "WS_wEMB11031966"),
                              srcs = NULL,
                              mapType = 'scRNA',
                              notes = "First whole embryo")
)
#Fill the default maptype and expand to correct length
for(nom in names(tt)){
  if(is.null(tt[[nom]]$mapType))
    tt[[nom]]$mapType = rep('scRNA',length(tt[[nom]]$IDs))
  if(length(tt[[nom]]$mapType)==1)
    tt[[nom]]$mapType = rep(tt[[nom]]$mapType,length(tt[[nom]]$IDs))
}
#Building on Chenqu Suo and Emma's annotation.  Starting point is sample metadata PAN.A01.v01.20210429.sample_metadata.csv, irods sources Pan_Fetal_sample_locations.txt, and cell annotations PAN.A01.v01.entire_data_normalised_log.20210429.full_obs.annotated.clean.csv.  All are saved in Data/
mDat = read.table('Data/PAN.A01.v01.20210429.sample_metadata.csv',sep=',',header=TRUE)
#Fill in sex
mDat$sex[mDat$donor=='F72']='female'
mDat$sex[mDat$donor=='F78']='male'
#This one is wrong, so fix it
mDat$sex[mDat$donor=='F64']='male'
#Fill in irods path
tmp = read.table('Data/Pan_Fetal_sample_locations.txt',sep=',',header=FALSE)
mDat$irodsPath = unlist(lapply(mDat$file,function(e) grep(e,tmp[,1],value=TRUE)))
#Drop all non-females
mDat = mDat[mDat$sex=='female',]
#Add or update
for(nom in mDat$donor){
  tmp = mDat[mDat$donor==nom,]
  ids = tmp$X
  srcs = c('/archive/HCA/10X',dirname(tmp$irodsPath[grepl('Human_colon_',tmp$X)]))
  if(nom %in% names(tt)){
    tt[[nom]]$srcs = unique(c(tt[[nom]]$srcs,srcs))
    idsToAdd = setdiff(ids,tt[[nom]]$IDs)
    tt[[nom]]$IDs = c(tt[[nom]]$IDs,idsToAdd)
    #Update maptype
    tt[[nom]]$mapType = c(tt[[nom]]$mapType,rep('scRNA',length(idsToAdd)))
    #tt[[nom]]$IDs = unique(c(tt[[nom]]$IDs,ids))
  }else{
    tt[[nom]] = list(IDs=ids,srcs=srcs,mapType=rep('scRNA',length(ids)))
  }
  #Make sure /archive is at the end of srcs
  srcs = tt[[nom]]$srcs
  if('/archive/HCA/10X' %in% srcs)
    srcs = c(srcs[which(srcs!='/archive/HCA/10X')],'/archive/HCA/10X')
  tt[[nom]]$srcs = srcs
}
#Some cleanup of duplicate/redundant/unwanted samples.  These will be deleted as created.
toDel = file.path(dstBase,c(
                  'GOSH028/cellranger302_count_31667_CG_SB_NB8368298_GRCh38-1_2_0',
                  'GOSH028/cellranger302_count_31667_CG_SB_NB8368299_GRCh38-1_2_0',
                  'GOSH029/cellranger302_count_31667_CG_SB_NB8368300_GRCh38-1_2_0',
                  'GOSH029/cellranger302_count_31667_CG_SB_NB8368301_GRCh38-1_2_0',
                  'IMP672/cellranger302_count_31387_CG_SB_NB8350521_GRCh38-3_0_0_and_H19',
                  'IMP672/cellranger302_count_31387_CG_SB_NB8350522_GRCh38-3_0_0_and_H19',
                  'IMP672/cellranger302_count_31387_CG_SB_NB8350523_GRCh38-3_0_0_and_H19',
                  'IMP672/cellranger302_count_31667_CG_SB_NB8368876_GRCh38-1_2_0',
                  'IMP672/cellranger302_count_31667_CG_SB_NB8368876_GRCh38-3_0_0_and_H19',
                  'IMP672/cellranger302_count_31667_CG_SB_NB8368877_GRCh38-1_2_0',
                  'IMP672/cellranger302_count_31667_CG_SB_NB8368877_GRCh38-3_0_0_and_H19',
                  'GOSH023/cellranger302_count_28881_4602STDY7843576_GRCh38-3_0_0',
                  'GOSH023/cellranger302_count_28881_4602STDY7843577_GRCh38-3_0_0',
                  'GOSH023/cellranger302_count_28881_4602STDY7843578_GRCh38-3_0_0',
                  'NodularPredominantHodgkins/cellranger131_count_4710STDY7018928_GRCh38',
                  'NodularPredominantHodgkins/cellranger131_count_4710STDY7018929_GRCh38',
                  'NodularPredominantHodgkins/cellranger131_count_4710STDY7018930_GRCh38',
                  'NodularPredominantHodgkins/cellranger131_count_4710STDY7018931_GRCh38',
                  'fAdr_19wk/cellranger302_count_32644_WSSS_F_Adr8710632_GRCh38-3_0_0_and_H19',
                  'fAdr_19wk/cellranger302_count_32644_WSSS_F_Adr8710633_GRCh38-3_0_0_and_H19',
                  'fAdr_19wk/cellranger302_count_32644_WSSS_F_Adr8710634_GRCh38-3_0_0_and_H19',
                  'fAdr_19wk/cellranger302_count_32644_WSSS_F_Adr8710635_GRCh38-3_0_0_and_H19',
                  'GOSH012/cellranger302_count_26952_4602STDY7654777_GRCh38-3_0_0_and_H19',
                  'GOSH012/cellranger302_count_26952_4602STDY7654778_GRCh38-3_0_0_and_H19',
                  'GOSH012/cellranger302_count_26952_4602STDY7654779_GRCh38-3_0_0_and_H19',
                  'GOSH012/cellranger302_count_26952_4602STDY7654780_GRCh38-3_0_0_and_H19',
                  'GOSH012/cellranger302_count_26952_4602STDY7654781_GRCh38-3_0_0_and_H19',
                  'GOSH012/cellranger302_count_26952_4602STDY7654782_GRCh38-3_0_0_and_H19',
                  'GOSH016/cellranger202_count_4602STDY7685335_GRCh38-1_2_0',
                  'GOSH016/cellranger202_count_4602STDY7685336_GRCh38-1_2_0',
                  'GOSH016/cellranger202_count_4602STDY7685337_GRCh38-1_2_0',
                  'GOSH016/cellranger202_count_4602STDY7685338_GRCh38-1_2_0',
                  'GOSH016/cellranger202_count_4602STDY7685339_GRCh38-1_2_0',
                  'GOSH016/cellranger302_count4602STDY7685335_GRCh38-3_0_0_and_H19',
                  'GOSH016/cellranger302_count4602STDY7685336_GRCh38-3_0_0_and_H19',
                  'GOSH016/cellranger302_count4602STDY7685337_GRCh38-3_0_0_and_H19',
                  'IMP582/cellranger302_count_27838_4602STDY7733090_GRCh38-3_0_0_and_H19',
                  'IMP582/cellranger302_count_27838_4602STDY7733091_GRCh38-3_0_0_and_H19',
                  'IMP582/cellranger302_count_27838_4602STDY7733092_GRCh38-3_0_0_and_H19',
                  'F72/cellranger210_count_31512_Human_colon_16S8159182_GRCh38-1_2_0',
                  'F72/cellranger210_count_31512_Human_colon_16S8159183_GRCh38-1_2_0',
                  'F72/cellranger210_count_31512_Human_colon_16S8159184_GRCh38-1_2_0',
                  'F72/cellranger210_count_31512_Human_colon_16S8159185_GRCh38-1_2_0',
                  'F72/cellranger210_count_31512_Human_colon_16S8159186_GRCh38-1_2_0',
                  'F73/cellranger210_count_31512_Human_colon_16S8159187_GRCh38-1_2_0',
                  'F73/cellranger210_count_31512_Human_colon_16S8159188_GRCh38-1_2_0',
                  'F73/cellranger210_count_31512_Human_colon_16S8159189_GRCh38-1_2_0',
                  'F73/cellranger210_count_31512_Human_colon_16S8159190_GRCh38-1_2_0',
                  'P13/spaceranger110_count_34914_WS_PLA_S9101764_GRCh38-3_0_0_premrna',
                  'P13/spaceranger110_count_34914_WS_PLA_S9101765_GRCh38-3_0_0_premrna',
                  'P13/spaceranger110_count_34914_WS_PLA_S9101766_GRCh38-3_0_0_premrna',
                  'P13/spaceranger110_count_34914_WS_PLA_S9101767_GRCh38-3_0_0_premrna',
                  'P13/spaceranger110_count_36173_Pla_Camb9518737_GRCh38-3_0_0_premrna',
                  'P13/spaceranger110_count_36645_Pla_Camb9779195_GRCh38-3_0_0_premrna'
                  ))
#Get data
for(nom in names(tt)){
  if(localMapping){
    message(sprintf("Getting data for sample %s",nom))
    tgtDir=file.path(dstBase,nom)
    tt[[nom]]$lustrePath = tgtDir
    ids = tt[[nom]]$IDs
    fetchDataFromIrods(ids,tgtDir,srcs=tt[[nom]]$srcs,getAll=TRUE,skipIfExists=!forceRefetch)
    #Delete in this highly inefficient way
    unlink(toDel,recursive=TRUE)
    #Get the paths for each ID
    x = list.files(tgtDir,full.names=TRUE)
    dds = lapply(ids,function(e) grep(e,x,value=TRUE))
    dds = lapply(dds,function(e){
                   if(file.exists(file.path(e,'outs','possorted_genome_bam.bam'))){
                     return(file.path(e,'outs'))
                   }else{
                     return(e)
                   }})
    names(dds) = ids
    if(!all(lengths(dds))==1)
      stop("Duplicate entries for the same ID")
    #Should now be one per thing, so flatten
    #x = unlist(dds)
    #names(x) = make.unique(rep(names(dds),lengths(dds)))
    tt[[nom]]$lustrePathByID = unlist(dds)
  }
}
#Now record exactly where the DNA BAMs are, they don't need any remapping
projBase='/nfs/cancer_ref01/nst_links/live/'
for(nom in names(tt)){
  for(bamNom in c('normBAM','mumBAM','dadBAM')){
    nb = tt[[nom]][[bamNom]]
    if(!is.null(nb)){
      if(file.exists(nb)){
        src = nb
      }else{
        src = grep('sample.dupmarked.bam$',list.files(file.path(projBase,nb),full.names = TRUE),value=TRUE)
        if(length(src)!=1)
          warning(sprintf("Could not find %s file %s in expected location.",bamNom,nb))
        #dst = file.path(dstBase,nom,paste0(bamNom,'.bam'))
        #suppressWarnings(file.symlink(src,dst))
        #suppressWarnings(file.symlink(paste0(src,'.bai'),paste0(dst,'.bai')))
        #src = normalizePath(src)
      }
      tt[[nom]][[paste0(bamNom,'_lustre')]] = src
    }
  }
}
dataSources = tt


if(liveServer)
  stop("Staging finished on farm.  Swap to processing server")


