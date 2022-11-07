#' Maps data and prepares it for model fit 
#'
#' Gets data from whereever, makes sure it's mapped to cellranger, and records metadata about everything in a variable \code{dataSources}
#'
#' External variables:
#' @param dataSources Defined by running staging.R 
#' @return scDat, a data.frame with the successfully mapped things and their locations.  dataSources, a more verbose version that contains full metadata that can be dug into.


#############
# Libraries #
#############

##########
# Params #
##########

dstBase = 'Data/mappedData'
stageDir = 'Data/singleCellRNA' #Final directory everything gets symlinked to
refGenome = normalizePath('~/src/refdata-gex-GRCh38-2020-A')
dcpUpdate = FALSE #If the DCP has re-indexed everything the manifests need to be redownloaded.  Which is a pain.

#############
# Functions #
#############

#' Format HCA metadata
#'
#' HCA metadata is provided in a somewhat standard format.  Check the assumptions we're making and try and auto-generate cellranger compatible file names.  The default assumption here is that the "sequencing input", given by sequencing_input.biomaterial_core.biomaterial_id, is the lowest level of grouping that signifies a collection of files that should be processed together.
#'
#' Note that for almost all metadata sequencing_input.biomaterial_core.biomaterial_id and cell_suspension.biomaterial_core.biomaterial_id are the same thing.  But there are a subset where they are not and the sequencing_input seems the more dependable (and logically is the lower level) of the two.
#'
#' @param mDat data.frame containing HCA metadata.
#' @param dataName Unique name for this data.
#' @param useDonorName Use the \code{donorName} string to name donors.  If false, name numerically.
#' @param channelName What to use to name the channels.
#' @param donorName What to use to name the channels.
#' @return Updated version of mDat
formatMetadataHCA = function(mDat,dataName,useDonorName=FALSE,channelName = mDat$sequencing_input.biomaterial_core.biomaterial_id,donorName = mDat$donor_organism.biomaterial_core.biomaterial_id){
  #Decide on namings 
  if(useDonorName){
    humanDonorName = sprintf('%s_%s',dataName,donorName)
  }else{
    humanDonorName = sprintf('%s%03d',dataName,match(donorName,unique(donorName)))
  }
  mDat$humanDonorName = humanDonorName
  mDat$groupName = channelName
  #We only want the fastq files
  mDat = mDat[mDat$file_type == 'sequence_file' &
              mDat$file_format %in% c('fastq','fastq.gz','fq.gz') &
              mDat$read_index %in% c('read1','read2','index1'),]
  message('Found library prep methods')
  print(table(mDat$library_preparation_protocol.library_construction_approach))
  message('Derived from')
  print(table(mDat$library_preparation_protocol.nucleic_acid_source))
  message("Sequencing input types")
  print(table(mDat$sequencing_input_type))
  #The basic expectation I'm making is that per "sequencing input", there is one donor, specimen, sample, library prep protocol, cell suspension
  x = split(mDat,mDat$groupName)
  #First check whatever we're using for channel/donor is unique
  if(!all(lengths(lapply(split(mDat$groupName,mDat$sequencing_input.biomaterial_core.biomaterial_id),unique))==1))
    stop("Specified channelName is not unique by channel!")
  if(!all(lengths(lapply(split(mDat$humanDonorName,mDat$sequencing_input.biomaterial_core.biomaterial_id),unique))==1))
    stop("Specified donorName is not unique by channel!")
  if(!all(sapply(x,function(e) length(unique(e$donor_organism.biomaterial_core.biomaterial_id)))==1))
    stop("Some sequencing inputs have more than one donor")
  if(!all(sapply(x,function(e) length(unique(e$specimen_from_organism.provenance.document_id)))==1))
    stop("Some sequencing inputs have more than one specimen")
  if(!all(sapply(x,function(e) length(unique(e$sample.biomaterial_core.biomaterial_id)))==1))
    stop("Some sequencing inputs have more than one sample")
  if(!all(sapply(x,function(e) length(unique(e$library_preparation_protocol.library_construction_approach)))==1))
    stop("Some sequencing inputs have more than one library prep protocol")
  if(!all(sapply(x,function(e) length(unique(e$cell_suspension.biomaterial_core.biomaterial_id)))==1))
    stop("Some sequencing inputs have more than one cell suspension")
  #Now check if we can uniquely resolve file name for all
  #First assumption is that within each, sequencing_process.provenance.document_id, there's one set of reads.
  if(!all(sapply(split(mDat$read_index,mDat$sequencing_process.provenance.document_id),function(e) all(table(e)==1)))){
    warning("More than one group of reads per sequencing process.  You'll have to manually build file names.")
    return(mDat)
  }else{
    message("Auto constructing cellranger compatible filenames")
    #Construct name if we can
    for(nom in names(x)){
      tgts = x[[nom]]
      tgts$fileName = sprintf('%s_S001_L%03d_%s_001.fastq.gz',
                              nom,
                              match(tgts$sequencing_process.provenance.document_id,unique(tgts$sequencing_process.provenance.document_id)),
                              c(read1='R1',read2='R2',index1='I1')[tgts$read_index]
                              )
      x[[nom]] = tgts
    }
    mDat = do.call(rbind,x)
  }
  return(mDat)
}

#' HCA download
#'
#' Download fastq files from HCA using curl.
#'
#' @param url Usually the column named file_url
#' @param fileName The derived file name that will work with cellranger
#' @param tmpDir Directory where we will download to
#' @param sha256 Checksums for files to test.
#' @param retries How many times to retry if checksum fails?
downloadHCA = function(urls,fileName,tmpDir,sha256=NULL,retries=3){
  dir.create(tmpDir,showWarnings=FALSE)
  #Fix dcp urls
  urls = gsub('\\?catalog\\=dcp11\\&','?catalog=dcp12&',urls)
  base='--create-dirs
  --compressed
  --location
  --globoff
  --fail
  --fail-early
  --continue-at -
  --write-out "Downloading to: %{filename_effective}"\n\n'
  tFile = tempfile()
  write(paste0(base,paste(sprintf('url="%s"\n',urls),sprintf('output="%s"\n',fileName),collapse='')),tFile)
  system(sprintf('cd %s && cat %s|curl --config -',tmpDir,tFile),ignore.stdout=TRUE,ignore.stderr=TRUE)
  #At least check all the files exist!
  outFiles = file.path(tmpDir,fileName)
  checkPass=all(file.exists(outFiles))
  #Also check sha256 if given
  if(checkPass && !is.null(sha256)){
    check = system(sprintf('sha256sum %s',paste(outFiles,collapse=' ')),intern=TRUE)
    check = data.frame(t(do.call(cbind,strsplit(check,'  ',fixed=TRUE))))
    checkPass = checkPass && all(check[,1]==sha256)
  }
  #Retry if failed
  if(!checkPass){
    if(retries<1)
      stop('Checksum mismatch or missing files despite retrying.')
    #Beautiful recursion
    downloadHCA(urls,fileName,tmpDir,sha256=sha256,retries=retries-1)
  }
  #All done, cleanup
  message("Files downloaded and integrity checked.")
  unlink(tFile)
}

#' Run spaceranger
#'
#' Runs spaceranger on temporary fastq files, then move the results and deletes the temporary files.
#'
#' 
runSpaceranger = function(tmpBase,useIntrons=FALSE,nParallel=12,mem=64){
  cmd = sprintf('cd %s && cellranger count --id=map_%s --transcriptome=%s --fastqs=%s --expect-cells=3000 --localcores=%d --localmem=%d',tmpBase,grp,refGenome,grp,nParallel,mem)
}

#' Run cellranger
#'
#' Runs cellranger on temporary fastq files, then moves the results and deletes the temporary files.
#'
#' @param tmpBase Directory with folder containing fastq files.
#' @param grp Name of folder containing the fastqs in the expected format.
#' @param tgtDir Where to move things to when done.
#' @param useIntrons Should we include introns in the count?  Usually only for snRNA
#' @param nParallel How many cores to use.
#' @param mem How much memory to use.
runCellranger = function(tmpBase,grp,tgtDir,useIntrons=FALSE,nParallel=12,mem=64){
  cmd = sprintf('cd %s && cellranger count --id=map_%s --transcriptome=%s --fastqs=%s --expect-cells=3000 --localcores=%d --localmem=%d',tmpBase,grp,refGenome,grp,nParallel,mem)
  if(useIntrons)
    cmd = paste0(cmd,' --include-introns')
  system(cmd)
  #Check it worked
  tmp = file.path(tmpBase,paste0('map_',grp))
  w = grep('Pipestance completed successfully!',readLines(file.path(tmp,'_log')))
  if(length(w)==0 || !file.exists(file.path(tmp,'outs','possorted_genome_bam.bam'))){
    #Make it a warning instead of a stop
    warning(sprintf("Something went wrong mapping %s with cellranger.",grp))
  }else{
    #Move results if it worked.
    dir.create(tgtDir,showWarnings=FALSE)
    system(sprintf('mv %s/* %s',file.path(tmpBase,paste0('map_',grp)),tgtDir),ignore.stderr=TRUE,ignore.stdout=TRUE)
  }
  #Delete old files
  unlink(file.path(tmpBase,grp),recursive=TRUE)
}

#' Process HCA
#'
#' Run the standard HCA processing pipeline.  If the column 'processFail' exists in the metadata, any sample marked as TRUE here will not be processed.
#'
#' @param mDat data.frame of HCA metadata with extra columns fileName, humanDonorName, and groupName
#' @param ... Passed to runCellranger. 
processHCA = function(mDat,...){
  tmpBase=tempdir()
  #Run cellranger on everything
  for(sampName in unique(mDat$humanDonorName)){
    mDatDonor = mDat[mDat$humanDonorName==sampName,]
    #Drop anything that will fail
    if(!is.null(mDat$processFail))
      mDatDonor = mDatDonor[!mDatDonor$processFail,]
    for(grp in unique(mDatDonor$groupName)){
      tmpDir = file.path(tmpBase,grp)
      tgtDir = file.path(normalizePath(dstBase),sampName,grp)
      dir.create(file.path(dstBase,sampName),showWarnings=FALSE)
      lockFile=file.path(tgtDir,'00_lock')
      if(!file.exists(lockFile) && !dir.exists(file.path(tgtDir,'outs'))){
        #Create the lock file to show we're doing this now
        dir.create(tgtDir)
        system(sprintf('touch %s',lockFile))
        tgts = mDatDonor[mDatDonor$groupName==grp,]
        message(sprintf("Downloading %s",grp))
        #Download
        downloadHCA(tgts$file_url,tgts$fileName,tmpDir,tgts$file_sha256)
        #Map and copy
        message(sprintf("Mapping %s",grp))
        runCellranger(tmpBase,grp,tgtDir,nParallel=nParallel,mem=localMem,...)
        #We're done, remove the lock
        unlink(lockFile)
      }
    }
  }
}

#' Reprocess from local files
#'
#' Remap files that have been downloaded onto lustre from irods.  The original irods files will be removed as they can always be re-downloaded from irods if needed.  The starting point for this is a cellranger BAM that will be converted to fastq and remapped.
#'
#' There's a lockfile system to make sure we don't try and process the same file twice simultaneously.
#'
#' @param dat A dataSources entry, with a "lustrePathByID" entry.
#' @param ... Passed to runCellranger.
#' @return The modified dataSource with lustrePathByID updated.
reprocessLocal = function(dat,...){
  #Check we have the expected files
  if(!(!is.null(dat$lustrePathByID) && all(dir.exists(dat$lustrePathByID))))
    stop("Input data not in expected format")
  tmpBase = tempdir()
  for(sampName in names(dat$lustrePathByID)){
    srcDir = dat$lustrePathByID[sampName]
    #Store the existing files to delete them on success
    trashPile = list.files(srcDir,full.names=TRUE)
    fqDir = file.path(tmpBase,sampName)
    #Check if someone else has done it and we need to update and go on
    if(file.exists(file.path(srcDir,'outs','possorted_genome_bam.bam'))){
      dat$lustrePathByID[sampName] = file.path(srcDir,'outs')
      next
    }
    #Skip if already done.  When done, the path will end in 'outs', whereas the unremapped ones don't have the outs directory
    if(basename(srcDir)!='outs'){
      #Check for a lock file that suggests it's being processed.
      lockFile=file.path(srcDir,'00_lock')
      if(file.exists(lockFile)){
        message(sprintf("Sample %s is being processed elsewhere.  Skipping...",sampName))
        next
      }
      #Make the lock file
      system(sprintf('touch %s',lockFile))
      trashPile = c(trashPile,lockFile)
      #First dump it to local file
      message(sprintf("Extracting fastqs for %s",sampName))
      if(any(grepl('\\.cram$',trashPile))){
        dir.create(fqDir)
        #We're reprocessing stuff from irods
        crams = grep('\\.cram$',trashPile,value=TRUE)
        for(cram in crams)
          system(sprintf('cp %s %s',cram,fqDir))
        system(sprintf('parallel -j %d bash %s ::: %s',nParallel,normalizePath('Code/cram2fastq.sh'),paste(file.path(fqDir,basename(crams)),collapse=' ')))
        #Do the renaming
        lanes = sapply(strsplit(basename(crams),'#'),`[`,1)
        lCount=1
        for(lane in unique(lanes)){
          sCount=1
          for(a in basename(crams[lanes==lane])){
            for(b in c('I1','R1','R2')){
              file.rename(file.path(fqDir,sprintf('%s_%s_001.fastq.gz',a,b)),file.path(fqDir,sprintf('%s_S%d_L%03d_%s_001.fastq.gz',sampName,sCount,lCount,b)))
            }
            sCount = sCount+1
          }
          lCount=lCount+1
        }
        unlink(file.path(fqDir,basename(crams)))
      }else{
        cmd = sprintf('bamtofastq --nthreads %d %s %s',min(nParallel,8),file.path(srcDir,'possorted_genome_bam.bam'),fqDir)
        system(cmd)
      }
      #Now run cellranger from this set of fastqs 
      message(sprintf("Mapping %s",sampName))
      useIntrons = dat$mapType[match(sampName,names(dat$lustrePathByID))]=='snRNA'
      if(dat$mapType[match(sampName,names(dat$lustrePathByID))]=='spatial'){
        runSpaceranger()
      }else{
        runCellranger(tmpBase,sampName,srcDir,nParallel=nParallel,mem=localMem,useIntrons=useIntrons,...)
      }
      #Do a final check that everything has copied happily
      if(!file.exists(file.path(srcDir,'outs','possorted_genome_bam.bam')))
        stop("Something went wrong copying files.")
      #Delete the source material
      unlink(trashPile,recursive=TRUE)
      #Update lustrePath
      dat$lustrePathByID[sampName] = file.path(srcDir,'outs')
    }
  }
  return(dat)
}


#' Standardise metadata
#'
#' Convert metadata into standard list format with minimum information.  Each list entry should be the donor name, with at a minimum \code{IDs}, the groupName, \code{lustrePath}, and lustrePathByID, which is the path to the 'outs' folder for each channel, grouped by ID.
#'
#' @param mDat The metadata used in the mapping.
#' @param nMultiplexed A vector of length matching \code{mDat} giving the number of donors multiplexed in each row.  Or single number if all samples the same.  Default is no multiplexing.
#' @param extras A list containing function to extract each field.  Each function should take one argument, the \code{mDat} object.
#' @param srcs Optional parameter storing the source.
standardiseMetadata = function(mDat,nMultiplexed=1,extras=list(),srcs=NULL){
  if(length(nMultiplexed)==1)
    nMultiplexed = rep(nMultiplexed,nrow(mDat))
  if(length(nMultiplexed)!=nrow(mDat))
    stop("nMultiplexed of invalid length.")
  tt = list()
  for(nom in unique(mDat$humanDonorName)){
    tt[[nom]] = list(IDs = unique(mDat$groupName[mDat$humanDonorName==nom]),
                     srcs = srcs
                     )
    if(localMapping)
      tt[[nom]]$lustrePathByID = setNames(file.path(dstBase,nom,unique(mDat$groupName[mDat$humanDonorName==nom]),'outs'),unique(mDat$groupName[mDat$humanDonorName==nom]))
    #Add multiplexing info
    w = which(mDat$humanDonorName==nom)
    nMult = unique(nMultiplexed[w])
    if(length(nMult)!=1)
      stop("Invalid multiplexing argument.")
    tt[[nom]]$nMultiplexed=nMult
    for(eNom in names(extras)){
      tt[[nom]][[eNom]] = extras[[eNom]](mDat[mDat$humanDonorName==nom,])
    }
  }
  return(tt)
}

#################################
# Remap previously mapped files #
#################################

if(localMapping){
  for(nom in names(dataSources)){
    message(sprintf("Reprocessing samples from %s",nom))
    #Skip reprocessing for these, they've been mapped well enough by default.
    if(nom %in% c('WholeEmbryo1'))
      next
    dataSources[[nom]] = reprocessLocal(dataSources[[nom]])
  }
}



#################
# Data from HCA #
#################

#Define and download all manifests
maniSrc=c(AIDApilotdata='https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp12&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%22f0f89c14-7460-4bab-9d42-22228a91f185%22%5D%7D%7D&objectKey=manifests%2F49bd7c23-a0b5-5be4-8594-536303bdcbc6.tsv',
          nasalMucosaLifespan = 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp12&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%228d566d35-d8d3-4975-a351-be5e25e9b2ea%22%5D%7D%7D&objectKey=manifests%2F02f9fa2b-d3de-5be1-ba3b-054386a375e0.tsv',
          HnsccImmuneLandscape = 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp12&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%223089d311-f9ed-44dd-bb10-397059bad4dc%22%5D%7D%7D&objectKey=manifests%2F963cf757-b950-5129-acc3-d94ca3d347a0.tsv',
          PulmonaryFibrosis = 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp12&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%22c1a9a93d-d9de-4e65-9619-a9cec1052eaa%22%5D%7D%7D&objectKey=manifests%2F252e5ee5-42e7-59ba-b9b6-de41fd04695f.tsv',
          Microglia = 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp12&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%22b51f49b4-0d2e-4cbd-bbd5-04cd171fc2fa%22%5D%7D%7D&objectKey=manifests%2F11ff84a1-290d-517d-8745-15d5a9ff78c5.tsv',
          oralMucosaAtlas = 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp12&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%2250151324-f3ed-4358-98af-ec352a940a61%22%5D%7D%7D&objectKey=manifests%2F7ffc5daa-4a77-5daf-a73e-7ff34aa0bd12.tsv',
          vagWall = 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp12&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%2271eb5f6d-cee0-4297-b503-b1125909b8c7%22%5D%7D%7D&objectKey=manifests%2F32fcb5c7-d5b6-56d0-ad0c-48d8e4007b40.tsv',
          IleumChrons = 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp12&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%22504e0cee-1688-40fa-b936-361c4a831f87%22%5D%7D%7D&objectKey=manifests%2F5dfe58eb-8f4a-52a3-bb0c-a5c37d0f6544.tsv',
          ColonMesenchymeIBD = 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp12&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%22f8aa201c-4ff1-45a4-890e-840d63459ca2%22%5D%7D%7D&objectKey=manifests%2F454e43e1-f774-53d3-a558-d6d3a8ef25b0.tsv',
          ImmuneCensus = 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp12&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%22cc95ff89-2e68-4a08-a234-480eca21ce79%22%5D%7D%7D&objectKey=manifests%2F02240c4d-cc47-514d-a778-9c2237bd8691.tsv',
          TissueStability = 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp12&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%22c4077b3c-5c98-4d26-a614-246d12c2e5d7%22%5D%7D%7D&objectKey=manifests%2F4cdf46a3-8ac7-563b-bb59-a5d6134f173b.tsv',
          COPD = 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp12&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%22ad04c8e7-9b7d-4cce-b8e9-01e31da10b94%22%5D%7D%7D&objectKey=manifests%2F2b25e9e7-82ba-5209-812e-feaaef060058.tsv',
          Urine = 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp12&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%224af795f7-3e1d-4341-b867-4ac0982b9efd%22%5D%7D%7D&objectKey=manifests%2F5a43e4d0-8c81-5add-ab7e-46782d20da17.tsv',
          breast = 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp12&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%22a004b150-1c36-4af6-9bbd-070c06dbc17d%22%5D%7D%7D&objectKey=manifests%2F40a3c2fe-3862-5ff2-a075-44c8d7381e08.tsv',
          BreastCancer = 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp12&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%227c75f07c-608d-4c4a-a1b7-b13d11c0ad31%22%5D%7D%7D&objectKey=manifests%2F0c26c0c0-fe95-5b0e-bb3b-f732fc87b7bf.tsv',
          HumanFirstTrimesterPlacentaDecidua = 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp12&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%221cd1f41f-f81a-486b-a05b-66ec60f81dcf%22%5D%7D%7D&objectKey=manifests%2Ffac3109a-861d-5aba-9c19-52c002849029.tsv',
          Aorta = 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp12&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%2207073c12-8006-4710-a00b-23abdb814904%22%5D%7D%7D&objectKey=manifests%2Fbe86a11b-7b1e-5326-88ff-532954a2e3a2.tsv',
          LungEndothelium = 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp12&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%22d7b7beae-652b-4fc0-9bf2-bcda7c7115af%22%5D%7D%7D&objectKey=manifests%2F83d69949-c2ce-595a-8e2e-7705abb53c24.tsv',
          EpithelialIBD = 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp12&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%22c893cb57-5c9f-4f26-9312-21b85be84313%22%5D%7D%7D&objectKey=manifests%2Fa03ba448-f11b-51ec-b510-c6fc6fa0d342.tsv',
          DiabeticNephropathy = 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp12&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%22577c946d-6de5-4b55-a854-cd3fde40bff2%22%5D%7D%7D&objectKey=manifests%2F2834c309-40fd-5ac7-86f2-1ef728bbd690.tsv',
          decidua = 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp12&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%223cfcdff5-dee1-4a7b-a591-c09c6e850b11%22%5D%7D%7D&objectKey=manifests%2Fd7de94cb-d941-5b22-9808-31bd8b3f1882.tsv',
          Ovary = 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp12&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%22faeedcb0-e046-4be7-b1ad-80a3eeabb066%22%5D%7D%7D&objectKey=manifests%2F138dcbf4-2822-52a8-80db-acaff946f5cb.tsv',
          Teeth = 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp12&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%22d3446f0c-30f3-4a12-b7c3-6af877c7bb2d%22%5D%7D%7D&objectKey=manifests%2Fc9d72a75-a433-5a3e-8ddc-d2487b08319e.tsv',
          FemaleGonads = 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp12&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%2256e73ccb-7ae9-4fae-a738-acfb69936d7a%22%5D%7D%7D&objectKey=manifests%2F54db5602-693e-55b5-8167-e9adec268cc1.tsv',
          kidneyCortex = 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp12&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%222af52a13-65cb-4973-b513-39be38f2df3f%22%5D%7D%7D&objectKey=manifests%2Fd67b3cac-febe-54a9-ba88-55c355155248.tsv',
          retina = 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp12&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%228bd2e5f6-9453-4b9b-9c56-59e3a40dc87e%22%5D%7D%7D&objectKey=manifests%2F92ceaf48-7c2a-54cb-85be-86aaaec04809.tsv',
          fallopianTubes = 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp12&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%2221ea8ddb-525f-4f1f-a820-31f0360399a2%22%5D%7D%7D&objectKey=manifests%2Fcd905ec5-650f-5339-9fcb-bbd1df84f388.tsv',
          HCA_cfLungMetadata = 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp12&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%22e526d91d-cf3a-44cb-80c5-fd7676b55a1d%22%5D%7D%7D&objectKey=manifests%2F4056eb68-c426-50d2-a7d1-e57b47e645be.tsv',
          HCA_kidneyMetadata = 'https://service.azul.data.humancellatlas.org/manifest/files?catalog=dcp12&format=compact&filters=%7B%22projectId%22%3A+%7B%22is%22%3A+%5B%2294023a08-611d-4f22-a8c9-90956e091b2e%22%5D%7D%7D&objectKey=manifests%2F34a503c9-a804-5196-893e-1aeb18a6b36e.tsv'
)
if(dcpUpdate){
  #Save old ones in a dated folder
  oldManiDir = file.path('Data',paste0('HCA_Manifests_',format(Sys.time(),'%Y_%m_%d_%H_%M_%S')))
  dir.create(oldManiDir)
  #Check all files exist
  srcs = file.path('Data',paste0(names(maniSrc),'.tsv'))
  if(!all(file.exists(srcs)))
    stop("Manifest missing, check manifest definition.")
  #Make a copy of existing ones
  dsts  = file.path(oldManiDir,basename(srcs))
  for(i in seq_along(srcs))
    file.copy(srcs[i],dsts[i])
  if(!all(file.exists(dsts)))
    stop("Backup failed")
  #OK, now download the new ones.
  for(i in seq_along(srcs))
    system(sprintf('wget -O %s "%s"',normalizePath(srcs[i]),maniSrc[i]))
  #Compare old and new
  for(i in seq_along(srcs)){
    a = read.delim(dsts[i],sep='\t',header=TRUE)
    b = read.delim(srcs[i],sep='\t',header=TRUE)
    message(sprintf("For manifest %s, old one had %d rows and %d samples, new has %d rows and %d samples",names(maniSrc[i]),nrow(a),length(unique(a$donor_organism.biomaterial_core.biomaterial_id)),nrow(b),length(unique(b$donor_organism.biomaterial_core.biomaterial_id))))
  }
}


#############
# AIDA pilot
#Project here https://data.humancellatlas.org/explore/projects/f0f89c14-7460-4bab-9d42-22228a91f185.  Downloaded metadata 
#The update on 2022_01_14 seems to essentially drop a bunch of files from the old metadata table.  Unfortunately this breaks the creation of the human readable names.
mDat = read.delim('Data/AIDApilotdata.tsv',sep='\t',header=TRUE)
mDat = formatMetadataHCA(mDat,'HCA_AIDA',useDonorName=FALSE)
#Use this to ensure consistency in human name with the old naming convention
tmp = readRDS('Data/AIDAdonorNameMapper.RDS')
mDat$humanDonorName = tmp[mDat$donor_organism.biomaterial_core.biomaterial_id]
mDat = mDat[mDat$donor_organism.sex=='female',]
write.table(mDat,file.path(dstBase,'HCA_AIDA_toMap.tsv'),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
dataSources$HCA_AIDA = mDat
if(localMapping)
  processHCA(mDat)


###############
# Nasal mucosa
#Project https://data.humancellatlas.org/explore/projects/8d566d35-d8d3-4975-a351-be5e25e9b2ea
#NOTE: Samples are multiplexed
mDat = read.delim('Data/nasalMucosaLifespan.tsv',sep='\t',header=TRUE)
mDat = formatMetadataHCA(mDat,'HCA_NasalMucosa',useDonorName=FALSE)
x = split(mDat,mDat$groupName)
x = lapply(x,function(e) {
             tmp = gsub('_.*','',e$file_name)
             e$fileName = sprintf('%s_S1_L%03d_%s_001.fastq.gz',
                                  e$groupName,
                                  match(tmp,unique(tmp)),
                                  c(read1='R1',read2='R2',index1='I1')[e$read_index]
                                  )
             return(e)})
mDat = do.call(rbind,x)
write.table(mDat,file.path(dstBase,'HCA_NasalMucosa_toMap.tsv'),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
dataSources$HCA_NasalMucosa = mDat
if(localMapping)
  processHCA(mDat)



####################################
# Head and neck cancer immune cells
#Projects https://data.humancellatlas.org/explore/projects/3089d311-f9ed-44dd-bb10-397059bad4dc
mDat = read.delim('Data/HnsccImmuneLandscape.tsv',sep='\t',header=TRUE)
mDat = formatMetadataHCA(mDat,'HCA_headNeckCancer',useDonorName=TRUE)
#Drop useless crappy males
mDat = mDat[mDat$donor_organism.sex=='female',]
write.table(mDat,file.path(dstBase,'HCA_headNeckCancer_toMap.tsv'),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
dataSources$HCA_headNeckCancer = mDat
#Some samples just fail, so drop these
mDat = mDat[mDat$groupName!='HD_PBMC_1_cells',]
if(localMapping)
  processHCA(mDat)


#####################
# Pulmonary Fibrosis
#Project https://data.humancellatlas.org/explore/projects/c1a9a93d-d9de-4e65-9619-a9cec1052eaa
mDat = read.delim('Data/PulmonaryFibrosis.tsv',sep='\t',header=TRUE)
mDat = formatMetadataHCA(mDat,'HCA_PF',useDonorName=TRUE)
#Create filenames
tt = split(mDat,mDat$groupName)
for(grp in names(tt)){
  tgts = tt[[grp]]
  #Don't know if this should increment S or L and don't know if it matters
  tmp = gsub('_[RI][12].*','',tgts$file_name)
  tgts$fileName = sprintf('%s_S001_L%03d_%s_001.fastq.gz',
                            grp,
                            match(tmp,unique(tmp)),
                            c(read1='R1',read2='R2',index1='I1')[tgts$read_index]
                            )
  tt[[grp]] = tgts
}
mDat = do.call(rbind,tt)
write.table(mDat,file.path(dstBase,'HCA_PF_toMap.tsv'),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
dataSources$HCA_PF = mDat
#Drop broken samples
mDat$processFail =FALSE
mDat$processFail[mDat$groupName %in% c('THD0005','TILD010','VUHD65','VUHD66','VUHD70','VUILD61-1')]=TRUE
if(localMapping)
  processHCA(mDat)


############
# Microglia
#Project https://data.humancellatlas.org/explore/projects/b51f49b4-0d2e-4cbd-bbd5-04cd171fc2fa
#Mostly bulk and SS2
mDat = read.delim('Data/Microglia.tsv',sep='\t',header=TRUE)
mDat = formatMetadataHCA(mDat,'HCA_MicroGlia',FALSE)
mDat$fileName = paste0(mDat$groupName,'_',gsub('.*_S','S',mDat$file_name))
mDat = mDat[mDat$library_preparation_protocol.library_construction_approach=="10X 3' v2 sequencing",]
write.table(mDat,file.path(dstBase,'HCA_MicroGlia_toMap.tsv'),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
dataSources$HCA_MicroGlia = mDat
if(localMapping)
  processHCA(mDat)


##############
# Oral mucosa
#Project https://data.humancellatlas.org/explore/projects/50151324-f3ed-4358-98af-ec352a940a61/project-metadata
#NOTE: Not confident in naming scheme, but should at least run.
mDat = read.delim('Data/oralMucosaAtlas.tsv',sep='\t',header=TRUE)
mDat = formatMetadataHCA(mDat,'HCA_OralMucosa',TRUE)
#Construct names
tt = split(mDat,mDat$groupName)
for(grp in names(tt)){
  tgts = tt[[grp]]
  #Don't know if this should increment S or L and don't know if it matters
  tmp = gsub('_[RI][12].*','',tgts$file_name)
  tgts$fileName = sprintf('%s_S001_L%03d_%s_001.fastq.gz',
                            grp,
                            match(tmp,unique(tmp)),
                            c(read1='R1',read2='R2',index1='I1')[tgts$read_index]
                            )
  tt[[grp]] = tgts
}
mDat = do.call(rbind,tt)
#Drop useless crappy males
mDat = mDat[mDat$donor_organism.sex=='female',]
write.table(mDat,file.path(dstBase,'HCA_OralMucosa_toMap.tsv'),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
dataSources$HCA_OralMucosa = mDat
if(localMapping)
  processHCA(mDat)
 

##########
# VagWall
#Project https://data.humancellatlas.org/explore/projects/71eb5f6d-cee0-4297-b503-b1125909b8c7
#Looks pretty standard, should work without issue.  Some bulk there too
mDat = read.delim('Data/vagWall.tsv',sep='\t',header=TRUE)
mDat = formatMetadataHCA(mDat,'HCA_VagWall',TRUE)
mDat = mDat[mDat$library_preparation_protocol.library_construction_approach=="10X 3' v2 sequencing",]
write.table(mDat,file.path(dstBase,'HCA_VagWall_toMap.tsv'),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
dataSources$HCA_VagWall = mDat
if(localMapping)
  processHCA(mDat)



#########
# Chrons
#Project https://data.humancellatlas.org/explore/projects/504e0cee-1688-40fa-b936-361c4a831f87
mDat = read.delim('Data/IleumChrons.tsv',sep='\t',header=TRUE)
mDat = formatMetadataHCA(mDat,'HCA_Chrons',TRUE)
#Drop those that make no sense (duplicated R1)
tmp = gsub('_[0-9].fastq.gz','',mDat$file_name)
tt = table(tmp,mDat$read_index)
mDat = mDat[tmp %in% rownames(tt)[apply(tt,1,max)==1],]
#Construct names
tt = split(mDat,mDat$groupName)
for(grp in names(tt)){
  tgts = tt[[grp]]
  #Don't know if this should increment S or L and don't know if it matters
  tmp = gsub('_[0-9].fastq.gz','',tgts$file_name)
  tgts$fileName = sprintf('%s_S001_L%03d_%s_001.fastq.gz',
                            grp,
                            match(tmp,unique(tmp)),
                            c(read1='R1',read2='R2',index1='I1')[tgts$read_index]
                            )
  tt[[grp]] = tgts
}
mDat = do.call(rbind,tt)
write.table(mDat,file.path(dstBase,'HCA_Chrons_toMap.tsv'),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
dataSources$HCA_Chrons = mDat
#Drop this one, cellranger can't detect the chemistry
mDat$processFail = mDat$groupName=='SRX9053453'
#No sex information
if(localMapping)
  processHCA(mDat)


###################
# Colon mesenchyme
#Project https://data.humancellatlas.org/explore/projects/f8aa201c-4ff1-45a4-890e-840d63459ca2
#Pretty standard naming, should be fine.  Bunch of SS2 as well
mDat = read.delim('Data/ColonMesenchymeIBD.tsv',sep='\t',header=TRUE)
mDat = formatMetadataHCA(mDat,'HCA_ColMes',TRUE)
#Basically named sensibly already
mDat$fileName = paste0(mDat$groupName,gsub('.*_S','_S',mDat$file_name))
mDat = mDat[mDat$library_preparation_protocol.library_construction_approach=="10X 3' v2 sequencing",]
mDat = mDat[mDat$donor_organism.genus_species == 'Homo sapiens',]
#Drop useless crappy males
mDat = mDat[mDat$donor_organism.sex=='female',]
write.table(mDat,file.path(dstBase,'HCA_ColMes_toMap.tsv'),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
dataSources$HCA_ColMes = mDat
if(localMapping)
  processHCA(mDat)


################
# Immune census
#Project https://data.humancellatlas.org/explore/projects/cc95ff89-2e68-4a08-a234-480eca21ce79
mDat = read.delim('Data/ImmuneCensus.tsv',sep='\t',header=TRUE)
mDat = formatMetadataHCA(mDat,'HCA_ImmCen',TRUE)
#Basically named sensibly already
mDat$fileName = paste0(mDat$groupName,gsub('.*_S','_S',mDat$file_name))
#Drop useless crappy males
mDat = mDat[mDat$donor_organism.sex=='female',]
write.table(mDat,file.path(dstBase,'HCA_ImmCen_toMap.tsv'),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
dataSources$HCA_ImmCen = mDat
if(localMapping)
  processHCA(mDat)


###################
# Ischaemic tissue
#Project https://data.humancellatlas.org/explore/projects/c4077b3c-5c98-4d26-a614-246d12c2e5d7
#Bunch of bulk here too
mDat = read.delim('Data/TissueStability.tsv',sep='\t',header=TRUE)
mDat = formatMetadataHCA(mDat,'HCA_Isc',TRUE)
mDat = mDat[mDat$library_preparation_protocol.library_construction_approach=="10X v2 sequencing",]
#Drop useless crappy males
mDat = mDat[mDat$donor_organism.sex=='female',]
write.table(mDat,file.path(dstBase,'HCA_Isc_toMap.tsv'),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
dataSources$HCA_Isc = mDat
#And bad samples
mDat$processFail = mDat$groupName %in% c('A15-OES-0-TL2-1_cells','A15-OES-0-TL3-1_cells','A15-OES-0-TL4-1_cells')
if(localMapping)
  processHCA(mDat)


#######
# COPD
#Project https://data.humancellatlas.org/explore/projects/ad04c8e7-9b7d-4cce-b8e9-01e31da10b94
mDat = read.delim('Data/COPD.tsv',sep='\t',header=TRUE)
mDat = formatMetadataHCA(mDat,'HCA_COPD',TRUE)
#Drop useless crappy males
mDat = mDat[mDat$donor_organism.sex=='female',]
write.table(mDat,file.path(dstBase,'HCA_COPD_toMap.tsv'),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
dataSources$HCA_COPD = mDat
if(localMapping)
  processHCA(mDat)


########
# Urine
#Project https://data.humancellatlas.org/explore/projects/4af795f7-3e1d-4341-b867-4ac0982b9efd
#NOTE: 12 multiplexed samples
mDat = read.delim('Data/Urine.tsv',sep='\t',header=TRUE)
mDat = formatMetadataHCA(mDat,'HCA_Urine',FALSE)
tt = split(mDat,mDat$groupName)
for(grp in names(tt)){
  tgts = tt[[grp]]
  #Don't know if this should increment S or L and don't know if it matters
  tmp = gsub('_R[0-9].fastq.gz','',tgts$file_name)
  tgts$fileName = sprintf('%s_S001_L%03d_%s_001.fastq.gz',
                            grp,
                            match(tmp,unique(tmp)),
                            c(read1='R1',read2='R2',index1='I1')[tgts$read_index]
                            )
  tt[[grp]] = tgts
}
mDat = do.call(rbind,tt)
write.table(mDat,file.path(dstBase,'HCA_Urine_toMap.tsv'),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
dataSources$HCA_Urine = mDat
if(localMapping)
  processHCA(mDat)



#########
# Breast
#Project https://data.humancellatlas.org/explore/projects/a004b150-1c36-4af6-9bbd-070c06dbc17d
mDat = read.delim('Data/breast.tsv',sep='\t',header=TRUE)
mDat = formatMetadataHCA(mDat,'HCA_Breast')
write.table(mDat,file.path(dstBase,'HCA_Breast_toMap.tsv'),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
dataSources$HCA_Breast = mDat
if(localMapping)
  processHCA(mDat)


################
# Breast cancer
#Project https://data.humancellatlas.org/explore/projects/7c75f07c-608d-4c4a-a1b7-b13d11c0ad31
#NOTE: Data looks interest (matched tumour/normal tissue), but metadata is a mess.  Need to look at original paper to work it out I think.  This is my best guess from the metadata alone
mDat = read.delim('Data/BreastCancer.tsv',sep='\t',header=TRUE)
#mDat = formatMetadataHCA(mDat,'HCA_BreastCancer',TRUE)
mDat = mDat[!grepl('\\.tar$',mDat$file_name),]
#This seems as close to consistent as I can get
mDat$humanDonorName = paste0('HCA_BreastCancer_',mDat$donor_organism.biomaterial_core.biomaterial_id)
mDat$groupName = mDat$cell_suspension.biomaterial_core.biomaterial_id
#These bastards even fucked up the read index
mDat$read_index_guess = gsub('.*[\\._](R[0-9])(_[0-9]+)?\\.fastq\\.gz$','\\1',mDat$file_name)
#I think these are correct, if you assign all R3s to R2s
mDat$read_index_fixed = ifelse(mDat$read_index_guess=='R1','read1','read2')
tt = split(mDat,mDat$groupName)
for(grp in names(tt)){
  tgts = tt[[grp]]
  #This pairs things up, it will run but is probably nonsense...
  tmp = gsub('[\\._]R[0-9](_[0-9]+)?\\.fastq\\.gz$','\\1',tgts$file_name)
  tgts$fileName = sprintf('%s_S001_L%03d_%s_001.fastq.gz',
                            grp,
                            match(tmp,unique(tmp)),
                            c(read1='R1',read2='R2',index1='I1')[tgts$read_index_fixed]
                            )
  tt[[grp]] = tgts
}
mDat = do.call(rbind,tt)
#Drop the guess/fixed columns for consistency
mDat$read_index_guess=mDat$read_index_fixed=NULL
write.table(mDat,file.path(dstBase,'HCA_BreastCancer_toMap.tsv'),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
dataSources$HCA_BreastCancer = mDat
mDat$processFail = mDat$humanDonorName  %in% sprintf('HCA_BreastCancer_BC0%d',c(1:8))
if(localMapping)
  processHCA(mDat)

##############
#### Placenta
#Project https://data.humancellatlas.org/explore/projects/1cd1f41f-f81a-486b-a05b-66ec60f81dcf
#NOTE: There's some usable data here, but it's almost all drop-seq so not worth the bother (3 samples total)
mDat = read.delim('Data/HumanFirstTrimesterPlacentaDecidua.tsv',sep='\t',header=TRUE)
mDat = formatMetadataHCA(mDat,'HCA_Placenta',TRUE)
#Drop useless crappy males
mDat = mDat[mDat$donor_organism.sex=='female',]
write.table(mDat,file.path(dstBase,'HCA_Placenta_toMap.tsv'),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
dataSources$HCA_Placenta = mDat
#If I'm doing it, need to drop drop-seq
mDat$processFail = mDat$library_preparation_protocol.library_construction_approach!="10X 3' v2 sequencing"
if(localMapping)
  processHCA(mDat)

#########
# Aorata
#Project https://data.humancellatlas.org/explore/projects/07073c12-8006-4710-a00b-23abdb814904
mDat = read.delim('Data/Aorta.tsv',sep='\t',header=TRUE)
mDat = formatMetadataHCA(mDat,'HCA_Aorta',TRUE)
#Drop useless crappy males
mDat = mDat[mDat$donor_organism.sex=='female',]
write.table(mDat,file.path(dstBase,'HCA_Aorta_toMap.tsv'),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
dataSources$HCA_Aorta = mDat
if(localMapping)
  processHCA(mDat)


###################
# Lung Endothelium
#Project https://data.humancellatlas.org/explore/projects/d7b7beae-652b-4fc0-9bf2-bcda7c7115af
mDat = read.delim('Data/LungEndothelium.tsv',sep='\t',header=TRUE)
mDat = formatMetadataHCA(mDat,'HCA_LungEndo')
mDat$fileName = paste0(mDat$groupName,gsub('.*(_S[0-9]*_L[0-9]*)','\\1',mDat$file_name))
write.table(mDat,file.path(dstBase,'HCA_LungEndo_toMap.tsv'),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
dataSources$HCA_LungEndo = mDat
if(localMapping)
  processHCA(mDat)


######
# IBD
#Project https://data.humancellatlas.org/explore/projects/c893cb57-5c9f-4f26-9312-21b85be84313
mDat = read.delim('Data/EpithelialIBD.tsv',sep='\t',header=TRUE)
mDat = formatMetadataHCA(mDat,'HCA_EpiIBD')
mDat$fileName = paste0(mDat$groupName,gsub('.*bamtofastq_S1_L001_([RI][0-9])_([0-9]+)','_S1_L\\2_\\1_001',mDat$file_name))
write.table(mDat,file.path(dstBase,'HCA_EpiIBD_toMap.tsv'),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
dataSources$HCA_EpiIBD = mDat
if(localMapping)
  processHCA(mDat)

#######################
# Diabetic Nephropathy
#Project https://data.humancellatlas.org/explore/projects/577c946d-6de5-4b55-a854-cd3fde40bff2
#Is single nuclei, so need to work out how we're going to handle it
mDat = read.delim('Data/DiabeticNephropathy.tsv',sep='\t',header=TRUE)
mDat = formatMetadataHCA(mDat,'HCA_DiabeticKid')
mDat$fileName = paste0(mDat$groupName,'_',gsub('.*_S','S',mDat$file_name))
#Drop useless crappy males
mDat = mDat[mDat$donor_organism.sex=='female',]
write.table(mDat,file.path(dstBase,'HCA_DiabeticKid_toMap.tsv'),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
dataSources$HCA_DiabeticKid = mDat
if(localMapping)
  processHCA(mDat,useIntrons=TRUE)

##########
# Decidua
#Project https://data.humancellatlas.org/explore/projects/3cfcdff5-dee1-4a7b-a591-c09c6e850b11
#NOTE: Multiplexed...
mDat = read.delim('Data/decidua.tsv',sep='\t',header=TRUE)
mDat = formatMetadataHCA(mDat,'HCA_Decidua')
write.table(mDat,file.path(dstBase,'HCA_Decidua_toMap.tsv'),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
dataSources$HCA_Decidua = mDat
mDat$processFail = mDat$groupName %in% c('SRX9806677','SRX9806676')
if(localMapping)
  processHCA(mDat)

##########
# Ovaries
#Project https://data.humancellatlas.org/explore/projects/faeedcb0-e046-4be7-b1ad-80a3eeabb066
#Actually quite a lot of data here
mDat = read.delim('Data/Ovary.tsv',sep='\t',header=TRUE)
mDat = formatMetadataHCA(mDat,'HCA_Ovaries')
mDat$fileName = paste0(mDat$groupName,'_S1_L',gsub('.*_([0-9]+)\\.fastq.gz','\\1',mDat$file_name),'_',c(read1='R1',read2='R2',index1='I1')[mDat$read_index],'_001.fastq.gz')
write.table(mDat,file.path(dstBase,'HCA_Ovaries_toMap.tsv'),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
dataSources$HCA_Ovaries = mDat
mDat$processFail = mDat$groupName %in% c('SAMN09765859','SAMN09765846','SAMN09765854')
if(localMapping)
  processHCA(mDat)


########
# Teeth
#Project https://data.humancellatlas.org/explore/projects/d3446f0c-30f3-4a12-b7c3-6af877c7bb2d
#NOTE: This has not been exported to be mapped.
mDat = read.delim('Data/Teeth.tsv',sep='\t',header=TRUE)
mDat = formatMetadataHCA(mDat,'HCA_Teeth')
x = split(mDat,mDat$groupName)
x = lapply(x,function(e) {
             tmp = gsub('_.*','',e$file_name)
             e$fileName = sprintf('%s_S1_L%03d_%s_001.fastq.gz',
                                  e$groupName,
                                  match(tmp,unique(tmp)),
                                  c(read1='R1',read2='R2',index1='I1')[e$read_index]
                                  )
             return(e)})
mDat = do.call(rbind,x)
mDat$processFail = mDat$groupName %in% c('SRX9770748','SRX9770749','SRX9770750','SRX9770751','SRX9770752')
if(localMapping)
  processHCA(mDat)

################
# Female gonads
#Project https://data.humancellatlas.org/explore/projects/56e73ccb-7ae9-4fae-a738-acfb69936d7a
mDat = read.delim('Data/FemaleGonads.tsv',sep='\t',header=TRUE)
mDat = formatMetadataHCA(mDat,'HCA_FemGonads')
x = split(mDat,mDat$groupName)
x = lapply(x,function(e) {
             tmp = gsub('_.*','',e$file_name)
             e$fileName = sprintf('%s_S1_L%03d_%s_001.fastq.gz',
                                  e$groupName,
                                  match(tmp,unique(tmp)),
                                  c(read1='R1',read2='R2',index1='I1')[e$read_index]
                                  )
             return(e)})
mDat = do.call(rbind,x)
write.table(mDat,file.path(dstBase,'HCA_FemGonads_toMap.tsv'),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
dataSources$HCA_FemGonads = mDat
if(localMapping)
  processHCA(mDat)

################
# Kidney Cortex
#Project https://data.humancellatlas.org/explore/projects/2af52a13-65cb-4973-b513-39be38f2df3f
#NOTE: Literally just one sample...
mDat = read.delim('Data/kidneyCortex.tsv',sep='\t',header=TRUE)
mDat = formatMetadataHCA(mDat,'HCA_KidCortex')
mDat = mDat[grep('10x 5',mDat$library_preparation_protocol.library_construction_approach),]
#Drop useless crappy males
mDat = mDat[mDat$donor_organism.sex=='female',]
write.table(mDat,file.path(dstBase,'HCA_KidCortex_toMap.tsv'),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
dataSources$HCA_KidCortex = mDat
if(localMapping)
  processHCA(mDat)


#########
# Retina
#Project https://data.humancellatlas.org/explore/projects/8bd2e5f6-9453-4b9b-9c56-59e3a40dc87e/project-metadata
mDat = read.delim('Data/retina.tsv',sep='\t',header=TRUE)
mDat = formatMetadataHCA(mDat,'HCA_Retina')
write.table(mDat,file.path(dstBase,'HCA_Retina_toMap.tsv'),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
dataSources$HCA_Retina = mDat
if(localMapping)
  processHCA(mDat)


##################
# Fallopian tubes
#Project https://data.humancellatlas.org/explore/projects/21ea8ddb-525f-4f1f-a820-31f0360399a2
mDat = read.delim('Data/fallopianTubes.tsv',sep='\t',header=TRUE)
mDat = formatMetadataHCA(mDat,'HCA_Fallopian')
write.table(mDat,file.path(dstBase,'HCA_Fallopian_toMap.tsv'),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
dataSources$HCA_Fallopian = mDat
if(localMapping)
  processHCA(mDat)


##########
# CF Lung
#Project https://data.humancellatlas.org/explore/projects/e526d91d-cf3a-44cb-80c5-fd7676b55a1d
mDat = read.delim('Data/HCA_cfLungMetadata.tsv',sep='\t',header=TRUE)
#The read_index seems to be inconsistent in a subset of them.  Specifically, it says it's from a different chemistry from the read groups.  Easiest solution just to drop all indicies 
mDat = mDat[mDat$read_index!='index1',]
mDat = formatMetadataHCA(mDat,'HCA_CF',useDonorName=TRUE)
write.table(mDat,file.path(dstBase,'HCA_CF_toMap.tsv'),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
dataSources$HCA_CF = mDat
mDat$processFail = mDat$library_preparation_protocol.library_construction_approach=='Drop-seq'
#We could switch to this, but the naming convention is different in a 1-to-1 mapping.
if(localMapping)
  processHCA(mDat)


#########
# Kidney
#Project https://data.humancellatlas.org/explore/projects/94023a08-611d-4f22-a8c9-90956e091b2e/project-metadata
mDat = read.delim('Data/HCA_kidneyMetadata.tsv',sep='\t',header=TRUE)
#Super straightforward, but again has a different naming convention
mDat = formatMetadataHCA(mDat,'HCA_kidney')
write.table(mDat,file.path(dstBase,'HCA_kidney_toMap.tsv'),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
dataSources$HCA_kidney = mDat
if(localMapping)
  processHCA(mDat)

##################################
# Fix up and homogenise metadata #
##################################

#Sanger samples first
tt = dataSources[!grepl("^(HCA|EXT)_",names(dataSources))]
dd = lengths(lapply(tt,function(e) e$IDs))
dd = data.frame(ID = unlist(lapply(tt,function(e) e$IDs)),
                donor = rep(names(tt),dd),
                org = 'sanger',
                dataType = unlist(lapply(names(tt),function(e) if(length(tt[[e]]$mapType)==dd[e]){tt[[e]]$mapType}else{rep(tt[[e]]$mapType,dd[e])}))
)
dd$dataSrc = NA
dd$mappedDataPath = NA
dd$notes = NA
dd$multiplexedIDs = NA
for(nom in names(tt)){
  if(!is.null(tt[[nom]]$srcs))
    dd$dataSrc[dd$donor==nom]=paste(tt[[nom]]$srcs,collapse=', ')
  if(!is.null(tt[[nom]]$notes))
    dd$notes[dd$donor==nom] = tt[[nom]]$notes
  if(!is.null(tt[[nom]]$lustrePathByID)){
    m = match(paste0(nom,'___',names(tt[[nom]]$lustrePathByID)),paste0(dd$donor,'___',dd$ID))
    if(any(is.na(m)))
      stop("Failed to match IDs")
    dd$mappedDataPath[m] = normalizePath(tt[[nom]]$lustrePathByID)
  }
}
#Save out here if exporting
if(!localMapping)
  write.table(dd,file.path(dstBase,'Sanger_toMap.tsv'),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
#Add in the DNA
for(nom in names(tt)){
  for(bamNom in c('normBAM','mumBAM','dadBAM')){
    if(!is.null(tt[[nom]][[bamNom]])){
      dd = rbind(dd,
                 data.frame(ID = basename(tt[[nom]][[bamNom]]),
                            donor = nom,
                            org = 'sanger',
                            dataType = gsub('BAM','DNA',bamNom),
                            dataSrc = tt[[nom]][[bamNom]],
                            mappedDataPath = tt[[nom]][[paste0(bamNom,'_lustre')]],
                            notes = NA,
                            multiplexedIDs=NA)
                 )
    }
  }
}
datMap = dd
#Now HCA ones
tt = dataSources[grepl("^HCA_",names(dataSources))]
dd = do.call(rbind,tt)
dd$datasetName = rep(names(tt),sapply(tt,nrow))
#Save all the info if exporting
if(!localMapping)
  write.table(dd,file.path(dstBase,'HCA_toMap.tsv'),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
#Smush down into standard format
tmp = strsplit(dd$donor_organism.biomaterial_core.biomaterial_id,' || ',fixed=TRUE)
dd = data.frame(ID = dd$groupName,
                donor = dd$humanDonorName,
                org = 'HCA',
                dataType = c(`single cell`='scRNA',`single nucleus`='snRNA')[dd$library_preparation_protocol.nucleic_acid_source],
                dataSrc = dd$datasetName,
                mappedDataPath = NA,
                notes=NA,
                multiplexedIDs=ifelse(lengths(tmp)==1,NA,sapply(tmp,paste,collapse=', ')))
#Now need to stick in the paths if they exist
if(localMapping){
  dd$mappedDataPath = file.path(dstBase,dd$donor,dd$ID,'outs')
  dd$mappedDataPath[!dir.exists(dd$mappedDataPath)]=NA
}
dd = unique(dd)
datMap = rbind(dd,datMap)
#Add placenta
tt = dataSources[grepl('^EXT_',names(dataSources))]
dd = data.frame(ID = unlist(lapply(tt,function(e) e$IDs)),
                donor = unlist(lapply(tt,function(e) e$donor)),
                org = 'external',
                dataType = unlist(lapply(tt,function(e) e$mapType)),
                dataSrc = unlist(lapply(tt,function(e) e$src)))
dd$mappedDataPath = NA
dd$notes = NA
dd$multiplexedIDs = NA
if(!localMapping)
  write.table(dd,file.path(dstBase,'EXT_toMap.tsv'),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
datMap = rbind(dd,datMap)


if(!localMapping){
  message("The assumption is that at this point all the _toMap.tsv files will be sent somewhere and mapped.  If this has not been done, everything will fail.")
  #Now need to get the mapped paths
  srcs = list.files(mapDir,full.names=TRUE)
  srcs = srcs[dir.exists(srcs)]
  srcs = setNames(lapply(srcs,function(e) {tmp = list.files(e,full.names=TRUE); tmp[dir.exists(tmp)]}),basename(srcs))
  if(length(srcs)==0)
    stop("No mapped output files found.")
  srcs = data.frame(group = rep(names(srcs),lengths(srcs)),
                    path = unlist(srcs),
                    ID = basename(unlist(srcs)))
  rownames(srcs)=NULL
  #Deal with duplicates
  pres = c('OUT_custom','OUT_extras','OUT_remaining_Sanger_10x','OUT_sanger_batch2','OUT_sanger_batch1','OUT_placenta')
  srcs = srcs[order(srcs$ID,factor(srcs$group,levels=c(pres,setdiff(srcs$group,pres)))),]
  srcs = srcs[!duplicated(srcs$ID),]
  w = which(is.na(datMap$mappedDataPath))
  m = match(datMap$ID[w],srcs$ID)
  ww = !is.na(m)
  datMap$mappedDataPath[w[ww]] = srcs$path[m[ww]]
}

#############################################
# Stage to folder organised as we'd like it #
#############################################

#Only want to proceed with ceratin data and can only proceed with things that have been successfully mapped
dd = datMap[datMap$dataType %in% c('scRNA','snRNA'),]
dd = dd[!is.na(dd$mappedDataPath),]
#Drop those with no bam file
dd = dd[sapply(dd$mappedDataPath,function(e) any(grepl('\\.bam$',list.files(e)))),]
dd$stagedPath = NA
dd$stagedPathBAM = NA
dd$stagedPathCnts = NA
message("Staging mapped files to final directory.")
for(nom in unique(dd$donor)){
  w = which(dd$donor==nom)
  tgtDir = file.path(stageDir,nom)
  tgtFiles = file.path(tgtDir,dd$ID[w])
  unlink(tgtDir,recursive=TRUE)
  dir.create(tgtDir,showWarnings=FALSE)
  suppressWarnings(file.symlink(normalizePath(dd$mappedDataPath[w]),tgtFiles))
  #Now point at the things
  dd$stagedPath[w] = tgtFiles
  is10X = !file.exists(file.path(tgtFiles,'Aligned.sortedByCoord.out.bam'))
  dd$stagedPathBAM[w] = file.path(tgtFiles,ifelse(is10X,'possorted_genome_bam.bam','Aligned.sortedByCoord.out.bam'))
  dd$stagedPathCnts[w] = file.path(tgtFiles,ifelse(is10X,'filtered_feature_bc_matrix','output/Gene/filtered/'))
}
#This is the final validated object that everything else can work off of
scDat = dd

#######################
# Define used samples #
#######################

#These are the ones that actually get used.
areUsed = list(goldStandard = c('Wilms2','Wilms3','GOSH025'),
               #cancerValidation = c("GOSH028", "GOSH029", "GOSH032", "IMP672", "VHL_RCC", "RCC2"),
               cancerValidation = c("GOSH029", "GOSH032", "IMP672", "VHL_RCC", "RCC2"),
               bigIndividual = 'F45',
               placenta = c("P13", "Hrv98", "Hrv99", "Hrv100", "F37", "F40"),
               #kidney= c("HCA_KidCortex005", "HCA_kidney002", "HCA_kidney007", "HCA_kidney008", "HCA_kidney009", "HCA_kidney010", "HCA_kidney012", "HCA_kidney013", "HCA_kidney014", "HCA_kidney015", "HCA_kidney017", "HCA_kidney018", "RCC2", "VHL_RCC", "Wilms2", "Wilms3", "KidTransplant", "F35", "F41", "F45"),
               kidney= c("HCA_KidCortex005", "HCA_kidney002", "HCA_kidney007", "HCA_kidney008", "HCA_kidney009", "HCA_kidney010", "HCA_kidney012", "HCA_kidney013", "HCA_kidney014", "HCA_kidney015", "HCA_kidney017", "HCA_kidney018", "RCC2", "VHL_RCC", "Wilms2", "Wilms3", "KidTransplant"),
               oralMucosa=c("HCA_OralMucosa_PD170", "HCA_OralMucosa_BM150", "HCA_OralMucosa_BM152", "HCA_OralMucosa_BM156", "HCA_OralMucosa_BM157", "HCA_OralMucosa_BM158", "HCA_OralMucosa_GM136", "HCA_OralMucosa_GM143", "HCA_OralMucosa_GM144", "HCA_OralMucosa_GM147", "HCA_OralMucosa_GM148", "HCA_OralMucosa_GM183", "HCA_OralMucosa_GM184", "HCA_OralMucosa_GM238", "HCA_OralMucosa_GM241", "HCA_OralMucosa_GM242", "HCA_OralMucosa_PD153", "HCA_OralMucosa_PD161"),
               lung = c("HCA_PF_nonfibrotic_control_9", "HCA_PF_donor_IPF_11", "HCA_PF_donor_IPF_5", "HCA_PF_donor_IPF_8", "HCA_PF_nonfibrotic_control_8", "HCA_PF_nonfibrotic_control_6", "HCA_PF_donor_IPF_10", "HCA_PF_donor_sarcoidosis_2", "HCA_PF_donor_nsip_1", "HCA_PF_donor_IPF_1", "HCA_CF_donor_BG9_7", "HCA_CF_donor_BG9_18", "HCA_CF_donor_BG9_19", "HCA_CF_donor_BG9_20", "HCA_CF_donor_BG9_21", "HCA_CF_donor_BG9_22", "HCA_CF_donor_BG9_23", "HCA_CF_donor_BG9_24", "HCA_CF_donor_BG9_25", "HCA_CF_donor_BG9_30", "HCA_CF_donor_BG9_31", "HCA_CF_donor_BG9_34", "HCA_CF_donor_BG9_35", "HCA_CF_donor_BG9_36", "HCA_CF_donor_BG9_37", "HCA_CF_donor_cc05p", "HCA_CF_donor_dd09p", "HCA_CF_donor_dd49p", "HCA_CF_donor_h18", "HCA_CF_donor_cf019", "HCA_CF_donor_KKCFFT32N", "HCA_CF_donor_KKD019N", "HCA_CF_donor_KKD020N", "HCA_CF_donor_KKD024N", "HCA_CF_donor_ND15989", "HCA_CF_donor_CF301018", "HCA_CF_donor_ND15656", "HCA_CF_donor_ND15989_ALI", "HCA_CF_donor_ND15656_ALI"),
               vagina = c("HCA_VagWall_Nor1", "HCA_VagWall_Nor2", "HCA_VagWall_Nor3", "HCA_VagWall_Nor4", "HCA_VagWall_Nor5", "HCA_VagWall_POP1", "HCA_VagWall_POP10", "HCA_VagWall_POP11", "HCA_VagWall_POP12", "HCA_VagWall_POP13", "HCA_VagWall_POP14", "HCA_VagWall_POP15", "HCA_VagWall_POP16", "HCA_VagWall_POP2", "HCA_VagWall_POP3", "HCA_VagWall_POP4", "HCA_VagWall_POP5", "HCA_VagWall_POP6", "HCA_VagWall_POP7", "HCA_VagWall_POP8", "HCA_VagWall_POP9"),
               pbmc = c("HCA_AIDA006", "HCA_AIDA135", "HCA_AIDA089", "HCA_AIDA171", "HCA_AIDA091", "HCA_AIDA051", "HCA_AIDA295", "HCA_AIDA234", "HCA_AIDA009", "HCA_AIDA107", "HCA_AIDA125", "HCA_AIDA057", "HCA_AIDA326", "HCA_AIDA016", "HCA_AIDA247", "HCA_AIDA335", "HCA_AIDA298", "HCA_AIDA336", "HCA_AIDA140", "HCA_AIDA332", "HCA_AIDA267", "HCA_AIDA170", "HCA_AIDA331", "HCA_AIDA028", "HCA_AIDA138", "HCA_AIDA219", "HCA_AIDA221", "HCA_AIDA026", "HCA_AIDA229", "HCA_AIDA110", "HCA_AIDA094", "HCA_AIDA300", "HCA_AIDA206", "HCA_AIDA106", "HCA_AIDA253", "HCA_AIDA314", "HCA_AIDA197", "HCA_AIDA008", "HCA_AIDA120", "HCA_AIDA052", "HCA_AIDA238", "HCA_AIDA175", "HCA_AIDA306", "HCA_AIDA320", "HCA_AIDA313", "HCA_AIDA291", "HCA_AIDA265", "HCA_AIDA304", "HCA_AIDA157", "HCA_AIDA224", "HCA_AIDA075", "HCA_AIDA119", "HCA_AIDA053", "HCA_AIDA079", "HCA_AIDA329", "HCA_AIDA309", "HCA_AIDA022", "HCA_AIDA142", "HCA_AIDA191", "HCA_AIDA337", "HCA_AIDA252", "HCA_AIDA216", "HCA_AIDA273", "HCA_AIDA281", "HCA_AIDA280", "HCA_AIDA330", "HCA_AIDA214", "HCA_AIDA080", "HCA_AIDA199", "HCA_AIDA014", "HCA_AIDA065", "HCA_AIDA064", "HCA_AIDA018", "HCA_AIDA254", "HCA_AIDA196", "HCA_AIDA308", "HCA_AIDA045", "HCA_AIDA285", "HCA_AIDA210", "HCA_AIDA147", "HCA_AIDA277", "HCA_AIDA299", "HCA_AIDA155", "HCA_AIDA015", "HCA_AIDA132", "HCA_AIDA236", "HCA_AIDA162", "HCA_AIDA104", "HCA_AIDA180", "HCA_AIDA312", "HCA_AIDA041", "HCA_AIDA227", "HCA_AIDA217", "HCA_AIDA189", "HCA_AIDA040", "HCA_AIDA255", "HCA_AIDA174", "HCA_AIDA270", "HCA_AIDA039", "HCA_AIDA283", "HCA_AIDA186", "HCA_AIDA067", "HCA_AIDA188", "HCA_AIDA128", "HCA_AIDA081", "HCA_AIDA077", "HCA_AIDA334", "HCA_AIDA136", "HCA_AIDA237", "HCA_AIDA302", "HCA_AIDA084", "HCA_AIDA123", "HCA_AIDA121", "HCA_AIDA031", "HCA_AIDA061", "HCA_AIDA097", "HCA_AIDA141", "HCA_AIDA074", "HCA_AIDA168", "HCA_AIDA056", "HCA_AIDA262", "HCA_AIDA220", "HCA_AIDA187", "HCA_AIDA190", "HCA_AIDA193", "HCA_AIDA090", "HCA_AIDA163", "HCA_AIDA007", "HCA_AIDA047", "HCA_AIDA059", "HCA_AIDA069", "HCA_AIDA266", "HCA_AIDA185", "HCA_AIDA145", "HCA_AIDA099", "HCA_AIDA111", "HCA_AIDA251", "HCA_AIDA020", "HCA_AIDA093", "HCA_AIDA208", "HCA_AIDA152", "HCA_AIDA042", "HCA_AIDA264", "HCA_AIDA293", "HCA_AIDA181", "HCA_AIDA070", "HCA_AIDA122", "HCA_AIDA134", "HCA_AIDA278", "HCA_AIDA002", "HCA_AIDA143", "HCA_AIDA019", "HCA_AIDA017", "HCA_AIDA256", "HCA_AIDA092", "HCA_AIDA013", "HCA_AIDA192",  "HCA_AIDA023", "HCA_AIDA318", "HCA_AIDA115", "HCA_AIDA085", "HCA_AIDA276", "HCA_AIDA164", "HCA_AIDA102", "HCA_AIDA066", "HCA_AIDA301", "HCA_AIDA315", "HCA_AIDA296", "HCA_AIDA046", "HCA_AIDA228", "HCA_AIDA200", "HCA_AIDA034", "HCA_AIDA159", "HCA_AIDA310", "HCA_AIDA261", "HCA_AIDA169", "HCA_AIDA303",  "HCA_AIDA194", "HCA_AIDA204", "HCA_AIDA275", "HCA_AIDA209", "HCA_AIDA257", "HCA_AIDA060", "HCA_AIDA258", "HCA_AIDA153", "HCA_AIDA235", "HCA_AIDA271", "HCA_AIDA083", "HCA_AIDA030", "HCA_AIDA072", "HCA_AIDA095", "HCA_AIDA151", "HCA_AIDA050", "HCA_AIDA211", "HCA_AIDA279", "HCA_AIDA156", "HCA_AIDA198", "HCA_AIDA048", "HCA_AIDA100"))

areUsedAnalysis = areUsed
areUsedAnalysis$kidney = grep('^F[0-9]+',areUsedAnalysis$kidney,invert=TRUE,value=TRUE)
t1 = setNames(as.numeric(gsub(' year','',dataSources$HCA_AIDA$donor_organism.organism_age)),dataSources$HCA_AIDA$humanDonorName)
areUsedAnalysis$pbmc = unique(names(t1)[t1<35])

###########
# Finish! #
###########

