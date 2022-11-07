#' Various utility functions

#############
# Libraries #
#############

library(Seurat)
library(Matrix)

#############
# Functions #
#############

#' Progress monitoring functions
#'
#' Get a list of those directories that are currently being mapped, based on the presence of the 00_lock file.
#'
#' @return A list of donors and samples are currently running
getRunning = function(){
  x = list.files('Data/mappedData',full.names=TRUE)
  x = x[dir.exists(x)]
  xx = lapply(x,list.files,full.names=TRUE)
  xx = lapply(xx,function(e) e[dir.exists(e)])
  names(xx) = basename(x)
  xxx = lapply(xx,function(e) setNames(file.exists(file.path(e,'00_lock')),basename(e)))
  w = which(sapply(xxx,any))
  return(xxx[w])
}

#' Which ones are yet to be remapped
#'
#' Looks across all folders and indentifies which are yet to be remapped based on the absence of an "outs" folder.
#'
#' @param dataSources The usual dataSources object.  If given, exclude anything from consideration that is not also in dataSources.
#' @return A list of unmapped donors/samples.
getUnmapped = function(dataSources=NULL){
  x = list.files('Data/mappedData',full.names=TRUE)
  x = x[dir.exists(x)]
  xx = lapply(x,list.files,full.names=TRUE)
  xx = lapply(xx,function(e) e[dir.exists(e)])
  names(xx) = basename(x)
  xxx = lapply(xx,function(e) setNames(!file.exists(file.path(e,'outs')),basename(e)))
  w = which(sapply(xxx,any))
  out = xxx[w]
  if(!is.null(dataSources)){
    out = out[names(out) %in% names(dataSources)]
    out = setNames(lapply(names(out),function(e) out[[e]][intersect(names(out[[e]]),dataSources[[e]]$IDs)]),names(out))
    out = out[sapply(out,any)]
  }
  return(out)
}

#' Extract chemistry form web summary
#'
#' Scrapes the web summary html to work out which chemistry was used.  By default, only considers those folders where remapping has finished.
#'
#' @param tgt A 10X output directory (containing a web_summary.html file) to get the chemistry for.
#' @return A data.frame with the chemistries.
getChemistries = function(tgt=NULL){
  if(is.null(tgt)){
    x = list.files('Data/mappedData',full.names=TRUE)
    x = x[dir.exists(x)]
    xx = lapply(x,list.files,full.names=TRUE)
    xx = lapply(xx,function(e) e[dir.exists(e)])
    names(xx) = basename(x)
    xxx = lapply(xx,function(e) setNames(file.exists(file.path(e,'outs')),basename(e)))
    w = which(sapply(xxx,any))
    y=xxx[w]
    y = data.frame(donorID = rep(names(y),lengths(y)),
                   sampleID  = unlist(lapply(y,names)),
                   isMapped = unlist(y),
                   row.names=NULL)
    y$hasHTML = file.exists(file.path('Data','mappedData',y$donorID,y$sampleID,'outs','web_summary.html'))
    y = y[y$isMapped & y$hasHTML,]
    y$baseDir = file.path('Data','mappedData',y$donorID,y$sampleID,'outs')
  }else{
    y = data.frame(donorID=tgt,sampleID='.',baseDir=tgt)
  }
  #Extract chemistry from each
  y$chem = NA
  for(i in seq(nrow(y))){
    message(sprintf("[%d of %d] Getting chemistry for %s lane %s",i,nrow(y),y$donorID[i],y$sampleID[i]))
    tmp = readLines(file.path(y$baseDir[i],'web_summary.html'),n=10000) 
    w = grep("Chemistry",tmp)
    #The old html format
    if(nchar(tmp[w])<200){
      chem = gsub('.*<td>(.*)<.*','\\1',tmp[w+1])
    }else{
      tmp = strsplit(gsub('.*..Chemistry','"Chemistry',tmp[w]),']')[[1]][1]
      chem = strsplit(gsub('"','',tmp),',')[[1]][2]
    }
    y$chem[i]=chem
  }
  return(y)
}

#' Run celltypist
#'
#' Export a count matrix via a file to python, run celltypist, then return the results.
#'
#' @param cnts The count matrix (in sparse matrix format) that we want to call things on.
#' @param ... Passed to the celltypist.annotate call.
#' @return A list containing the raw logit probabilities (logitMat) and the label predictions (labMat)
runCelltypist = function(cnts,...){
  #Create three temp files and write out matrix
  base = tempfile()
  mtxFile = paste0(base,'.mtx')
  gnsFile = paste0(base,'_gns.tsv')
  clsFile = paste0(base,'_cls.tsv')
  writeMM(cnts,mtxFile)
  write.table(rownames(cnts),gnsFile,sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)
  write.table(colnames(cnts),clsFile,sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)
  #Extra params
  theDots = list(...)
  #Default for transpose_input
  if(!'transpose_input' %in% names(theDots))
    theDots$transpose_input = TRUE
  if(length(theDots)==0){
    extraParams=''
  }else{
    if(any(lengths(theDots))!=1)
      stop("Any extra characters must be of length 1.")
    extraParams = paste0(names(theDots),'=',sprintf(ifelse(sapply(theDots,is.character),"'%s'",'%s'),as.character(unlist(theDots))))
    #Fix up bools into python format
    extraParams = gsub('FALSE','False',gsub('TRUE','True',extraParams))
    #Wrap it up
    extraParams = paste0(',',paste(extraParams,collapse=','))
  }
  #The python code to run
  tmpdir=tempdir()
  cmd = sprintf("import celltypist\npreds=celltypist.annotate('%s',gene_file='%s',cell_file='%s'%s)\npreds.to_table('%s')",mtxFile,gnsFile,clsFile,extraParams,tmpdir)
  system(paste0('python -c "',cmd,'"'))
  #Load output 
  logitMat = read.table(file.path(tmpdir,"decision_matrix.csv"),header=TRUE,sep=',')
  rownames(logitMat) = logitMat$X
  logitMat$X=NULL
  labMat = read.table(file.path(tmpdir,"predicted_labels.csv"),header=TRUE,sep=',')
  unlink(c(mtxFile,gnsFile,clsFile))
  return(list(logitMat=logitMat,labMat=labMat))
}

#' Train a logistic regression model with celltypist
#'
#' Passes (using files) the counts and cell type labels to celltypist.  There are many parameters to the celltypist.train function, which can be set using the ...
#'
#' @param cnts Sparse matrix to use to train the model.
#' @param outPath Where to save the model.  Should end in pkl.
#' @param labels The cell type labels to use for model training.
#' @param ... Extra parameters to pass to celltypist.train.
#' @return Nothing, just saves the model where we should.
trainCelltypistModel = function(cnts,outPath,labels=colnames(cnts),...){
  #Create three temp files and write out matrix
  base = tempfile()
  mtxFile = paste0(base,'.mtx')
  gnsFile = paste0(base,'_gns.tsv')
  labFile = paste0(base,'_cls.tsv')
  writeMM(cnts,mtxFile)
  write.table(rownames(cnts),gnsFile,sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)
  write.table(labels,labFile,sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)
  #Extra params
  theDots = list(...)
  #Default for transpose_input
  if(!'transpose_input' %in% names(theDots))
    theDots$transpose_input = TRUE
  if(length(theDots)==0){
    extraParams=''
  }else{
    if(any(lengths(theDots))!=1)
      stop("Any extra characters must be of length 1.")
    extraParams = paste0(names(theDots),'=',sprintf(ifelse(sapply(theDots,is.character),"'%s'",'%s'),as.character(unlist(theDots))))
    #Fix up bools into python format
    extraParams = gsub('FALSE','False',gsub('TRUE','True',extraParams))
    #Wrap it up
    extraParams = paste0(',',paste(extraParams,collapse=','))
  }
  #The python code to run
  tmpdir=tempdir()
  cmd = sprintf("import celltypist\nfit=celltypist.train('%s',labels='%s',genes='%s'%s)\nfit.write('%s')",mtxFile,labFile,gnsFile,extraParams,outPath)
  system(paste0('python -c "',cmd,'"'))
  unlink(c(mtxFile,gnsFile,labFile))
}


#' Work out plot dimensions
#'
#' If you want to plot N things in a grid, work out how many rows and columns to use, guided by a target aspect ratio.
#'
#' @param N Number of things to plot.
#' @param ar Target aspect ratio.
#' @return Vector of two numbers giving the number of rows and columns.
gridDimensions = function(N,ar=4/3){
  #Work out plot dimensions, using aspect ratio
  x=sqrt(ar*N)
  #Which of the two options with rounding gives the closest match to the aspect ratio
  ar1 = floor(x)/ceiling(N/floor(x))
  ar2 = ceiling(x)/ceiling(N/ceiling(x))
  if(abs(ar1-ar) >= abs(ar2-ar)){
    x = floor(x)
    y = ceiling(N/x)
  }else{
    x = ceiling(x)
    y = ceiling(N/x)
  }
  return(c(x,y))
}

#' Construct layout for grid of N plots
#'
#' Given N plots that will be plotted in an x by y grid, construct an appropriate layout call.  This call will leave empty any entries that are redundant, for example when you want to plot 13 things in a 3 by 5 grid, the last 2 entries should be empty.  
#'
#' The size of the plots on the left and bottom edges are also expanded slightly to leave room for axis labels.
#'
#' @param N Number of things to plot.
#' @param x Number of columns.
#' @param y Number of rows.
#' @param widthExtra How much bigger to make the left column.
#' @param heightExtra How much bigger to make the bottom row.  Assumes bottom labels are bigger.
#' @return Positions in the grid that should 
gridLayout = function(N,x,y,widthExtra=0.2,heightExtra=0.5){
  #Create matrix with plots, then add one at left/bottom for margins
  mat = rep(0,x*y)
  mat[seq(N)] = seq(N)
  mat = matrix(mat,nrow=y,byrow=TRUE)
  #Add column to left
  mat = cbind(rep(0,y),mat)
  mat = rbind(mat,rep(0,x+1))
  layout(mat,
         widths=c(widthExtra,rep(1,x)),
         heights=c(rep(1,y),heightExtra))
  #layout(mat,
  #       widths = c(1+widthExtra,rep(1,(x-1))),
  #       heights= c(rep(1,(y-1)),1+heightExtra))
  #Record the idx of plots at different boundaries
  labIdxs = list(left=setdiff(mat[,2],0),
                 right = setdiff(mat[,x+1],0),
                 bottom = setdiff(mat[y,],0),
                 top = setdiff(mat[1,],0))
  #labIdxs = which(mat==0,arr.ind=TRUE)
  ##Drop those that are in the boundary plots
  #labIdxs = labIdxs[labIdxs[,1]!=y+1 & labIdxs[,2]!=1,,drop=FALSE]
  #if(length(labIdxs)!=0){
  #  labIdxs[,1]=labIdxs[,1]-1
  #}
  #tmp = which(mat!=0,arr.ind=TRUE)
  #tmp = tmp[tmp[,1]==y,]
  #labIdxs = rbind(labIdxs,tmp)
  return(labIdxs)
}


#' Sets margins for each plot.
#'
#' Given N things to be plotted in x by y grid, work out the margin for entry i and j.  labIdxs should contain the coordinates of the entries that should have a label at the bottom.  The margin for ordinary plots is set by \code{baseMar}, which is then modified by adding \code{bottomMar} and \code{leftMar} as appropriate.
#'
#' Note that any plots where the label on the bottom is to be drawn have xpd=NA set, so that the bottom labels won't be truncated.  This can be reversed by restoring xpd=FALSE, which is the default.
#'
#' @param x Number of columns.
#' @param y Number of rows.
#' @param i Particular column.
#' @param j Particular row.
#' @param labIdxs List indicating which plots fall on the margin of the plotting area.
#' @param baseMar The margin of the default plots.
#' @param bottomMar Extra margin at bottom.
#' @param leftMar Extra margin at left.
#' @return List with xaxt, yaxt, ii
setGridMargins = function(x,y,i,j,labIdxs,baseMar = c(0,0,2,0.5),bottomMar=8,leftMar=3){
  ii = j + x*(i-1)
  labBottom = ii %in% labIdxs$bottom
  labLeft = ii %in% labIdxs$left
  xaxt=ifelse(labBottom,'s','n')
  yaxt=ifelse(labLeft,'s','n')
  par(mar=baseMar)
  ####labBottom = paste0(i,'_',j) %in% paste0(labIdxs[,1],'_',labIdxs[,2])
  ####labLeft = j==1
  ###xaxt='n'
  ###yaxt='n'
  ###if(labBottom & labLeft){
  ###  #par(mar=baseMar +c(bottomMar,leftMar,0,0))
  ###  yaxt=xaxt='s'
  ###}else if(labBottom && i==y){
  ###  #par(mar=baseMar+c(bottomMar,0,0,0))
  ###  xaxt='s'
  ###}else if(labLeft){
  ###  #par(mar=baseMar+c(0,leftMar,0,0))
  ###  yaxt='s'
  ###}else{
  ###  #par(mar=baseMar)
  ###}
  ###ii = j + x*(i-1)
  ###if(labBottom){
  ###  xaxt='s'
  ###  #par(xpd=NA)
  ###}
  return(list(xaxt=xaxt,yaxt=yaxt,ii=ii))
}



