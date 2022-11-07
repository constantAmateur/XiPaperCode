#' Complete end to end pipeline for analysis
#'
#' Requires a handful of other things as well as the bits of code called.
#' irodsCramFetch.sh in Code/ used in staging.R
#' cram2fastq.sh in Code/ used in prepData.R

setwd('~/trueHome/Projects/clonalXX')

##########
# Params #
##########

nParallel=12
forceRefetch=FALSE
localMem = 96 #How much memory to use for cellranger in Gb
localMapping = FALSE #If you're getting someone else to map things from the _toMap.tsv files, set this to false.  Otherwise, set to TRUE and the mapping will be done as part of these scripts.
mapDir = 'Data/cellGenMapped' #Only used if localMapping is false
nSamps=1000 #For bayesian analysis.  How many samples to take.

#############
# Functions #
#############

source('Code/utils.R')
source('Code/analysisFunctions.R')
source('Code/standardPlots.R')

##########################
# Stage files from irods #
##########################

#NOTE: This part needs to be run once on an irods connected machine if doing the mapping ourselves
source('Code/staging.R')

##############################
# Download and prep all data #
##############################

source('Code/prepData.R') 

############################
# Demultiplex where needed #
############################

source('Code/demultiplex.R')

########################
# Check data integrity #
########################

#source('Code/checkDataIntegrity.R')

##############
# Fit models #
##############

source('Code/fitModels.R')

############################
# Tissue specific analysis #
############################

source('Code/kidney.R')
source('Code/lung.R')
source('Code/oralMucosa.R')
source('Code/pbmc.R')
source('Code/placenta.R')
source('Code/vag.R')

###################
# Validate method #
###################

source('Code/validateMethod.R')

#######################
# Population analysis #
#######################

source('Code/popSizeComp.R')

###############
# Make tables #
###############

source('Code/makeTables.R')
