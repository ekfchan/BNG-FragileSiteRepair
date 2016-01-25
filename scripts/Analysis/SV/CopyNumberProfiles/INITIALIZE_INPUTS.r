# FILE: INITIALIZE_INPUTS.r
# AUTHOR: Zeljko Jovan Dzakula
# DATE: 7/21/2014

################################################################################
################################################################################
################################################################################
################################################################################
# The following code assumes that the script INPUT_PARAMETER_CHOICES.r has 
# been sourced.
#
# The following arguments are mandatory:
# ..............................................................................
# GENERAL:
#     scriptFolder: No default
outputFolder <- NA; # No default
outputFilePrefix <- NA; # No default
#
# NORMALIZE COUNTS:
#
# COMPARE WITH SVs
#     inputSvDir: No default
#
# ..............................................................................
# The following arguments are optional:
#
# ALIGN MOLECULES:
performAlignment <- TRUE; # Default=TRUE; Possible values: TRUE, FALSE
inputMoleculesFile <- NA; # No default, must be BNX 
reference <- referenceChoices[[ 1 ]]; # Default=hg19; Possible values: hg19, hg38, BioNano
minLength <- minLengthChoices[[ 1 ]]; # Default=150; Possible values: 100, 125, 150, 180
tValue <- tValueChoices[[ 1 ]]; # Default=9; Possible values: 7, 9, 11
refAligner <- "/home/users/wandrews/execs/2940/RefAligner";
#
# NORMALIZE COUNTS:
binSize <- binSizeChoices[[ 1 ]]; # Possible values: 50kbp, 200kbp
#     binParametersFile: Default=NULL (in which case the actual value is constructed from minLength and tValue)
#
performScaling <- TRUE;
#
# DETECT ANEUPLOIDY:
detectAneuploidy <- TRUE;
#     euploidChrRepresentationsFile: Default=50 kbp bins, hg19 (median, MAD)
#
# DETECT SUBCHROMOSOMAL FEATURES:
detectSubchromosomalFeatures <- TRUE;
#
# COMPARE WITH SVs:
compareWithSVs <- TRUE;
breakPointTolerance <- 2.5e4; # Default=50 kbp
svFolder <- NA;
#
