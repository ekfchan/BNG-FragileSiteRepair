# FILE: 000_COPY_NUMBER_WRAPPER.r
# AUTHOR: Zeljko Jovan Dzakula
# DATE: 7/2/2014
#
# Example command line: 
#R CMD BATCH --no-save --no-restore "--args performAlignment=FALSE inputMoleculesFile=NULL reference='hg19' minLength=150 tValue=9 refAligner='.../RefAligner' binSize='50kbp' performScaling=TRUE detectAneuploidy=TRUE detectSubchromosomalFeatures=TRUE compareWithSVs=TRUE breakPointTolerance=2.5e4 scriptFolder='/home/.../Pipeline/Analysis/RUtil/' copyNumberPackageFolder='/home/.../Pipeline/Analysis/SV/CopyNumberProfiles/' outputFolder='/home/.../TEST/' outputFilePrefix='TEST' svFolder='.../contigs/exp_refineFinal1_sv/merged_smaps/' mapMoleculeFile='/home/.../contigs/alignmolvref/merge/alignmolvref_3_merge_r.cmap'" '/home/.../Pipeline/Analysis/SV/CopyNumberProfiles/000_COPY_NUMBER_WRAPPER.r' '/home/.../TEST/TEST.log' &


################################################################################
################################################################################
################################################################################
################################################################################

paste0 <- function( ... ) { paste( ..., sep="" ) }

#
# Collect input arguments in a list of character vectors
commandLineArguments <- commandArgs( trailingOnly = TRUE );
#
# Check if arguments are passed
if ( length( commandLineArguments ) == 0 ) {
    stop( "No arguments supplied\n" );
} else {
    # Cycle through each element of the list and evaluate the expressions
    for ( item in seq_along( commandLineArguments ) ) {
        if ( item == 1 ) {
            cat( "\n\nInput parameters:\n" );
        } # if item
        eval( parse( text=commandLineArguments[[ item ]] ) );
        cat( paste0( commandLineArguments[[ item ]], "\n" ) );
    } # for item
} # if length else
cat( "\n" );


pkgTest <- function(x)
  {
    if (!require(x, character.only = TRUE, lib.loc=copyNumberPackageFolder))
    {
      install.packages(paste0(copyNumberPackageFolder,"/BioNanoCN_1.0.tar.gz"), lib=copyNumberPackageFolder, dep=TRUE, repos=NULL, type="source")
      if(!require(x, character.only = TRUE, lib.loc=copyNumberPackageFolder)) stop("Package not found")
    }
  }
  
  
#pkgTest("BioNanoCN");
################################################################################
################################################################################
################################################################################
################################################################################
#library(BioNanoCN, lib.loc=copyNumberPackageFolder)

cat( "GROM Copy Number Version 1.0\n" );

# Define the locations of external code
# Import all necessary code from external scripts:
#................................................
#
source( paste0( scriptFolder, "/readmaps.R" ) );
source( paste0( scriptFolder, "/overlap.R" ) );
source( paste0( copyNumberPackageFolder, "/copyNumberCluster.r" ) ); 
source( paste0( copyNumberPackageFolder, "/UTILS.r" ) );
source( paste0( copyNumberPackageFolder, "/001_WORKFLOW.r" ) );
source( paste0( copyNumberPackageFolder, "/002_ALIGN_MOLECULES.r" ) );
source( paste0( copyNumberPackageFolder, "/003_NORMALIZE_COUNTS.r" ) );
source( paste0( copyNumberPackageFolder, "/004_PROCESS_COPY_NUMBERS.r" ) );
source( paste0( copyNumberPackageFolder, "/005_DETECT_ANEUPLOIDY.r" ) );
source( paste0( copyNumberPackageFolder, "/006_DETECT_SUBCHROMOSOMAL_FEATURES.r" ) );
source( paste0( copyNumberPackageFolder, "/007_COMPARE_WITH_SVs.r" ) );
source( paste0( copyNumberPackageFolder, "/008_GENERATE_REPORT.r" ) );

################################################################################
################################################################################
################################################################################
################################################################################

autosomalCopyNumbers <- read.table( paste0( copyNumberPackageFolder,
    "/DATA/hg19/50kbp/EUPLOID_AUTOSOMAL_COPY_NUMBERS.txt" ),
    header=TRUE );

chrXCopyNumbers <- read.table( paste0( copyNumberPackageFolder,
    "/DATA/hg19/50kbp/CHR_X_COPY_NUMBERS.txt" ),
    header=TRUE );

chrYCopyNumbers <- read.table( paste0( copyNumberPackageFolder,
    "/DATA/hg19/50kbp/CHR_Y_COPY_NUMBERS.txt" ),
    header=TRUE );
    
load( file=paste0( copyNumberPackageFolder, 
    "/DATA/hg19/50kbp/chrSelector_hg19_50kbp.RData" ) );

# Define possible choices for input parameters:
source( paste0( copyNumberPackageFolder, "/INPUT_PARAMETER_CHOICES.r" ) );
#
################################################################################
################################################################################
################################################################################
################################################################################

runCopyNumber(
    performAlignment=performAlignment,
    inputMoleculesFile=inputMoleculesFile, # No default, must be BNX 
    reference=reference, # Default=hg19; Possible values: hg19, hg38, BioNano
    minLength=minLength, # Default=150; Possible values: 100, 125, 150, 180
    tValue=tValue, # Default=9; Possible values: 7, 9, 11
    refAligner=refAligner,
    binSize=binSize, # Possible values: 50kbp, 200kbp
    performScaling=performScaling, 
    detectAneuploidy=detectAneuploidy,
    detectSubchromosomalFeatures=detectSubchromosomalFeatures,
    compareWithSVs=compareWithSVs,
    breakPointTolerance=breakPointTolerance,
    scriptFolder=scriptFolder,
    copyNumberPackageFolder=copyNumberPackageFolder,
    outputFolder=outputFolder,
    outputFilePrefix=outputFilePrefix,
    svFolder=svFolder,
    mapMoleculeFile=mapMoleculeFile );


