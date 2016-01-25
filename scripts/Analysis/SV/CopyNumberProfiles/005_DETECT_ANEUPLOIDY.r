# FILE: 005_DETECT_ANEUPLOIDY.r
# AUTHOR: Zeljko Jovan Dzakula
# DATE: 7/2/2014

# Caution: 
# the hardcoded path DATA/hg19/50kbp/ needs to change when we add 
# options (hg38, BioNano reference, various bin sizes):

################################################################################
################################################################################
################################################################################
################################################################################
trimCopyNumberProfile <- function( copyNumberProfile ) {
    tryResult <- try( {
        selector <- which( copyNumberProfile < 0 );
        copyNumberProfile[ selector ] <- 0;
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in trimCopyNumberProfile" );
        save.image( file=paste0( outputFolder, "005_DETECT_ANEUPLOIDY-1.RData" ) );
    } else {
        cat( "Execution of trimCopyNumberProfile successfully completed\n" );
    } # if error else
    return( copyNumberProfile );
} # trimCopyNumberProfile

################################################################################
################################################################################
################################################################################
################################################################################
estimateScalingFactor <- function( meanCopyNumber ) {
    tryResult <- try( {
        tolerance <- 0.25;
        for ( scalingFactor in 1:4 ) {
            meanValue <- scalingFactor * meanCopyNumber;
            testValue <- abs( ( meanValue - round( meanValue ) ) / tolerance );
            if ( testValue < 1 ) {
                return( scalingFactor );
            } # if pValue
        } # for scalingFactor
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in estimateScalingFactor" );
        save.image( file=paste0( outputFolder, "005_DETECT_ANEUPLOIDY-2.RData" ) );
    } else {
        cat( "Execution of estimateScalingFactor successfully completed\n" );
    } # if error else
    return( NA );
} # estimateScalingFactor

################################################################################
################################################################################
################################################################################
################################################################################
evaluateChrCopyNumber <- function( copyNumberProfile, mask, chrSelector ) {
    chrCopyNumber <- NULL;
    tryResult <- try( {
        dummy <- rep( NA, 24 );
        chrCopyNumber <- data.frame(
            "Median"=dummy,
            "MAD"=dummy,
            "Mean"=dummy,
            "SD"=dummy,
            "Length"=dummy,
            "ScalingFactor"=dummy );
        for ( chr in 1:24 ) {
            chrProfile <- ( copyNumberProfile * mask )[ chrSelector[[ chr ]][ 1 ]:chrSelector[[ chr ]][ 2 ] ];
            chrCopyNumber[ chr, "Median" ] <- median( chrProfile, na.rm=TRUE );
            chrCopyNumber[ chr, "MAD" ] <- mad( chrProfile, na.rm=TRUE );
            chrCopyNumber[ chr, "Mean" ] <- mean( chrProfile, na.rm=TRUE );
            chrCopyNumber[ chr, "SD" ] <- sd( chrProfile, na.rm=TRUE );
            chrCopyNumber[ chr, "Length" ] <- chrSelector[[ chr ]][ 2 ] - chrSelector[[ chr ]][ 1 ] + 1;
            chrCopyNumber[ chr, "ScalingFactor" ] <- estimateScalingFactor( 
                chrCopyNumber[ chr, "Mean" ] );
        } # for chr
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in evaluateChrCopyNumber" );
        save.image( file=paste0( outputFolder, "005_DETECT_ANEUPLOIDY-3.RData" ) );
    } else {
        cat( "Execution of evaluateChrCopyNumber successfully completed\n" );
    } # if error else
    return( chrCopyNumber );
} # evaluateChrCopyNumber

################################################################################
################################################################################
################################################################################
################################################################################
scaleCopyNumberProfile <- function( copyNumberProfile, mask, chrSelector ) {
    scalingFactor <- 1;
    tryResult <- try( {
        copyNumberProfile <- trimCopyNumberProfile( copyNumberProfile );
        chrCopyNumber <- evaluateChrCopyNumber( copyNumberProfile, mask, chrSelector );
        scalingFactor <- max( chrCopyNumber$ScalingFactor, na.rm=TRUE );
        if ( is.na( scalingFactor ) ) {
            warning( "Scaling factor is undetermined. No scaling will be done." );
            scalingFactor <- 1;
        } # if is.na
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in scaleCopyNumberProfile" );
        save.image( file=paste0( outputFolder, "005_DETECT_ANEUPLOIDY-4.RData" ) );
    } else {
        cat( "Execution of scaleCopyNumberProfile successfully completed\n" );
    } # if error else
    return( scalingFactor * copyNumberProfile );
} # scaleCopyNumberProfile

################################################################################
################################################################################
################################################################################
################################################################################
classifyAutosomalAneuploidy <- function( copyNumberProfile, mask, chrSelector, 
    populationData=autosomalCopyNumbers, cutoff=10 ) {
    result <- NULL;
    tryResult <- try( {
        chrCopyNumber <- evaluateChrCopyNumber( 
            copyNumberProfile$SmoothCopyNumber, mask, chrSelector );
        autosomes <- 1:22;
        #Z <- ( chrCopyNumber$Median[ autosomes ] - populationData$Median[ autosomes ] ) / 
        #    sqrt( chrCopyNumber$MAD[ autosomes ]^2 + populationData$MAD[ autosomes ]^2 );
        Z <- ( chrCopyNumber$Mean[ autosomes ] - populationData$Mean[ autosomes ] ) / 
            sqrt( chrCopyNumber$SD[ autosomes ]^2 / chrCopyNumber$Length[ autosomes ] + 
                  populationData$SD[ autosomes ]^2 );
        names( Z ) <- paste0( "Chr", autosomes  );
        monosomy <- Z[ which( Z < -cutoff ) ];
        polysomy <- Z[ which( Z > cutoff ) ];
        result <- list( "Monosomy"=monosomy, "Polysomy"=polysomy );
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in classifyAutosomalAneuploidy" );
        save.image( file=paste0( outputFolder, "005_DETECT_ANEUPLOIDY-5.RData" ) );
    } else {
        cat( "Execution of classifyAutosomalAneuploidy successfully completed\n" );
    } # if error else
    return( result );
} # classifyAutosomalAneuploidy

################################################################################
################################################################################
################################################################################
################################################################################
predictGender <- function( ) {
    tryResult <- try( {
        1 # TODO!
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in predictGender" );
        save.image( file=paste0( outputFolder, "005_DETECT_ANEUPLOIDY-6.RData" ) );
    } else {
        cat( "Execution of predictGender successfully completed\n" );
    } # if error else
} # predictGender

################################################################################
################################################################################
################################################################################
################################################################################
classifySexChromosomeAneuploidy <- function() {
    tryResult <- try( {
        1 # TODO!
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in classifySexChromosomeAneuploidy" );
        save.image( file=paste0( outputFolder, "005_DETECT_ANEUPLOIDY-7.RData" ) );
    } else {
        cat( "Execution of classifySexChromosomeAneuploidy successfully completed\n" );
    } # if error else
} # classifySexChromosomeAneuploidy

################################################################################
################################################################################
################################################################################
################################################################################
