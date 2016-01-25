# FILE: 004_PROCESS_COPY_NUMBERS.r
# AUTHOR: Zeljko Jovan Dzakula
# DATE: 7/22/2014

################################################################################
################################################################################
################################################################################
################################################################################
prepareSummary <- function( inputProfile, autosomalBinSelection, chrSelector ) {
    result <- NULL;
    tryResult <- try( {
        sdAutosomal <- sd( inputProfile[ autosomalBinSelection ], na.rm=TRUE );
        madAutosomal <- mad( inputProfile[ autosomalBinSelection ], na.rm=TRUE );
        perChr <- matrix( nrow=4, ncol=24 );
        rownames( perChr ) <- c( "SDPerChr", "MADPerChr", "MeanPerChr", "MedianPerChr" );
        colnames( perChr ) <- paste0( "Chr", c( 1:22, "X", "Y" ) );
        for ( chr in 1:24 ) {
            tmp <- inputProfile[ ( chrSelector[[ chr ]][ 1 ] ):( chrSelector[[ chr ]][ 2 ] ) ];
            perChr[ "SDPerChr", chr ] <- sd( tmp, na.rm=TRUE );
            perChr[ "MADPerChr", chr ] <- mad( tmp, na.rm=TRUE );
            perChr[ "MeanPerChr", chr ] <- mean( tmp, na.rm=TRUE );
            perChr[ "MedianPerChr", chr ] <- median( tmp, na.rm=TRUE );
        } # for chr
        result <- list(
            "SDAutosomal"=sdAutosomal,
            "MADAutosomal"=madAutosomal,
            "StatsPerChr"=perChr ); 
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in prepareSummary" );
        save.image( file=paste0( outputFolder, "004_PROCESS_COPY_NUMBERS-1.RData" ) );
    } else {
        cat( "Execution of prepareSummary successfully completed\n" );
    } # if error else
    return( result );
} # prepareSummary

################################################################################
profileSummary <- function( copyNumberProfile, mask, chrSelector ) {
    result <- NULL;
    tryResult <- try( {
        autosomalBinSelection <- 1:( chrSelector[[ 22 ]][ 2 ] );
        #
        summaryScaled <- prepareSummary( copyNumberProfile$RawCounts * mask / copyNumberProfile$TotalAutosomalCounts * length( autosomalBinSelection ), 
            autosomalBinSelection, chrSelector );
        summaryCopyNumber <- prepareSummary( copyNumberProfile$CopyNumber * mask, autosomalBinSelection, chrSelector );
        summarySmooth <- prepareSummary( copyNumberProfile$SmoothCopyNumber * mask, autosomalBinSelection, chrSelector );
        #
        autosomalBins <- length( which( !is.na( mask[ autosomalBinSelection ] ) ) );
        nAutosomalBins <- length( autosomalBins );
        countsPerBin <- sum( copyNumberProfile$RawCounts[ autosomalBins ], na.rm=TRUE ) / nAutosomalBins;
        expectedRelativeError <- 1 / sqrt( countsPerBin );
        #
        result <- list(
            "Scaled"=summaryScaled,
            "CopyNumber"=summaryCopyNumber,
            "Smooth"=summarySmooth,
            "CountsPerBin"=countsPerBin,
            "ExpectedRelativeError"=expectedRelativeError
            );
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in profileSummary" );
        save.image( file=paste0( outputFolder, "004_PROCESS_COPY_NUMBERS-2.RData" ) );
    } else {
        cat( "Execution of profileSummary successfully completed\n" );
    } # if error else
    return( result );
} # profileSummary

