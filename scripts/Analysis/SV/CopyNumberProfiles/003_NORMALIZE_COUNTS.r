# FILE: 003_NORMALIZE_COUNTS.r
# AUTHOR: Zeljko Jovan Dzakula
# DATE: 7/21/2014

################################################################################
################################################################################
################################################################################
################################################################################
evaluateRawCountsSlow <- function( cMap, bins ) {
    rawCounts <- NULL;
    tryResult <- try( {
        binNames <- rownames( bins );
        nBins <- length( binNames );
        rawCounts <- rep( NA, nBins );
        names( rawCounts ) <- binNames;
        for ( chr in 1:24 ) {
            chrBins <- which( bins$Chr == chr );
            for ( bin in chrBins ) {
                binStart <- bins$Start[[ bin ]];
                binEnd <- bins$End[[ bin ]];
                labels <- which( ( cMap$CMapId == chr ) & 
                                 ( cMap$Position >= binStart ) &
                                 ( cMap$Position < binEnd ) );
                nLabels <- length( labels );
                if ( nLabels > 0 ) {
                    rawCounts[[ bin ]] <- 
                        sum( cMap$Coverage[ labels ], na.rm=TRUE ) / nLabels;
                } else {
                    #rawCounts[[ bin ]] <- 0;
                } # if nLabels
            } # for bin
        } # for chr
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in evaluateRawCountsSlow" );
        save.image( file=paste0( outputFolder, "003_NORMALIZE_COUNTS-1.RData" ) );
    } else {
        cat( "Execution of evaluateRawCountsSlow successfully completed\n" );
    } # if error else
    return( rawCounts );
} # evaluateRawCountsSlow

evaluateRawCounts <- function( cMap, bins, binLabels ) {
    rawCounts <- NULL;
    tryResult <- try( {
        binNames <- rownames( bins );
        nBins <- length( binNames );
        rawCounts <- rep( NA, nBins );
        names( rawCounts ) <- binNames;
        for ( bin in seq_along( binLabels$Labels ) ) {
            labels <- binLabels$Labels[[ bin ]];
            nLabels <- binLabels$NLabels[[ bin ]];
            if ( nLabels > 0 ) {
                rawCounts[[ bin ]] <- 
                    sum( cMap$Coverage[ labels ], na.rm=TRUE ) / nLabels;
            } else {
                #rawCounts[[ bin ]] <- 0;
            } # if nLabels
        } # for chr
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in evaluateRawCounts" );
        save.image( file=paste0( outputFolder, "003_NORMALIZE_COUNTS-2.RData" ) );
    } else {
        cat( "Execution of evaluateRawCounts successfully completed\n" );
    } # if error else
    return( rawCounts );
} # evaluateRawCounts

################################################################################
evaluateLabelDensityBias <- function( bins, counts ) {
    tryResult <- try( {
        x <- bins[ , "NLabels" ] - median( bins[ , "NLabels" ] );
        linearModel <- lm( counts ~ x );
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in evaluateLabelDensityBias" );
        save.image( file=paste0( outputFolder, "003_NORMALIZE_COUNTS-3.RData" ) );
        return( NULL );
    } else {
        cat( "Execution of evaluateLabelDensityBias successfully completed\n" );
        return( linearModel[[ "coefficients" ]][[ "x" ]] );
    } # if error else
} # evaluateLabelDensityBias

################################################################################
normalizeCounts <- function( labelDensityBias, scaledRawCounts, binParameters, 
                             chrSelector, autosomalBinSelection ) {
    normalizedProfile <- NULL;
    tryResult <- try( {
        normalizedProfile <- scaledRawCounts - labelDensityBias * binParameters[ , "Slope" ];
        normalizedProfile <- normalizedProfile / binParameters[ , "Intercept" ];
        normalizedProfile[ autosomalBinSelection ] <- 2 * normalizedProfile[ autosomalBinSelection ];
        selector <- chrSelector[[ 23 ]][ 1 ]:chrSelector[[ 23 ]][ 2 ];
        normalizedProfile[ selector ] <- 2 * normalizedProfile[ selector ];
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in normalizeCounts" );
        save.image( file=paste0( outputFolder, "003_NORMALIZE_COUNTS-4.RData" ) );
    } else {
        cat( "Execution of normalizeCounts successfully completed\n" );
    } # if error else
    return( normalizedProfile );
} # normalizeCounts

################################################################################
slidingWindowSmooth <- function( countProfile, chrSelector, kernel ) {
    smoothProfile <- NULL;
    tryResult <- try( {
        windowSize <- length( kernel );
        smoothProfile <- NA * countProfile;
        flankSize <- floor( windowSize / 2 );
        for ( chr in 1:24 ) {
            binSelector <- ( chrSelector[[ chr ]][ 1 ] ):( chrSelector[[ chr ]][ 2 ] );
            nBins <- length( binSelector );
            for ( i in 1:windowSize ) {
                smoothProfile[ binSelector[ 1 ] + i - 1 ] <- NA;
            } # for i
            for ( i in ( windowSize + 1 ):( nBins - windowSize ) ) {
                selector <- ( i - flankSize ):( i + flankSize );
                kernelSelector <- which( !is.na( countProfile[ binSelector[ 1 ] + selector - 1 ] ) );
                localKernel <- kernel[ kernelSelector ] / sum( kernel[ kernelSelector ] );
                smoothProfile[ binSelector[ 1 ] + i - 1 ] <- sum( 
                    countProfile[ binSelector[ 1 ] + selector[ kernelSelector ] - 1 ] * localKernel, 
                    na.rm=TRUE );
            } # for i
            for ( i in ( nBins - windowSize + 1 ):nBins ) {
                smoothProfile[ binSelector[ 1 ] + i - 1 ] <- NA;
            } # for i
        } # for chr
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in slidingWindowSmooth" );
        save.image( file=paste0( outputFolder, "003_NORMALIZE_COUNTS-5.RData" ) );
    } else {
        cat( "Execution of slidingWindowSmooth successfully completed\n" );
    } # if error else
    return( smoothProfile );
} # slidingWindowSmooth

################################################################################
processCounts <- function( mapMoleculeFile, bins, binLabels, binParameters, mask, chrSelector ) {
    copyNumberProfile <- NULL;
    tryResult <- try( {
        if ( !testFile( mapMoleculeFile ) ) {
            warning( paste0( "Illegal input: ", mapMoleculeFile ) );
            return( NULL );
        } # if !file
        cMap <- readCMap( mapMoleculeFile ); 
        binNames <- rownames( bins );
        nBins <- length( binNames )
        #rawCounts <- evaluateRawCountsSlow( cMap, bins );
        rawCounts <- evaluateRawCounts( cMap, bins, binLabels );
        names( rawCounts ) <- binNames;
        autosomalBinSelection <- 1:( chrSelector[[ 22 ]][ 2 ] );
        scaledCountProfile <- rawCounts;
        totalAutosomalCounts <- sum( rawCounts[ autosomalBinSelection ], na.rm=TRUE );
        scaledCountProfile <- scaledCountProfile / totalAutosomalCounts * length( autosomalBinSelection );    
        labelDensityBias <- evaluateLabelDensityBias( bins, scaledCountProfile );
        normalizedProfile <- normalizeCounts(
                labelDensityBias, 
                scaledCountProfile,
                binParameters, chrSelector, autosomalBinSelection );
        normalizedProfile[ normalizedProfile < 0 ] <- 0;
        kernel <- c( 1, 3, 5, 3, 1 );
        smoothProfile <- slidingWindowSmooth( normalizedProfile * mask, chrSelector, kernel );
        copyNumberProfile <- list(
            "TotalAutosomalCounts"=totalAutosomalCounts,
            "LabelDensityBias"=labelDensityBias,
            "RawCounts"=rawCounts,
            "CopyNumber"=normalizedProfile,
            "SmoothCopyNumber"=smoothProfile
            );
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in processCounts" );
        save.image( file=paste0( outputFolder, "003_NORMALIZE_COUNTS-6.RData" ) );
    } else {
        cat( "Execution of processCounts successfully completed\n" );
    } # if error else
    return( copyNumberProfile );
} # processCounts
