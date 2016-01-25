# FILE: 006_DETECT_SUBCHROMOSOMAL_FEATURES.r
# AUTHOR: Zeljko Jovan Dzakula
# DATE: 7/29/2014

################################################################################
################################################################################
################################################################################
################################################################################

################################################################################
detectBreakpoints <- function( copyNumberProfile, chrSelector, mask, scaleFactor=1, performClustering=FALSE ) {
    breakpoints <- list();
    tryResult <- try( {
        for ( chr in 1:24 ) {
            chrBinSelection <- ( chrSelector[[ chr ]][ 1 ] ):( chrSelector[[ chr ]][ 2 ] );
            chrProfile <- copyNumberProfile[ chrBinSelection ] * mask[ chrBinSelection ];
            chrProfile <- scaleFactor * chrProfile[ !is.na( chrProfile ) ];
            chrProfile[ chrProfile < 0 ] <- 0;
            x <- seq( 1:length( chrProfile ) ) * 5e4 / 1e3;
            x <- x[ !is.na( chrProfile ) ];
            if ( performClustering ) {
                nClusters <- length( chrProfile ) / 20;
                copyNumberCluster <- evaluateClusterCopyNumber( inputProfile=chrProfile, nClusters=nClusters, 
                    weight=nClusters, maxNIterations=100, tolerance=1e-6,
                    trueProfile=NULL, plotFlag=TRUE );
            } else {
                copyNumberCluster <- NULL;
            } # if performClustering else
            copyNumberRunningMedian <- runningMedian( inputProfile=chrProfile, firstWindow=11, secondWindow=21 );
            #
    #        pdf( file=paste0( projectFolder, "RESULTS/", id, "/Chr", chr, ".PDF" ),
    #             width=12, height=8 );
    #        par( mfrow=c( 1, 1 ) );
    #        plot( x, chrProfile, type="l", main=paste0( id, ": Chr", chr ),
    #              xlab="Chromosome Position (kbp)", ylab="GROM Copy Number" );
    #        lines( x, copyNumberRunningMedian[[ "Elevation" ]], col="red", lwd=2 );
    #        #lines( x, copyNumberCluster[[ "Elevation" ]], col="green", lwd=2 );
    #        dev.off();
            #
            breakpoints[[ chr ]] <- list(
                "Cluster"=copyNumberCluster,
                "RunningMedian"=copyNumberRunningMedian );
        } # for chr
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in detectBreakpoints" );
        save.image( file=paste0( outputFolder, "006_DETECT_SUBCHROMOSOMAL_FEATURES-1.RData" ) );
    } else {
        cat( "Execution of detectBreakpoints successfully completed\n" );
    } # if error else
    return( breakpoints );
} # detectBreakpoints

################################################################################
getNeighbors <- function( chrProfile, breakpoints, pointIndex ) {
    neighbors <- NULL;
    tryResult <- try( {
        if ( ( pointIndex <= 0 ) | ( pointIndex > length( breakpoints ) ) ) {
            return( neighbors );
        } # if breakPoint
        if ( pointIndex == 1 ) {
            leftPoint <- 1;
        } else {
            leftPoint <- breakpoints[ pointIndex - 1 ];
        } # if breakPoint
        if ( pointIndex == length( breakpoints ) ) {
            rightPoint <- length( chrProfile );
        } else {
            rightPoint <- breakpoints[ pointIndex + 1 ];
        } # if breakPoint
        centralPoint <- breakpoints[ pointIndex ];
        leftSet <- chrProfile[ leftPoint:( centralPoint - 1 ) ];
        rightSet <- chrProfile[ centralPoint:rightPoint ];
        neighbors <- list(
            "LeftPoint"=leftPoint,
            "CentralPoint"=centralPoint,
            "RightPoint"=rightPoint,
            "LeftSet"=leftSet,
            "RightSet"=rightSet
            );
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in getNeighbors" );
        save.image( file=paste0( outputFolder, "006_DETECT_SUBCHROMOSOMAL_FEATURES-2.RData" ) );
    } else {
        #cat( "Execution of getNeighbors successfully completed\n" );
    } # if error else
    return( neighbors ); 
} # getNeighbors

################################################################################
testSignificance <- function( chrProfile, breakpoints, pointIndex ) {
    pValue <- NA;
    tryResult <- try( {
        neighbors <- getNeighbors( chrProfile, breakpoints, pointIndex );
        if ( is.null( neighbors ) ) {
            return( pValue );
        } # if NULL
        if ( is.null( neighbors$LeftSet ) | is.null( neighbors$RightSet ) ) {
            return( pValue );
        } # if NULL
        leftMedian <- round( median( neighbors$LeftSet, na.rm=TRUE ) );
        rightMedian <- round( median( neighbors$RightSet, na.rm=TRUE ) );
        if ( !is.na( leftMedian ) & !is.na( rightMedian ) ) {
            if ( leftMedian == rightMedian ) {
                pValue <- 1.0;
            } else {
                testObject <- wilcox.test( neighbors$LeftSet, neighbors$RightSet, 
                    alternative="two.sided" );
                pValue <- testObject[[ "p.value" ]];
            } # if == 
        } # if !NA
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in testSignificance" );
        save.image( file=paste0( outputFolder, "006_DETECT_SUBCHROMOSOMAL_FEATURES-3.RData" ) );
    } else {
        #cat( "Execution of testSignificance successfully completed\n" );
    } # if error else
    return( pValue );
} # testSignificance

################################################################################
assignPValues <- function( copyNumberProfile, breakpoints, chrSelector, mask ) {
    tryResult <- try( {
        for ( chr in 1:24 ) {
            chrBinSelection <- ( chrSelector[[ chr ]][ 1 ] ):( chrSelector[[ chr ]][ 2 ] );
            chrProfile <- copyNumberProfile[ chrBinSelection ] * mask[ chrBinSelection ];
            chrProfile <- chrProfile[ !is.na( chrProfile ) ] 
            chrProfile[ chrProfile < 0 ] <- 0;
            for ( method in "RunningMedian" ) { # names( breakpoints[[ chr ]] ) ) {
                if ( !is.null( breakpoints[[ chr ]][[ as.character( method ) ]] ) ) {
                    nPoints <- length( breakpoints[[ chr ]][[ as.character( method ) ]]$BreakPoints );
                    if ( nPoints == 0 ) {
                        pValues <- NA;
                    } else {
                        pValues <- rep( NA, nPoints );
                        for ( pointIndex in 1:nPoints ) {
                            pValues[[ pointIndex ]] <- testSignificance( chrProfile, 
                                breakpoints[[ chr ]][[ as.character( method ) ]]$BreakPoints, 
                                pointIndex );
                        } # for pointIndex
                    } # if nPoints else
                    breakpoints[[ chr ]][[ as.character( method ) ]]$PValues <- pValues;
                } # if !NULL
            } # for method
        } # for chr
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in assignPValues" );
        save.image( file=paste0( outputFolder, "006_DETECT_SUBCHROMOSOMAL_FEATURES-4.RData" ) );
    } else {
        #cat( "Execution of assignPValues successfully completed\n" );
    } # if error else
    return( breakpoints );
} # assignPValues

################################################################################
# This is correct but useless:
#refineBreakpoint <- function( chrProfile, breakpoints, pointIndex ) {
#    refinedBreakpoint <- NA;
#    neighbors <- getNeighbors( chrProfile, breakpoints, pointIndex );
#    if ( is.null( neighbors ) ) {
#        return( refinedBreakpoint );
#    } # if NULL
#    integral <- sum( neighbors$RightSet, na.rm=TRUE ) + 
#                sum( neighbors$LeftSet, na.rm=TRUE );
#    totalWidth <- neighbors$RightPoint - neighbors$LeftPoint + 1;
#    leftCopyNumber <- mean( neighbors$LeftSet, na.rm=TRUE );
#    rightCopyNumber <- mean( neighbors$RightSet, na.rm=TRUE );
#    if ( leftCopyNumber == rightCopyNumber ) {
#        # Do nothing - refinedBreakpoint is already NA
#    } else {
#        refinedBreakpoint <- ( integral - totalWidth * rightCopyNumber ) / ( leftCopyNumber - rightCopyNumber );
#        refinedBreakpoint <- neighbors$LeftPoint + refinedBreakpoint;
#        names( refinedBreakpoint ) <- names( breakpoints )[ pointIndex ];
#    } # if == else
#    return( refinedBreakpoint );
#} # refineBreakpoints

################################################################################
unwrapGenomicCoordinates <- function( breakpoints, binSize=5e4 ) {
    tryResult <- try( {
        for ( chr in 1:24 ) {
            for ( method in "RunningMedian" ) { # names( breakpoints[[ chr ]] ) ) {
                if ( !is.null( breakpoints[[ chr ]][[ as.character( method ) ]] ) ) {
                    condensedIndices <- breakpoints[[ chr ]][[ as.character( method ) ]]$BreakPoints;
                    if ( !is.null( condensedIndices ) ) {
                        breakpoints[[ chr ]][[ as.character( method ) ]]$CondensedBreakPointIndices <- condensedIndices;
                            breakpoints[[ chr ]][[ as.character( method ) ]]$BreakPoints;
                        binPositions <- sapply( names( breakpoints[[ chr ]][[ as.character( method ) ]]$BreakPoints ),
                            function( x ) { 
                                x <- as.character( x );
                                x <- as.integer( as.character( unlist( strsplit( x, "_" ) )[[ 2 ]] ) );
                                return( x ) } );
                        genomicCoordinates <- binPositions * binSize;
                        breakpoints[[ chr ]][[ as.character( method ) ]]$BreakPointIndices <- binPositions;
                        breakpoints[[ chr ]][[ as.character( method ) ]]$BreakPoints <- genomicCoordinates;
                    } else {
			breakpoints[[ chr ]][[ as.character( method ) ]]$CondensedBreakPointIndices <- NA;
                        breakpoints[[ chr ]][[ as.character( method ) ]]$BreakPointIndices <- NA;
                    } # if !NULL else
                } # if !NULL
            } # for method
        } # for chr
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in unwrapGenomicCoordinates" );
        save.image( file=paste0( outputFolder, "006_DETECT_SUBCHROMOSOMAL_FEATURES-5.RData" ) );
    } else {
        cat( "Execution of unwrapGenomicCoordinates successfully completed\n" );
    } # if error else
    return( breakpoints );
} # unwrapGenomicCoordinates

################################################################################
unwrapElevations <- function( copyNumberProfile, chrSelector, breakpoints ) {
    tryResult <- try( {
        for ( chr in 1:24 ) {
            for ( method in "RunningMedian" ) { # names( breakpoints[[ chr ]] ) ) {
                if ( !is.null( breakpoints[[ chr ]][[ as.character( method ) ]] ) ) {
                    elevations <- rep( NA, chrSelector[[ chr ]][[ 2 ]] - chrSelector[[ chr ]][[ 1 ]] + 1 );
                    binNames <- names( copyNumberProfile$CopyNumber )[ chrSelector[[ chr ]][[ 1 ]]:chrSelector[[ chr ]][[ 2 ]] ];
                    names( elevations ) <- binNames;
                    selector <- which( !is.na( match( as.character( binNames ), 
                        as.character( names( breakpoints[[ chr ]]$RunningMedian$Elevation ) ) ) ) );
                    elevations[ selector ] <- breakpoints[[ chr ]]$RunningMedian$Elevation;
                    breakpoints[[ chr ]][[ as.character( method ) ]]$Elevation <- elevations;
                } # if !NULL
            } # for method
        } # for chr
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in unwrapElevations" );
        save.image( file=paste0( outputFolder, "006_DETECT_SUBCHROMOSOMAL_FEATURES-6.RData" ) );
    } else {
        cat( "Execution of unwrapElevations successfully completed\n" );
    } # if error else
    return( breakpoints );
} # unwrapElevations

################################################################################
defineBlocks <- function( breakpoints, binSize=5e4 ) {
    tryResult <- try( {
        for ( chr in 1:24 ) {
            for ( method in "RunningMedian" ) { # names( breakpoints[[ chr ]] ) ) {
                if ( !is.null( breakpoints[[ chr ]][[ as.character( method ) ]] ) ) {
                    elevations <- breakpoints[[ chr ]][[ as.character( method ) ]]$Elevation;
                    breakpointIndices <- breakpoints[[ chr ]][[ as.character( method ) ]]$BreakPointIndices;
                    notNA <- which( !is.na( elevations ) );
                    if ( is.na( breakpointIndices ) ) {
                        blocks <- data.frame(
                            "Start"=NA,
                            "End"=NA,
                            "StartIndex"=NA,
                            "EndIndex"=NA,
                            "CopyNumber"=NA,
                            "Length"=NA );
                        blocks$StartIndex[[ 1 ]] <- min( which( !is.na( elevations ) ) ); # The first non-NA bin
                        blocks$EndIndex[[ 1 ]] <- max( which( !is.na( elevations ) ) ); # The last non-NA bin
                        blocks$CopyNumber[[ 1 ]] <- unique( elevations[ 
                            blocks$StartIndex[[ 1 ]]:blocks$EndIndex[[ 1 ]] ] )[[ 1 ]];
                        
                    } else {
                        nBlocks <- length( breakpointIndices ) + 1;
                        dummy <- rep( NA, nBlocks );
                        blocks <- data.frame(
                            "Start"=dummy,
                            "End"=dummy,
                            "StartIndex"=dummy,
                            "EndIndex"=dummy,
                            "CopyNumber"=dummy,
                            "Length"=dummy );
                        blocks$StartIndex[[ 1 ]] <- min( which( !is.na( elevations ) ) ); # The first non-NA bin
                        blocks$EndIndex[[ 1 ]] <- breakpointIndices[[ 1 ]];
                        blocks$CopyNumber[[ 1 ]] <- unique( elevations[ 
                            blocks$StartIndex[[ 1 ]]:blocks$EndIndex[[ 1 ]] ] )[[ 1 ]];
                        notNA <- which( !is.na( elevations ) );
                        if ( nBlocks > 2 ) {
                            for ( block in 2:( nBlocks - 1 ) ) {
                                blocks$StartIndex[[ block ]] <- breakpointIndices[[ block - 1 ]];
                                blocks$EndIndex[[ block ]] <- notNA[ max( which( notNA < breakpointIndices[[ block ]] ) ) ];
                                blocks$CopyNumber[[ block ]] <- unique( elevations[ 
                                    blocks$StartIndex[[ block ]]:blocks$EndIndex[[ block ]] ] )[[ 1 ]];
                            } # for block
                        } # if nBlocks
                        blocks$StartIndex[[ nBlocks ]] <- breakpointIndices[[ nBlocks - 1 ]];
                        blocks$EndIndex[[ nBlocks ]] <- max( notNA );
                        blocks$CopyNumber[[ nBlocks ]] <- unique( elevations[ 
                            blocks$StartIndex[[ nBlocks ]]:blocks$EndIndex[[ nBlocks ]] ] )[[ 1 ]];
                    } # if is.na else
                    blocks$Start <- ( blocks$StartIndex - 1 ) * binSize + 1;
                    blocks$End <- blocks$EndIndex * binSize;
                    blocks$Length <- blocks$End - blocks$Start + 1;
                    breakpoints[[ chr ]][[ as.character( method ) ]]$Blocks <- blocks;
                } # if !NULL
            } # for method
        } # for chr
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in unwrapElevations" );
        save.image( file=paste0( outputFolder, "006_DETECT_SUBCHROMOSOMAL_FEATURES-6.RData" ) );
    } else {
        cat( "Execution of unwrapElevations successfully completed\n" );
    } # if error else
    return( breakpoints );
} # defineBlocks

