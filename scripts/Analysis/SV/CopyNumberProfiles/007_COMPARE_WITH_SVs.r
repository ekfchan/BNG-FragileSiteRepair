# FILE: 007_COMPARE_WITH_SVs.r
# AUTHOR: Zeljko Jovan Dzakula
# DATE: 8/6/2014

################################################################################
################################################################################
################################################################################
################################################################################
compareBreakpointsSMap <- function( breakpointsVsSV, breakpoints, sMap, breakPointTolerance ) {
    tryResult <- try( {
        count <- length( breakpointsVsSV );
        for ( chr in 1:24 ) {
            copyNumberBreakpoints <- breakpoints[[ chr ]]$RunningMedian$BreakPoints;
            selector <- which( ( sMap$RefcontigID1 == chr ) | ( sMap$RefcontigID2 == chr ) );
            if ( length( selector ) > 0 ) {
                for ( breakpoint in copyNumberBreakpoints ) {
                    overlappingSVs <- c( 
                        which( sMap$RefcontigID1[ selector ] == chr & sapply( 
                            sMap$RefStartPos[ selector ], function( x ) overlap( c(
                            breakpoint - breakPointTolerance, 
                            breakpoint + breakPointTolerance ), c(
                            x - breakPointTolerance, 
                            x + breakPointTolerance ) ) ) ),
                        which( sMap$RefcontigID2[ selector ] == chr & sapply( 
                            sMap$RefEndPos[ selector ], function( x ) overlap( c(
                            breakpoint - breakPointTolerance, 
                            breakpoint + breakPointTolerance ), c(
                            x - breakPointTolerance, 
                            x + breakPointTolerance ) ) ) ) );
                    overlappingSVs <- unique( overlappingSVs );
                    if ( length( overlappingSVs ) > 0 ) {
                        for ( index in overlappingSVs ) {
                            count <- count + 1;
                            overlapItem <- list(
                                "Chr"=chr,
                                "CopyNumberBreakpoint"=breakpoint,
                                "SVType"=sMap$Type[[ selector[[ index ]] ]],
                                "SVChr1"=sMap$RefcontigID1[[ selector[[ index ]] ]],
                                "SVChr2"=sMap$RefcontigID2[[ selector[[ index ]] ]],
                                "RefStartPos"=sMap$RefStartPos[[ selector[[ index ]] ]],
                                "RefEndPos"=sMap$RefEndPos[[ selector[[ index ]] ]],
                                "QryContigID"=sMap$QryContigID[[ selector[[ index ]] ]],
                                "QryStartPos"=sMap$QryStartPos[[ selector[[ index ]] ]],
                                "QryEndPos"=sMap$QryEndPos[[ selector[[ index ]] ]] );
                            breakpointsVsSV[[ count ]] <- overlapItem;
                        } # for index 
                    } # if length
                } # for breakpoint            
            } # if length
        } # for chr
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in compareBreakpointsSMap" );
        save.image( file=paste0( outputFolder, "007_COMPARE_WITH_SVs-1.RData" ) );
    } else {
        cat( "Execution of compareBreakpointsSMap successfully completed\n" );
    } # if error else
    return( breakpointsVsSV );
} # compareBreakpointsSMap
                
################################################################################
convertToDataFrame <- function( breakpointsVsSV ) {
    result <- NULL;
    tryResult <- try( {
        nRows <- length( breakpointsVsSV );
        dummy <- rep( NA, nRows );
        result <- data.frame(
            "Chr"=dummy,
            "CopyNumberBreakpoint"=dummy,
            "SVType"=dummy,
            "SVChr1"=dummy,
            "SVChr2"=dummy,
            "RefStartPos"=dummy,
            "RefEndPos"=dummy,
            "QryContigID"=dummy,
            "QryStartPos"=dummy,
            "QryEndPos"=dummy );
        if ( nRows > 0 ) {
            for ( index in 1:nRows ) {
                result$Chr[[ index ]] <- breakpointsVsSV[[ index ]]$Chr;
                result$CopyNumberBreakpoint[[ index ]] <- breakpointsVsSV[[ index ]]$CopyNumberBreakpoint;
                result$SVType[[ index ]] <- as.character( breakpointsVsSV[[ index ]]$SVType );
                result$SVChr1[[ index ]] <- breakpointsVsSV[[ index ]]$SVChr1;
                result$SVChr2[[ index ]] <- breakpointsVsSV[[ index ]]$SVChr2;
                result$RefStartPos[[ index ]] <- breakpointsVsSV[[ index ]]$RefStartPos;
                result$RefEndPos[[ index ]] <- breakpointsVsSV[[ index ]]$RefEndPos;
                result$QryContigID[[ index ]] <- breakpointsVsSV[[ index ]]$QryContigID;
                result$QryStartPos[[ index ]] <- breakpointsVsSV[[ index ]]$QryStartPos;
                result$QryEndPos[[ index ]] <- breakpointsVsSV[[ index ]]$QryEndPos;
            } # for index
        } # if nRows
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in convertToDataFrame" );
        save.image( file=paste0( outputFolder, "007_COMPARE_WITH_SVs-2.RData" ) );
    } else {
        cat( "Execution of convertToDataFrame successfully completed\n" );
    } # if error else
    return( result );    
} # convertToDataFrame

################################################################################
compareCopyNumberSV <- function( breakpoints, svFolder, breakPointTolerance ) {
    breakpointsVsSV <- NULL;
    tryResult <- try( {
        outputFolder <- addSlash( outputFolder );
        if ( !testPath( svFolder ) ) {
            return( breakpointsVsSV );
        } # if !testPath
        breakpointsVsSV <- list();
        svFolder <- addSlash( svFolder );
        sMapFiles <- dir( svFolder, pattern=".smap" );
        count <- 0;
        for ( sMapFile in sMapFiles ) {
            inputFile <- paste0( svFolder, sMapFile );
            if ( testFile( inputFile ) ) {
                sMap <- readSMap( inputFile );
                breakpointsVsSV <- compareBreakpointsSMap( 
                    breakpointsVsSV, breakpoints, sMap, breakPointTolerance );
            } # if !testFile
        } # for sMapFile
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in compareCopyNumberSV" );
        save.image( file=paste0( outputFolder, "007_COMPARE_WITH_SVs-3.RData" ) );
    } else {
        cat( "Execution of compareCopyNumberSV successfully completed\n" );
    } # if error else
    return( breakpointsVsSV );
} # compareCopyNumberSV

