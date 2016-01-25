# FILE: 008_GENERATE_REPORT.r
# AUTHOR: Zeljko Jovan Dzakula
# DATE: 7/2/2014

################################################################################
exportCnFile <- function( bins, copyNumberProfile, breakpoints, mask, outputFolder, fileNamePrefix ) {
    tryResult <- try( {
        outputFolder <- addSlash( outputFolder );
        if ( !testPath( outputFolder ) ) {
            warning( paste0( "Illegal argument: outputFolder ", outputFolder, " does not exist." ) );
            return;
        } # if !testPath
        if ( !isNotNULL( fileNamePrefix ) | !isNotNA( fileNamePrefix ) ) {
            warning( paste0( "Illegal argument: fileNamePrefix=", fileNamePrefix ) );
            return;
        } # if !fileNamePrefix
        elevations <- rep( NA, nrow( bins ) );
        names( elevations ) <- as.character( rownames( bins ) );
        for ( chr in 1:24 ) {
            binNames <- names( breakpoints[[ chr ]]$RunningMedian$Elevation );
            elevations[ as.character( binNames ) ] <- breakpoints[[ chr ]]$RunningMedian$Elevation;
        } # for chr
        cnExtension <- ".cn";
        tmp <- cbind( "SNP"=as.character( rownames( bins ) ),
                      "Chromosome"=as.integer( bins$Chr ), 
                      "PhysicalPosition"=as.integer( 0.5 * ( bins$Start + bins$End ) ), 
                      "RawCounts"=as.numeric( copyNumberProfile$RawCounts ), 
                      "CopyNumber"=as.numeric( copyNumberProfile$CopyNumber ), 
                      "SmoothCopyNumber"=as.numeric( copyNumberProfile$SmoothCopyNumber ) * mask, 
                      "BlockElevation"=as.numeric( elevations ),
                      "Mask"=as.integer( mask ) );
        rownames( tmp ) <- NULL;
        selector <- which( is.na( mask ) );
        tmp <- tmp[ -selector, ];
        write.table( tmp, file=paste0( outputFolder, fileNamePrefix, cnExtension ), 
                     sep="\t", quote=FALSE, row.names=FALSE );
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in exportCnFile" );
        save.image( file=paste0( outputFolder, "008_GENERATE_REPORT-1.RData" ) );
    } else {
        cat( "Execution of exportCnFile successfully completed\n" );
    } # if error else
} # exportCnFile

################################################################################
printPerChrSummary <- function( copyNumberSummary, section ) {
    tryResult <- try( {
        tmp <- copyNumberSummary[[ as.character( section ) ]];
        cat( "..............................................\n" );
        cat( paste0( section, " Summary:\n" ) );
        cat( paste0( "SDAutosomal=", round( tmp$SDAutosomal, 4 ), "\n" ) );
        cat( paste0( "MADAutosomal=", round( tmp$MADAutosomal, 4 ), "\n" ) );
        cat( paste0( "Chr\tSDPerChr", "\t", "MADPerChr", "\t", "MeanPerChr", "\t", "MedianPerChr", "\t", "CV\n" ) );
        tmp <- tmp$StatsPerChr;
        for ( chr in colnames( tmp ) ) {
            cat( paste0( as.character( chr ), "\t", 
                         round( tmp[ "SDPerChr", as.character( chr ) ], 4 ), "\t", 
                         round( tmp[ "MADPerChr", as.character( chr ) ], 4 ), "\t", 
                         round( tmp[ "MeanPerChr", as.character( chr ) ], 4 ), "\t", 
                         round( tmp[ "MedianPerChr", as.character( chr ) ], 4 ), "\t",
                         round( tmp[ "SDPerChr", as.character( chr ) ] / tmp[ "MeanPerChr", as.character( chr ) ], 4 ), "\n" ) );
        } # for chr
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in printPerChrSummary" );
        save.image( file=paste0( outputFolder, "008_GENERATE_REPORT-2.RData" ) );
    } else {
        #cat( "Execution of printPerChrSummary successfully completed\n" );
    } # if error else
} # printPerChrSummary 

################################################################################
exportCopyNumberSummary <- function( copyNumberProfile, copyNumberSummary, outputFolder, fileNamePrefix ) {
    tryResult <- try( {
        outputFolder <- addSlash( outputFolder );
        if ( !testPath( outputFolder ) ) {
            warning( paste0( "Illegal argument: outputFolder ", outputFolder, " does not exist." ) );
            return;
        } # if !testPath
        if ( !isNotNULL( fileNamePrefix ) | !isNotNA( fileNamePrefix ) ) {
            warning( paste0( "Illegal argument: fileNamePrefix=", fileNamePrefix ) );
            return;
        } # if !fileNamePrefix
        sink( paste0( outputFolder, fileNamePrefix, ".summary" ) );
        cat( paste0( "TotalAutosomalCounts=", round( copyNumberProfile$TotalAutosomalCounts, 4 ), "\n" ) ); 
        cat( paste0( "LabelDensityBias=", round( copyNumberProfile$LabelDensityBias, 40 ), "\n" ) );
        cat( paste0( "CountsPerBin=", round( copyNumberSummary$CountsPerBin, 4 ), "\n" ) );          
        cat( paste0( "ExpectedRelativeError=", round( copyNumberSummary$ExpectedRelativeError, 4 ), "\n" ) );
        printPerChrSummary( copyNumberSummary, "Scaled" );
        printPerChrSummary( copyNumberSummary, "CopyNumber" );
        printPerChrSummary( copyNumberSummary, "Smooth" );
    } ); # try
    sink();
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in exportCopyNumberSummary" );
        save.image( file=paste0( outputFolder, "008_GENERATE_REPORT-3.RData" ) );
    } else {
        cat( "Execution of exportCopyNumberSummary successfully completed\n" );
    } # if error else
} # exportCopyNumberSummary

################################################################################
exportCnSvBreakpoints <- function( breakpointsVsSV, outputFolder, fileNamePrefix ) {
    tryResult <- try( {
        outputFolder <- addSlash( outputFolder );
        if ( !testPath( outputFolder ) ) {
            warning( paste0( "Illegal argument: outputFolder ", outputFolder, " does not exist." ) );
            return;
        } # if !testPath
        if ( !isNotNULL( fileNamePrefix ) | !isNotNA( fileNamePrefix ) ) {
            warning( paste0( "Illegal argument: fileNamePrefix=", fileNamePrefix ) );
            return;
        } # if !fileNamePrefix
        if ( is.null( breakpointsVsSV ) ) {
            warning( "Illegal argument: breakpointsVsSV is NULL" );
            return;
        } # if is.null
        if ( as.character( class( breakpointsVsSV ) ) != "data.frame" ) {
            warning( "Illegal argument: breakpointsVsSV has wrong type" );
            return;
        } # if !data.frame
        write.table( breakpointsVsSV, file=paste0( outputFolder, fileNamePrefix, "_COPY_NUMBER_vs_SV.txt" ), 
                     sep="\t", quote=FALSE, row.names=FALSE );
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in exportCnSvBreakpoints" );
        save.image( file=paste0( outputFolder, "008_GENERATE_REPORT-4.RData" ) );
    } else {
        cat( "Execution of exportCnSvBreakpoints successfully completed\n" );
    } # if error else
} # exportCnSvBreakpoints

################################################################################
exportAutosomalAneuploidy <- function( autosomalAneuploidy, outputFolder, fileNamePrefix ) {
    tryResult <- try( {
        outputFolder <- addSlash( outputFolder );
        if ( !testPath( outputFolder ) ) {
            warning( paste0( "Illegal argument: outputFolder ", outputFolder, " does not exist." ) );
            return;
        } # if !testPath
        if ( !isNotNULL( fileNamePrefix ) | !isNotNA( fileNamePrefix ) ) {
            warning( paste0( "Illegal argument: fileNamePrefix=", fileNamePrefix ) );
            return;
        } # if !fileNamePrefix
        if ( !is.null( autosomalAneuploidy ) ) {
            sink( paste0( outputFolder, fileNamePrefix, "_AUTOSOMAL_ANEUPLOIDY.txt" ) );
            cat( "Type\tChromosome\tZValue\n" );
            if ( !is.null( autosomalAneuploidy$Monosomy ) ) {
            if ( length( autosomalAneuploidy$Monosomy ) > 0 ) {
                for ( count in 1:length( autosomalAneuploidy$Monosomy ) ) {
                    chr <- names( autosomalAneuploidy$Monosomy )[[ count ]];
                    nChar <- nchar( chr );
                    if ( nChar > 3 ) {
                        chr <- substr( chr, 4, nChar );
                        cat( paste0( "Monosomy", "\t", chr, "\t", 
                            autosomalAneuploidy$Monosomy[[ count ]], "\n" ) );
                    } # if nChr
                } # for count
            } # if length
            } # if !NULL
            if ( !is.null( autosomalAneuploidy$Polysomy ) ) {
            if ( length( autosomalAneuploidy$Polysomy ) > 0 ) {
                for ( count in 1:length( autosomalAneuploidy$Polysomy ) ) {
                    chr <- names( autosomalAneuploidy$Polysomy )[[ count ]];
                    nChar <- nchar( chr );
                    if ( nChar > 3 ) {
                        chr <- substr( chr, 4, nChar );
                        cat( paste0( "Polysomy", "\t", chr, "\t", 
                            autosomalAneuploidy$Polysomy[[ count ]], "\n" ) );
                    } # if nChar
                } # for count
            } # if length
            } # if !NULL
        } # !null
    } ); # try
    sink();
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in exportAutosomalAneuploidy" );
        save.image( file=paste0( outputFolder, "008_GENERATE_REPORT-5.RData" ) );
    } else {
        cat( "Execution of exportAutosomalAneuploidy successfully completed\n" );
    } # if error else
} # exportAutosomalAneuploidy



