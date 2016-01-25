#http://www.biostars.org/p/9335/

matcher <- function(pattern, x) {
    return(unlist(gregexpr(pattern, x)))
}

doOne <- function(inputCharacter, cigar) {
    pattern <- inputCharacter
    as.numeric(matcher(pattern, cigar))
}

parseCigar <- function(cigar, chars=c("M","D","I")) {
    temp <- sapply (chars, doOne, cigar)
    temp <- unlist(temp)
    temp <- temp[temp > 0]
    temp <- sort(temp)
    
    tempChars <- substr(as.character(names(temp)), 1, 1)
    nItems <- length(temp)
    result <- data.frame("letter"=rep(NA, nItems),
                         "length"=rep(0, nItems))
    
    previousIndex <- 1
    for (index in 1:nItems) {
        result[index, "letter"] <- as.character(tempChars[index])
        currentIndex <- temp[index] - 1
        howMany <- substr(cigar, previousIndex, currentIndex)
        result[index, "length"] <- as.numeric(howMany)
        previousIndex <- currentIndex + 2
    }
    return(result)
} # parseCigar

enumerateLabels <- function( singleMoleculeCMap, contigCMap, xMap, colname="ContigSiteID" ) {
    result <- rep( NA, nrow( singleMoleculeCMap ) );
    moleculeIDs <- unique( singleMoleculeCMap[ , "CMapId" ] );
    for ( molID in moleculeIDs ) {
        #print( molID );
        xMapIndex <- which( xMap[ , "QryContigID" ] == molID );
        if ( !is.null( xMapIndex ) & ( !length( xMapIndex ) < 1 ) ) {
            queryStart <- xMap[ xMapIndex, "QryStartPos" ];
            queryEnd <- xMap[ xMapIndex, "QryEndPos" ];
            contigStart <- xMap[ xMapIndex, "RefStartPos" ];
            contigEnd <- xMap[ xMapIndex, "RefEndPos" ];
            orientation <- xMap[ xMapIndex, "Orientation" ];
            confidence <- xMap[ xMapIndex, "Confidence" ];
            cigar <- xMap[ xMapIndex, "HitEnum" ];
            cigarItems <- parseCigar( cigar );
            #
            labelIndices <- which( singleMoleculeCMap[ , "CMapId" ] == molID );
            labelPositions <- singleMoleculeCMap[ labelIndices, "Position" ];
            contigPositions <- contigCMap[ , "Position" ];
            nContigLabels <- length( contigPositions );
            #
            firstContigLabel <- which( contigPositions >= contigStart )[ 1 ];
            lastContigLabel <- which( contigPositions >= contigEnd )[ 1 ] - 1;
            currentContigLabel <- firstContigLabel;
            if ( as.character( orientation ) == "+" ) {
                firstQueryLabel <- labelIndices[ which( labelPositions >= queryStart )[ 1 ] ];
                lastQueryLabel <- labelIndices[ which( labelPositions >= queryEnd )[ 1 ] - 1 ];
                currentLabel <- firstQueryLabel;
                for ( cigarItem in 1:nrow( cigarItems ) ) {
                    code <- cigarItems[ cigarItem, "letter" ];
                    if ( as.character( code ) == "M" ) {
                        for ( count in 1:cigarItems[ cigarItem, "length" ] ) {
                            result[ currentLabel ] <- currentContigLabel;
                            currentLabel <- currentLabel + 1;
                            currentContigLabel <- currentContigLabel + 1;
                        } # for count
                    } else if ( as.character( code ) == "D" ) {
                        for ( count in 1:cigarItems[ cigarItem, "length" ] ) {
                            currentContigLabel <- currentContigLabel + 1;
                        } # for count
                    } else if ( as.character( code ) == "I" ) {
                        for ( count in 1:cigarItems[ cigarItem, "length" ] ) {
                            currentLabel <- currentLabel + 1;
                        } # for count
                    } else {
                        warning( paste0( "Unrecognized cigar code: ", code, "." ) );
                    } # if code else
                } # for cigarItem
            } else if ( as.character( orientation ) == "-" ) {
                firstQueryLabel <- labelIndices[ which( labelPositions >= queryEnd )[ 1 ] ];
                lastQueryLabel <- labelIndices[ which( labelPositions >= queryStart )[ 1 ] ];
                currentLabel <- lastQueryLabel;
                for ( cigarItem in 1:nrow( cigarItems ) ) {
                    code <- cigarItems[ cigarItem, "letter" ];
                    if ( as.character( code ) == "M" ) {
                        for ( count in 1:cigarItems[ cigarItem, "length" ] ) {
                            result[ currentLabel ] <- currentContigLabel
                            currentLabel <- currentLabel - 1;
                            currentContigLabel <- currentContigLabel + 1;
                        } # for count
                    } else if ( as.character( code ) == "D" ) {
                        for ( count in 1:cigarItems[ cigarItem, "length" ] ) {
                            currentContigLabel <- currentContigLabel + 1;
                        } # for count
                    } else if ( as.character( code ) == "I" ) {
                        for ( count in 1:cigarItems[ cigarItem, "length" ] ) {
                            currentLabel <- currentLabel - 1;
                        } # for count
                    } else {
                        warning( paste0( "Unrecognized cigar code: ", code, "." ) );
                    } # if code else
                } # for cigarItem
            } else {
                warning( paste0( "Erroneous input: orientation can only be + or -." ) );
            } # if orientation
        } # if !null
    } # for molID
    result <- cbind( singleMoleculeCMap, result );
    #colnames( result )[ ncol( result ) ] <- "ContigSiteID";
    colnames( result )[ ncol( result ) ] <- colname; #EL - 05122014
    return( result );
} # enumerateLabels

