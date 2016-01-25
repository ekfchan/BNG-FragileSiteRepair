# FILE: UTILS.r
# AUTHOR: Zeljko Jovan Dzakula
# DATE: 7/21/2014

################################################################################
################################################################################
################################################################################
################################################################################
trim <- function( x ) gsub( "^\\s+|\\s+$", "", x )

################################################################################
isNotNULL <- function( input ) {
    success <- FALSE;
    tryResult <- try( {
        success <- TRUE;
        if ( is.null( input ) ) {
            warning( "Illegal argument: input value is NULL." );
            return( !success );
        } # if NULL
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in isNotNULL" );
        save.image( file=paste0( outputFolder, "UTILS-1.RData" ) );
    } else {
        cat( "Execution of isNotNULL successfully completed\n" );
    } # if error else
    return( success );
} # isNotNULL

################################################################################
isNotNA <- function( input ) {
    success <- FALSE;
    if ( is.null( input ) ) {
        return( success );
    } # if NULL
    tryResult <- try( {
        success <- TRUE;
        if ( is.na( input ) ) {
            warning( "Illegal argument: input value is NA." );
            return( !success );
        } # if NA
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in isNotNA" );
        save.image( file=paste0( outputFolder, "UTILS-2.RData" ) );
    } else {
        cat( "Execution of isNotNA successfully completed\n" );
    } # if error else
    return( success );
} # isNotNA

################################################################################
removeTrailingSlash <- function( path ) {
    tryResult <- try( {
        nChar <- nchar( path );
        if ( substr( path, nChar, nChar ) == "/" ) {
            path <- substr( path, 1, nChar - 1 );
        } # if /
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in removeTrailingSlash" );
        save.image( file=paste0( outputFolder, "UTILS-3.RData" ) );
    } else {
        cat( "Execution of removeTrailingSlash successfully completed\n" );
    } # if error else
    return( path );
} # removeTrailingSlash

################################################################################
removeTrailingSlashes <- function( path ) {
    tryResult <- try( {
        if ( !isNotNULL( path ) ) {
            return( path );
        } # if NULL
        if ( !isNotNA( path ) ) {
            return( path );
        } # if NA
        nChar <- nchar( path );
        while( substr( path, nChar, nChar ) == "/" ) {
            path <- removeTrailingSlash( path );
            nChar <- nchar( path );
        } # while /
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in removeTrailingSlashes" );
        save.image( file=paste0( outputFolder, "UTILS-4.RData" ) );
    } else {
        cat( "Execution of removeTrailingSlashes successfully completed\n" );
    } # if error else
    return( path );
} # removeTrailingSlashes

################################################################################
convertTilda <- function( path ) {
    tryResult <- try( {
        if ( !isNotNULL( path ) | !isNotNA( path ) ) {
            warning( "Illegal argument: path is NULL or NA" );
            return( path );
        } # if 
        nChar <- nchar( path );
        if ( nChar < 1 ) {
            warning( "Illegal argument: path is an empty string" );
            return( path );
        } # if nChar
        path <- file.path( path.expand( path ), fsep="/" );
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( paste0( "Exception thrown in convertTilda. path=", path, "\n" ) );
        save.image( file=paste0( outputFolder, "UTILS-5.RData" ) );
    } else {
        cat( "Execution of convertTilda successfully completed\n" );
    } # if error else
    return( path );
} # convertTilda

################################################################################
convertBackslashes <- function( path ) {
    tryResult <- try( {
        path <- gsub( "\\\\", "/", path );
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in convertBackslashes" );
        save.image( file=paste0( outputFolder, "UTILS-6.RData" ) );
    } else {
        cat( "Execution of convertBackslashes successfully completed\n" );
    } # if error else
    return( path );
} # convertBackslashes

################################################################################
getPathElements <- function( path ) {
    pathElements <- c();
    tryResult <- try( {
        if ( !isNotNULL( path ) | !isNotNA( path ) ) {
            warning( "Illegal argument: path is NULL or NA" );
            return( pathElements );
        } # if 
        if ( trim( as.character( path ) ) == "" ) {
            warning( "Illegal argument: path is an empty string" );
            return( pathElements );
        } # if path
        pathElements <- unlist( strsplit( path, "/" ) );
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in getPathElements" );
        save.image( file=paste0( outputFolder, "UTILS-7.RData" ) );
    } else {
        cat( "Execution of getPathElements successfully completed\n" );
    } # if error else
    return( pathElements );
} # getPathElements

################################################################################
convertDots <- function( path ) {
    tryResult <- try( {
        pathElements <- getPathElements( path );
        nPathElements <- length( pathElements );
        if ( nPathElements > 0 ) {
            currentDir <- getwd();
            if ( as.character( pathElements[[ 1 ]] ) == "." ) {
                path <- paste( currentDir, pathElements[ -1 ], sep="", collapse="/" );
            } else if ( as.character( pathElements[[ 1 ]] ) == ".." ) {
                count <- 0;
                for ( index in seq_along( pathElements ) ) {
                    if ( as.character( pathElements[[ index ]] ) == ".." ) {
                       count <- count + 1;
                    } else {
                       break;
                    } # if .. else
                } # for index
                currentPathElements <- getPathElements( currentDir );
                nCurrentElements <- length( currentPathElements );
                if ( nCurrentElements > count ) {
                    path <- paste( paste( currentPathElements[ 1:( nCurrentElements - count ) ], sep="", collapse="/" ), "/",
                        paste( pathElements[ -( 1:count ) ], sep="", collapse="/" ), "/", sep="", collapse="/" );
                } else {
                    path <- paste( pathElements[ -( 1:count ) ], sep="", collapse="/" );
                } # if nElements else
            } # if pathElement else
        } # if nPathElements 
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in convertDots" );
        save.image( file=paste0( outputFolder, "UTILS-8.RData" ) );
    } else {
        cat( "Execution of convertDots successfully completed\n" );
    } # if error else
    return( path );
} # convertDots 

################################################################################
makeDir <- function( path ) {
    tryResult <- try( {
        pathElements <- getPathElements( path );
        isWindowsPath <- ( nchar( pathElements[[ 1 ]] ) == 2 ) &
                         ( substr( pathElements[[ 1 ]], 2, 2 ) == ":" );
        if ( isWindowsPath ) {
            currentPath <- "";
        } else {
            currentPath <- "/";
        } # if isWindowsPath else
        for ( pathElement in pathElements ) {
            currentPath <- paste0( currentPath, pathElement, "/" );
            if ( !file.exists( currentPath ) ) {
                dir.create( currentPath );
            } # if !exists
        } # for pathElement
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in makeDir" );
        save.image( file=paste0( outputFolder, "UTILS-9.RData" ) );
    } else {
        cat( "Execution of makeDir successfully completed\n" );
    } # if error else
} # makeDir

################################################################################
fixPath <- function( inputPath ) {
    success <- isNotNULL( inputPath ) & isNotNA( inputPath );
    if ( !success ) {
        return( inputPath );
    } # if !success
    if ( trim( as.character( inputPath ) ) == "" ) {
        warning( "Illegal argument: inputPath is an empty string" );
        return( inputPath );
    } # if inputPath
    path <- trim( inputPath );
    tryResult <- try( {
        if ( is.null( path ) ) {
            return( path );
        } # if NULL
        path <- convertTilda( inputPath );
        path <- convertBackslashes( path );
        path <- convertDots( path );
        #path <- removeTrailingSlashes( path );
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in fixPath" );
        save.image( file=paste0( outputFolder, "UTILS-10.RData" ) );
    } else {
        cat( "Execution of fixPath successfully completed\n" );
    } # if error else
    return( path );
} # fixPath

################################################################################
testPath <- function( inputPath ) {
    success <- FALSE;
    tryResult <- try( {
        success <- isNotNULL( inputPath ) & isNotNA( inputPath );
        if ( !success ) {
            return( success );
        } # if !success
        path <- fixPath( inputPath );
        path <- removeTrailingSlashes( path );
        if ( !file.exists( path ) ) {
            warning( "Illegal argument: input path ", inputPath, " does not exist." );
            makeDir( path );
        } # if NA
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in testPath" );
        save.image( file=paste0( outputFolder, "UTILS-11.RData" ) );
    } else {
        cat( "Execution of testPath successfully completed\n" );
    } # if error else
    return( success );
} # testPath

################################################################################
addSlash <- function( word ) {
    tryResult <- try( {
        success <- isNotNULL( word ) & isNotNA( word );
        if ( success ) {
            nChar <- nchar( word );
            lastChar <- substr( word, nChar, nChar );
            if ( as.character( lastChar ) != "/" ) {
                word <- paste0( word, "/" );
            } # if !=/
        } # if !NULL
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in addSlash" );
        save.image( file=paste0( outputFolder, "UTILS-12.RData" ) );
    } else {
        cat( "Execution of addSlash successfully completed\n" );
    } # if error else
    return( word );
} # addSlash

################################################################################
testFile <- function( fileName ) {
    success <- FALSE;
    tryResult <- try( {
        success <- isNotNULL( fileName ) & isNotNA( fileName );
        if ( !file.exists( fileName ) ) {
            warning( "Illegal argument: input file ", fileName, " does not exist." );
            return( !success );
        } # if NA
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in testFile" );
        save.image( file=paste0( outputFolder, "UTILS-13.RData" ) );
    } else {
        cat( "Execution of testFile successfully completed\n" );
    } # if error else
    return( success );
} # testFile

################################################################################
testRange <- function( number, from, to, inclusive=TRUE ) {
    success <- FALSE;
    tryResult <- try( {
        success <- isNotNULL( number ) & isNotNA( number ) &
                   isNotNULL( from ) & isNotNA( from ) &
                   isNotNULL( to ) & isNotNA( to );
        if ( !success ) {
            return( success );
        } # if !success
        segment <- sort( c( from, to ) );
        if ( inclusive ) {
            success <- ( number >= segment[[ 1 ]] ) & ( number <= segment[[ 2 ]] );
        } else {
            success <- ( number > segment[[ 1 ]] ) & ( number < segment[[ 2 ]] );
        } # if inclusive else
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in testRange" );
        save.image( file=paste0( outputFolder, "UTILS-14.RData" ) );
    } else {
        cat( "Execution of testRange successfully completed\n" );
    } # if error else
    return( success );
} # testRange

################################################################################
testValues <- function( value, values ) {
    success <- FALSE;
    tryResult <- try( {
        success <- isNotNULL( value ) & isNotNA( value ) & 
                   isNotNULL( values ) & isNotNA( values );
        if ( !success ) {
            return( success );
        } # if !success
        success <- length( which( values == value ) ) > 0;
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in testValues" );
        save.image( file=paste0( outputFolder, "UTILS-15.RData" ) );
    } else {
        cat( "Execution of testValues successfully completed\n" );
    } # if error else
    return( success );
} # testValues
