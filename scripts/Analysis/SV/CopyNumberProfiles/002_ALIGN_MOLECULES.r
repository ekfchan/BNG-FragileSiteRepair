# FILE: 002_ALIGN_MOLECULES.r
# AUTHOR: Zeljko Jovan Dzakula
# DATE: 7/21/2014

################################################################################
################################################################################
################################################################################
################################################################################
# The following script assumes that SCRIPTS/UTILS.r was sourced
alignMolsVsReference <- function( 
    inputBnx,
    outputFolder,
    outputFilePrefix,
    refAligner="/home/users/wandrews/execs/2940/RefAligner", 
    reference="/home/users/wandrews/analysis/human/HG_b19_chromosome_bspq.cmap",
    minLength=150,
    tValue=9,
    flags=paste0( "-sf .25 -sd .11 -minsites 9 -biaswt 0 -S -1000 ", 
                  "-res 3.3 -resSD 0.75 -outlier 1e-4 -M 5 -f -mres 2.0 ",
                  "-hashgen 5 3 2.4 1.4 0.05 5.0 1 1 -hash -hashdelta 10 ", 
                  "-Hash_Bits 23 -MHash_Bits 25 -HHash_Bits 18 ", 
                  "-insertThreads 32 -nosplit 2 ", 
                  "-queryThreads 32 -extend 1 -maxthreads 64 -BestRef 1 ", 
                  "-stdout -stderr" )
    ) {
    success <- FALSE;
    tryResult <- try( {
        outputFolder <- addSlash( outputFolder );
        success <- testFile( inputBnx ) & testPath( outputFolder );
        if ( !success ) {
            return( success );
        } # if !success
        output <- paste0( outputFolder, "/", outputFilePrefix );
        actualTValue <- 10^( -tValue );
        commandLine=paste0( refAligner, " ",
            "-ref ", reference, " ",
            "-i ", inputBnx, " ",
            "-o ", output, " ", 
            "-minlen ", minLength, " ",
            "-T ", actualTValue, " ",
            flags );
        system( commandLine );
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in alignMolsVsReference" );
        save.image( file=paste0( outputFolder, "002_ALIGN_MOLECULES-1.RData" ) );
    } else {
        cat( "Execution of alignMolsVsReference successfully completed\n" );
    } # if error else
    return( success );
} # alignMolsVsReference

 
