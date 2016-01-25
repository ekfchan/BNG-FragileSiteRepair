# FILE: 001_WORKFLOW.r
# AUTHOR: Zeljko Jovan Dzakula
# DATE: 7/2/2014
#
################################################################################
################################################################################
################################################################################
################################################################################
runCopyNumber <- function(
    performAlignment=FALSE, 
    inputMoleculesFile=NA, # No default, must be BNX 
    reference=referenceChoices[[ 1 ]], # Default=hg19; Possible values: hg19, hg38, BioNano
    minLength=minLengthChoices[[ 1 ]], # Default=150; Possible values: 100, 125, 150, 180
    tValue=tValueChoices[[ 1 ]], # Default=9; Possible values: 7, 9, 11
    refAligner=NULL,
    binSize=binSizeChoices[[ 1 ]], # Possible values: 50kbp, 200kbp
    performScaling=TRUE, 
    detectAneuploidy=TRUE,
    detectSubchromosomalFeatures=TRUE,
    compareWithSVs=TRUE,
    breakPointTolerance=2.5e4, 
    scriptFolder=NULL,
    copyNumberPackageFolder=NULL,
    outputFolder,
    outputFilePrefix,
    svFolder,
    mapMoleculeFile
    ) { 
    tryResult <- try( {
        ############################################################################
        # Validate input arguments
        scriptFolder <- fixPath( scriptFolder );
        copyNumberPackageFolder <- fixPath( copyNumberPackageFolder );
        outputFolder <- fixPath( outputFolder );
        svFolder <- fixPath( svFolder );
        inputMoleculesFile <- fixPath( inputMoleculesFile );
        mapMoleculeFile <- fixPath( mapMoleculeFile );
        
        ############################################################################
        # Process input arguments
        source( paste0( copyNumberPackageFolder, "/ERROR_MESSAGES.r" ) );
        #
        if ( as.character( reference ) == "hg19" ) {
            if ( as.character( binSize ) == "50kbp" ) {
                # Do nothing
            } else if ( as.character( binSize ) == "200kbp" ) {
                stop( binSizeError );
            } else {
                stop( binSizeError );
            } # if binSize else
        } else if ( as.character( reference ) == "hg38" ) {
            stop( referenceError );
        } else if ( as.character( reference ) == "BioNano" ) {
            stop( referenceError );
        } else {
            stop( referenceError );
        } # if referenceChoice else
        binFile <- paste0( copyNumberPackageFolder, "/DATA/", reference, "/", binSize, "/bins_", 
                           reference, "_", binSize, ".RData" ); 
        binLabelsFile <- paste0( copyNumberPackageFolder, "/DATA/", reference, "/", binSize, "/binLabels_", 
                           reference, "_", binSize, ".RData" ); 
        binParametersFile <- paste0( copyNumberPackageFolder, 
            "/DATA/", reference, "/", binSize, "/binParameters_", reference, "_", binSize, 
            "_MinLength_", minLength, "_T_", tValue, ".RData" ); 
        maskFile <- paste0( copyNumberPackageFolder, 
            "/DATA/", reference, "/", binSize, "/mask_", reference, "_", binSize, 
            "_MinLength_", minLength, "_T_", tValue, ".RData" ); 
        referenceFile <- paste0( copyNumberPackageFolder, "/DATA/", reference, "/", reference, "_chromosome_bspq.cmap" );
        chrSelectorFile <- paste0( copyNumberPackageFolder, 
            "/DATA/", reference, "/", binSize, "/chrSelector_", reference, "_", binSize, ".RData" );
        #
        load( binFile );
        load( binLabelsFile );
        load( binParametersFile );
        load( maskFile );
        load( chrSelectorFile );
        #
        ############################################################################
        # Start optimistically:
        success <- TRUE;
        
        ############################################################################
        # ALIGN MOLECULES:
        if ( performAlignment ) {
            success <- success & alignMolsVsReference( 
                inputBnx=inputMoleculesFile,
                outputFolder=outputFolder,
                outputFilePrefix=outputFilePrefix,
                refAligner=refAligner, 
                reference=referenceFile,
                minLength=minLength,
                tValue=tValue );
            if ( !success ) {
                warning( "Mapping of single molecules to reference failed." );
                return( success );
            } # if !success
            mapMoleculeFile <- paste0( outputFolder, "/", outputFilePrefix, "_r.cmap" );
        } # if performAlignment
        if ( !testFile( mapMoleculeFile ) ) {
            stop( paste0( "Fatal error: mapMoleculeFile=", mapMoleculeFile ) );
        } # if !mapMoleculeFile
        #    
        ############################################################################
        # NORMALIZE COUNTS:
        copyNumberProfile <- processCounts( mapMoleculeFile, bins, binLabels, binParameters, mask, chrSelector );
        #
        ############################################################################
        # Determine scaling factor (integer chr elevation)
        if ( performScaling ) {
            copyNumberProfile$CopyNumber <- 
                scaleCopyNumberProfile( copyNumberProfile$CopyNumber, mask, chrSelector );
            copyNumberProfile$SmoothCopyNumber <- 
                scaleCopyNumberProfile( copyNumberProfile$SmoothCopyNumber, mask, chrSelector );
        } # if performScaling
        #
        ############################################################################
        # PROCESS COPY NUMBERS:
        copyNumberSummary <- profileSummary( copyNumberProfile, mask, chrSelector );
        #
        ############################################################################
        # DETECT ANEUPLOIDY:
        if ( detectAneuploidy ) {
            autosomalAneuploidy <- classifyAutosomalAneuploidy( copyNumberProfile, mask, chrSelector, cutoff=10 );
            exportAutosomalAneuploidy( autosomalAneuploidy, outputFolder, fileNamePrefix=outputFilePrefix );
        } # if detectAneuploidy
        #
        ############################################################################
        # DETECT SUBCHROMOSOMAL FEATURES:
        # Always detect breakpoints, regardless of the flag detectSubchromosomalFeatures
        # The flag detectSubchromosomalFeatures is deprecated
        # Reason: exportCnFile requires breakpoints so IrysView can draw copy number blocks
        if ( TRUE ) { # detectSubchromosomalFeatures ) {
            breakpoints <- detectBreakpoints( copyNumberProfile$SmoothCopyNumber, 
                chrSelector, mask, scaleFactor=1, performClustering=FALSE );
            breakpoints <- assignPValues( 
                copyNumberProfile$SmoothCopyNumber, breakpoints, chrSelector, mask );
            # This is correct but useless:
            #breakpoints <- refineBreakpoints( 
            #    copyNumberProfile$SmoothCopyNumber, breakpoints, chrSelector );
            breakpoints <- unwrapGenomicCoordinates( breakpoints, binSize=5e4 );
            breakpoints <- unwrapElevations( copyNumberProfile, chrSelector, breakpoints );
            breakpoints <- defineBlocks( breakpoints, binSize=5e4 );
            save( list="breakpoints", file=paste0( addSlash( outputFolder ), "breakpoints.RData" ) );
            #
            ############################################################################
            # COMPARE WITH SVs:
            if ( compareWithSVs ) {
                breakpointsVsSV <- compareCopyNumberSV( breakpoints, svFolder, breakPointTolerance );
                breakpointsVsSV <- convertToDataFrame( breakpointsVsSV );
                exportCnSvBreakpoints( breakpointsVsSV, outputFolder, fileNamePrefix=outputFilePrefix );
            } else {
                # Do nothing
            } # if compareWithSVs
        } # if detectSubchromosomalFeatures
        #
        ############################################################################
        # GENERATE REPORT:
        exportCnFile( bins, copyNumberProfile, breakpoints, mask, outputFolder, fileNamePrefix=outputFilePrefix );
        exportCopyNumberSummary( copyNumberProfile, copyNumberSummary, outputFolder, fileNamePrefix=outputFilePrefix );
        #
    } ); # try
    if ( class( tryResult ) == "try-error" ) {
        warning( "Exception thrown in runCopyNumber" );
        save.image( file=paste0( outputFolder, "001_WORKFLOW-1.RData" ) );
    } else {
        warnings();
        cat( "Execution of runCopyNumber successfully completed\n" );
    } # if error else
} # runCopyNumber

