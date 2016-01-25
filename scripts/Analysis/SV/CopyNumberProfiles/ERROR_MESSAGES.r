# FILE: ERROR_MESSAGES.r
# AUTHOR: Zeljko Jovan Dzakula
# DATE: 7/21/2014

################################################################################
################################################################################
################################################################################
################################################################################
referenceError <- paste0( "Illegal argument - unknown reference: ", reference, ". ",
                         "Currently, the only supported valuse include", referenceChoices, "." );
                         #"Valid values are hg19, hg38, and BioNano" );
binSizeError <- paste0( "Illegal argument - unknown bin size choice: ", binSize, ". ",
                    "Currently, the only supported values include ", binSizeChoices, "." );
                    #"Valid values are hg19, hg38, and BioNano" );

