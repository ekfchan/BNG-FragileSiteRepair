# PATH: C:\BIONANO_PROJECTS\SCRIPTS\R\ANALYSIS\
# FILE: copyNumberCluster.r
# AUTHOR: Zeljko Jovan Dzakula
# DATE: 3/1/2014

################################################################################
################################################################################
################################################################################
################################################################################
evaluateRoughClusterProfile <- 
    function( inputProfile, nClusters, weight=nClusters, 
              maxNIterations=100, tolerance=1e-6, trueProfile=NULL, plotFlag=FALSE ) {
    nPoints <- length( inputProfile );
    centroidX <- sort( round( runif( n=nClusters, min=1, max=nPoints ) ) );
    centroidY <- round( runif( n=nClusters, min=0, max=ceiling( max( inputProfile ) ) ) );
    yRange <- diff( range( inputProfile ) );
    previousPenalty <- 6.22e23;
    for ( iteration in 1:maxNIterations ) {
        penalty <- 0;
        # Assign points to centroids
        distances <- matrix( nrow=nPoints, ncol=nClusters );
        assignment <- rep( NA, nPoints );
        for ( point in 1:nPoints ) {
            for ( centroid in 1:nClusters ) {
                deltaX <- weight * ( centroidX[ centroid ] - point ) / nPoints;
                deltaY <- ( centroidY[ centroid ] - inputProfile[ point ] ) / yRange;
                distances[ point, centroid ] <- deltaX * deltaX + deltaY * deltaY;
            } # for centroid
            bestCentroid <- which.min( distances[ point, ] );
            index <- runif( n=1, min=1, max=length( bestCentroid ) );
            bestCentroid <- bestCentroid[ index ];
            assignment[ point ] <- bestCentroid;
            penalty <- penalty + distances[ point, bestCentroid ];
        } # for point
        # Reevaluate centroids
        for ( centroid in 1:nClusters ) {
            currentAssignment <- which( assignment == centroid );
            centroidX[ centroid ] <- median( currentAssignment, na.rm=TRUE );
            centroidY[ centroid ] <- round( median( inputProfile[ currentAssignment ], na.rm=TRUE ) );
        } # for centroid
        reset <- which( is.na( centroidX ) );
        if ( length( reset ) > 0 ) {
            centroidX[ reset ] <- round( runif( n=length( reset ), min=1, max=nPoints ) );
            centroidY[ reset ] <- round( runif( n=length( reset ), min=0, max=12 ) );
            indices <- order( centroidX );
            centroidX <- centroidX[ indices ];
            centroidY <- centroidY[ indices ];
        } # if eliminate
        cat( paste0( "Iteration=", iteration,
                     ", Penalty=", penalty,
                     ", Clusters=", nClusters, "\n" ) );
        if ( abs( penalty - previousPenalty ) < tolerance ) {
            break();
        } # if penalty
        previousPenalty <- penalty;
        #
        if ( plotFlag == TRUE ) {
            par( mfrow=c( 4, 1 ) );
            plot( centroidY[ assignment ], pch=16, ylim=c( 0, 7 ) );
            abline( h=0:12, col="green" );
            if ( !is.null( trueProfile ) ) {
                lines( trueProfile, col="red", lwd=3 );
            } # if !null
            #
            plot( inputProfile, type="l", ylim=c( 0, 7 ) );
            abline( h=0:12, col="green" );
            if ( !is.null( trueProfile ) ) {
                lines( trueProfile, col="red", lwd=3 );
            } # if !null
            #
            plot( inputProfile - centroidY[ assignment ], type="l" );
            abline( h=-10:12, col="green" );
            #
            deviation <- inputProfile - centroidY[ assignment ];
            z <- ( deviation - median( deviation ) ) / mad( deviation );
            pValues <- 1 - pnorm( z );
            logP <- -log10( pValues );
            colorValues <- rep( rgb( 0.8, 0.8, 0.8, 0.5 ), nPoints );
            cutoff <- median( logP ) + 6 * mad( logP );
            colorValues[ logP > cutoff ] <- "red";
            plot( logP, col=colorValues, pch=16 );
            abline( h=cutoff, col="green" );
        } # if plotFlag
    } # for iteration
    #
    elevation <- centroidY[ assignment ];
    breakPoints <- findBreakPoints( elevation );
    clusterProfile <- list(
        "ClusterPosition"=centroidX,
        "ClusterElevation"=centroidY,
        "ClusterAssignment"=assignment,
        "Elevation"=elevation,
        "BreakPoints"=breakPoints
        );
    plotClusterCopyNumber( inputProfile, clusterProfile, trueProfile );
    return( clusterProfile );
} # evaluateRoughClusterProfile

################################################################################
################################################################################
################################################################################
################################################################################
medianClusterAssignment <- function( inputAssignment, windowSize=3 ) {
    nPoints <- length( inputAssignment );
    assignment <- rep( NA, nPoints );
    flankSize <- round( ( windowSize - 1 ) / 2 );
    for ( point in 1:nPoints ) {
        lowerBound <- max( 1, point - flankSize );
        upperBound <- min( nPoints, point + flankSize );
        assignment[ point ] <- round( median( inputAssignment[ lowerBound:upperBound ], na.rm=TRUE ) );
    } # for point
    outliers <- which( is.na( match( unique( inputAssignment ), assignment ) ) );
    if ( length( outliers ) > 0 ) {
        cat( paste( paste( "Outliers: ", outliers, collapse=" " ), "\n" ) );
    } # if length
    return( assignment );
} # medianClusterAssignment

################################################################################
################################################################################
################################################################################
################################################################################
fixEdges <- function( inputAssignment ) {
    assignment <- inputAssignment;
    indices <- unique( assignment );
    nPoints <- length( indices );
    for ( index in 1:( nPoints - 1 ) ) {
        maxPosition1 <- max( which( assignment == index ) );
        minPosition2 <- min( which( assignment == ( index + 1 ) ) );
        if ( minPosition2 < maxPosition1 ) {
            center <- mean( c( minPosition2, maxPosition1 ) );
            newMaxPosition1 <- floor( center );
            newMinPosition2 <- ceiling( center );
            if ( newMaxPosition1 == newMinPosition2 ) {
                if ( rbinom( n=1, size=1, p=0.5 ) == 1 ) {
                    newMinPosition2 <- newMinPosition2 + 1;
                } else {
                    newMaxPosition1 <- newMaxPosition1 - 1;
                } # if 1 else 
            } # if ==
            assignment[ minPosition2:newMaxPosition1 ] <- index;
            assignment[ newMinPosition2:maxPosition1 ] <- index + 1;
        } # if min2 < max1
    } # for index
    return( assignment );
} # fixEdges

################################################################################
################################################################################
################################################################################
################################################################################
findBreakPoints <- function( elevation ) {
    nPoints <- length( elevation );
    breakPoints <- c();
    for ( point in 2:( nPoints - 1 ) ) {
        elevation[ point ] <- median( elevation[ point + -1:1 ] );
        if ( elevation[ point ] != elevation[ point - 1 ] ) {
            breakPoints <- c( breakPoints, point );
        } # if elevation
    } # for point
    return( sort( breakPoints ) );
} # findBreakPoints

################################################################################
################################################################################
################################################################################
################################################################################
evaluateClusterCopyNumber <- function( 
    inputProfile, nClusters, weight=nClusters, maxNIterations=100, tolerance=1e-6,
    trueProfile=NULL, plotFlag=FALSE ) {
    clusterProfile <- evaluateRoughClusterProfile( 
        inputProfile, nClusters, weight=weight, maxNIterations=maxNIterations, tolerance=tolerance,
        trueProfile=trueProfile, plotFlag=plotFlag );
    #
    clusterAssignment <- medianClusterAssignment( clusterProfile[[ "ClusterAssignment" ]], windowSize=3 );
    clusterAssignment <- fixEdges( clusterAssignment );
    names( clusterAssignment ) <- names( inputProfile );
    clusterProfile[[ "ClusterAssignment" ]] <- clusterAssignment;
    #
    elevation <- clusterProfile[[ "ClusterElevation" ]][ clusterProfile[[ "ClusterAssignment" ]] ];
    names( elevation ) <- names( inputProfile );
    clusterProfile[[ "Elevation" ]] <- elevation;
    #
    breakPoints <- findBreakPoints( clusterProfile[[ "Elevation" ]] );
    if ( !is.null( breakPoints ) ) {
        names( breakPoints ) <- names( inputProfile )[ breakPoints ];
    } # if !null
    clusterProfile[[ "BreakPoints" ]] <- breakPoints;
    return( clusterProfile );
} # evaluateClusterCopyNumber

################################################################################
################################################################################
################################################################################
################################################################################
plotClusterCopyNumber <- function( inputProfile, clusterProfile, trueProfile=NULL ) {
    par( mfrow=c( 1, 1 ) );
    plot( inputProfile );
    abline( h=0:100, col="green" );
    lines( clusterProfile[[ "Elevation" ]], col="blue" );
    if ( !is.null( trueProfile ) ) {
        lines( trueProfile, col="red" );
    } # if !null
    abline( v=clusterProfile[[ "BreakPoints" ]], col="blue", lty=2 );
} # plotClusterCopyNumber

################################################################################
################################################################################
################################################################################
################################################################################
# Histeresis:
#for(i in 2:length(z)) { if(abs(w[i-1]-y[i])<0.7) { w[i]<-w[i-1] } else { w[i]<-round(y[i]) } }

################################################################################
################################################################################
################################################################################
################################################################################
runningMedian <- function( inputProfile, firstWindow=11, secondWindow=21 ) {
    noNaProfile <- inputProfile[ !is.na( inputProfile ) ];
    elevation <- runmed( round( runmed( noNaProfile, firstWindow ) ), secondWindow );
    names( elevation ) <- names( noNaProfile );
    breakPoints <- findBreakPoints( elevation );
    if ( !is.null( breakPoints ) ) {
        names( breakPoints ) <- names( inputProfile )[ breakPoints ];
    } # if !null
    result <- list(
        "Elevation"=elevation,
        "BreakPoints"=breakPoints );
    return( result );
} # runningMedian

