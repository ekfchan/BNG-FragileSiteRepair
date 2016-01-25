overlap <- function( block1, block2 ) {
    block1 <- range( block1 );
    block2 <- range( block2 );
    return( ( block1[[ 1 ]] < block2[[ 2 ]] ) &
            ( block1[[ 2 ]] > block2[[ 1 ]] ) );
} # overlap

overlaps <- function( block1DataFrame, block2 ) {
    nBlocks <- ncol( block1DataFrame );
    result <- rep( FALSE, nBlocks );
    for ( block1 in 1:nBlocks ) {
        result[[ block1 ]] <- overlap( block1DataFrame[ , block1 ], block2 );
    } # for block1
    return( result );
} # overlaps



