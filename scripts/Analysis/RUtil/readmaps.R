
#utility fns, mostly for reading BioNano file formats (cmap, xmap, smap)


paste0 <- function( ... ) { paste( ..., sep="" ) }


addTrailingSlash <- function(instr) {
  if( substr(instr, nchar(instr), nchar(instr) ) != "/" )
    return( paste0( instr, "/") )
  else
    return( instr )
}


readCMap_OLD <- function(filename) {
    cmap <- NULL
    try( {
        
        #message(filename)
        cmap <- read.table(filename, comment.char="#", header=FALSE, stringsAsFactors=FALSE, nrows=50, fill = TRUE) #adjust nrows as needed
        cmap <- as.data.frame(scan(filename, what=cmap, comment.char="#", fill = TRUE, quiet=TRUE))
        
        if (ncol(cmap) == 9) {
            colnames(cmap) <- c(  "CMapId", "ContigLength", "NumSites",
                                  "SiteID", "LabelChannel", "Position", 
                                  "StdDev", "Coverage", "Occurrence"
            )
        } else if (ncol(cmap) == 11) {
            colnames(cmap) <- c(  "CMapId", "ContigLength", "NumSites",
                                  "SiteID", "LabelChannel", "Position", 
                                  "StdDev", "Coverage", "Occurrence",
                                  "GmeanSNR", "lnSNRsd" #new columns
            )
        } else if (ncol(cmap) == 12) {
        colnames(cmap) <- c(  "CMapId", "ContigLength", "NumSites",
                              "SiteID", "LabelChannel", "Position", 
                              "StdDev", "Coverage", "Occurrence",
                              "GmeanSNR", "lnSNRsd", #new columns
                              "SNR" #new columns
            )
        }
    } ) # try
    return(cmap)
} #readCMap

readCMap <- function(filename) {
    return(load_map(filename))
}

readXMapOld <- function(filename) {
    xmap <- NULL
    try( {
        
      ##message(filename)
      colnames <- c( "XmapEntryID", "QryContigID", "RefContigID", 
                    "QryStartPos", "QryEndPos", "RefStartPos", "RefEndPos", 
                    "Orientation", "Confidence", "HitEnum" )
      xmap <- read.table(filename, comment.char="#", header=FALSE, stringsAsFactors=FALSE, col.names=colnames)
      xmap <- as.data.frame(scan(filename, what=xmap, comment.char="#", quiet=TRUE))
        
        #colnames(xmap) <- c( "XmapEntryID", "QryContigID", "RefcontigID", 
        #                     "QryStartPos", "QryEndPos", "RefStartPos", "RefEndPos", 
        #                     "Orientation", "Confidence", "HitEnum" )
    } ) #try
    return(xmap)
} #readXMap


readXMap <- function(filename) {
    return(load_map(filename))
}


load_map <- function(filename, sep="") {
  #header0 <- read.table(filename, header=FALSE, stringsAsFactors=FALSE, comment.char="", nrows=1000, fill=TRUE)
  header0 <- readLines(filename, n=100)
  if( all(grepl("#",header0)) ) #if all lines are comments, return NULL
    return(NULL)
  header1 <- header0[regexpr("^#h", header0)>=0]
  #print(header0)
  #print(header1)
  #header<-read.table(pipe(paste("grep -E '^#h'", filename)), header=FALSE, stringsAsFactors=FALSE, comment.char="")
  header<-read.table(textConnection(header1[1]), header=FALSE, stringsAsFactors=FALSE, comment.char="")
  #print(header); cat(ncol(header),"\n")

  #A <- read.table(filename, comment.char="#", header=FALSE, stringsAsFactors=FALSE, nrows=5000, sep=sep, col.names=header[2:length(header)]) #will die if number of header columns is wrong
  A <- read.table(filename, comment.char="#", header=FALSE, stringsAsFactors=FALSE, nrows=5000, sep=sep)
  #print(A[1:10,])
  stopifnot(ncol(header)-1 >= ncol(A)) #first entry in header is '#h' which is not a column; be sure to have enough columns
  #names(A) <- unlist(header[2:length(header)])
  names(A) <- header[2:(ncol(A)+1)] #both need to increase by 1
  #print(names(A))

  if(dim(A)[1]>=5000) {
    B <- as.data.frame(scan(filename, what=A, comment.char="#", sep=sep, quiet=TRUE))
    A <- B
  }
  return(A)
}


readSMapOLD <- function( filename ) {
    smap <- NULL;
    try ( {
        
        #message(filename)
        smap <- read.table(filename, comment.char="#", header=FALSE, stringsAsFactors=FALSE, nrows=100)
        
        if (ncol(smap) == 12) {
            colnames( smap ) <- c( "SmapEntryID", "QryContigID", "RefcontigID1", "RefcontigID2",
                                   "QryStartPos", "QryEndPos", "RefStartPos", "RefEndPos",
                                   "Confidence", #check for column for orientation
                                   "Type", "XmapID1", "XmapID2" 
            )
        } else if (ncol(smap) == 13) {
            colnames( smap ) <- c( "SmapEntryID", "QryContigID", "RefcontigID1", "RefcontigID2",
                                   "QryStartPos", "QryEndPos", "RefStartPos", "RefEndPos",
                                   "Orientation", "Confidence", #check for column for orientation
                                   "Type", "XmapID1", "XmapID2" 
            )
        }
    },
    silent = TRUE #added 03202014 to suppres read smap error there are no lines
    ) #try
    return(smap)
} #readSMap

readSMap <- function( filename ) {
    return(load_map(filename))
}
