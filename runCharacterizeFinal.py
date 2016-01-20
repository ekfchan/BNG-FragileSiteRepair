
description = """Standalone script for running the Characterize Module of the Pipeline,
ie, aligning genome maps (ie, bioNano contigs as .cmap) to sequence contigs or a reference (.cmap or .spots).
Output will be in a sub-dir of the dir which contains the query maps called 'alignref'.
Alternatively, load an xmap file and its corresponding _r.cmap and _q.cmap files and report summary statistics."""

import os, sys
import argparse


def getArgs() :    
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-t', dest='RefAligner', help='Path to RefAligner (required unless xmap is specified (-x))')
    parser.add_argument('-r', dest='referenceMap', help='Path to reference maps (.cmap or .spots), 1 file only (required unless xmap specified (-x) and _r.cmap is present in same dir as xmap)', default="")
    parser.add_argument('-q', dest='queryMap', help='Path to query maps (.cmap), 1 file only (required--if xmap specified (-x), this should be input (-i argument) for that command)', default="")
    parser.add_argument('-x', dest='xmap', help='Path to .xmap, 1 file only (optional, if specified, no alignment is done, if not specified, -t, -r, and -q must be specified)')
    parser.add_argument('-p', dest='pipelineDir', help='Pipeline dir (optional, defaults to current directory)')
    parser.add_argument('-a', dest='optArguments', help='Path to optArguments.xml (optional, default in Pipeline dir if found, otherwise required)')
    parser.add_argument('-n', dest='numThreads', help='Number of threads (cores) to use (optional, default 4)', default=4, type=int)
    parser.add_argument('-v', dest='pvalue', help='Pvalue (-T) used for alignment', default="1e-12")
    result = parser.parse_args()

    #check all Pipeline dependencies
    if result.pipelineDir :
        cwd = result.pipelineDir
    else :
        cwd = os.getcwd()

    #this is the only one imported here and in runCharacterize
    if not os.path.isfile(os.path.join(cwd,"utilities.py")):
        print "utilities.py missing in dir", cwd, "check -p argument, or run this script in Pipeline dir"
        sys.exit(1)
    import utilities as util

    #xmap -- optional
    runaligns = True #default is to run the alignment
    xmappath = None
    if result.xmap :
        xmappath = result.xmap
        if not util.checkFile(xmappath, ".xmap") :
            print "Xmap path is supplied ("+xmappath+") but not found or doesn't end in .xmap."
            sys.exit(1)
        runaligns = False

    #RefAligner -- only required if xmap not specified
    rabin = result.RefAligner
    if not xmappath and not util.checkExecutable(rabin):
        print "RefAligner not found at", rabin, "\nPlease supply RefAligner full path as -t arg."
        sys.exit(1)

    #reference maps -- only required if xmap not specified
    refcmap = result.referenceMap
    if runaligns and not util.checkFile(refcmap, ".cmap") and not util.checkFile(refcmap, ".spots") :
        print "Reference map file ("+refcmap+") not found or does not end in .cmap or .spots. Check -r argument."
        sys.exit(1)

    #query maps -- only required if xmap not specified
    qrypath = result.queryMap
    #if runaligns and not util.checkFile(qrypath, ".cmap") :
    if not util.checkFile(qrypath, ".cmap") : #always required
        print "Query map file ("+qrypath+") not found or does not end in .cmap or .spots. Check -q argument."
        sys.exit(1)
    #if runaligns :
    contigdir  = os.path.split(qrypath)[0] #dir of query maps
    contigbase = os.path.split(qrypath)[1] #filename
    #else :
    #    contigdir  = os.path.split(xmappath)[0]
    #    contigbase = os.path.split(xmappath)[1] #filename
    contigbase = contigbase[:contigbase.find(".")] #remove suffix

    #optargs file
    optargs = None
    if result.optArguments : #supplied on command line
        optargs = result.optArguments
        if not util.checkFile(optargs, ".xml") :
            print "optArguments path is supplied ("+optargs+") but not found or doesn't end in .xml, check -a argument."
            sys.exit(1)
    elif runaligns : #load from Pipeline dir if running alignments
        optafile = "optArguments_human.xml"
        optargs = os.path.join(cwd, optafile)
        if not util.checkFile(optargs):
            print "%s missing in Pipeline directory (%s). Try supplying path explicitly using -a." % (optafile, cwd)
            sys.exit(1)

    #nthreads
    nthreads = result.numThreads
    if nthreads <= 0 :
        print "Number of threads value invalid (must be >= 0): "+nthreads
        sys.exit(1)

    #pvalue
    if result.pvalue : #supplied on command line
        pvalue = result.pvalue
    else :
        pvalue = "1e-12"    

    #yes, this is messy...but I don't want another class (besides varsPipeline) and they just go to runCharacterize
    return cwd, rabin, refcmap, contigdir, contigbase, runaligns, xmappath, optargs, nthreads, pvalue



def runCharacterize(cwd, rabin, refcmap, contigdir, contigbase, runaligns, xmappath, optargs, nthreads, pvalue):
    '''Load Pipeline files from first arg; configure CharacterizeModule; run alignments if runaligns;
    report on those alignments or the xmap provided as xmappath.
    '''

    printargs = True

    if not os.path.isfile(os.path.join(cwd,"utilities.py")):
        print "utilities.py missing in dir", cwd, "check -p argument, or run this script in Pipeline dir"
        sys.exit(1)
    import utilities as util

    if not util.checkFile(os.path.join(cwd,"Pipeline.py")):
        print "Pipeline.py missing in dir", cwd, "check -p argument, or run this script in Pipeline dir"
        sys.exit(1)
    import Pipeline

    if not util.checkFile(os.path.join(cwd,"CharacterizeModule.py")):
        print "CharacterizeModule.py missing in dir", cwd, "check -p argument, or run this script in Pipeline dir"
        sys.exit(1)
    import CharacterizeModule as cm

    if not util.checkFile(os.path.join(cwd,"MapClassesRev.py")):
        print "MapClassesRev.py missing in dir", cwd, "check -p argument, or run this script in Pipeline dir"
        sys.exit(1)
    import MapClassesRev

    #use Pipeline objects

    varsP = Pipeline.varsPipeline()

    varsP.optArgumentsFileIn   = optargs
    varsP.RefAlignerBin        = rabin
    varsP.latestMergedCmap     = os.path.join(contigdir, contigbase+".cmap") #file suffix required to be .cmap
    varsP.contigFolder         = contigdir
    varsP.nThreads             = nthreads #necessary otherwise job won't start
    varsP.ref                  = refcmap
    varsP.stdoutlog            = True #enable -stdout -stderr args to RefAligner
    varsP.curCharacterizeCmaps = [varsP.latestMergedCmap]

    if runaligns :
        varsP.contigAlignTarget = contigdir+"/alignref_final" #this is output dir
        varsP.runSV = False
        varsP.groupContigs = False
        varsP.stageComplete = contigbase
        varsP.outputContigFolder = contigdir
        varsP.memoryLogpath  = os.path.join(contigdir, "memory_log.txt")
        varsP.stdoutlog = True
        varsP.pipeReportFile = os.path.join(contigdir, "pipeReport.txt")
        varsP.parseArguments() #parses optArgumentsFile
        varsP.replaceParam("characterizeFinal", "-T", pvalue)
        if printargs :
            print "\nRunning Characterization with arguments:\n" + " ".join(varsP.argsListed('characterizeFinal')) + '\n'
        if hasattr(util, "InitStatus") : #if old version, skip
            util.InitStatus(os.path.join(contigdir, "status.xml")) #needed otherwise call to status_log fails
        charmod = cm.Characterize(varsP, 1) #create Characterize object from CharacterizeModule -- this also calls generateJobList
        xmappath = charmod.xmapTarget #set in Characterize.generateJobList
        charmod.runJobs()
    else :
        #varsP.contigAlignTarget = contigdir #this is dir in which _q and _r cmaps must be located -- contigdir is from cmap; this should be from xmap
        varsP.contigAlignTarget = os.path.split(xmappath)[0]
        print "Loading alignments from\n" + xmappath + "\n"

#no longer using this in Pipeline
    #print MapClassesRev.TopLevelCharacterization(varsP, [os.path.join(varsP.contigAlignTarget, contigbase)])

    print cm.characterizeContigs(varsP, xmappath) #this is redundant with above

#end runCharacterize


if __name__ == "__main__" :
    args = getArgs()
    runCharacterize(*args)
