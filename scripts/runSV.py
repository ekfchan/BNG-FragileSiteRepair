
description = """Standalone script for running the SV Module of the Pipeline,
ie, aligning genome maps (ie, bioNano contigs as .cmap) to sequence contigs or a reference (.cmap)
for the purpose of structural variation (SV) detection."""


import os, sys
import argparse


def getArgs() :    
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-t', dest='RefAligner', help='Path to RefAligner or dir containing it (required)', type=str) 
    parser.add_argument('-r', dest='referenceMap', help='Path to reference maps (.cmap), 1 file only (required)', type=str)
    parser.add_argument('-q', dest='queryDir', help='Path to dir containing query maps (.cmaps) (required)', type=str)
    #parser.add_argument('-x', dest='xmap', help='Path to .xmap, 1 file only (optional, if specified, no alignment is done, if not specified, -t, -r, and -q must be specified)') #not supported
    parser.add_argument('-o', dest='outputDir', help='output dir (optional, defaults to input map dir with suffix "_sv")', default="", type=str)
    parser.add_argument('-p', dest='pipelineDir', help='Pipeline dir (optional, defaults to script dir, or current directory)', default="", type=str)
    parser.add_argument('-a', dest='optArguments', help='Path to optArguments.xml (optional, default optArguments_human.xml in Pipeline dir if found, otherwise required)', default="", type=str)
    parser.add_argument('-T', dest='numThreads', help='Total number of threads (cores) to use (optional, default 4)', default=4, type=int)
    parser.add_argument('-j', dest='maxthreads', help='Threads per Job, -maxthreads (non-cluster only;optional, default 4)', default=4, type=int)
    parser.add_argument('-b', dest='bedFile', help='.bed file with gaps in reference for flagging SVs which overlap N-base gaps (optional)', default="", type=str)
    parser.add_argument('-e', dest='errFile', help='.err file to use for noise parameters--will supersede noise parameters in the optArgument supplied (but that file must still be supplied for non-noise parameters)', default="", type=str)
    parser.add_argument('-E', dest='errbinFile', help='.errbin file to use for noise parameters--will supersede noise parameters in the optArgument supplied (but that file must still be supplied for non-noise parameters)', default="", type=str)
    parser.add_argument('-C', help='Run on cluster, read XML file for submission arguments (optional--will not use cluster submission if absent)', dest='cxml', default=None)
    parser.add_argument('-s', help='Disable grouping of SV jobs (default grouped; optional)', dest='groupsv', action='store_false')
    result = parser.parse_args()

    #check all Pipeline dependencies
    if result.pipelineDir :
        cwd = result.pipelineDir
    else :
        cwd = os.path.split(os.path.realpath(__file__))[0] #this is path of this script
        if not os.path.isfile(os.path.join(cwd,"utilities.py")) : #if still not here, last try is actual cwd
            cwd = os.getcwd() #still check this below

    #this is the only one imported here and in runCharacterize
    if not os.path.isfile(os.path.join(cwd,"utilities.py")):
        print "utilities.py missing in dir", cwd, "check -p argument, or run this script in Pipeline dir"
        sys.exit(1)
    import utilities as util

    #xmap -- don't use this
    runaligns = True #default is to run the alignment
    xmappath = None
    #if result.xmap :
    #    xmappath = result.xmap
    #    if not util.checkFile(xmappath, ".xmap") :
    #        print "Xmap path is supplied ("+xmappath+") but not found or doesn't end in .xmap."
    #        sys.exit(1)
    #    runaligns = False

    #RefAligner -- check for either path to RefAligner, or dir containing it, depending on cluster args
    rabin = result.RefAligner
    #replicate Pipeline behavior: RefAligner is always required
    if os.path.isdir(rabin) :
        rabin = os.path.join(rabin, "RefAligner")
    if not util.checkExecutable(rabin):
        print "RefAligner not found or not executable at", rabin, "\nPlease supply RefAligner dir or full path as -t arg."
        sys.exit(1)

    #reference maps -- only required if xmap not specified
    refcmap = os.path.realpath(result.referenceMap)
    if runaligns and not util.checkFile(refcmap, ".cmap") : #and not util.checkFile(refcmap, ".spots") :
        print "Reference map file ("+refcmap+") not found or does not end in .cmap or .spots. Check -r argument."
        sys.exit(1)

    #query maps
    qrypath = os.path.realpath(result.queryDir)
    #if runaligns and not util.checkFile(qrypath, ".cmap") :
    #    print "Query map file ("+qrypath+") not found or does not end in .cmap or .spots. Check -q argument."
    #    sys.exit(1)
    if not util.checkDir(qrypath, checkWritable=False, makeIfNotExist=False) : #does NOT have to be writeable
        print "Query dir ("+qrypath+") not found or not a dir. Check -q argument."
        sys.exit(1)
    if runaligns :
        contigdir  = qrypath #os.path.split(qrypath)[0] #dir of query maps
        contigbase = os.path.split(qrypath)[1] #filename
    else :
        contigdir  = os.path.split(xmappath)[0]
        contigbase = os.path.split(xmappath)[1] #filename
    #contigbase = contigbase[:contigbase.find(".")] #remove suffix

    #optargs file
    optargs = None
    if result.optArguments : #supplied on command line
        optargs = result.optArguments
        if not util.checkFile(optargs, ".xml") :
            print "optArguments path is supplied ("+optargs+") but not found or doesn't end in .xml, check -a argument."
            sys.exit(1)
    elif runaligns : #load from Pipeline dir if running alignments
        optargs = os.path.join(cwd,"optArguments_human.xml")
        if not util.checkFile(optargs):
            print "optArguments.xml missing in Pipeline directory ("+cwd+"). Try supplying path explicitly using -a."
            sys.exit(1)

    #cluster args
    clustargs = None
    if result.cxml :
        clustargs = os.path.realpath(result.cxml)
        if not util.checkFile(clustargs, ".xml") :
            print "clusterArguments path is supplied ("+clustargs+") but not found or doesn't end in .xml, check -C argument."
            sys.exit(1)

    #nthreads
    nthreads = result.numThreads
    if nthreads <= 0 :
        print "Number of threads value invalid (must be > 0): %i" % nthreads
        sys.exit(1)

    #maxthreads
    maxthreads = result.maxthreads
    if maxthreads <= 0 :
        print "Max threads value invalid (must be > 0): %i" % maxthreads
        sys.exit(1)

    #bed file
    bedfile = result.bedFile #must make local for return statement below
    if bedfile : #must check for empty string BEFORE you do realpath, or it returns cwd
        bedfile = os.path.realpath(result.bedFile)
        if not util.checkFile(bedfile, ".bed") :
            print "bed file supplied but not found or incorrect suffix:", bedfile
            sys.exit(1)

    #.errbin file
    errbinfile = result.errbinFile
    if errbinfile :
        errbinfile = os.path.realpath(result.errbinFile)
        if not util.checkFile(errbinfile, ".errbin") :
            print "errbin file supplied but not found or incorrect suffix:", errbinfile
            sys.exit(1)

    #.err file
    errfile = result.errFile
    if errfile and errbinfile :
        print "Warning: .err and .errbin arguments supplied; ignoring .err file"
        errfile = ""
    elif errfile :
        errfile = os.path.realpath(result.errFile)
        if not util.checkFile(errfile, ".err") :
            print "err file supplied but not found or incorrect suffix:", errfile
            sys.exit(1)

    outdir = os.path.realpath(result.outputDir)

    #yes, this is messy...but I don't want another class (besides varsPipeline) and they just go to runCharacterize
    return cwd, rabin, refcmap, contigdir, contigbase, runaligns, xmappath, optargs, nthreads, maxthreads, bedfile, errfile, outdir, errbinfile, clustargs, result.groupsv


def runSV(cwd, rabin, refcmap, contigdir, contigbase, runaligns, xmappath, optargs, nthreads, maxthreads, bedfile, errfile, outdir, errbinfile, clustargs, groupsv):
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

    if not util.checkFile(os.path.join(cwd,"SVModule.py")):
        print "SVModule.py missing in dir", cwd, "check -p argument, or run this script in Pipeline dir"
        sys.exit(1)
    import SVModule as svm

    if errfile and not util.checkFile(os.path.join(cwd,"SampleCharModule.py")):
        print "SampleCharModule.py missing in dir", cwd, "check -p argument, or run this script in Pipeline dir"
        sys.exit(1)
    elif errfile :
        import SampleCharModule as scm

    #use Pipeline objects

    varsP = Pipeline.varsPipeline()

    varsP.optArgumentsFileIn   = optargs
    varsP.RefAlignerBin        = rabin
    varsP.latestMergedCmap     = os.path.join(contigdir, contigbase+".cmap") #file suffix required to be .cmap
    varsP.contigFolder         = os.path.split(contigdir)[0]
    varsP.nThreads             = nthreads #necessary otherwise job won't start -- max threads per node
    varsP.maxthreads           = maxthreads #threads per job
    varsP.ref                  = refcmap
    varsP.stdoutlog            = True #enable -stdout -stderr args to RefAligner
    varsP.curCharacterizeCmaps = [varsP.latestMergedCmap]
    varsP.contigSubDirectories = True #needed for prepareContigIO
    varsP.doAlignMolvRef       = False #do not look for copy number
    varsP.groupSV              = groupsv #mimic Pipeline behavior: group or not 

    if runaligns :
        #varsP.contigAlignTarget = outdir
        varsP.runSV = False
        varsP.groupContigs = False
        varsP.stdoutlog    = True #use -stdout -stderr
        varsP.stageComplete = contigbase
        varsP.outputContigPrefix = getContigPrefix(util, contigdir) #if outdir is not supplied, this is used as dir prefix; also used as file pref for -o arg
        varsP.outputContigFolder = contigdir #cmaps are copied from here

        if not outdir :
            outdir = contigdir+"_sv" #this will be outdir of sv jobs
        if os.path.isdir(outdir) :
            if not util.checkDir(outdir) : #check writeable
                print "\nERROR: Output dir is not writeable:\n", outdir, "\n"                
                sys.exit(1)
            elif outdir == contigdir :
                print "\nERROR: Output dir cannot be same as input dir:\n", outdir, "\n"                
                sys.exit(1)                
            print "\nWARNING: Output dir already exists, results will be overwritten:\n", outdir, "\n"
        elif not util.checkDir(outdir) : #does not exist, make, if False, can't make or not writeable
            print "\nERROR: Output dir cannot be created or is not writeable:\n", outdir, "\n"
            sys.exit(1)

        if clustargs :
            os.putenv('SGE_ROOT', '/var/lib/gridengine') #do I want this???
            varsP.onCluster = True
            varsP.clusterLogDir = os.path.join(outdir, 'ClusterLogs')
            util.checkDir(varsP.clusterLogDir) #make it
            varsP.checkCluster()
            varsP.clusterArgumentsFileIn = clustargs #required for parseArguments
            varsP.parseArguments(readingClusterFile=True)
            if varsP.error :
                print varsP.message
                sys.exit(1)
            varsP.RefAlignerBin += "${BINARY_SUFFIX:=}" #copy from varsPipeline, handled by external script on phi host

        varsP.pipeReportFile = os.path.join(outdir, "sv_jobs_log.txt")
        varsP.infoReportFile = os.path.join(outdir, "sv_log.txt")
        varsP.memoryLogpath  = os.path.join(outdir, "memory_log.txt")
        if bedfile :
            varsP.bedFile = bedfile
        util.InitStatus( os.path.join(outdir, "status.xml") )
        varsP.parseArguments() #parses optArgumentsFile
        varsP.checkDependencies()
        varsP.RefAlignerBinOrig = rabin
        varsP.prerunLog() #general information in log -- needed for refaligner_version
        if printargs :
            print "\nRunning SV detection with arguments ("+os.path.split(optargs)[1]+"):\n" + " ".join(varsP.argsListed('svdetect')) + '\n'

        noisep = {}
        if errbinfile :
            noisep = {"readparameters": errbinfile}
            print "Using noise parameters from "+errbinfile+"\n"
        elif errfile :
            noisep = scm.readNoiseParameters(errfile.replace(".err",""))
            if noisep.has_key('readparameters') : #remove this because it's redundant, and it can cause problems with RefAligner compatibility
                del noisep['readparameters']
            if not noisep : #readNoiseParameters returns empty dict on failure
                print "ERROR reading noise parameters, check .err file:", errfile
                sys.exit(1)
            print "Using noise parameters from "+errfile+":\n" + " ".join(["-"+str(k)+" "+str(v) for k,v in noisep.iteritems()])+"\n"

        #make merged cmap to replace merged _q.cmap if not produced by RefAligner
        varsP.contigPathTxtFile = os.path.join(outdir, "contig_list.txt") #mergeIntoSingleCmap creates this file
        print "Creating merged cmap"
        varsP.mergeIntoSingleCmap(outdir)
        print "Merged cmap created:", varsP.latestMergedCmap, "\n"

        varsP.outputContigFolder = contigdir #cmaps are copied from here
        svmodule = svm.SVdetect(varsP, noisep, outdir, skipderes=True)
        #this got duplicated above
        #if hasattr(util, "InitStatus") : #if old version, skip -- do this after SVdetect.__init__ bc makes outdir
        #    util.InitStatus(os.path.join(outdir, "status.xml")) #needed otherwise call to status_log fails
        svmodule.runJobs()
        svmodule.checkResults()
        util.SummarizeErrors(varsP) 

    else :
        varsP.contigAlignTarget = contigdir #this is dir in which _q and _r cmaps must be located
        print "ERROR: feature not supported" #not implemented to not run jobs

    
#end runSV


#the SV module needs contig prefix for output dir and filenames--get any cmap and use its name
def getContigPrefix(util, contigdir) :
    cmaps = util.getListOfFilesFromDir(contigdir, suffix=".cmap")
    if not len(cmaps) : #this is an error--no maps to run on means we shouldn't proceed
        print "ERROR: no cmaps found in query dir, exiting:", contigdir
        sys.exit(1)
    #old behavior: always take first filename in cmaps
    #mapname = os.path.split(cmaps[0])[1]
    #if mapname.find("_contig") != -1 : #has this string
    #    return mapname[:mapname.find("_contig")] #name before this is prefix
    #new behavior: loop and skip any filenames with "_refined_", taking first without
    mapname = ""
    for cmap in cmaps :
        #print "cmap =", cmap #debug
        if cmap.rfind("_r.cmap") != -1 : #this is the _unrefined_ contig, so ignore these
            continue
        #this line is buggy in the case of, eg, hybrid scaffold which has _q.cmap not from final refine; instead, skip _all_ "_q.cmap"
        #if cmap.find("_refined_") != -1 and cmap.endswith("_q.cmap") : #must skip _refined_contigX_q.cmap files bc these are not the contigs
        if cmap.rfind("_q.cmap") != -1 :
            continue
        mapname = os.path.split(cmap)[1]
        #print "mapname =", mapname #debug
        if mapname.find("contig") != -1 : #has this string
            pref = mapname[:mapname.find("contig")] #name before this is prefix
            return pref.rstrip("_") #since _ was removed from find string, it's left here, but can't be for varsPipeline.findContigs
    if not mapname : #all contigs were _refined_
        print "ERROR: all cmaps found in query dir are _q.cmap files--invalid dir:", contigdir
        sys.exit(1)        
    return mapname.strip(".cmap") #just take this as prefix


if __name__ == "__main__" :
    args = getArgs()
    runSV(*args)
