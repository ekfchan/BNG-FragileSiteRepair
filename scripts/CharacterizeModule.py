
import os

import mapClasses
#import MapClassesRev
#import RefinementModule as rm
import Multithreading as mthread
import utilities as util

"""
@package CharacterizeModule 
Get general stats and mapping stats (if reference) for contigs

"""


util.setVersion("$Id: CharacterizeModule.py 4051 2015-08-17 19:15:56Z wandrews $")


class dummyCharacterize() :
    """For getting noise parameters in case of bypassing characterize."""
    def __init__(self, varsP) :
        self.curCharacterizeFileRoots = []
        self.varsP = varsP #bc Characterize uses this for totAssemblyLenMb
        #this is problematic for bypass (because mergeIntoSingleCmap isn't called)--don't need it
        #if not len(varsP.curCharacterizeCmaps) : #need this, set in mergeIntoSingleCmap
        #    return
        #ccc = varsP.curCharacterizeCmaps[0]
        #outFileName = os.path.split(ccc)[1].replace(".cmap", "")
        #outfile = os.path.join(varsP.contigAlignTarget,outFileName) #WRONG bc contigAlignTarget is wrong...try this

        outdir = os.path.join(varsP.outputContigFolder, self.varsP.characterizeDirName) #'alignref'
        if not util.checkDir(outdir, makeIfNotExist=False) : #if this doesn't exist, we can't get what we need
            return
        outfile = None
        for qfile in os.listdir(outdir) :
            if qfile.endswith(".err") : #just take first .err file
                outfile = qfile
                break
        if not outfile : #if no .err files found, give up
            return
        outfile = os.path.join(outdir, outfile.replace(".err",""))
        self.curCharacterizeFileRoots.append(outfile)
        #also want to get varsP.totAssemblyLenMb
        self.varsP.totAssemblyLenMb = mapClasses.multiCmap(varsP.latestMergedCmap, lengthonly=True).totalLength / 1e6
#end class dummyCharacterize


class Characterize(mthread.jobWrapper):
    """Align contigs to reference to characterize quality of assembly.
    If no reference, report basic contig stats only.
    
    """
    def __init__(self, varsP, argset=-1):
        '''argset is toggle between CharacterizeDefault and CharacterizeFinal argumets:
        -1 is default, 1 is final
        '''
        self.varsP = varsP
        self.argStr = ("Final" if argset == 1 else "Default") #!=1 and !=-1 is error in generateJobList, but not here
        self.stageName = 'Characterize'+self.argStr+' '+self.varsP.stageComplete
        util.LogStatus("progress", "stage_start", self.stageName)
        mthread.jobWrapper.__init__(self, varsP, self.stageName,clusterArgs=varsP.getClusterArgs('characterizeDefault'))
        self.xmapTarget = None
        self.curCharacterizeFileRoots = []
        outdir = self.varsP.characterizeDirName # = 'alignref'
        if argset == 1 : #this is final
            outdir += '_final'
        varsP.contigAlignTarget = os.path.join(varsP.outputContigFolder, outdir)
        if not(os.path.exists(varsP.contigAlignTarget)):
            os.mkdir(varsP.contigAlignTarget)
        self.generateJobList(argset)
        
    def generateJobList(self,argset=-1):
        if not self.varsP.ref : #no jobs if no ref
            return
        jobargs = [self.varsP.RefAlignerBin, '-ref', self.varsP.ref]
        if argset == -1 and self.varsP.argData.has_key('characterizeDefault') : # don't use nominal default
            opta = self.varsP.argsListed('characterizeDefault')
        elif argset == 1 and self.varsP.argData.has_key('characterizeFinal') : #extend (on default) -- make this default
            opta = self.varsP.argsListed('characterizeFinal')
        else : #this is an error
            self.varsP.updatePipeReport("ERROR in CharacterizeModule.generateJobList: invalid argset %s\n" % str(argset))
            return
        
        for i, curCharacterizeCmap in enumerate(self.varsP.curCharacterizeCmaps):
            if self.varsP.numCharacterizeJobs == 1:
                jobName = 'Char'+self.argStr+'_%s' % self.varsP.stageComplete
            else:
                jobName = 'Char'+self.argStr+'_%s_%d' % (self.varsP.stageComplete, i+1)
            outFileName = os.path.split(curCharacterizeCmap)[-1].replace(".cmap", "")
            outfile = os.path.join(self.varsP.contigAlignTarget,outFileName)
            self.curCharacterizeFileRoots.append(outfile)
            expectedResultFile = outfile+".xmap"
            self.xmapTarget = expectedResultFile
            currentArgs = jobargs + ["-i", curCharacterizeCmap, "-o", outfile]
            stdoutf = None
            if self.varsP.stdoutlog :
                currentArgs.extend( ['-stdout', '-stderr'] )
                stdoutf = outfile+".stdout"
            currentArgs += ['-maxthreads', str(self.varsP.nThreads)]
            currentArgs += ['-output-veto-filter', '_intervals.txt$']
            currentArgs += opta
            s1Job = mthread.singleJob(currentArgs, jobName, expectedResultFile, jobName.replace(' ',''),maxThreads=self.varsP.nThreads,clusterLogDir=self.varsP.clusterLogDir, expectedStdoutFile=stdoutf)
            self.addJob(s1Job)
            if i==0:
                self.logArguments()
        
    def runJobs(self):
        if not self.varsP.ref :
            return
        self.multiThreadRunJobs(self.varsP.nThreads, sleepTime = 0.2)
    
    def checkResults(self):
        #old heading says complete here and then summary after contig list; new says summary here
        outstr = 'Stage Summary: %s\n' % self.stageName
        if not self.varsP.ref : #still want contig stats
            infoReport = "Skipping Characterize because no reference (-r)\n"
            self.varsP.updatePipeReport(infoReport, printalso=False) #put this in pipereport just as an fyi
            infoReport += outstr
            #infoReport += 'Stage Complete: %s\n' % self.groupName #set in jobWrapper constructor
            #infoReport += MapClassesRev.ContigCharacterizationNoRef(self.varsP,self.groupName)
            infoReport += characterizeContigs(self.varsP)
            self.varsP.updateInfoReport(infoReport + '\n')
            return
        self.doAllPipeReport()
        #infoReport = 'Stage Complete: %s\n' % self.groupName #set in jobWrapper constructor
        #infoReport += MapClassesRev.TopLevelCharacterization(self.varsP,self.curCharacterizeFileRoots,self.groupName)
        #infoReport += 'OLD characterize\n' #debug
        infoReport = characterizeContigs(self.varsP, self.xmapTarget)
        self.varsP.updateInfoReport(outstr + infoReport + '\n')
        util.LogStatus("progress", "stage_complete", self.stageName)

#end class Characterize


class referenceProcess(mthread.jobWrapper) :
    """Pre-process the reference for SV jobs, applying -mres.
    """

    def __init__(self, varsP) :
        jobName = "reference_process"
        opta_section = "referenceSvdetect"
        default_mres = "2.9"
        mres = "-mres"
        self.varsP = varsP
        usedefault = False
        if self.varsP.argData.has_key(opta_section) : #check if in optargs
            opta = self.varsP.argsListed(opta_section)
            if not mres in opta : #must have mres
                self.varsP.updatePipeReport("Warning in referenceProcess: "+mres+" missing in optArguments section "+opta_section+"\n")
                usedefault = True
        else :
            self.varsP.updatePipeReport("Warning in referenceProcess: optArguments section "+opta_section+" missing\n")
            usedefault = True
        if usedefault :
            opta = [mres, default_mres]

        mresstr = opta[opta.index(mres)+1] #get string for mres value for output name
        mresstr = mresstr.replace(".","")

        if not util.checkDir(self.varsP.refFolder) :
            self.varsP.updatePipeReport( "ERROR in referenceProcess: could not make output dir %s\n" % self.varsP.refFolder )
            return None
        refpref = os.path.basename(self.varsP.ref[:self.varsP.ref.rfind(".")]) + "_res" + mresstr
        outarg = os.path.join(self.varsP.refFolder, refpref) #refFolder is new output folder for this job
        expectedResultFile = outarg+".cmap" #if ref is spots, is this spots?
        args = [self.varsP.RefAlignerBin, '-o', outarg, '-i', self.varsP.ref, '-f', '-merge'] + opta
        stdoutf = None
        if self.varsP.stdoutlog :
            args.extend( ['-stdout', '-stderr'] )
            stdoutf = outarg+".stdout"
        args += ['-maxthreads', str(self.varsP.nThreads)]

        super(referenceProcess, self).__init__(self.varsP, jobName, clusterArgs=self.varsP.getClusterArgs("assembly"))

        job = mthread.singleJob(args, jobName, expectedResultFile, jobName, maxThreads=self.varsP.nThreads, clusterLogDir=self.varsP.clusterLogDir, expectedStdoutFile=stdoutf)
        self.addJob(job)

        util.LogStatus("progress", "stage_start", jobName)
        self.varsP.runJobs(self, "referenceProcess")
        self.doAllPipeReport()
        if not self.allResultsFound() : #this is an error, but we'll continue processing without SV detect
            err = "ERROR in referenceProcess: job failed, disabling SV detect"
            self.varsP.updatePipeReport( err+"\n" )
            util.LogError("error", err)
            #self.varsP.runSV = False #no need since this class is used in SVModule
        else :
            self.varsP.refDeresed = expectedResultFile #store good result for SV detect
            self.varsP.updatePipeReport( "referenceProcess: using reference %s for svdetect\n" % self.varsP.refDeresed )
        util.LogStatus("progress", "stage_complete", jobName)            
    #end __init__
            
    def runJobs(self):
        self.multiThreadRunJobs(self.varsP.nThreads, sleepTime = 0.2)

#end class referenceProcess


#deprecated
#after calling alignContigRef with strictness > 0, call this to make dict of contig id : aligned len
#def makeStrictAlignDict(contigpaths, aligndir):
#    """Depricated: not in use 
#    """
#    maxcontigs = -1
#    #nmultxmap = 0 #; nstrict = 0 #test
#    retdict = {}
#    for idx, cpath in enumerate(contigpaths) :
#        xmappath = aligndir+os.path.split(cpath)[-1].replace(".cmap", ".xmap") #need cmap file name
#        #make sure xmappath exists; if not, bail silently
#        if not os.path.exists( xmappath ) :
#            continue
#        xmapobj = mapClasses.xmap(xmappath)
#        if len(xmapobj.xmapLookup) == 0 :
#            continue
#        #lenq = getMappedQryLen(xmapobj) #currently not used
#        lenr = getMappedRefLen(xmapobj)
#
#        #contigQry is equivalent to the contigID--this is the key
#        #print "adding", xmapobj.xmapLookup[1].contigQry, ":", lenq #debug
#        #value is just len of aligned portion--ref for compare with totalignlen
#        retdict[xmapobj.xmapLookup[1].contigQry] = lenr 
#
#        #if len(xmapobj.xmapLookup) > 1  : #test--turn off
#            #print "strict multi-xmap:", xmapobj.xmapLookup[1].contigQry #test--turn off
#            #nmultxmap += 1 #else 0)
#
#        if maxcontigs != -1 and idx >= maxcontigs :
#            break
#
#    return retdict


#deprecated
#Put a dictionary for strings for cmap entries for the merged reference
#def refIndexToChrStr(index) :
#    refchrdict = {} 
#    if refchrdict.has_key( index ) :
#        return refchrdict[index]
#    else :
#        return "NA"


#get the mapped length and the alignment parameters from the align to ref call
#mapped, params = getMappedStats(aligndir, cpath)
#return the alignParams object so that it can be used in characterizeContigs
def getMappedErrStats(aligndir, cpath):
    """Only used in characterizeContigs if listcontigs.
    """
    errfile = aligndir+os.path.split(cpath)[-1].replace(".cmap", ".err") #path checking exists in alignparams constructor
    return alignParams(errfile)


    
#if not listcontigs, do not print list of all contigs, ie, summary only
def characterizeContigs(varsP, xmappath=None, listcontigs=False) :
    """Log simple contigs stats, and optionally align stats from xmappath.
    """
    #print header of table
    unitscale = 1e-6
    dorefalign = bool(xmappath) #i'm never actually calling refaligner here--this is just using xmappath
    #dorefidchr = False
    #dorefcid = False
    printrange = False
    #printsegrms = False
    #dochrstr = False
    iscluster = True
    haveref = bool(varsP.ref)

    #refcmap = mapClasses.multiCmap() #not used
    aligndir = varsP.contigAlignTarget

    try :
        #refcmap = mapClasses.multiCmap(varsP.ref)
        #reflen = refcmap.totalLength #note: total length of _all_ contigs
        reflen = mapClasses.multiCmap(varsP.ref, lengthonly=True).totalLength
        #in summary table, this is a denominator--make sure it's non-zero, don't bail (still get contig summary)
        if reflen <= 0 :
            #print "Warning in CharacterizeModule.characterizeContigs: bad reflen", reflen, "defaulting to 1" #not necessary
            reflen = 1.
    except:
        reflen = 1.

    outstr = "" #Contig Characterization:\n"

    if listcontigs and dorefalign :
        outstr += "cID  len"
        outstr += "  Cov"
        #if dorefidchr or dorefcid :
        #    outstr += "  rID" #ref index for either of these
        #if dorefidchr :
        #    outstr += "  rpos"
        outstr += "  alignlen  alignlen/len"
        if printrange :
            outstr += "  Qry1  Qry2  Ref1  Ref2"
        #outstr += ("  segRMS" if printsegrms else "")
        outstr += "  Conf  Conf/lenkb"
        outstr += "  FP  FN  sf  sd  bpp" #"  res" #--ditch res (not bpp)
        #if dochrstr :
        #    outstr += "  Chr"
        outstr += "\n"

    totcontiglen = 0; totalignlen = 0; nmapcontigs = 0; defalignlen = 0; totalignqlen = 0
    contiglens = [] #lens of all contigs in bases
    uniqueseg = {} #the argument to util.uniqueRange--stores all the map lengths which go into totalignlen--now each of these is a value, and the keys are the reference contig id or chromosome depending on dorefidchr
    avgfplist = []; avgfnlist = [] #average FP/FN rates
    #if dorefidchr :
    #    chrsum = refcmap.makeChrSummaryDict() #see mapClasses.multiCmap
    for citr, cpath in enumerate([varsP.latestMergedCmap]) : #always use contigpaths
        mapi = mapClasses.multiCmap(cpath) 
        totcontiglen += mapi.totalLength
        contiglens += mapi.getAllMapLengths() #getAllMapLengths is list of all map lengths

        #store a list of the contig ids in this multiCmap, then remove them if they're in the xmap
        # if they're not, print at the end
        mapids = mapi.getAllMapIds() #this is once per cmap, not once per characterizeModule call--becuase it's keys, it's already a copy, no need to copy explicitly
        ncontigs = len(mapids) #this is ncontigs in this file, ie, in mapi (see below)

        xmapobj = mapClasses.xmap() #empty map to fix xmapobj scope
        if dorefalign : #get xmap object
            #xmappath = aligndir+os.path.split(cpath)[-1].replace(".cmap", ".xmap") #need cmap file name
            #xmappath = self.xmapTarget
            #if xmappath exists isn't a file, nothing will be loaded
            if os.path.isfile( xmappath ) : #was if not isfile : continue
                xmapobj = mapClasses.xmap(xmappath)

        for xitr, xmapentry in enumerate(xmapobj.xmapLookup.values()) :
            #print the contig id from the xmap
            #this is sorted by ref position; could sort the list this loop is over by the contigQry data member,
            # _but_, I think I like reference-oriented better because you can see gap spanning

            #get map length from multicmap.getMapLength--returns 0 for any exception
            contiglen = mapi.getMapLength(xmapentry.contigQry)
            if contiglen <= 0 : #this strikes me as clumsy...but I don't want to return non-zero from multiCmap.getMapLength
                contiglen = 1.
            contigcov = mapi.getMapAvgCoverage(xmapentry.contigQry)

            if listcontigs :
                outstr += "%5i" % xmapentry.contigQry
                outstr += "  %9.1f  %2i" % (contiglen, contigcov)

            #don't print lenr for each contig--just total them
            lenr = xmapentry.getMappedRefLen()
            lenq = xmapentry.getMappedQryLen()
            refid = xmapentry.contigRef #int
            #if dorefidchr : #this is the encoding of ref contig id to chromosome and start position
            #    chrpos = mapClasses.intToChrPos(refid, verbose=False) #this returns a tuple (shouldn't fail, but verbose false)
            #    refidstr = "  %2s  %6i" % chrpos
            #    chrs = chrpos[0] #just for readability
            #    if chrsum.has_key(chrs) : #the method that fills chrsum (makeChrSummaryDict) also uses intToChrPos
            #        chrsum[chrs][0] += lenr #values are list, first ele is total aligned length, second is ref length (fixed)
            #elif dorefcid :
            #    refidstr = "  %3i" % refid #refid is int
            #else : #nothing for neither, but still need empty string
            refidstr = ""

            conf = xmapentry.Confidence #confidence from the xmap, and ratio of it to length in kb
            if listcontigs :
                alignpars = getMappedErrStats(aligndir, cpath) #an empty err file is produced for case of no align
                avgfplist.append( alignpars.fp ) 
                avgfnlist.append( alignpars.fn ) 
                outstr += "%s  %9.1f  %.3f" % (refidstr, lenq, lenq/contiglen) #length for refidstr set above
                if printrange :
                    outstr += "  %5.0f  %5.0f  %5.0f  %5.0f" % (xmapentry.QryStart/pn, xmapentry.QryStop/pn, xmapentry.RefStart/pn, xmapentry.RefStop/pn)
                #outstr += ("  %5.0f" % 0 if printsegrms else "") #don't print anything
                outstr += "  %3.0f  %5.3f" % (conf, conf*1000./lenq) #1000 is for kb
                outstr += "  " + alignpars.getParamString()

            totalignlen  += lenr
            totalignqlen += lenq

            #uniqueseg is now a dict to take into account which chromosome the query contig is on
            #note need refid bc need to separate different contigs on the _same_ chromosome
            if not uniqueseg.has_key(refid) : #if first contig on chromosome, need to init new list
                uniqueseg[refid] = []
            uniqueseg[refid].append( [xmapentry.RefStart, xmapentry.RefStop] )

            #process mapids--remove contig id (contigQry) from mapids if they're in the xmap so non-aligning contigs can be printed
            if xmapentry.contigQry in mapids :
                mapids.remove(xmapentry.contigQry)

            #note: the feature of multiple alignments (strict vs default) is no longer implemented
            defalignlen  += lenr #currently, just default and strict

            #if listcontigs and dochrstr :
            #    outstr += "  " + refIndexToChrStr( xmapentry.contigRef )
            #    outstr += "\n"
            
        #end loop on xmap entries

        #now that all xmap entries are processed, all contigs with an alignment are removed from mapids,
        # so we can get n contigs align using this and ncontigs
        nmapcontigs += ncontigs - len(mapids) #sum multiple cmaps

        #and print the data for the contigs which don't align--just id, length, and coverage
        #these lines are kind of redundant, but I guess that's ok
        if listcontigs :
            for ids in mapids :
                outstr += "%5i" % ids
                #get map length from multicmap.getMapLength--returns 0 for any exception
                contiglen = mapi.getMapLength(ids) #it's ok if it's 0 bc it's never a denominator here
                contigcov = mapi.getMapAvgCoverage(ids)
                outstr += "  %9.1f  %2i\n" % (contiglen, contigcov)

    #end loop on contigs

    varsP.totAssemblyLenMb = totcontiglen*unitscale
    ncontigs = len(contiglens) #contigpaths is just files--contiglens is all contigs
    avgcontiglen = (float(totcontiglen)/ncontigs if ncontigs > 0 else 0)

    #print averages
    if listcontigs and not iscluster : #only do avg if not merged, otherwise just one noise parameter
        avgfp    = sum(avgfplist)/len(avgfplist)
        avgfn    = sum(avgfnlist)/len(avgfnlist)
        outstr += "AVG    %9.1f           %9.1f                     %5.3f  %5.3f\n" % (avgcontiglen, totalignqlen/nmapcontigs, avgfp, avgfn)

    if unitscale > 1e-6 : #if not megabases
        fstr = "%9.0f"
    else : #megabases
        fstr = "%8.3f" 

    outstr += "N Genome Maps: %i\n" % ncontigs
    outstr += ("Total Genome Map Len (Mb): "+fstr+"\n") % (totcontiglen*unitscale)
    outstr += ("Avg. Genome Map Len  (Mb): "+fstr+"\n") % (avgcontiglen*unitscale)
    outstr += ("Median Genome Map Len(Mb): "+fstr+"\n") % (util.getMedian(contiglens)*unitscale)
    outstr += ("Genome Map n50       (Mb): "+fstr+"\n") % (util.getn50(contiglens)*unitscale)

    if haveref :
        outstr += ("Total Ref Len   (Mb): "+fstr+"\n") % (reflen*unitscale)
        outstr += ("Total Genome Map Len / Ref Len : "+fstr+"\n") % (totcontiglen/reflen)
    if dorefalign :
        #print the chromosome summary before the strict/default/total align stats
        #if dorefidchr :
        #    outstr += "Chromosome Summary:\n"
        #    outstr += "Chromosome  align len  ref len  (ratio):\n"
        #    for chrs, align in chrsum.iteritems() :
        #        outstr += "%3s  %9.0f  %9.0f  (%5.3f)\n" % (chrs, align[0], align[1], align[0]/align[1])

        ratio = (float(nmapcontigs)/ncontigs if ncontigs > 0 else 0)
        outstr += ("N Genome Maps total align      : %i (%.2f)\n") % (nmapcontigs, ratio)
        outstr += ("Total Aligned Len (Mb)            : "+fstr+"\n") % (totalignlen*unitscale)
        outstr += ("Total Aligned Len / Ref Len       : "+fstr+"\n") % (totalignlen/reflen)
        uniquelen = 0
        for segs in uniqueseg.values() : # need to sum on dict entries
            util.uniqueRange(segs) #this modifies list in place
            uniquelen += util.totalLengthFromRanges( segs )
        outstr += ("Total Unique Aligned Len (Mb)     : "+fstr+"\n") % (uniquelen*unitscale)
        outstr += ("Total Unique Aligned Len / Ref Len: "+fstr+"\n") % (uniquelen/reflen)

    return outstr

#end characterizecontigs
