import os
import shutil
from collections import defaultdict
from math import *
from string import Template

import Multithreading as mthread
import mapClasses as mc
import CharacterizeModule as cm

"""@package SVModule 
   Run structural variation (SV) detection on a given set of cmaps.
   Copy cmaps to new dir (for convenience of having input and output in same dir).
   Run SV detection using RefAligner to produce smaps.
   Report on results in informaticsReport. Merge smaps into bed files.

"""


import utilities as util
util.setVersion("$Id: SVModule.py 4302 2015-11-23 21:44:36Z wandrews $")



class SVdetect(mthread.jobWrapper):
    """SVdetect: main class for SVModule
    This should be called immediately after the stage for which you want SVs detected.
    """
    def __init__(self, varsP, noisep={}, fixedoutputdir="", skipderes=False): #targetStage, 
        """SVdetect.__init__: instantiating this class will also copy the cmaps to run detect on and generate its jobs.
        Noisep are noise parameters returned by SampleCharModule.readNoiseParameters.
        """
        #the stage to run on is obtained from self.varsP.outputContigPrefix
        #self.targetStage = targetStage #not necessary
        self.varsP = varsP
        self.moduleName = "svdetect" #used as arg to varsP.argsListed, varsP.getClusterArgs -- this shouldn't change
        self.varsP.stageName=self.moduleName
        self.verbose = 1
        self.output_reference_end = "contig1" #this is the one reference map which is kept--change all xmaps to refer to it
        self.group_jobs = varsP.groupSV #True #False #False should revert to old behavior (no grouping)
        self.combine_output = True #combine xmaps, smaps, and _q.cmaps
        self.ws = "\t" #Irysview requires tab separated
        self.noisep = noisep #process in generateJobList

        #configs for smapToBed
        self.smapTypeIdx = -1 #this is the column of the 'type' string in the smap--will change, so get it in getTypeIndex
        self.smapLinkIdx = -1 #if -1, there is no LinkID field; if >= 0, this is column number of it (should be last column)
        self.smapXmapIdx = -1 #analogous to LinkIdx
        #self.smapExcludeTypes = ["complex", "end", "inversion"] #the types to exclude keep changing--do not use
        self.smapIncludeTypes = ["insertion", "deletion"] #only used for bed files
        self.smapFilterTypes  = ["complex", "duplicate", "tiny", "end", "repeat", "overlap", "error"] #for filtered smap
        self.svcolors = {}
        self.svcolors["insertion"] = "0,128,0" #green
        self.svcolors["deletion" ] = "255,0,0" #red
        self.svcolors["default"  ] = "0,0,128" #blue, presumably--use a default to avoid crashes if names change/are old

        #initialize necessary data members
        self.mergedsmap = ""
        self.refoutpath = ""
        self.merged_qcmap = ""
        self.merged_rcmap = ""
        self.outxmap_list = [] #output of mergeSmaps

        if skipderes :
            self.ref = self.varsP.ref
        else :
            cm.referenceProcess(self.varsP) #make condensed reference, result stored here:
            if not self.varsP.refDeresed : #this means referenceProcess failed, so skip SV
                return None #it seems that despite this argument, the constructor still returns self

            self.ref = self.varsP.refDeresed

        #relies on outputContigPrefix: this is why you need to call this immediately after the stage for which you want SVs detected.
        newContigPrefix = self.varsP.outputContigPrefix + "_sv" 
        self.stageName = self.moduleName+"_"+newContigPrefix

        #setup for prepareContigIO:
        #specifically, the line: "self.inputContigPrefix = self.outputContigPrefix"
        #should just reset self.outputContigPrefix--store it here; see runJobs
        self.oldOutputContigPrefix = self.varsP.outputContigPrefix

        #these are the cmaps that need to be copied to the target dir, ie, varsP.outputContigFolder
        self.sourceContigFolder = self.varsP.outputContigFolder 

        #self.varsP.prepareContigIO(newContigPrefix) #for runSV.py, make this call conditional on fixedoutputdir
        if fixedoutputdir :
            self.varsP.outputContigFolder = fixedoutputdir
        else :
            self.varsP.prepareContigIO(newContigPrefix, self.varsP.stageName)

        super(SVdetect, self).__init__(self.varsP, self.stageName, clusterArgs=self.varsP.getClusterArgs(self.moduleName))

        #set dir for merged smap/bed files; will create at smap merge stage--see self.mergeSmaps()
        self.mergesmapdir = os.path.join(self.varsP.outputContigFolder, "merged_smaps")

        util.LogStatus("progress", "stage_start", self.stageName)
        self.copyTargetCmaps()
        self.generateJobList()
    #end __init__


    #using data members made in constructor (ie, i/o dirs),
    # find the cmaps which should be copied, and copy them
    # note: part of this is very similar to utilities.getListOfCmapsFromDir, but it needs more functionality
    #also, make list of contigs to run on to save reading from disk again in generateJobList
    def copyTargetCmaps(self) :
        """SVdetect.copyTargetCmaps: copy the cmaps which are compared to the reference for 
        SV detection to a new folder.
        """
        if self.verbose :
            outstr =  "SVdetect: Copying cmaps from: %s\n" % self.sourceContigFolder
            outstr += "  with prefix %s\n" % self.oldOutputContigPrefix
            outstr += "SVdetect: To output folder: %s" % self.varsP.outputContigFolder #this should be where the contigs will go
            self.varsP.updatePipeReport(outstr+"\n\n")

        self.targetCmapList = []
        for ifile in os.listdir(self.sourceContigFolder) :
            #define which files to skip:
            #any file which doesn't end in ".cmap"
            if not ifile.endswith(".cmap") :
                continue
            #any map for which the prefix is in all caps--these are the combined
            pref = ifile[:ifile.find(".")]
            if pref.upper() == pref :
                continue
            #any file (cmap) which ends with "_q" or "_r"
            if ifile.endswith("_q.cmap") or ifile.endswith("_r.cmap") :
                continue

            if ifile.find(self.oldOutputContigPrefix) == -1 : #for when last stage was merge, take only last merge
                continue

            #check for empty cmaps (zero size file), of so, skip
            if os.stat(os.path.join(self.sourceContigFolder,ifile)).st_size < 1 :
                continue

            #this should be a cmap to be processed--copy
            shutil.copy(os.path.join(self.sourceContigFolder,ifile), self.varsP.outputContigFolder)
            #store in list for further processing
            self.targetCmapList.append( os.path.join(self.varsP.outputContigFolder, ifile) )
    #end copyTargetCmaps


    #use addJob (inherited from Multithreading.jobWrapper) for each contig in self.targetCmapList
    def generateJobList(self) :
        """SVdetect.generateJobList: create RefAligner jobs for each contig to detect SVs on.
        """
        self.smapList = [] #store all smap paths for later processing in checkResults -- create before possible return
        if not len(self.targetCmapList) : #there are no jobs to run
            return

        group = self.group_jobs
        ncontig_per_group = 5 #number of contigs per group, currently fixed

        endargs = self.varsP.argsListed(self.moduleName)
        baseargs = [self.varsP.RefAlignerBin]
        if '-RefSplit' in endargs : #for refsplit, need -ref input
            baseargs += ['-ref', self.ref]
        else :
            #-sv is pairsplit which is only implemented for pairwise, so even the reference must be a -i argument
            #validity of reference should have been checked in varsP.checkInputs() at pipeline startup
            baseargs += ['-i', self.ref]
        if self.varsP.bedFile :
            baseargs += ['-bed', self.varsP.bedFile]
        for k,v in self.noisep.iteritems(): #this is empty dict if could not read .err file from characterize
            endargs.extend(["-"+k, str(v)])
        #if self.noisep : #debug
        #    print "Args = ", " ".join(["-"+k+" "+str(v) for k,v in self.noisep.iteritems()])

	wSum, wList = util.contigComputeWeights(self.targetCmapList)
	wLimit=wSum*1.0/len(self.targetCmapList)
        #print "wLimit:", wLimit, "wSum:", wSum, "wList:", wList #debug
	mapOrder=sorted(range(0, len(self.targetCmapList)), key=lambda i: -wList[i])

        #print ("mapOrder (%i)"%len(mapOrder)), mapOrder #debug
        if group : #make the groups: simplest rule is top contigs are combined with bottom (ncontig_per_group-1)
            single_frac = 0.25 #fraction of contigs (longest) to leave as one per job
            nsingle = max(int(len(mapOrder) * single_frac), 1) #if len(mapOrder) < 4, then this is 0--use at least 1
            #ngroup = int(len(mapOrder) / float(ncontig_per_group)) + 1 #no nsingle
            ngroup = int((len(mapOrder)-nsingle) / float(ncontig_per_group)) + 1 if len(mapOrder) > 1 else 0 #need if for single contig
            #groupIds = [[mapOrder[i]] for i in range(ngroup)] #0 to ngroup-1; no single
            #print "nsingle=", nsingle, "ngroup=", ngroup, "len mapOrder=", len(mapOrder) #debug
            groupIds =  [[mapOrder[i]] for i in range(nsingle)] #single contig ids
            if ngroup > 0 :
                groupIds += [[mapOrder[i]] for i in range(nsingle,nsingle+ngroup)] #0 to ngroup-1 -- move to groupIds
            #print "ncontig", len(mapOrder), "nsingle", nsingle, "ngroup", ngroup, "groupIds", groupIds #debug
            #gi = 0 #no single
            gi = len(groupIds)-1 #start from end
            #for i in reversed(range(ngroup,len(mapOrder))) : #ngroup to end; no single
            for i in reversed(range(ngroup+nsingle,len(mapOrder))) : #ngroup to end
                groupIds[gi].append(mapOrder[i]) #ngroup to end
                #print gi, i, mapOrder[i], groupIds[gi] #ngroup to end -- DEBUG
                if len(groupIds[gi]) == ncontig_per_group : #this group is done
                    #gi += 1 #going forward
                    gi -= 1
            #print "ids", groupIds #terse -- DEBUG
            manifest = "SVdetect grouped job manifest\n" #log the 'manifest', ie, which contigs are in which group
            manifest += "\n".join( ["group "+str(i+1)+": "+" ".join([getContigStr(self.targetCmapList[x]) for x in l]) for i,l in enumerate(groupIds)] )
            self.varsP.updatePipeReport(manifest+"\n\n")

        #check that self.output_reference_end is present (ie, contig1), if not, find another one to use instead
        if not any(map(lambda x: x.find(self.output_reference_end+".cmap") != -1, self.targetCmapList)) : #!= -1 means found
            cmappath = (self.targetCmapList[groupIds[0][0]] if group else sorted(self.targetCmapList)[0])
            self.output_reference_end = os.path.basename(cmappath).replace(".cmap","") #still has prefix
            if self.output_reference_end.find("contig") != -1 : #if have this string, use everything after it; if not, null it
                self.output_reference_end = self.output_reference_end[self.output_reference_end.find("contig"):]
            else : #note this means that all of the _r.cmaps will be kept, which wastes disk space, but at least it's not missing anything
                self.output_reference_end = ""

        minthreads=self.varsP.getClusterArgs(self.moduleName, category="MinThreads")
        if minthreads:
            minthreads=Template(minthreads).substitute(maxthreads=self.varsP.maxthreads)
        else:
            minthreads=self.varsP.maxthreads

        looplist = (groupIds if group else mapOrder) #pick list to loop on based on group
        weight = 1.0
        for i,idx in enumerate(looplist) :
            if group : #idx is list of indices in targetCmapList
                #jobargs = zip(["-i"]*len(idx), [self.targetCmapList[j] for j in idx]) #broken...
                jobargs = []
                for j in idx :
                    jobargs += ["-i", self.targetCmapList[j]]
                #print jobargs #debug
                cmappath = self.targetCmapList[idx[0]]
            else : #idx is just idx in targetCmapList
                cmappath=self.targetCmapList[idx]
                jobargs = ["-i", cmappath]
                weight=wList[idx]
	    
	    threadBoost = ceil((weight*1.0/wLimit)**2) if wLimit > 0 else 1
	    if threadBoost<1:
		    threadBoost=1
		    
	    nthreads=float(minthreads)
	    nthreads=int(round(nthreads*threadBoost))
            #change for group: if group, always use maxthreads
	    if nthreads>self.varsP.maxthreads or group:
		    nthreads=self.varsP.maxthreads
	    
            cbase, cname = os.path.split(cmappath)
            cname = cname[:cname.find(".")] #strip suffix
            if group : #get job name
                jname = "group"+str(i+1)
                #cname = jname #for outputfile: it is called the group. Careful. This may break merging...actually, I guess smapList takes care of this?
                cname = cname[:cname.find("_contig")]+"_"+jname #instead of above, take same prefix and replace 'contigX' with 'groupY'
            else :
                jname = getContigStr(cname)
                #if cname.find("_contig") != -1 : #has this string
                #    jname = cname[cname.find("_contig")+1:] #get 'contigXXX' string from path
                #else :
                #    jname = cname[-15:] #just take last N chars
            #jobName = "SV_detect__"+cname #old
            jobName = jname+"__SV_detect" #new
            outputfile = os.path.join(cbase,cname)
            jobargs += ["-o", outputfile ]
            stdoutf = None
            if self.varsP.stdoutlog :
                jobargs.extend( ['-f', '-stdout', '-stderr'] )
                stdoutf = outputfile+".stdout"
            vetostr = '.align$' #if contig 1, just veto .align
            if cmappath.find(self.output_reference_end+".cmap") == -1 : #keep _r for contig1 only
                #jobargs.extend( ['-output-veto-filter', '_r.cmap$'] ) #do not output the _r.cmap
                vetostr = "(_r.cmap|.align)$"
                if self.varsP.onCluster : #if cluster, need to enclose in quotes
                    vetostr = "\'"+vetostr+"\'"
            jobargs.extend( ['-output-veto-filter', vetostr] ) #if not contig 1, veto both
            #print jobargs #debug
            expectedResultFile = outputfile+".smap" #smap is SV map (could also just check xmap)
            self.smapList.append(expectedResultFile)
            #expectedOutputString = cname #afaik, this is only used for cluster args -- also used in report: set same as jobName
            s1Job = mthread.singleJob(baseargs+['-maxthreads', str(nthreads)]+jobargs+endargs, 
                                      jobName, 
                                      expectedResultFile, 
                                      jobName, #expectedOutputString,
                                      maxThreads=nthreads,
                                      clusterLogDir=self.varsP.clusterLogDir,
                                      expectedStdoutFile=stdoutf
                                      )
            self.addJob(s1Job)
        #end for idx in mapOrder 

        outfile = os.path.split(self.smapList[0])[1]
        outfile = outfile[:outfile.rfind("_")]
        self.merge_map_prefix = os.path.join(self.mergesmapdir, outfile+"_merged") #no trailing "."

        self.logArguments()
    #end generateJobList


    #for big bear compatibility, this is only allowed to call multiThreadRunJobs
    def runJobs(self):
        self.multiThreadRunJobs(self.varsP.nThreads, sleepTime = 0.2)


    def checkResults(self):
        """SVdetect: record results of sv calls in pipelineReport and informaticsReport
        Also merge smaps and create merged bed file.
        """        
        self.doAllPipeReport() #see Multithreading.jobWrapper -- this is the pipelineReport
        self.checkSmapResults(not self.group_jobs) #add SV info to informatics report -- if grouped, skip detailed report (not implemented)
        #merge _q.cmap files, need to do first to set data member, return false only if can't make output dir
        #if _r.cmap is not found, the merged files won't point to the right place (but make them anyway)
        if self.combine_output and self.mergeQCmaps() :
            self.copyRCmap()
            #xmaplist = self.mergeSmaps(True) #returns list of xmap paths for xmaps which contain SVs
            #self.mergeSmaps(False) #merge all xmaps
            #self.mergeSmaps(False, xmaplist=xmaplist) #merge only xmaps with SVs, empty xmap if no SVs
            #self.mergeSmaps(True, filtersv=True) #filtered smap -- last smap produced will be used in mergeInversion

            #above old; below new

            xmaplist = self.mergeSmaps_new() #returns bool now for success/fail

            if xmaplist : #if empty, error, do not do these
                self.smapToBed() #don't do this if merge fails, or no SVs

        #fix the varsP data members needed for prepareContigIO for next stage
        self.varsP.outputContigPrefix = self.oldOutputContigPrefix
        util.LogStatus("progress", "stage_complete", self.stageName)

        # Call new PERL methods - Edited by XZ
        self.mergeInversion()
        # Call R scripts
        #self.callIndelConfidence()
        self.callCopyNumber() #R scripts 
    #end checkResults


    #open each smap, analyze results for SVs
    #info to get: same as in run_sv:
    #dict of sv type to number of occurrances
    #sv rate: number of svs over number of maps
    #detailedreport will also print the length of the maps, the total aligned length, and the number of SVs
    # _for each map_. This can be large for human.
    #note: if jobs are grouped, skip detailed report (not implemented)
    def checkSmapResults(self, detailedreport=True) :
        """SVdetect: loop on all maps produced by all sv calls, collect stats
        """
        self.useSmapList = [] #will contain only those smaps found on disk -- make before potential return
        if not len(self.smapList) : #no smaps, nothing to do
            return

        mapunit = 1e6 #print map lens in Mb
        svtypeslen = "24" #for printing
        svtype_dict = defaultdict(int) #don't make any assumptions about sv types--use a default dict

        #this is a warning only
        # if happens, proceed using len(smapList) as number of contigs so as to not penalize the sv rate based on failed jobs
        if len(self.smapList) != len(self.jobList) :
            self.varsP.updatePipeReport("Warning in SVModule: number of jobs not equal to number of smaps\n\n")

        #check presence of each smap prior to loop below because if it doesn't exist, there's not much point
        # and it's a bit misleading to have a 0 entry in the summary when the job failed becuase 0 looks
        # like no alignment
        # do this inside loop below

        #get full path of _r.cmap output based on output_reference_end and targetCmapList
        cbase, cname = os.path.split(self.targetCmapList[0])
        if cname.rfind("contig") != -1 :
            cname = cname[:cname.rfind("contig")] #strip contig suffix
        elif cname.rfind(".") != -1 :
            cname = cname[:cname.rfind(".")] #strip file suffix
        if self.output_reference_end :
            self.refoutpath = os.path.join(cbase,cname+self.output_reference_end+"_r.cmap")
        else :
            self.refoutpath = os.path.join(cbase,cname+"_r.cmap")
        #above defines what the output reference is expected to be, but must check it is actually there -- see copyRCmap
        if not util.checkFile(self.refoutpath) :
            self.refoutpath = ""

        outstr = "SV detect: %s\n" % self.stageName
        outstr += "map filename : len (Mb) : align len (Mb) (ratio) : N SV\n" #header for below
        maxpathlen = 0
        nsvtot = nsvfilt = 0
        self.sfilt = lambda x : not any([x.find(s)!=-1 for s in self.smapFilterTypes]) #add ones which pass this
        for smappath in self.smapList :
            if self.smapTypeIdx == -1 :
                self.getTypeIndex(smappath)
            svdict = self.checksmap(smappath) #new fn returns dict of svtype : number
            if svdict == None : #empty dict is no SV, None is no file
                continue
            self.useSmapList.append(smappath)
            util.addDictToDict(svtype_dict, svdict)
            nsv = (sum(svdict.values()) if svdict else 0)
            nsvfilt += (sum([v for k,v in svdict.iteritems() if self.sfilt(k)]) if svdict else 0)
            nsvtot += nsv
            if detailedreport :
                smapfile = os.path.split(smappath)[1]
                if len(smapfile) > maxpathlen :
                    maxpathlen = len(smapfile)
                cmaplen = mc.cmap(smappath.replace(".smap",".cmap"), lengthonly=True).length/mapunit
                alignlen = self.checkxmap(smappath.replace(".smap",".xmap"), mapunit) #just for kicks
                #print in Mb (mapunit) so use 3 decimal places (6.3)
                arat = (alignlen/cmaplen if cmaplen > 0 else 0)
                outstr += ("%-"+str(maxpathlen)+"s : %6.3f : %6.3f (%5.2f) : %i\n") % (smapfile, cmaplen, alignlen, arat, nsv)
        #end loop on smaps

        #for group, job failure is not checked; for non-group, each smap is checked above, so most failures are caught
        ncontigs = (len(self.targetCmapList) if self.group_jobs else len(self.smapList))
        rate     = (nsvtot/float(ncontigs) if ncontigs > 0 else 0)
        outstr += ("%"+svtypeslen+"s :   N : N/tot\n") % "SV type"
        for key in sorted(svtype_dict.keys()) : #for uniformity of output
            rat = (svtype_dict[key]/float(nsvtot) if nsvtot else 0)
            outstr += ("%"+svtypeslen+"s : %4i : %.4f\n") % (key, svtype_dict[key], rat)
        outstr += "N contigs  : %6i\n" % ncontigs
        outstr += "N sv total : %6i  (%.4f per contig)\n" % (nsvtot, rate)
        outstr += "N sv filter: %6i  (%.4f per contig)\n" % (nsvfilt, (nsvfilt/float(ncontigs) if ncontigs > 0 else 0))

        self.varsP.updateInfoReport(outstr + '\n')
    #end checkSmapResults


    def getTypeIndex(self, usefile):
        """Instead of hardcoding the index of the smap Type column, get it from the header.
        """
        if not util.checkFile(usefile, ".smap") :
            self.varsP.updatePipeReport("Warning in checksmap: File missing or not smap: "+usefile+"\n")
            return

        svdict = defaultdict(int) #int will default to 0, conveniently
        f1 = open(usefile)
        for line in f1 :
            if line.find("SmapEntryID") == -1 : #this is the line with the column descriptions
                continue
            tok = line.split()
            if "Type" in tok :
                self.smapTypeIdx = tok.index("Type") - 1 #always have "#h" start the line, so indices are shifted by 1
            if "LinkID" in tok :
                self.smapLinkIdx = tok.index("LinkID") - 1 #always have "#h" start the line, so indices are shifted by 1
            if "XmapID1" in tok :
                self.smapXmapIdx = tok.index("XmapID1") - 1
            break #done
        f1.close()
        #print "\ngetTypeIndex: index =", self.smapTypeIdx, "\n"
        #self.varsP.updatePipeReport("getTypeIndex: LinkID index = %i\n" % self.smapLinkIdx) #debug
    #end getTypeIndex

    #check if svs found: 
    #return dict of sv type to number of SVs, where type is read from the smap
    # empty dict if no SVs; return None if no file found to avoid further errors
    def checksmap(self, usefile) :
        """SVdetect.checksmap: open smap produced by sv call, get stats on SVs detected
        """
        if not util.checkFile(usefile, ".smap") :
            self.varsP.updatePipeReport("Warning in checksmap: File missing or not smap: "+usefile+"\n")
            return None #None means no file found; empty dict means no SVs

        svdict = defaultdict(int) #int will default to 0, conveniently
        f1 = open(usefile)
        for line in f1 :
            if line[0] == "#" :
                continue
            thissvtype = line.split()[self.smapTypeIdx] #see getTypeIndex
            svdict[thissvtype] += 1
        f1.close()
        return svdict
    #end checksmap

    #check if any alignment, if so, return total length
    #treat mapunit as denominator, ie, 1e6 for Mb
    def checkxmap(self, usefile, mapunit=1. ) :
        """SVdetect.checkxmap: open xmap produced by sv call, get total aligned length
        """
        if not util.checkFile(usefile, ".xmap") : #os.path.isfile(usefile) :
            outstr = "Warning in checkxmap: File missing: " + usefile
            self.varsP.updatePipeReport(outstr+"\n")
            return 0

        #return mc.xmap(usefile).getSumMappedQryLen(unitscale=mapunit) #old version without header replacement
        xmap = mc.xmap(usefile)
        if self.refoutpath : #be sure this was set, otherwise no change to header
            xmap.editHeaderMaps(self.refoutpath, query=False)
            xmap.writeToFile(usefile)
        return xmap.getSumMappedQryLen(unitscale=mapunit) 
    #end checkxmap


    def mergeQCmaps(self) :
        """Loop on all smaps (self.useSmapList), find their corresponding _q.cmap,
        and merge them into one cmap in mergesmapdir.
        """
        if not len(self.useSmapList) :
            return False #fail

        if not util.checkDir(self.mergesmapdir) : #set in self.__init__; makes this dir
            self.varsP.updatePipeReport("ERROR in SVModule.mergeQCmaps: unable to create merged smap dir %s\n" % self.mergesmapdir)
            return False #False for fail

        qcmap = "_q.cmap"
        self.merged_qcmap = self.merge_map_prefix+qcmap
        outfile = ""
        firstmap = True
        nmap = 0
        for smappath in self.useSmapList :
            qpath = smappath.replace(".smap", qcmap)
            if not util.checkFile(qpath) :
                #failed jobs are warned elsewhere, no need to repeat, plus, Pipeline may not produce _q.cmap
                #self.varsP.updatePipeReport("SVdetect.mergeQCmaps: missing file %s\n" % qpath)
                continue
            nmap += 1
            qfile = open(qpath,"rb")
            for line in qfile :
                if firstmap or line[0] != "#" :
                    outfile += line
            qfile.close()
            firstmap = False
        #end for self.useSmapList
        #self.varsP.updatePipeReport("nmap = %i\n" % nmap) #debug
        if nmap > 0 and outfile : #there is data to write
            if nmap != len(self.useSmapList) : #warn that some jobs have missing _q.cmaps
                self.varsP.updatePipeReport("SVdetect.mergeQCmaps: expected %i _q.cmap files, found %i\n" % (len(self.useSmapList), nmap))
            qfile = open(self.merged_qcmap,"wb")
            qfile.write(outfile)
            qfile.close()
        elif self.varsP.latestMergedCmap :
            self.varsP.updatePipeReport("SVdetect.mergeQCmaps: no _q.cmap files found; using merged cmap %s\n" % self.varsP.latestMergedCmap)
            shutil.copy(self.varsP.latestMergedCmap, self.merged_qcmap)
        else :
            self.varsP.updatePipeReport("SVdetect.mergeQCmaps: no _q.cmap files found, merged cmap is invalid\n")
            return False
        return True #True for success
    #end mergeQCmaps


    def copyRCmap(self) :
        """Copy the _r.cmap to the merge dir for easier processing.
        """
        self.merged_rcmap = self.merged_qcmap.replace("_q.cmap","_r.cmap")
        try : #just to be safe...
            shutil.copy(self.refoutpath, self.merged_rcmap)
            return True
        except IOError : #if source isn't present
            if self.refoutpath : #if this was nulled in checkSmapResults, it's not really a warning
                err_msg = "Warning in SVModule.copyRCmap: reference map not found: %s" % self.refoutpath
                self.varsP.updatePipeReport(err_msg+"\n")
                util.LogError("warning", err_msg)
            #old behavior was to return false here; new is to assume that the input reference is also the output reference (should occur only when abs(bpp - 500.0) < 1e-12), and copy that instead below
        #return if no exception above
        try :
            shutil.copy(self.ref, self.merged_rcmap) #input reference
        except IOError : #this should never happen
            self.merged_rcmap = "" #null in order to not change xmap headers
            err_msg = "Warning in SVModule.copyRCmap: reference map not found: %s" % self.ref
            self.varsP.updatePipeReport(err_msg+"\n")
            util.LogError("warning", err_msg)
            return False
        return True #no exception


    #this opens all the smaps again even though they were all opened already in checkSmapResults
    def mergeSmaps(self, dosmap, xmaplist=None, filtersv=False) :
        """SVdetect.mergeSmaps: create one file with all smaps merged.
        If dosmap is False, merge xmaps instead.
        If xmaplist is None, merge all xmaps. If it is a list of full paths, merge only those xmaps,
        if an empty list (no SVs), generate empty _sv.xmap.
        If filtersv, remove SVs (smap entries) based on self.smapFilterTypes.
        """

        assert not (xmaplist and dosmap), "Cannot do smap if xmaplist supplied" #if xmaplist, must not have dosmap

        looplist = (xmaplist if xmaplist else self.useSmapList) #contains only smaps on disk; see checkSmapResults
        if not dosmap and not xmaplist : #need to get xmap paths from smaps
            looplist = [x.replace(".smap",".xmap") for x in self.useSmapList]
        mapstr = ("smap" if dosmap else "xmap")
        if not looplist : #if empty, no smaps, nothing to merge
            self.varsP.updatePipeReport("SVdetect: no "+mapstr+"s found, skipping merge\n")
            return [] #return empty list for fail

        if dosmap and not filtersv : #only print once
            self.varsP.updatePipeReport("SVdetect: merging smaps to %s\n" % self.mergesmapdir)
        # now, find and merge all smaps in usedir
        nmaps = 0
        nsv = 0
        idx = 1 #index for combine_output
        firstmap = True #copy header of first map
        mergedoutput = "" #this is the smap
        outpathlist = [] #if dosmap, populate this with xmap paths iff an sv is present
        headeronly = (not dosmap and xmaplist != None and not xmaplist) #empty list
        for smappath in looplist : 
            # open smap, read it
            try :
                usefile = open(smappath)
            except IOError:
                err_msg = "ERROR in SVModule.mergeSmaps: unable to open file %s, skipping in merged smap\n" % smappath
                self.varsP.updatePipeReport(err_msg)
                util.LogError("error", err_msg)
                continue
            nmaps += 1
            thissv = False
            filter_partial = False #if there is an inversion event which is filtered, the following line may be a partial, in which case it is filtered also
            for line in usefile :
                if firstmap and line[0] == "#" :
                    #also modify header of smap
                    mergedoutput += self.modifyXmapSmapHeader(line) 
                    #if line.find("Query Maps") != -1 and self.merged_qcmap :
                    #    mergedoutput += line.split(":")[0] + ":" + self.ws + self.merged_qcmap + "\n"
                    #elif line.find("Reference Maps") != -1 and self.merged_rcmap :
                    #    mergedoutput += line.split(":")[0] + ":" + self.ws + self.merged_rcmap + "\n"
                    #elif line.find("Xmap Entries") != -1 :
                    #    mergedoutput += line.split(":")[0] + ":" + self.ws + self.merge_map_prefix+".xmap" + "\n"
                    #elif line.find("Smap Entries") != -1 :
                    #    mergedoutput += line.split(":")[0] + ":" + self.ws + self.merge_map_prefix+".smap" + "\n"
                    #else :
                    #    mergedoutput += line
                elif firstmap and headeronly :
                    break #skip remainder of file
                elif line[0] != "#" :
                    entries = line.split()
                    if dosmap and filtersv and not self.sfilt(entries[self.smapTypeIdx]) : #see checkSmapResults
                        if line.find("inversion") != -1 : #flag for inversion_partial--next entry should be removed
                            filter_partial = True
                        continue
                    elif filter_partial : #check if this is a partial if this flag is set; if it is, skip
                        filter_partial = False #will not be a partial for inversion_duplicate, but flag should still be false
                        if line.find("partial") != -1 : 
                            continue
                    nsv += 1
                    thissv = True
                    if self.combine_output :
                        #need to check if LinkID is present in smap (smapLinkIdx) and its value is not -1
                        if dosmap and self.smapLinkIdx > -1 and entries[self.smapLinkIdx] != "-1" :
                            #LinkID should always be one more or less than the index
                            linkid = entries[self.smapLinkIdx]
                            #self.varsP.updatePipeReport("LinkID entry %s, value %s, compare %s\n" % (entries[0], linkid, entries[0])) #debug
                            if linkid == str(int(entries[0])-1) : #link to previous entry
                                mergedoutput += str(idx) + self.ws + self.ws.join(entries[1:self.smapLinkIdx]) + self.ws + str(idx-1) + self.ws + self.ws.join(entries[self.smapLinkIdx+1:]) + "\n"
                            elif linkid == str(int(entries[0])+1) : #link to next entry
                                mergedoutput += str(idx) + self.ws + self.ws.join(entries[1:self.smapLinkIdx]) + self.ws + str(idx+1) + self.ws + self.ws.join(entries[self.smapLinkIdx+1:]) + "\n"
                            #the general case does not yet need to be handeled; print warning, add the line with LinkID -1
                            else : 
                                err_msg = "Warning in SVModule.mergeSmaps: unexpected LinkID in smap %s for line:\n%s" % (smappath, line)
                                self.varsP.updatePipeReport(err_msg+"\n")
                                util.LogError("warning", err_msg)
                                mergedoutput += str(idx) + self.ws + self.ws.join(entries[1:-1]) + self.ws + "-1\n"
                        else : #no LinkID, or it's -1 so it doesn't need to change
                            mergedoutput += str(idx) + self.ws + self.ws.join(entries[1:]) + "\n"
                        idx += 1
                    else :
                        mergedoutput += line
            # end for line in smap
            usefile.close()
            if dosmap and thissv :
                outpathlist.append( smappath.replace(".smap",".xmap") )
            if firstmap : #the first map has been processed
                firstmap = False
                if headeronly : #skip remaining maps
                    break
        # end for smaps in usedir
        if dosmap :
            self.varsP.updatePipeReport("SVdetect: merging smaps: %i Maps with %s SVs\n" % (nmaps, nsv))
        else :
            self.varsP.updatePipeReport("SVdetect: merging xmaps: %i Maps\n" % (nmaps))
        suf = ("_sv." if xmaplist or headeronly else ".")
        suf = ("_filter." if dosmap and filtersv else suf)
        outfile = self.merge_map_prefix + suf + mapstr
        if mapstr == "xmap" :
            self.outxmap_list.append( outfile )
        self.varsP.updatePipeReport(("SVdetect: Writing merged "+mapstr+": %s\n") % outfile)
        outf = open(outfile, "w+")
        outf.write(mergedoutput)
        outf.close()
        if dosmap :
            self.mergedsmap = outfile #for use in smapToBed and elsewhere -- note this will be the last call to this method
            return outpathlist
    #end mergeSmaps


    #mergeSmaps does not handle XmapIDs in the merged smap. Easier if both are made simultaneously.
    def mergeSmaps_new(self) :
        """Create merged smap, filtered merged smap, and merged xmap."""

        if not self.useSmapList :
            self.varsP.updatePipeReport("SVdetect: empty smap list, skipping merge\n")
            return False #return False for fail

        self.varsP.updatePipeReport("SVdetect: merging smaps to %s\n" % self.mergesmapdir)
        #now, find and merge all smaps
        nmaps = 0 #actually n files (grouped)
        #nsv = 0 #redundant with indices
        sidx = 1 #index for combine_output
        sfidx = 1 #index for filtered smap
        xidx = 1 #index in smap
        firstmap = True #copy header of first map
        mergedxmap = mergedoutput = filtoutput = "" #output
        outpathlist = [] #if dosmap, populate this with xmap paths iff an sv is present

        for smappath in self.useSmapList : 
            # open smap, read it
            try :
                usesmap = open(smappath)
            except IOError:
                err_msg = "ERROR in SVModule.mergeSmaps: unable to open file %s, skipping in merged smap\n" % smappath
                self.varsP.updatePipeReport(err_msg)
                util.LogError("error", err_msg)
                continue
            #same for xmap
            xmappath = smappath.replace(".smap",".xmap")
            try :
                usexmap = open(xmappath)
            except IOError:
                usesmap.close() #don't leave above file handle open
                err_msg = "ERROR in SVModule.mergeSmaps: unable to open file %s, skipping in merged xmap\n" % xmappath
                self.varsP.updatePipeReport(err_msg)
                util.LogError("error", err_msg)
                continue
            nmaps += 1

            filter_partial = False #if there is an inversion event which is filtered, the following line may be a partial, in which case it is filtered also

            #two loops: the smap should come first bc it uses N xmap previous
            for line in usesmap :
                thisoutput = "" #be sure to clear -- need this bc use line for both mergedoutput and filtoutput
                thifoutput = "" #filtered
                if firstmap and line[0] == "#" :
                    thisoutput = self.modifyXmapSmapHeader(line)
                    mergedoutput += thisoutput
                    filtoutput += thisoutput
                elif line[0] != "#" :
                    entries = line.split()
                    #fix the xmapIDs
                    if self.smapXmapIdx > -1 :
                        #print "xmapID:", 
                        #xidx starts at 1, but first smap file's entries should not increase
                        entries[self.smapXmapIdx  ] = str(int(entries[self.smapXmapIdx  ]) + xidx-1) #XmapID1
                        if entries[self.smapXmapIdx+1] != "-1" : #for eg inversion_partial, this must remain at -1
                            entries[self.smapXmapIdx+1] = str(int(entries[self.smapXmapIdx+1]) + xidx-1) #XmapID2
                    #need to check if LinkID is present in smap (smapLinkIdx) and its value is not -1
                    if self.smapLinkIdx > -1 and entries[self.smapLinkIdx] != "-1" :
                        #LinkID should always be one more or less than the index
                        linkid = entries[self.smapLinkIdx]
                        #self.varsP.updatePipeReport("LinkID entry %s, value %s, compare %s\n" % (entries[0], linkid, entries[0])) #debug
                        if linkid == str(int(entries[0])-1) : #link to previous entry
                            thisoutput = self.ws + self.ws.join(entries[1:self.smapLinkIdx]) + self.ws + str(sidx-1) + self.ws + self.ws.join(entries[self.smapLinkIdx+1:]) + "\n" #remove idx
                            thifoutput = self.ws + self.ws.join(entries[1:self.smapLinkIdx]) + self.ws + str(sfidx-1) + self.ws + self.ws.join(entries[self.smapLinkIdx+1:]) + "\n" #use filtered index
                        elif linkid == str(int(entries[0])+1) : #link to next entry
                            thisoutput = self.ws + self.ws.join(entries[1:self.smapLinkIdx]) + self.ws + str(sidx+1) + self.ws + self.ws.join(entries[self.smapLinkIdx+1:]) + "\n" #remove idx
                            thifoutput = self.ws + self.ws.join(entries[1:self.smapLinkIdx]) + self.ws + str(sfidx+1) + self.ws + self.ws.join(entries[self.smapLinkIdx+1:]) + "\n" #use filtered index
                        #the general case does not yet need to be handeled; print warning, add the line with LinkID -1
                        else : 
                            err_msg = "Warning in SVModule.mergeSmaps: unexpected LinkID in smap %s for line:\n%s" % (smappath, line)
                            self.varsP.updatePipeReport(err_msg+"\n")
                            util.LogError("warning", err_msg)
                            #thisoutput = str(sidx) + self.ws + self.ws.join(entries[1:-1]) + self.ws + "-1\n"
                            entries[self.smapLinkIdx] = "-1" #remove LinkID
                            thisoutput = self.ws + self.ws.join(entries) + "\n" 
                    else : #no LinkID, or it's -1 so it doesn't need to change
                        thisoutput = self.ws + self.ws.join(entries[1:]) + "\n"
                    mergedoutput += str(sidx) + thisoutput #put idx back
                    sidx += 1
                    #now filtered smap
                    if not self.sfilt(entries[self.smapTypeIdx]) : #see checkSmapResults
                        if line.find("inversion") != -1 : #flag for inversion_partial--next entry should be removed
                            filter_partial = True
                        continue
                    elif filter_partial : #check if this is a partial if this flag is set; if it is, skip
                        filter_partial = False #will not be a partial for inversion_duplicate, but flag should still be false
                        if line.find("partial") != -1 : 
                            continue
                    filtoutput += str(sfidx) + (thifoutput if thifoutput else thisoutput)
                    sfidx += 1
            # end for line in smap
            usesmap.close()

            #now process xmap
            for line in usexmap :
                if firstmap and line[0] == "#" :
                    mergedxmap += self.modifyXmapSmapHeader(line)
                elif line[0] != "#" :
                    entries = line.split()
                    mergedxmap += str(xidx) + self.ws + self.ws.join(entries[1:]) + "\n"
                    xidx += 1
            #end for line in xmap
            usexmap.close()

            if firstmap : #the first map has been processed -- this is used in both smap and xmap loops
                firstmap = False

        # end for smaps in usedir

        #two log lines per file seems a bit redundant, but it's always been like that
        self.varsP.updatePipeReport("SVdetect: merging smaps: %i files with %s SVs\n" % (nmaps, sidx-1)) #-1 bc +=1 at end
        outfile = self.merge_map_prefix + ".smap"
        self.varsP.updatePipeReport(("SVdetect: Writing merged smap: %s\n") % outfile)
        outf = open(outfile, "w+")
        outf.write(mergedoutput)
        outf.close()

        #now make filtered smap
        self.varsP.updatePipeReport("SVdetect: merging filtered smaps: %i files with %s SVs\n" % (nmaps, sfidx-1))
        outfile = self.merge_map_prefix + "_filter.smap"
        self.varsP.updatePipeReport(("SVdetect: Writing merged smap: %s\n") % outfile)
        outf = open(outfile, "w+")
        outf.write(filtoutput)
        outf.close()
        self.mergedsmap = outfile #for use in smapToBed and elsewhere -- note this will be the last call to this method

        #now do xmap
        #suf = ("_sv." if xmaplist or headeronly else ".") #no longer make this one
        self.varsP.updatePipeReport("SVdetect: merging xmaps: %i files with %i alignments\n" % (nmaps, xidx-1))
        outfile = self.merge_map_prefix + ".xmap"
        self.varsP.updatePipeReport(("SVdetect: Writing merged xmap: %s\n") % outfile)
        outf = open(outfile, "w+")
        outf.write(mergedxmap)
        outf.close()
        self.outxmap_list.append( outfile )

        return True
    #end mergeSmaps_new


    def modifyXmapSmapHeader(self, line) :
        """Utility to fix xmap and smap headers."""
        mergedoutput = ""
        if line.find("Query Maps") != -1 and self.merged_qcmap :
            mergedoutput += line.split(":")[0] + ":" + self.ws + self.merged_qcmap + "\n"
        elif line.find("Reference Maps") != -1 and self.merged_rcmap :
            mergedoutput += line.split(":")[0] + ":" + self.ws + self.merged_rcmap + "\n"
        elif line.find("Xmap Entries") != -1 :
            mergedoutput += line.split(":")[0] + ":" + self.ws + self.merge_map_prefix+".xmap" + "\n"
        elif line.find("Smap Entries") != -1 :
            mergedoutput += line.split(":")[0] + ":" + self.ws + self.merge_map_prefix+".smap" + "\n"
        else :
            mergedoutput += line
        return mergedoutput
    #end modifyXmapSmapHeader


    # How to convert an smap to a bed:
    #  First column is chromosome number--this is same as ref contig id. 
    #  Second and third are start and stop position--also same as ref start and ref stop. 
    #  Third is type.
    #  Next four are always ["-1", "+", "0", "0"] for score, top/bottom, and highlighted start/stop
    #  Last is rgb color.
    # includetypes is the types of svs to include--see __init__
    # output is put in same path as smappath, with .smap replaced by .bed
    # also produce query bed file--only difference is contig id is query contig id, 
    #  and positions are query positions
    #  and filename has _query
    def smapToBed(self, makequerybed=True) :
        smappath = self.mergedsmap #convenience
        outpath = smappath.replace(".smap", ".bed")
        #if os.path.exists(outpath): #forget this--just overwrite it
        #    self.varsP.updatePipeReport("ERROR in SVdetect: Output bed file already exists: %s\n" % outpath)
        #    return
        if makequerybed :
            qoutpath = smappath.replace(".smap", "_query.bed")
            #assert not os.path.exists(qoutpath), "Output bed file already exists:"+outpath #forget existence check
        ws = self.ws #Irysview requires tab separated

        self.varsP.updatePipeReport("SVdetect: Converting smap %s to bed\n" % smappath)
        f1 = open(smappath)
        bedstr = "" #this will be the bed file
        qbedstr = "" #query bed (stays empty if not makequerybed)
        for line in f1 :
            if line[0] == "#" :
                continue
            tokens = line.split()
            svtype = tokens[self.smapTypeIdx] #see __init__
            #if svtype in self.smapExcludeTypes : #ditch exclude types for include types
            if not svtype in self.smapIncludeTypes :
                continue
            # 0 : SmapEntryID #1 : QryContigID
            ref = tokens[2] #RefcontigID1
            # 3 : RefcontigID2 #4-5 : QryStartPos - QryEndPos
            refstart = int(float(tokens[6])) #RefStartPos
            refstop  = int(float(tokens[7])) #RefEndPos
            # prevent extra space -- no numerals for % : instead of 0 0 for 2nd and 3rd to last columns, repeat coords
            bedstr += ("%-s"+ws+"%i"+ws+"%i"+ws+"%s"+ws+"-1"+ws+"+"+ws+"%i"+ws+"%i"+ws+"%s\n") % (ref, refstart, refstop, svtype, refstart, refstop, self.getSVcolor(svtype))
            if makequerybed :
                qry = tokens[1]
                qrystart = int(float(tokens[4])) #QryStartPos
                qrystop  = int(float(tokens[5])) #QryEndPos
                qbedstr += ("%-s"+ws+"%i"+ws+"%i"+ws+"%s"+ws+"-1"+ws+"+"+ws+"%i"+ws+"%i"+ws+"%s\n") % (qry, qrystart, qrystop, svtype, qrystart, qrystop, self.getSVcolor(svtype))
        f1.close()

        # write output file
        self.varsP.updatePipeReport("SVdetect: writing bed file %s\n" % outpath)
        outfile = open(outpath, "w+")
        outfile.write(bedstr)
        outfile.close

        if makequerybed :
            self.varsP.updatePipeReport("SVdetect: writing bed file %s\n" % qoutpath)
            outfile = open(qoutpath, "w+")
            outfile.write(qbedstr)
            outfile.close
        self.varsP.updatePipeReport("\n")
    #end smapToBed

    def getSVcolor(self, svtype) :
        return self.svcolors[svtype] if self.svcolors.has_key(svtype) else self.svcolors["default"]


    def copyMergeSmap(self, outpath) :
        '''If mergeInversion cannot be run, or if it fails, copy the previous merged smap so the
        required output file is still present.'''
        if util.checkFile(self.mergedsmap) :
            shutil.copy(self.mergedsmap, outpath) #input reference


    def copyMergeSmapWithError(self, outpath, err_msg) :
        '''Same as copyMergeSmap, but also log error.'''
        self.varsP.updatePipeReport(err_msg+"\n")
        util.LogError("warning", err_msg)
        self.copyMergeSmap(outpath)


    def mergeInversion(self) :
        '''Call Andy's Merge_inversions.pl script.'''

        miout = self.varsP.scriptsMergeInv_outsuf #this is the suffix of the smap which merge inversions makes
        miout = os.path.join(self.mergesmapdir, os.path.basename(self.mergedsmap).replace(".smap",miout))
        #need to change line in two xmaps to point to this file--it will be copied even if this stage is not run or fails
        for xmappath in self.outxmap_list :
            if not util.checkFile(xmappath, ".xmap") :
                self.varsP.updatePipeReport("merged xmap output file missing: %s\n" % xmappath)
                continue
            xmap = mc.xmap(xmappath)
            xmap.editHeaderMaps(miout, query=False, smap=True)
            xmap.writeToFile(xmappath)

        if not self.varsP.perl_path : #if this is None, there is no perl, so cannot run script
            err_msg = "Warning in SVModule.mergeInversions: perl is not found"
            return self.copyMergeSmapWithError(miout, err_msg)

        # see varsPipeline.getScriptPaths
        checkpath = self.varsP.scriptsMergeInv[1] if len(self.varsP.scriptsMergeInv) > 1 else self.varsP.scriptsMergeInv[0]
        if not util.checkFile(checkpath) : #script does not exist--this is an error
            err_msg = "Warning in SVModule.mergeInversions: script or executable missing: %s" % checkpath
            return self.copyMergeSmapWithError(miout, err_msg)

        #check self.mergedsmap just to be safe, otherwise no point invoking script
        if not self.mergedsmap or not util.checkFile(self.mergedsmap,".smap") :
            err_msg = "Warning in SVModule.mergeInversions: merged smap missing: %s" % self.mergedsmap
            return # self.copyMergeSmapWithError(miout, err_msg) #if the file isn't there, you can't copy it

        #the current version of Merge_Inversions.pl requires RefAligner to be at least 3567; skip this stage if that condition is not met
        if self.varsP.refaligner_version < self.varsP.scriptsMergeInv_ver :
            err_msg = "Warning in SVModule.mergeInversions: RefAligner version incompatible: %i; required version is >= %i" % (self.varsP.refaligner_version, self.varsP.scriptsMergeInv_ver)
            return self.copyMergeSmapWithError(miout, err_msg)

        # Default parameters
        Overlap = "30000"
        #Flank_Align_Conf_Cutoff = "12" #no longer used
        #Inv_Align_Conf_Cutoff = "9" #no longer used
        
        self.varsP.updatePipeReport("Starting SV merge-inversions stage\n", printalso=True)
	util.LogStatus("progress", "stage_start", "mergeInversions")

        # <merged_smap_dir> <merged_smap_file> <Overlap> <OUTPUT_DIR>
        perl_command =  self.varsP.scriptsMergeInv + list(os.path.split(self.mergedsmap))
        perl_command += [#self.varsP.outputContigFolder, #no longer used
                         #self.oldOutputContigPrefix+"_contig", #the _contig suffix is missing from this argument (added in findContigs) -- no longer used
                         Overlap, #Flank_Align_Conf_Cutoff, Inv_Align_Conf_Cutoff, #these two no longer used
                         self.mergesmapdir] #now output in same merged_smap dir as rest of results
        #print " ".join(perl_command)

        return_code, stdout, stderr = util.runJob(perl_command, returnstdout=True)

        if not util.checkFile(miout, ".smap") :
            err_msg = ("Merge-inversions output file missing: %s; copying merged smap" % miout)
            self.copyMergeSmapWithError(miout, err_msg)
            return_code = 1

        if return_code == 0:
            self.varsP.updatePipeReport(stdout+"\n")
            self.varsP.updatePipeReport("SV merge-inversions stage successfully finished\n\n")
        else:
            err_msg = "SV merge-inversions stage failed"
            self.varsP.updatePipeReport("%s:\ncommand: %s\nstdout: %s\nstderr: %s\n" % (err_msg, " ".join(perl_command), stdout, stderr))
            util.LogError("error", err_msg)
	util.LogStatus("progress", "stage_complete", "mergeInversions")
    #end mergeInversion


    def callIndelConfidence(self) :
        '''Call Zeljko's R script for indel confidence.'''

        if not self.varsP.R_path : #if no R, cannot run
            return

        stagename = "indelconfidence" #name in xml, also output dir
        outdir = os.path.join(self.varsP.outputContigFolder, stagename) #subdir of sv dir
        util.checkDir(outdir)
        svtypes = ['insertion', 'deletion'] #used for svType arg to R script, also for log filename
        self.varsP.updatePipeReport("Starting indel (SV) confidence stage\n")
        try : #argsListed will raise KeyError if the stage is not present in the xml
            con_args = self.varsP.argsListed(stagename)
        except KeyError : #set defaults here
            con_args = ["flankLength", "5e4", "alpha", "0.10", "beta", "0.20",
                        "quantileAcross", "0.95", "quantileFlanks", "0.05", "minNMol", "10"]
            self.varsP.updatePipeReport("Warning: indelconfidence args not found in xml, using default args:\n"+" ".join(con_args)+"\n\n")
        #print con_args
        con_args = [con_args[i]+"="+con_args[i+1] for i in range(0,len(con_args),2)]
        con_args[-1] = con_args[-1]+"'" #last ele here is end of args list which is in single quote
        #print con_args
        args = [self.varsP.R_path, "CMD", "BATCH", "--no-save", "--no-restore",
                "'--args", #R seems to be very particular for some reason about quoting--all args in single quote, each in double
                "scriptDir=\""+self.varsP.rutil_dir+"\"",
                "svDir=\""+self.varsP.outputContigFolder+"\"",
                "molVsContigDir=\""+self.varsP.alignMolDir+"\"", #should check that this is not empty, ie, in case AlignModule not run
                "molVsRefDir=\""+self.varsP.alignMolvrefDir+"\"", 
                "outputDir=\""+outdir+"\"" ]

        #import shlex #this is a tool for using quotes in strings, which we apparently need to do to use R -- this is not needed, see callCopyNumber below
        script = self.varsP.scriptsIndelConf # script path
        testfile = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test.R")
        for svt in svtypes :
            usearg = args + ["svType=\""+svt.upper()+"\""] + con_args
            logpath = os.path.join(outdir, stagename+"_"+svt+".log")
            usearg += [script, logpath] #last arg is log
            #print " ".join(usearg)
            #retcode = util.runJob(usearg) #quotes not handled properly
            #splitarg = shlex.split(" ".join(usearg)) #do not use this, it's unnecessary
            retcode = util.runJob(splitarg) #because batch mode is used, there is no stdout
            if retcode : #non-zero is fail
                err_msg = "Indel confidence call failed for type %s, check log %s\n\n" % (svt, logpath)
                self.varsP.updatePipeReport(err_msg)
                util.LogError("error", err_msg)

            break #TEST

        self.varsP.updatePipeReport("Finished indel (SV) confidence stage\n\n") 


    def callCopyNumber(self) :
        """Call Zeljko's R scripts for copy number profiles
        """
        if not self.varsP.R_path : #if no R, cannot run
            return
        if not self.varsP.doAlignMolvRef :
            err_msg = "Warning in SVModule.callCopyNumber: skipping copy number because no reference alignments run"
            self.varsP.updatePipeReport(err_msg+"\n")
            util.LogError("warning", err_msg)
            return 
        if self.checkReference() : #warn inside this method on failure
            return

        #because AlignModule cannot store in varsP (because it is run in backgrou), need to recreate merge path
        tmp = os.path.join(self.varsP.contigFolder, self.varsP.alignMolvrefName) #alignmolvref
        tmp = os.path.join(tmp, self.varsP.alignMolvrefMergeName) #sub-dir: merge
        mergepath = os.path.join(tmp, self.varsP.alignMolvrefName+"_"+self.varsP.alignMolvrefMergeName+"_r.cmap")
        if not util.checkFile(mergepath) : 
            err_msg = "Warning in SVModule.callCopyNumber: skipping copy number because missing merged alignments file %s" % mergepath
            self.varsP.updatePipeReport(err_msg+"\n")
            util.LogError("warning", err_msg)
            return 

        #all requirements to run have passed, start preparing output
        alignTarget = os.path.join(self.varsP.contigFolder, self.varsP.alignMolvrefName) #this is where alignmolvref goes
        outdir = os.path.join(alignTarget, "copynumber")
        if not util.checkDir(outdir) :
            err_msg = "Warning in SVModule.callCopyNumber: skipping copy number because unable to make output dir: %s" % outdir
            self.varsP.updatePipeReport(err_msg+"\n")
            util.LogError("warning", err_msg)
            return             

        self.varsP.updatePipeReport("Starting copy number profiles script\n", printalso=True)
	util.LogStatus("progress", "stage_start", "copyNumber")

        #prepare all R script args
        fixedargs1 = '--args performAlignment=FALSE inputMoleculesFile=NULL reference="hg19" minLength=150 tValue=9 refAligner=NULL binSize="50kbp" performScaling=FALSE detectAneuploidy=TRUE detectSubchromosomalFeatures=TRUE compareWithSVs=TRUE breakPointTolerance=2.5e4'
        scriptfolder = 'scriptFolder="%s"' % self.varsP.rutil_dir
        cnfolder = 'copyNumberPackageFolder="%s"' % self.varsP.scriptsCopyNum_dir
        outfolderarg = 'outputFolder="%s"' % outdir
        prefarg = 'outputFilePrefix="%s"' % self.varsP.expID
        svfolderarg = 'svFolder="%s"' % self.mergesmapdir #needs to be merged dir
        mapmolf = 'mapMoleculeFile="%s"' % mergepath
        scriptargs = " ".join([fixedargs1, scriptfolder, cnfolder, outfolderarg, prefarg, svfolderarg, mapmolf])
        scriptname = os.path.join(self.varsP.scriptsCopyNum_dir, self.varsP.scriptsCopyNum_top) #USE ME
        logpath = os.path.join(outdir, self.varsP.expID+'.log')
        args = [self.varsP.R_path, "CMD", "BATCH", "--no-save", "--no-restore",
                scriptargs, scriptname, logpath]

        #print args
        #print " ".join(args)
        retcode = util.runJob(args) #because batch mode is used, there is no stdout
        if retcode : #non-zero is fail
            err_msg = "Copy number profiles call failed, check log %s\n\n" % (logpath)
            self.varsP.updatePipeReport(err_msg)
            util.LogError("warning", err_msg)
        else :
            self.varsP.updatePipeReport("Copy number profiles successfully completed\n")
	util.LogStatus("progress", "stage_complete", "copyNumber")
    #end callCopyNumber

    def checkReference(self) :
        """Copy number: right now, only hg19 is supported--check this is hg19 using
        chromosome lengths, hardcoded in this method.
        """

        if self.varsP.load_hg19 and self.varsP.ref : #in this case, it was loaded from the CopyNumberProfiles dir, so assume it's correct
            return 0 #success

        #this check is redundant with the one in varsPipeline.checkInputs, but do it to be safe
        if not self.varsP.ref or not util.checkFile(self.varsP.ref) :
            err_msg = "Warning in SVModule.checkReference: reference missing: %s" % self.varsP.ref
            self.varsP.updatePipeReport(err_msg+"\n")
            util.LogError("warning", err_msg)
            return 1 #return 1 for error -- stop further processing

        hg19_chr = {1 : 249250621.0 ,
                    2 : 243199373.0 ,
                    3 : 198022430.0 ,
                    4 : 191154276.0 ,
                    5 : 180915260.0 ,
                    6 : 171115067.0 ,
                    7 : 159138663.0 ,
                    8 : 146364022.0 ,
                    9 : 141213431.0 ,
                    10:  135534747.0,
                    11:  135006516.0,
                    12:  133851895.0,
                    13:  115169878.0,
                    14:  107349540.0,
                    15:  102531392.0,
                    16:  90354753.0 ,
                    17:  81195210.0 ,
                    18:  78077248.0 ,
                    19:  59128983.0 ,
                    20:  63025520.0 ,
                    21:  48129895.0 ,
                    22:  51304566.0 ,
                    23:  155270560.0,
                    24:  59373566.0}

        refmcmap = mc.multiCmap(self.varsP.ref)

        valid = True #here, assume valid
        if len(refmcmap.cmapdict) != len(hg19_chr) :
            valid = False
        for chrn in sorted(hg19_chr.keys()) :
            if not valid :
                break
            if not refmcmap.cmapdict.has_key(chrn) :
                valid = False
            if refmcmap.cmapdict[chrn].length != hg19_chr[chrn] :
                valid = False
        if not valid:
            err_msg = "Warning in SVModule.checkReference: reference is not hg19: skipping copy number profiles" #probably not necessary to print why
            self.varsP.updatePipeReport(err_msg+"\n\n")
            util.LogError("warning", err_msg)
            return 1 #return 1 for error -- stop further processing

        return 0 #no errors, ref is hg19, return false (0)
    #end checkReference


#end class SVdetect(mthread.jobWrapper):


#given string like
#"/home/users/wandrews/analysis/ecoli_golden/run17/contigs/exp_refineFinal1_sv/exp_refineFinal1_contig15.cmap"
#get "contig15"
def getContigStr(cname, endlen=15) :
    if cname.rfind(".") != -1 : #strip suffix
        cname = cname[:cname.rfind(".")]
    if cname.rfind("contig") != -1 : #has this string
        jname = cname[cname.rfind("contig"):] #get 'contigXXX' string from path
    else :
        jname = cname[-endlen:] #just take last N chars
    return jname
