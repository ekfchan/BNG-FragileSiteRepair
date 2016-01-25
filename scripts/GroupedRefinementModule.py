import subprocess
import time
import os
import shutil
import pdb
import utilities as util
from math import *
from string import Template

import Multithreading as mthread

"""@package GroupedRefinementModule
Package runs refinement phases: RefineA, RefineB, refineNGS, and refineFinal.
Grouped refinement splits refinement into two parts in order to increase CPU
efficiency. First part does alignment of molecule maps to contigs, and second
part does refinement.

"""


util.setVersion("$Id: GroupedRefinementModule.py 4120 2015-09-17 20:58:31Z wandrews $")

           
#see comment above init
class Refine(mthread.jobWrapper):
    """Replaces the old RefineA and RefineB classes by merging them
    
    Combines RefineA, RefineB, refineNGS, and refineFinal
    """
    #refineStage is a string specifying what stage of refinement you're in
    #must be 'refineA', 'refineB', 'refineNGS', or 'refineFinal' (later)
    def __init__(self, StageName, varsP):
        self.refineStage = StageName
        self.multigroup = True #if False, force single group (not normally good)
        self.varsP = varsP
        ContigPrefix = self.varsP.expID + "_" + StageName
        
        if StageName=="extension0":
		self.varsP.extensionCount += 1
		
        
	for case in util.switch(StageName):
		if case("refine(B0|B1|Final0|Final1)", regexp=True):
			self.bunching=12
			self.ref_arg="-reff"
			break
		if case("refineA"):
			self.bunching=12
			self.ref_arg="-ref"
			break
		if case("refineNGS"):
			self.bunching=1
			self.ref_arg="-ref"
			self.varsP.inputContigPrefix = self.varsP.ngsContigPrefix
			self.varsP.inputContigFolder = self.varsP.ngsInDir
			break
		if case("extension[01]", regexp=True):
			self.bunching=12
			self.ref_arg="-reff"
			ContigPrefix = self.varsP.expID + "_"+ StageName+'_%s' % self.varsP.extensionCount
			break;
		if case():
			#varsP.error += 1 #these don't do anything
			#varsP.message += '  Error: Refine stage name invalid: '+str(StageName)+'\n'
			self.varsP.updatePipeReport("Internal error: unknown stage %s" % StageName)
			return

        clusargs = varsP.getClusterArgs(StageName) #get arguments before changing StageName, then add suffix
        StageName += (("_%i" % self.varsP.extensionCount) if StageName.startswith("extension") else "") #for status.xml only
        self.varsP.stageName=StageName
        util.LogStatus("progress", "stage_start", StageName)
        #super is more pythonic than referring to the base class explicitly (only matters for multiple inheritance)
        super(Refine, self).__init__(varsP, StageName, clusterArgs=clusargs)
        #intermediateContigPrefix = self.varsP.expID + self.StageName.replace("refine", "_r")
	self.varsP.prepareContigIO(ContigPrefix, StageName)	
        #modify results of varsP.prepareContigIO for special case of refineNGS        
        self.generateJobList()
        
    def runJobs(self):
        self.multiThreadRunJobs(self.varsP.nThreads, sleepTime = 0.2)
        
    def writeIDFile(self, nJobs):
        f1 = open(self.varsP.idFile, 'wb')
        f1.write(str(nJobs))
        f1.close()
        
    def groupContigName(self, i) :
        return os.path.join(self.varsP.outputContigFolder, self.varsP.outputContigPrefix+"_group"+str(i))
        
    def findGroupedContigs(self):
	print self.varsP.inputContigFolder, self.varsP.inputContigPrefix
	L=getattr(self.varsP, "count_"+self.varsP.inputContigPrefix)

	#L=[]
	#for i in range(1, count+1):
	#	L.append((str(i), os.path.join(self.varsP.inputContigFolder, self.varsP.inputContigPrefix+"_group"+str(i))))
	
	return L
        
    def groupContigs(self):
        contigFiles, contigIDs = self.varsP.findContigs(self.varsP.inputContigFolder, self.varsP.inputContigPrefix)
        #self.writeIDFile(len(contigFiles))
        M=max(map(int, contigIDs))
        #print(M, len(contigFiles))
        if M<len(contigFiles):
		print("Failure in computing max contig ID: have %d for %d files" %(M, len(contigFiles)))
		raise Exception,("Failure in computing max contig ID: have %d for %d files" %(M, len(contigFiles)))
		
        self.writeIDFile(M) 
        
        if self.bunching==1:
		return zip(contigIDs, contigIDs, contigFiles, [1]*len(contigIDs))

	if self.refineStage == "refineA":
		# For grouping refine A it is important that contigs are in order
		contigIDs=range(1, len(contigFiles)+1)
		base=os.path.join(self.varsP.inputContigFolder, self.varsP.inputContigPrefix)
		contigFiles=[base+"_contig"+str(m)+".cmap" for m in contigIDs]

	# Decrease bunching if we have too few contigs to analyze
	if len(contigFiles)< self.bunching*300:
		self.bunching=len(contigFiles)/300
		if self.bunching<1:
			self.bunching=1

	# Increase bunching if we have too many contigs to avoid going past open files limit (4096)
	if self.bunching>1 and len(contigFiles)>self.bunching*1000:
		self.bunching=len(contigFiles)/1000
		if self.bunching<1:
			self.bunching=1
	
	# Compute contig weights
	wSum, wList = util.contigComputeWeights(contigFiles)
	wLimit=wSum*self.bunching*1.0/len(contigFiles)
	
	# Order contigs to be descending in weight
	if self.refineStage == "refineA":
		# RefineA cannot be sorted this way, however, the output from Assembler is already sorted
		contigOrder=range(0, len(contigFiles))
	else:
		contigOrder=sorted(range(0, len(contigFiles)), key=lambda i: -wList[i])
	
	# Do not group the last tail_len contigs so that the jobs finish faster
	tail_len=60
        k=0
        waccum=0.0
        gid=1
        L=[]
        manifest_file=open(self.groupContigName("_manifest"), "wb")
        manifest_file.write("# Manifest file for " + str(self.refineStage))
        cnt=0
        group_file=None
        for m in contigOrder:
		contigFile=contigFiles[m]
		contigID=contigIDs[m]
		contigWeight=wList[m]

		# Make sure that extra large contigs are in one group
		if group_file and ((not self.multigroup) or (k>=self.bunching) or (waccum>wLimit) or (cnt+tail_len>len(contigFiles))) :
			group_file.close()
			w=L[len(L)-1]
			L[len(L)-1]=(w[0], w[1], w[2], (waccum*1.0/wLimit)**2)
			k=0
			waccum=0.0
			
		
		if k==0:
			#fname=self.groupContigName(contigID)
			fname=self.groupContigName(gid)
			#L.append((str(contigID), fname))
			L.append((str(gid), contigID, fname, (contigWeight*1.0/wLimit)**2))
			gid+=1
			#print "Creating", fname
			manifest_file.write("\n"+str(contigID))
			group_file=open(fname, "wb")
		else:
			manifest_file.write(" "+str(contigID))			
		group_file.write(contigFile+"\n")
		k+=1
		waccum+=contigWeight
		cnt+=1
			
	if k>0:
		group_file.close()
		w=L[len(L)-1]
		L[len(L)-1]=(w[0], w[1], w[2], (waccum*1.0/wLimit)**2)
        manifest_file.write("\n")
	manifest_file.close()
	return L
	

    def generateJobList(self):
        baseArgs1 = self.varsP.argsListed(self.refineStage)
        
	for case in util.switch(self.refineStage):
		if case("refine(B1|Final1)", regexp=True):
			baseArgs1 += self.varsP.argsListed('noise0')
			ContigGroupList=self.findGroupedContigs()
			r1args = [self.varsP.RefAlignerBin]
			break
		if case("refine(B0|Final0)", regexp=True):
			baseArgs1 += self.varsP.argsListed('noise0')
			ContigGroupListFull=self.groupContigs()
			setattr(self.varsP, "count_"+self.varsP.outputContigPrefix, (ContigGroupListFull))
			#print self.varsP.outputContigPrefix, getattr(self.varsP, "count_"+self.varsP.outputContigPrefix)
			#r1args = [self.varsP.RefAlignerBin, '-i', self.varsP.bnxFile]
			#InputFileList=[self.varsP.bnxFile]
			r1args = [self.varsP.RefAlignerBin]
			ContigGroupList = zip(range(1,self.varsP.nPairwiseJobs + 1), range(1,self.varsP.nPairwiseJobs + 1), [self.varsP.bnxFile.replace(".bnx", "_%s_of_%s.bnx" %(x, self.varsP.nPairwiseJobs)) for x in range(1,self.varsP.nPairwiseJobs + 1)], [1]*self.varsP.nPairwiseJobs)
			break
		if case("refineA"):
			baseArgs1 += self.varsP.argsListed('noise0')
			ContigGroupList=self.groupContigs()
			#print("Found %d groups for refineA" % (len(ContigGroupList)))
			#r1args = [self.varsP.AssemblerBin, '-i', self.varsP.bnxFile.replace(".bnx", "_sorted.bnx")] #need this before -contigs -- can no longer use all_sorted.bnx due to scan scaling: must refer to varsP.sorted_file
			#r1args = [self.varsP.AssemblerBin, '-i', self.varsP.sorted_file+".bnx"] #need this before -contigs
			r1args = [self.varsP.AssemblerBin, '-if', self.varsP.bnxFileList] #need this before -contigs; use split files in case splitting changed (eg due to scan scaling producing labels at < 20 bp)
			r1args += ['-contigs', os.path.join(self.varsP.inputContigFolder, self.varsP.inputContigPrefix) + '.contigs']
			break
		if case("refineNGS"):
			r1args = [self.varsP.RefAlignerBin, '-i', self.varsP.bnxFile]
			ContigGroupList=self.groupContigs()
			break
		if case("extension0"):
			baseArgs1 += self.varsP.argsListed('noise0')
			ContigGroupList=self.groupContigs()
			setattr(self.varsP, "count_"+self.varsP.outputContigPrefix, (ContigGroupList))
			#print self.varsP.outputContigPrefix, getattr(self.varsP, "count_"+self.varsP.outputContigPrefix), self.varsP.inputContigFolder, self.varsP.inputContigPrefix
			#r1args = [self.varsP.RefAlignerBin, '-i', self.varsP.bnxFile]
			#InputFileList=[self.varsP.bnxFile]
			r1args = [self.varsP.RefAlignerBin]
			ContigGroupList = zip(range(1,self.varsP.nPairwiseJobs + 1), range(1,self.varsP.nPairwiseJobs + 1), [self.varsP.bnxFile.replace(".bnx", "_%s_of_%s.bnx" %(x, self.varsP.nPairwiseJobs)) for x in range(1,self.varsP.nPairwiseJobs + 1)], [1]*self.varsP.nPairwiseJobs)
			break;
		if case("extension1"):
			baseArgs1 += self.varsP.argsListed('noise0')
			ContigGroupList=self.findGroupedContigs()
			r1args = [self.varsP.RefAlignerBin]
			break;
		if case():
			varsP.error += 1
			varsP.message += '  Error: Refine stage name invalid: '+str(StageName)+'\n'
			return


        stdarg = []
        if self.varsP.stdoutlog : #this is the same for all cases below
            stdarg = ['-stdout', '-stderr'] 

	#contigFiles, contigIDs = self.varsP.findContigs(self.varsP.inputContigFolder, self.varsP.inputContigPrefix)
        #nJobs = len(contigFiles)
        output1String = os.path.join(self.varsP.outputContigFolder, self.varsP.outputContigPrefix)
        #for jobNum in range(1,nJobs + 1):
            #contigID = contigIDs[jobNum - 1]
        for m in range(0, len(ContigGroupList)):
	    contigID=ContigGroupList[m][0]
	    rawContigID=ContigGroupList[m][1]
	    contig=ContigGroupList[m][2]
	    
	    # Figure out desired number of threads to use
	    threadBoost=ceil(ContigGroupList[m][3])
	    if threadBoost<1:
		    threadBoost=1
	    minthreads=self.varsP.getClusterArgs(self.refineStage, category="MinThreads")
	    if minthreads:
		minthreads=Template(minthreads).substitute(maxthreads=self.varsP.maxthreads)
	    else:
		minthreads=self.varsP.maxthreads
	    nthreads=float(minthreads)
	    nthreads=int(round(nthreads*threadBoost))
	    if nthreads>self.varsP.maxthreads:
		    nthreads=self.varsP.maxthreads
#        for contigID, contig in ContigGroupList :		
            jobName = self.refineStage + ' %5s' % contigID
	    for case in util.switch(self.refineStage):
		    if case("refineA"):
			endId=int(rawContigID)+self.bunching-1
			if m+1<len(ContigGroupList) :
				endId=int(ContigGroupList[m+1][1])-1
			currentArgs = [str(rawContigID), str(endId), '-maxthreads', str(nthreads)] #this must come after r1args because it's actually an argument to -contigs
			#currentArgs = r1args + currentArgs + baseArgs1 + ['-id', str(contigID), '-i', contig+"_mapped.bnx", '-o', output1String]
			currentArgs = r1args + currentArgs + ['-o', output1String] + stdarg + baseArgs1
			expectedOutputString = self.varsP.outputContigPrefix + '_contig' + str(rawContigID)
			expectedResultFile = os.path.join(self.varsP.outputContigFolder, expectedOutputString + '.cmap') #refineB
			expectedStdoutFile = output1String + "_id"+str(rawContigID)+".stdout"
			break
			
		    #if case("refineB1|refineFinal1|extension1", regexp=True):
			## TODO: make thread number configurable from clusterArgs
			#currentArgs = ['-maxthreads', str(16), self.ref_arg, contig]
			#currentArgs = r1args + currentArgs + baseArgs1 + ['-id', str(contigID), '-i', contig+"_mapped.bnx", '-o', output1String]
			#expectedOutputString = self.varsP.outputContigPrefix + '_contig' + str(contigID)
			#expectedResultFile = os.path.join(self.varsP.outputContigFolder, expectedOutputString + '.cmap') #refineB
			#break

		    if case("refineB1|refineFinal1|extension1", regexp=True):
			Inputs=zip(["-i"]*self.varsP.nPairwiseJobs, [contig.replace("_group", "_group"+str(i)+"_mapped_group")+".bnx" for i in range(1,self.varsP.nPairwiseJobs + 1)])
			Inputs=[x for t in Inputs for x in t]
                        #-id must come before -o, otherwise expectedStdoutFile is wrong
			currentArgs = ['-maxthreads', str(nthreads), '-id', str(contigID), '-o', output1String, self.ref_arg, contig]
			currentArgs = r1args + currentArgs + stdarg + baseArgs1 + Inputs 
			expectedOutputString = self.varsP.outputContigPrefix + '_contig' + str(rawContigID)
			expectedResultFile = os.path.join(self.varsP.outputContigFolder, expectedOutputString + '.cmap') #refineB
			expectedStdoutFile = output1String + "_id"+str(contigID)+".stdout"
			break
			
		    #if case("refineB0|refineFinal0|extension0", regexp=True):
			#currentArgs = ['-maxthreads', str(self.varsP.maxthreads), self.ref_arg, contig]
			#currentArgs = r1args + currentArgs + baseArgs1 + ['-mapped-unsplit', '1', '-refine', '0', '-mapped', contig+"_mapped", "-o", "/dev/null"]
			#expectedOutputString =  self.refineStage + "contig"+str(contigID) + "_mapped.bnx"
			#expectedResultFile = contig + "_mapped.bnx" #refineB
			#break
			
		    if case("refineB0|refineFinal0|extension0", regexp=True):
			currentArgs = [ '-maxthreads', str(nthreads), "-ref", os.path.join(self.varsP.inputContigFolder, util.uniquifyContigName(self.varsP.inputContigPrefix)+".cmap")]
                        outputfile = os.path.join(self.varsP.outputContigFolder, self.varsP.outputContigPrefix+'_group'+str(contigID))
                        #-id must come before -o, otherwise expectedStdoutFile is wrong
			currentArgs = r1args + ['-i', contig, '-id', str(contigID), '-o', outputfile] + stdarg + currentArgs + baseArgs1
                        currentArgs += ['-refine', '0', '-grouped', os.path.join(self.varsP.outputContigFolder, self.varsP.outputContigPrefix+'_group_manifest'), '-mapped', os.path.join(self.varsP.outputContigFolder, self.varsP.outputContigPrefix+'_group'+str(contigID)+"_mapped"), '-output-filter', ".*.bnx"]
			expectedOutputString = self.varsP.outputContigPrefix+'_group'+str(contigID)+"_mapped.bnx"
			expectedResultFile = outputfile + "_mapped_group1.bnx" 
			expectedStdoutFile = outputfile + "_id"+str(contigID)+ ".stdout"
			break
			
		    if case():
			self.varsP.updatePipeReport("Internal error: cannot handle stage %s" % (self.refineStage))
			raise ValueError
                
            if self.varsP.bnxStatsFile!=None:
		currentArgs.extend(['-XmapStatRead', self.varsP.bnxStatsFile])
		    		    
            s1Job = mthread.singleJob(currentArgs, 
                                    jobName, 
                                    expectedResultFile, 
                                    expectedOutputString,
                                    maxThreads=nthreads,
                                    clusterLogDir=self.varsP.clusterLogDir,
                                    expectedStdoutFile=expectedStdoutFile,
                                    )
            self.addJob(s1Job)
        self.logArguments()
            
    def checkResults(self, stageSuffix=""):
        '''Call jobWrapper (self) .doAllPipeReport, and varsP.mergeIntoSingleCmap.
        stageSuffix, if supplied, is appended to varsP.stageComplete in order to
        fix the stage name reported by the CharacterizeModule in the informaticsReport.
        '''
        self.doAllPipeReport()
        self.varsP.stageComplete = self.refineStage + stageSuffix
        if self.refineStage not in ['refineB0', 'refineFinal0', 'extension0'] :
		self.varsP.mergeIntoSingleCmap()
        StageName = self.refineStage + ("_%i" % self.varsP.extensionCount if self.refineStage.startswith("extension") else "") #for status.xml only
        util.LogStatus("progress", "stage_complete", StageName)

    def endStage(self) :
        """Call this in place of checkResults when this stage is bypassed."""
        if self.refineStage not in ['refineB0', 'refineFinal0', 'extension0'] :
            self.varsP.mergeIntoSingleCmap()
        StageName = self.refineStage + ("_%i" % self.varsP.extensionCount if self.refineStage.startswith("extension") else "") #for status.xml only
        util.LogStatus("progress", "stage_complete", StageName)

