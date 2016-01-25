#import subprocess
#import time
import os
#import pdb

import Multithreading as mthread
import manageTargetLocation as mtl
import molecule
import utilities as util

"""@package ImageProcessingModule Extract BNX data from raw tiff data

Currently not functional for cluster (file transfer considerations)
Read tiff, run DM-static, Rescale for lambda, generate bnx file
"""


util.setVersion("$Id: ImageProcessingModule.py 3566 2015-02-13 20:30:19Z wandrews $")


#quality is to turn on QX lines in bnx
#forceonecolor is for debugging only -- use False
#now return 1 for any errors, return 0 (or None) for success
def performImageAnalysis(varsP, bypass=False, quality=True, forceonecolor=False):
    """Top level function for instrument scaling, image handling, bnx encoding
    
    """
    #print "bypass = ", bypass #debug
    
    remoteDataLocations = parseExperimentFile(varsP.imgFile)
    if not remoteDataLocations :
        varsP.updatePipeReport("ERROR in performImageAnalysis: no images found in paths in "+varsP.imgFile+"\n")
        return 1 #this is an error--new convention is to return 1 on error
    
    processImagesJobSet = mthread.jobWrapper(varsP, groupName = 'Image Transfer and Processing', throttle=8)
    expID = 1
    devices = []
    allJobs = []
    bnxFiles = [] #only used if bypass--see below
    for remoteDataLocation in remoteDataLocations:
        expTag = '%02d' % expID
        expID += 1
        localPath = os.path.join(varsP.localRoot, expTag + '/')
        #if expTag == "04" : #debug
        #print "data\n:", remoteDataLocation, "\n" #debug
        curDevice = Device(varsP, expTag, remoteDataLocation, localPath, bypass, quality, forceonecolor)
        for sJob in curDevice.getTargetJobs():
            if sJob : #empty list is returned if target mol file already exists
                processImagesJobSet.addJob(sJob)
        devices.append(curDevice)
        if bypass :
            bnxFiles.append(curDevice.bnxFile)
     
    if bypass:
        return #bnxFiles #no longer need to return anything--this is not an error

    processImagesJobSet.multiThreadRunJobs(varsP.nThreads, sleepTime =0.25)
    #pipeReport += processImagesJobSet.makeRunReport()
    #pipeReport += processImagesJobSet.makeParseReport()
    #varsP.updatePipeReport(pipeReport, printalso=False)
    processImagesJobSet.doAllPipeReport()
    
    if varsP.lambdaRef:
        mapLambdaJobSet = mthread.jobWrapper(varsP, groupName = 'Map Lambda')
    for device in devices:
        device.findSNRCutoff()
        if varsP.lambdaRef:
            for sJob in device.getLambdaMapJobs():
                if sJob :
                    mapLambdaJobSet.addJob(sJob)
                
    if varsP.lambdaRef:
        mapLambdaJobSet.multiThreadRunJobs(varsP.nThreads, sleepTime =0.25)
        pipeReport =  mapLambdaJobSet.makeRunReport()
        pipeReport += mapLambdaJobSet.makeParseReport()
        for device in devices:
            device.processLambdaMapResults()
        varsP.updatePipeReport(pipeReport, printalso=False)    

    #still need to writeCorrectedBnx
    targetLog = ''
    scanLog = ''
    deviceLog = ''
    bnxFiles = []
    for i,device in enumerate(devices):
        device.writeCorrectedBnx()
        if device.bnxFile : #this is nulled if above failed
            bnxFiles.append(device.bnxFile)
        if i == 0:
            #targetLog += device.getTargetReport(headerOnly=True) + '\n'
            deviceLog += device.getDeviceReport(headerOnly=True) + '\n'
        scanLog   += device.getScanReport() + '\n'
        #targetLog += device.getTargetReport() + '\n'
        deviceLog += device.getDeviceReport() + '\n'

    #remove targetLog here; put it in pipeReport also, 
    # and do it at the beginning, not the end, of processing
    #varsP.updateInfoReport(targetLog + '\n' + scanLog + '\n' + deviceLog + '\n')
    varsP.updateInfoReport(scanLog + '\n' + deviceLog + '\n')

    #return bnxFiles #instead of returning here, merge here, then return
    return joinBnxFiles(varsP, bnxFiles) #see return of this fn

#end performImageAnalysis    


#this is moved here from a method of varsPipeline--it fits better here
def joinBnxFiles(varsP, bnxFiles):
    """After image processing, merge results into all.bnx.
    """
    #the old way was to use this fn which simply copies lines
    # while this is fine most of the time, RefAligner is more sophisticated,
    # so it should be more robust to use RefAligner
    #molecule.joinBnxFiles(bnxFiles, self.bnxFile)            

    #this used to be called writeIntermediate
    varsP.writeListToFile(bnxFiles, varsP.bnxTarget)

    # args, jobName, expectedResultFile, uniqueString
    args = [varsP.RefAlignerBin, "-if", varsP.bnxTarget, "-merge", "-bnx", "-o", varsP.bnxFile.replace(".bnx",""), "-f"]
    if varsP.stdoutlog :
        args.extend( ['-stdout', '-stderr'] )
    #print "joinBnxFiles: args:", args

    jobwrapper = mthread.jobWrapper(varsP, "joinBnxFiles")
    jobwrapper.addJob( mthread.singleJob(args, "joinBnxFiles", varsP.bnxFile, "joinBnxFiles") )
    jobwrapper.multiThreadRunJobs(1)
    jobwrapper.doAllPipeReport()

    success = jobwrapper.allResultsFound()
    if not success :
        varsP.updatePipeReport("ERROR in performImageAnalysis: joinBnxFiles failed. Check: "+varsP.bnxTarget+"\n")

    # this is just putting the path of bnxFile in bnxTarget
    # on second thought, if I don't do this, then SampleCharModule will run on each bnx individually
    #if success :
    #    varsP.writeListToFile([varsP.bnxFile], varsP.bnxTarget)

    #sense of return of allResultsFound is opposite of return of performImageAnalysis:
    # allResultsFound is True for all jobs success, False for any jobs fail
    # performImage analysis return is 1 for failure
    return not success


class Scan():
    """Metadata structure for copying, reading, interpreting single Tiff file
    
    Runs detect, Maps lambda to reference, writes corrected bnx file
    """
    def __init__(self, varsP, curExp, molTag, remoteTiff, localTiff, quality=True):
        self.minMolLen = 20 #this used to be hardcoded at 100, now it's 20
        self.varsP = varsP
        self.curExp = curExp
        self.numLabelChannels = self.curExp.nColors - 1
        self.molTag = molTag
        self.remoteTiff = remoteTiff
        self.localTiff = localTiff
        self.detectJobs = None
        self.quality = quality
        
        # Logging Values
        self.Occupancy = 0
        self.Mb        = 0
        #init lambda data members
        self.lambdaErrFile = None
        self.lambdaNMapped = None
        self.lambdaBpp = 0

    #this is to facilitate logging
    def nameStr(self, short=True) :
        path = (os.path.split(self.remoteTiff)[1] if short else self.remoteTiff)
        return self.molTag + " " + path
        

    def getDetectJobs(self, contingentJob=None):
        self.molFile = self.localTiff + '.mol'
        #if the molFile is already there, the image processing is already done; do not repeat
        if os.path.exists(self.molFile): 
            return []
        #print "remoteTiff:", self.remoteTiff #debug
        #print "localTiff:", self.localTiff, os.path.exists(self.localTiff) #debug
        #print "expectedOverlap", self.curExp.ExpectedOverlap, "minOverlap", minOverlap #debug
        #there was an issue with self.curExp.ExpectedOverlap being incorrectly computed due to a bad value in the xml (see manageTargetLocation.py)
        #hopefully this will work
        if self.curExp.ExpectedOverlap > 100 :
            oldoverlap = self.curExp.ExpectedOverlap
            self.curExp.ExpectedOverlap = 15
            #print "Warning: calculated expectedOverlap", oldoverlap, "too large; defaulting to", self.curExp.ExpectedOverlap
            self.varsP.updateInfoReport("Warning: "+self.nameStr()+": calculated expectedOverlap %i too large; defaulting to %i\n" % (oldoverlap, self.curExp.ExpectedOverlap), printalso=True)
        expolap = (self.curExp.ExpectedOverlap - 10 if self.curExp.ExpectedOverlap >= 10 else 0) # ExpectedOverlap - 10
        minOverlap = '%d' % (expolap)
        maxOverlap = '%d' % (self.curExp.ExpectedOverlap + 10) # ExpectedOverlap + 10
        dmOverlapArgs = ['-o', minOverlap, '-O', maxOverlap]
        dmArgs = self.varsP.argsListed('imgDetection')
        dmArgs = util.argumentReplaceList(dmArgs, ['-x', str(self.curExp.ScanColumnCount)])
        dmArgs = util.argumentReplaceList(dmArgs, ['-y', str(self.curExp.ScanRowCount)])
        dmArgs = util.argumentReplaceList(dmArgs, ['-p', str(self.curExp.Pitch)])
        nchan = (self.curExp.nColors - 1 if self.curExp.nColors >= 2 else 1) #must be at least 1
        colorArgs = ['-n', str(nchan)]
        sJobCpName = 'cp ' + shorten(self.remoteTiff) + ' to ' + shorten(self.localTiff)
        #print "cp\n"+self.remoteTiff, "\n"+self.localTiff #debug
        sJobCp = mthread.singleJob(['cp', self.remoteTiff, self.localTiff], sJobCpName, self.localTiff, 'cpTiff', throttleClass = True)
        if contingentJob:
            sJobCp.addContingentJob(contingentJob)
        sJobDMName = 'Detect ' + shorten(self.localTiff)
        curArgs = [self.varsP.DMstaticBin] + dmOverlapArgs + dmArgs + colorArgs + [self.localTiff]
        argumentString = " ".join(curArgs) + '\n'
        print " ".join(curArgs) #debug
        sJobDM = mthread.singleJob(curArgs, sJobDMName, self.molFile, 'Detect')
        sJobDM.addContingentJob(sJobCp)
        sJobDM.bpp = self.curExp.basesPerPixel
        sJobDM.molTag = self.molTag
        sJobDM.numLabelChannels = self.numLabelChannels
        #inputMoleculesReport += '   ' + self.molTag + '  ' + self.remoteTiff + '\n'
        dorm = True #default True (False for debug)
        joblist = [sJobCp, sJobDM]
        if dorm :
            sJobRmImgName = 'Detect Complete, rm ' + shorten(self.localTiff)
            sJobRmImg = mthread.singleJob(['rm', self.localTiff], sJobRmImgName, '', 'rmFile')
            sJobRmImg.addContingentJob(sJobDM)
            joblist.append( sJobRmImg )
        return joblist
        
    def getLambdaMapJob(self, snrCutoff=0, verbose=False):
        #Note, verbose will print once per job, so use for debugging only
        # add lambda alignment band to this
        lambdaFilter = self.varsP.argData['lambdaFilter']
        lamMinLen = float(lambdaFilter[lambdaFilter.index('-minlen')  +1]) if '-minlen'   in lambdaFilter else 40.
        lamMaxLen = float(lambdaFilter[lambdaFilter.index('-maxlen')  +1]) if '-maxlen'   in lambdaFilter else 60.
        lamMinLab = float(lambdaFilter[lambdaFilter.index('-minsites')+1]) if '-minsites' in lambdaFilter else 6.
        lamMaxLab = float(lambdaFilter[lambdaFilter.index('-maxsites')+1]) if '-maxsites' in lambdaFilter else 10.
        #old format below (dict, not list)
        #lamMinLen = int(lambdaFilter['-minlen'  ]) # 40
        #lamMaxLen = int(lambdaFilter['-maxlen'  ]) # 60
        #lamMinLab = int(lambdaFilter['-minsites']) # 6
        #lamMaxLab = int(lambdaFilter['-maxsites']) # 10
        if verbose :
            self.varsP.updateInfoReport("lamMinLen = %.0f\n" % lamMinLen, printalso=True)
            self.varsP.updateInfoReport("lamMaxLen = %.0f\n" % lamMaxLen, printalso=True)
            self.varsP.updateInfoReport("lamMinLab = %.0f\n" % lamMinLab, printalso=True)
            self.varsP.updateInfoReport("lamMaxLab = %.0f\n" % lamMaxLab, printalso=True)
        
        #need mol file to do this; if doesn't exist, return with warning
        if not(os.path.exists(self.molFile)):
            print "Skipping map lambda job", self.molTag, "because mol file missing:", self.molFile
            self.lambdaErrFile = None
            return

        bnxFileLambda = '%s_lambda.bnx' % self.molTag
        bnxFileLambda = os.path.join(os.path.split(self.molFile)[0], bnxFileLambda)
        #if lambda bnx exists, skip the isolation step
        if os.path.exists(bnxFileLambda) :
            print "Using lambda bnx", bnxFileLambda
        else :
            print '  Isolating Lambda %s' % self.molTag
            lab2File = self.molFile.replace('.mol', '.0.lab')
            scanDset = molecule.moleculeDataset(self.curExp.basesPerPixel, molTag=int(self.molTag))
            scanDset.readMolFile(self.molFile)
            scanDset.annotateLabels(lab2File)
            # Introduce optArguments for Lambda Band
            scanDsetLambda = molecule.filteredSubset(scanDset,snrCutoff,lamMinLab,lamMaxLab,lamMinLen,lamMaxLen,True)
            scanDsetLambda.writeBnxFile(bnxFileLambda, quality=self.quality)

        self.lambdaBnx = bnxFileLambda
        baseArgs = self.varsP.argsListed('mapLambda')
        outputTarget = bnxFileLambda.replace('.bnx', '')
        curArgs = [self.varsP.RefAlignerBin, '-i', bnxFileLambda, '-o', outputTarget, '-ref', self.varsP.lambdaRef] + baseArgs
        if self.varsP.stdoutlog :
            curArgs.extend( ['-stdout', '-stderr'] )
        jobTag = self.molTag + '_lambda'
        self.lambdaErrFile = outputTarget + '.err'

        #if the err file exists, no need to process
        if os.path.exists(self.lambdaErrFile) :
            print "Skipping map lambda job ", jobTag, "because err file exists", self.lambdaErrFile
            return
        
        return mthread.singleJob(curArgs, jobTag, self.lambdaErrFile, jobTag)
        
    def getLambdaScaleFactor(self):
        self.lambdaBpp = 500
        bppList = []
        nMapsList = []
        nMappedMapsList = []
        #lambdaErrFile is set to None in getLambdaMapJob if mol file does not exist
        if not self.lambdaErrFile or not os.path.exists(self.lambdaErrFile) :
            return
        f1 = open(self.lambdaErrFile)
        errSkips = ['#', 'I']
        for line in f1 :
            if line[0] in errSkips :
                continue
            tokens = line.split()
            bppList.append(float(tokens[5]))
            nMapsList.append(int(tokens[7]))
            nMappedMapsList.append(int(tokens[9]))
        f1.close()
        self.lambdaBpp = bppList[-1]
        self.lambdaN = nMapsList[-2]
        self.lambdaNMapped = nMappedMapsList[-2]
        self.BppCorrection = self.lambdaBpp/self.curExp.basesPerPixel

    #option swapchannels is for testing only--do not use
    #verbose enables extra printing and logging in informatics report,
    # but this is redundant with logging in Device.findSNRCutoff and
    # Device.writeCorrectedBnx, which calls this, so False by default
    def getCorrectedDataset(self, bpp, snrCutoff, swapchannels=False, verbose=False):
        '''Load mol and lab files into molecule.moleculeDataset,
        call scaleTo500 and filteredSubset in order to apply SNR filter as previously calculated.'''
        scanDset = molecule.moleculeDataset(bpp, molTag=int(self.molTag),numLabelChannels=self.numLabelChannels)
        ret = scanDset.readMolFile(self.molFile)
        if ret != 0 :
            return
        suff = '.0.lab'
        if swapchannels :
            suff = '.1.lab'
        lab2File = self.molFile.replace('.mol', suff)
        scanDset.annotateLabels(lab2File)
        self.Occupancy = scanDset.CalculateOccupancy(1)
        self.Mb = scanDset.MegabasesDetected
        if verbose :
            self.varsP.updateInfoReport("   Prefilter : %s %s" % (self.molTag, scanDset.makeExperimentReport()), printalso=True)
        if self.numLabelChannels > 1:
            lab2File = self.molFile.replace('.mol', '.1.lab')
            scanDset.annotateLabels(lab2File, labelChannel=1)
        scanDset.scaleTo500()
        if verbose :
            self.varsP.updateInfoReport("   Postscale : %s %s" % (self.molTag, scanDset.makeExperimentReport()), printalso=True) #this one is ok
        filt = molecule.filteredSubset(scanDset,minLabSnr=snrCutoff,minLabels=5,minLen=self.minMolLen)
        if verbose :
            self.varsP.updateInfoReport("   Postfilter: %s %s" % (self.molTag, filt.makeExperimentReport()), printalso=True) 
        return filt
        
    def getScanReport(self, headerOnly=False):
        scanReport = '% 10s% 8s% 8s% 8s' % ('H','ID','totalMB','Occup.')
        if self.varsP.lambdaRef and self.lambdaNMapped :
            scanReport += '% 8s% 8s% 8s' % ('N_Lam', '%Map' ,'Bpp_Lam')
        if headerOnly:
            return scanReport
        scanReport = ''
        scanReport += '% 10s% 8s% 8.1f% 8.3f' % ('Scan',self.molTag,self.Mb,self.Occupancy)
        if self.varsP.lambdaRef and self.lambdaNMapped :
            lambdarat = (float(self.lambdaNMapped)/self.lambdaN if self.lambdaN > 0 else 0)
            scanReport += '% 8d% 8.3f% 8.2f' % (self.lambdaNMapped, lambdarat, self.lambdaBpp)
        return scanReport
        
        
class Device():
    """Metadata for handling device and instrument specific information
    
    Performs lambda scaling, single scan interpolation
    Class is author of device bnx file
    """
    
    def __init__(self, varsP, expTag, remotePath, localPath, bypass, quality=True, forceonecolor=False):
        self.varsP = varsP
        self.nScans = 0
        self.expTag = expTag
        self.remotePath = remotePath
        self.localPath = localPath
        self.bnxFile = os.path.join(localPath, expTag + '.bnx')
        self.jobList = []
        #self.lambdaDiffProfiles = []
        self.sizeProfile = None
        #self.lambdaBandProfiles = []
        self.lambdaBnxFiles = []
        self.bnxFiles = []
        self.messages = ''
        self.SnrCutoff = [0,0] #list of size 2--one per color channel
        self.bypass = bypass
        #self.scanReportHeader = ''
        #self.scanReport = ''
        #self.deviceReportHeader = ''
        #self.deviceReport = ''
        self.scans = []
        self.scanNames = []
        self.curExp = mtl.repositoryFolder(remotePath, localPath, leadTiffTag=expTag)
        if forceonecolor :
            self.curExp.nColors = 2 #2 is 1. See next line.
        self.numLabelChannels = (self.curExp.nColors - 1 if self.curExp.nColors > 0 else 1)
        self.quality = quality #must be before setupScans
        self.setupScans()
        self.bpp = [500] * self.nScans
        self.bppCorrection = [1.0] * self.nScans
        self.lambdaBpp = []
        self.lambdaSnrCutoff = None #to skip adding a map job in self.getLambdaMapJobs
        self.lambdaNMapped = None
        
        # Logging Values
        self.labPer100kb            = [0,0] #list of size 2 ??
        self.roughSizeProfileHeader = ""
        self.roughSizeProfile       = 0
        self.roughMassProfileHeader = ""
        self.roughMassProfile       = 0
        
    def setupScans(self):
        """Create Scan object for each tiff file in the path for this object.
        Log the expTag, which corresponds to the output directory and the tiff files
        which correspond to the indices in this directory in both the pipelineReport
        and informaticsReport.
        """
        self.varsP.updatePipeReport("Loading scans for experiment %s\n" % self.expTag) 
        self.varsP.updateInfoReport("Loading scans for experiment %s\n" % self.expTag) #printed already
        for i,imagePair in enumerate(zip(self.curExp.remoteTiffList, self.curExp.localTiffList)):
            scanTag = '%02d' % (i+1)
            molTag = self.expTag + scanTag
            remoteTiff = imagePair[0]
            localTiff = imagePair[1]
            curScan = Scan(self.varsP, self.curExp, molTag, remoteTiff, localTiff, quality=self.quality)
            self.scans.append(curScan)
            self.scanNames.append(molTag)
            self.nScans += 1
            self.varsP.updatePipeReport(curScan.nameStr(short=False)+"\n") #log which paths correspond to which tiff files
            self.varsP.updateInfoReport(curScan.nameStr(short=False)+"\n")

    def getLambdaMapJobs(self):
        lambdaMapJobs = []
        for scan in self.scans:
            if not self.lambdaSnrCutoff :
                print "Skipping scan because no lambdaSnrCutoff:", scan.molTag
                continue
            addjob = scan.getLambdaMapJob(self.lambdaSnrCutoff)
            if addjob : #put the warning message inside getLambdaMapJob
                lambdaMapJobs.append(addjob)
                
        return lambdaMapJobs

    def processLambdaMapResults(self):
        lambdaFilter = self.varsP.argData['lambdaFilter']
        minNMapped = float(lambdaFilter[lambdaFilter.index('-mincount')+1]) if '-mincount' in lambdaFilter else 10.
        #old format is dict (below), new is list (above)
        #minNMapped = int(lambdaFilter['-mincount']) # 25
        prevBpp = 500
        prevCorrection = 1.0
        for i,scan in enumerate(self.scans):
            scan.getLambdaScaleFactor()
            if scan.lambdaNMapped >= minNMapped: # criteria of previous map rate
                prevBpp = scan.lambdaBpp
                prevCorrection = scan.BppCorrection
            self.bpp[i] = prevBpp
            self.bppCorrection[i] = prevCorrection
        
    def writeCorrectedBnx(self):
        self.varsP.updateInfoReport("  Correcting %s\n" % self.expTag, printalso=True)
        deviceDset = molecule.moleculeDataset(self.curExp.basesPerPixel,numLabelChannels=self.numLabelChannels)
        self.varsP.updateInfoReport("               ID  "+deviceDset.makeExperimentHeader(), printalso=True)
        goodscans = 0
        for i,scan in enumerate(self.scans):
            #print scan.localTiff #debug
            #print self.bppCorrection[i] #debug
            #if self.bppCorrection != 1.0 :
            correctedBpp = self.bppCorrection[i]*self.curExp.basesPerPixel
            curDataset = scan.getCorrectedDataset(correctedBpp, self.SnrCutoff)
            if not curDataset :
                continue
            goodscans += 1
            #this is redundant with the info report at the end of performImageAnalysis
            # partially, since that comes from the scan object, this is from the moleculeDataset
            self.varsP.updateInfoReport('   Corrected : %s %s' % (scan.molTag, curDataset.makeExperimentReport()), printalso=True)
            #self.varsP.updateInfoReport('\n', printalso=True) #debug
            deviceDset.addDset(curDataset)
        if not goodscans :
            self.bnxFile = "" #if no good scans, null the bnxFile
            return
        self.labPer100kb = deviceDset.calculateLabPer100Kb()
        deviceDset.writeBnxFile(self.bnxFile, quality=self.quality)
        #self.varsP.updateInfoReport("\n", printalso=True)

    #how to treat re-running:
    #goal is to avoid repeating processing of already processed Devices
    #if there are new scans, you need to keep those in allJobs
    #if there are no new scans, you need to do nothing, so return empty list
    #I don't see why you'd ever want to remove and remake the dir, so disable that by default
    def getTargetJobs(self, dormdir=False):
        localDataLocation = os.path.join(self.varsP.localRoot, self.expTag + '/')
        #print "localDataLocation:", localDataLocation #debug
        if dormdir :
            sJobRmName = 'Pre-Remove Folder: ' + shorten(localDataLocation)
            sJobRm = mthread.singleJob(['rm', '-f', '-r', localDataLocation], sJobRmName, '', 'rmDir')
            sJobMkdirName = 'Make Folder: ' + shorten(localDataLocation)
            sJobMkdir = mthread.singleJob(['mkdir', localDataLocation], sJobMkdirName, localDataLocation, 'mkDir')
            sJobMkdir.addContingentJob(sJobRm)
            allJobs = [sJobRm, sJobMkdir]
            contingentjob = sJobMkdir
        else :
            util.checkDir(localDataLocation) #will make dir localDataLocation
            allJobs = []
            contingentjob = None
        for scan in self.scans:
            scanjobs = scan.getDetectJobs(contingentjob)
            if not scanjobs : #no scan jobs means the scan has already been processed--clear all jobs
                self.varsP.updatePipeReport("Device.getTargetJobs: skipping path "+scan.nameStr()+"\n") #localDataLocation
            else :
                allJobs += scanjobs
        return allJobs
        
    def findSNRCutoff(self):
        self.varsP.updateInfoReport( '  Finding SNR cutoff %s\n' % self.expTag, printalso=True ) #put this in info report
        usedScanCount = 0
        #print "\nnscans = ", len(self.scans), "\n" #debug
        deviceDset = molecule.moleculeDataset(self.curExp.basesPerPixel,numLabelChannels=self.numLabelChannels)
        self.varsP.updateInfoReport("               ID  "+deviceDset.makeExperimentHeader(), printalso=True)
        for i, scan in enumerate(self.scans):
            #if i == 0:
                #self.scanReportHeader += deviceDset.makeExperimentHeader() +'\n'
                #if self.varsP.lambdaRef:
                #    deviceDsetLambda = molecule.moleculeDataset(self.curExp.basesPerPixel,numLabelChannels=1)
            molFile = scan.molFile
            if not os.path.exists(molFile) :
                self.varsP.updatePipeReport("Device.findSNRCutoff: MISSING FILE: " + molFile + "\n") 
                continue
            usedScanCount += 1
            scanDset = molecule.moleculeDataset(self.curExp.basesPerPixel, molTag=int(scan.molTag),numLabelChannels=self.numLabelChannels)
            scanDset.readMolFile(molFile)
            lab2File = molFile.replace('.mol', '.0.lab')
            scanDset.annotateLabels(lab2File)
            if self.numLabelChannels > 1:
                lab2File = molFile.replace('.mol', '.1.lab')
                scanDset.annotateLabels(lab2File, labelChannel=1)
            self.varsP.updateInfoReport( "   SNR cutoff: " + scan.molTag + "  " + scanDset.makeExperimentReport(), printalso=True)
            deviceDset.addDset(scanDset)
        if not usedScanCount :
            return
        self.roughMassProfileHeader = deviceDset.getRoughMassProfile(headerOnly=True)
        self.roughMassProfile = deviceDset.getRoughMassProfile()
        deviceDset.getLogSnrCutoff()
        self.SnrCutoff = deviceDset.LogSnrSimpleThresh
        self.lambdaSnrCutoff = self.SnrCutoff
        self.varsP.updateInfoReport( "   SNR total :       " + deviceDset.makeExperimentReport(), printalso=True )
        snrstr = " ".join( map(lambda x: "%.3f"%x, self.SnrCutoff) )
        self.varsP.updateInfoReport( "   SNR Cutoff: " + snrstr + "\n", printalso=True )
        #self.varsP.updateInfoReport( "   SNR Cutoff: " + str(self.SnrCutoff) + "\n", printalso=True )
    
    def getTargetReport(self,headerOnly=False):
        targetReport = '% 8s% 10s' % ('ExpID', 'Location')
        if headerOnly:
            return targetReport
        targetReport = '% 8s% 100s' % (self.expTag, self.remotePath)
        return targetReport
        
    def getScanReport(self, header=True):
        scanReport = ''
        self.Mb = 0
        self.nLabels = 0
        for i,scan in enumerate(self.scans):
            if i==0 and header:
                scanReport += scan.getScanReport(headerOnly=True) + '\n'
            scanReport += scan.getScanReport() + '\n'
        return scanReport
        
    def getDeviceReport(self, headerOnly=False):
        if headerOnly:
            if self.numLabelChannels > 1:
                deviceReport =  '% 5s % 8s % 8s % 8s % 8s % 8s' % ('Exp', 'nScans', 'SNR-CutA', 'SNR-CutB', 'Lab/100A', 'Lab/100B')
            else:
                deviceReport =  '% 5s % 8s % 8s % 8s' % ('Exp', 'nScans', 'SNR-Cut', 'Lab/100')
            deviceReport += self.roughMassProfileHeader
            return deviceReport

        #deviceReport = self.deviceReport #this is no longer used
        labelDens = self.labPer100kb
        if self.numLabelChannels > 1:
            deviceReport = '  % 3s % 8s % 8.1f % 8.1f % 8.1f % 8.1f % 8s' % (self.expTag, self.nScans, self.Mb, self.SnrCutoff[0], self.SnrCutoff[1], labelDens[0], labelDens[1]) #, self.roughMassProfile)
        else:
            deviceReport = '  % 3s % 8s % 8.1f % 8.1f % 8s' % (self.expTag, self.nScans, self.SnrCutoff[0], labelDens[0], self.roughMassProfile)
        return deviceReport
        

def splitColorChannelsBnx(bnxFile, color = 1):
    colorTag = str(color)
    newBnxFile = bnxFile.replace('.bnx', '_%s.bnx' % colorTag)
    f1 = open(bnxFile)
    f2 = open(newBnxFile, 'w')
    lineTags = ['#','0']
    while(True):
        line = f1.readline()
        if line == '':
            f1.close()
            f2.close()
            break
        if lineTags.__contains__(line[0]):
            f2.write(line)
            continue
        if line[0] == colorTag:
            line = '1' + line[1:]
            f2.write(line)
    
                
def initializeInfoFiles(varsP):
    f1 = open(varsP.extraInfoFile1, 'w')
    f1.close()
    f1 = open(varsP.extraInfoFile2, 'w')
    f1.close()
    
def addDataToFile(target, lead, data):
    f1 = open(target, 'a')
    wdat = '%04d' % lead
    for val in data:
        wdat += '\t%d' % val
    wdat += '\n'
    f1.write(wdat)
    f1.close()
    
    
def parseExperimentFile(targetFile):
    targetLocations = []
    for line in open(targetFile) : #this should have been checked in varsPipeline.checkInputs()
        if line[0] == '#':
            continue
        targetLocation = line.strip().replace("\\", "")
        #print targetLocation #debug
        if os.path.exists(targetLocation):
            targetLocations.append(targetLocation)
        else:
            print 'Location Not Found: %s' % targetLocation
            print 'Omitting...'
    #check that list is unique
    if not sorted(targetLocations) == sorted(set(targetLocations)) :
        print "ERROR in ImageProcessingModule.parseExperimentFile: duplicate entries in image paths file:", targetFile
        return [] #this will cause exit
    return targetLocations
        

def shorten(inName, maxLen = 20):
    if len(inName) >= maxLen + 2:
        return inName[:maxLen/2] + '..' + inName[-maxLen/2:]
    return inName


