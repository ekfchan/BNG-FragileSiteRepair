import os
#import pdb
import sys
from random import choice
import math
from collections import defaultdict

"""@package molecule Author of BNX file; aggregates detection results

"""


import utilities
utilities.setVersion("$Id: molecule.py 2494 2014-02-15 00:18:25Z wandrews $")


tab = '\t'
newline = '\n'

def lambdaAnchors(bpp = 500):
    """In silico lambda bspqI reference (hard coded)
    """
    lambdaRef = [2396,6488,8701,10369,13285,24768,27233,34326,34799,47711]
    distPattern = []
    for i, val0 in enumerate(lambdaRef):
        for j, val1 in enumerate(lambdaRef):
            if i >= j:
                continue
            distPattern.append(abs(val0-val1))
    scaleFactor = 1./bpp
    return [scaleFactor*x for x in distPattern]
    
    
def subsampleBnx(varsP, sizeLowerLimit = 150., sortKeys = False, bypass = False):
    """Reduce size of bnx file by random subselection to prescribed content
    """
    suffix = '_sub%dMb' % varsP.subsampleMb
    newBnx = varsP.bnxFile.replace('.bnx', suffix + '.bnx')
    if bypass:
        varsP.bnxFile = newBnx
        return 0
    bpp = getBppFromBnx(varsP.bnxFile)
    dset_orig = moleculeDataset(bpp, PitchNm = 700)
    dset_orig.readBnxFile(varsP.bnxFile)
    if dset_orig.MegabasesDetected < varsP.subsampleMb:
        print '  Bnx File %s has less than targetMB % 6.1f' % (varsP.bnxFile, varsP.subsampleMb)
        return 1
    dset_subsample = moleculeDataset(bpp, PitchNm = 700)
    allMols = dset_orig.molDB.keys()
    usedMols = []
    totalContent = 0
    while True :
        molID = choice(allMols)
        mol = dset_orig.molDB[molID]
        if mol.lengthKb >= sizeLowerLimit :
            totalContent += 0.001 * mol.lengthKb
        dset_subsample.addMol(mol)
        usedMols.append(molID)
        allMols.remove(molID)
        if totalContent >= varsP.subsampleMb:
            print '  Subsampled %3.1fMb' % totalContent
            break
    dset_subsample.writeBnxFile(newBnx)
    varsP.bnxFile = newBnx
    return 0


def convertMolLabToBnx(molFile, labFile, bpp, labelSNR = 3.0, moltag = '') : #, sizeLowerLimit = 150., sizeUpperLimit = 100000.):
    dataset = moleculeDataset(bpp, molTag = moltag, PitchNm = 700)
    dataset.readMolFile(molFile)
    dataset.annotateLabels(labFile, labelSNR = labelSNR)
    dataset.scaleTo500()
    bnxFile = molFile.replace('.tiff.mol','.bnx')
    dataset.writeBnxFile(bnxFile)
    return bnxFile, dataset    
    

def convertDetectOutToBnxNew(molFile, labFile="", bpp=500., labelSNR = 3.0, moltag = 0, numLabelChannels = 1, sizeLowerLimit = 0., sizeUpperLimit = 0.):
    """Top level method to convert .mol + .lab(s) to .bnx file format
    """
    NanoStudio = True
    #moleculeDataset bpp,molTag = 0, NanoStudio = False, FovSizeMicrons = 80., PitchNm = 500., numLabelChannels = 1):
    scanDset = moleculeDataset(bpp, moltag, NanoStudio, numLabelChannels = numLabelChannels)
    scanDset.readMolFile(molFile)
    for i in range(numLabelChannels):
        if not labFile :
            labFile = molFile.replace('.mol', '.' + str(i) + '.lab')
        if os.path.exists(labFile):
            scanDset.annotateLabels(labFile, labelSNR = labelSNR, labelChannel = i)
        else :
            print "Missing lab file:", labFile
            return
    deviceDsetFilt = filteredSubset(scanDset,[labelSNR,0],2,200,sizeLowerLimit,sizeUpperLimit,False)
    bnxFile = molFile.replace('.mol','.bnx')
    deviceDsetFilt.writeBnxFile(bnxFile)
    print(' Writing %s' % bnxFile)
    
    
def convertDetectOutToBnx(molFile, bpp, labelSNR = 3.0, moltag = '', numLabelChannels = 1, sizeLowerLimit = 150., sizeUpperLimit = 100000.):
    dataset = moleculeDataset(bpp, molTag = moltag, PitchNm = 700, numLabelChannels = numLabelChannels)
    dataset.readMolFile(molFile)
    curLabChannel = 0
    for i in range(numLabelChannels):
        labFile = molFile.replace('.mol', '.' + str(i) + '.lab')
        if os.path.exists(labFile):
            dataset.annotateLabels(labFile, labelSNR = labelSNR, labelChannel = i)
    dataset.scaleTo500()
    if dataset.nLabelsDetected[0] > 1e5:
        dataset.getLogSnrCutoff()
        dataset.labelSNRcutoff = self.labelSNRcutoff
    else:
        dataset.labelSNRcutoff = 3.5
    bnxFile = molFile.replace('.tiff.mol','.bnx')
    dataset.writeBnxFile(bnxFile)
    return bnxFile, dataset
    
    
def joinBnxFiles(bnxList, newbnxTarget):
    """Aggregate BNX files, typically from a number of IRYS scans or devices
    """
    lineCount = 0
    bnxTemp = newbnxTarget + '_temp'
    fout = open(bnxTemp, 'w')
    searchSeq = 'N Fragments:'
    lnCount = 0
    ct = 0
    header = ''
    for bnxFile in bnxList:
        f1 = open(bnxFile)
        while(True):
            line = f1.readline()
            if line == '':
                f1.close()
                ct += 1
                break
            if line[0] == '#':
                if ct == 0:
                    header += line
                    continue
                continue
            lnCount += 1
            fout.write(line)
    fout.close()
    
    fin = open(bnxTemp, 'r')
    fout = open(newbnxTarget, 'w')
    fout.write(header)
    while(True):
        line = fin.readline()
        if line == '':
            fin.close()
            fout.close()
            break
        fout.write(line)


#maxLabels and maxLen has no effect if <= 0
def filteredSubset(dSet, minLabSnr=[-100,-100], minLabels=0, maxLabels=0, minLen=0, maxLen=0, rejectStitch=False):
    """Filter Molecule/Label set by: Length, Label Count, Label SNR, stitched
    """
    newDset = moleculeDataset(dSet.bpp, numLabelChannels=dSet.numLabelChannels)
    for mol in dSet.molDB.values():
        if mol.isStitched and rejectStitch:
            continue
        if mol.lengthKb < minLen or (maxLen > 0 and mol.lengthKb > maxLen) :
            continue
        nLabels = 0
        labList = []
        colorList = []
        for i,labSet in enumerate(mol.labels):
            for lab in labSet.values():
                if lab.snr >= minLabSnr[i]:
                    nLabels += 1
                    labList.append(lab)
                    colorList.append(i)
        if nLabels < minLabels or (maxLabels > 0 and nLabels > maxLabels):
            continue
        newMol = molecule(mol.ID, mol.lengthPx, mol.lengthKb,0,0,0,0,0, mol.intensity, mol.isStitched,nLabelChannels=dSet.numLabelChannels)
        for i,lab in enumerate(labList):
            newMol.add_label(lab, colorList[i])
        newDset.addMol(newMol)
    return newDset


class moleculeDataset():
    """This class contains the detection and alignment data from nanostudio"""
    #PitchNm is the pitch (distance between channels) in nm
    def __init__(self, bpp,molTag = 0, NanoStudio = False, FovSizeMicrons = 80., PitchNm = 700., numLabelChannels = 1):
        self.bpp = bpp
        self.originalBpp = bpp
        self.molDB = {}
        self.mappedMolecules = set([])
        self.nMolecules = 0
        self.MegabasesDetected = 0.
        self.MegabasesDetectedPreFilter = 0.
        self.PitchNm = PitchNm
        NumChannels = (FovSizeMicrons * 1000) / PitchNm
        MbPerChannel = (512. * self.bpp) / 1e6
        self.FOVCapacityMb = MbPerChannel * NumChannels
        self.Occupancy = 0.
        if numLabelChannels < 1 : #don't want zero size list
            numLabelChannels = 1
        self.nLabelsDetected = [0] * numLabelChannels
        self.nLabelsDetectedPreFilter = 0
        self.MegabasesMapInput = 0
        self.NanoStudio = NanoStudio
        self.curMolTag = molTag
        self.numLabelChannels = numLabelChannels
        self.labelSNRcutoff = 3.5
        self.LogSnrSimpleThresh = [0] * numLabelChannels
        self.messages = ''
        self.molAvgSnr = -1
    
    def calculateLabPer100Kb(self):
        dat = []
        for i in range(self.numLabelChannels):
            if self.MegabasesDetected:
                dat.append(self.nLabelsDetected[i] / (self.MegabasesDetected * 10))
            else:
                dat.append(0)
        return dat
    
    def diffAnalysis(self,binSize = 0.1):
        nBins = 500
        diffPattern = [0] * nBins
        scaleFactor = 1/binSize
        for mol in self.molDB.values():
            diffSet = mol.getDiffSet(0)
            for val in diffSet:
                bini = int(scaleFactor * val)
                if bini < nBins:
                    diffPattern[bini] += 1
        return diffPattern
        
    def addDset(self, dSet):
        """Add contents of another instance of this class
        """
        self.molDB.update(dSet.molDB)
        self.nMolecules += dSet.nMolecules
        self.MegabasesDetected += dSet.MegabasesDetected
        for i in range(min(self.numLabelChannels,dSet.numLabelChannels)):
                self.nLabelsDetected[i] += dSet.nLabelsDetected[i] 
            
    def getRoughMassProfile(self, headerOnly=False):
        """Creates course size distribution output for report
        """
        bins = [0,50,100,150,200,250]
        header = '% 6s' % 'Mb'
        for val in bins:
            header += '% 7dkb' % val
        if headerOnly:
            return header
        mass = [0.] * 6
        totalMass = 0.
        for mol in self.molDB.values():
            bin = int(mol.lengthKb) / 50
            bin = min(5,bin)
            mass[bin] += mol.lengthKb
            totalMass += mol.lengthKb
        normMass = [x/totalMass for x in mass]
        data = '% 8.1f' % (totalMass * 0.001)
        for val in normMass:
            data += '% 8.1f%%' % (100*val)
        return data
    
    def getSizeProfile(self, nbins=500):
        """Get fine molecular length distribution
        """
        bins = range(nbins)
        sprofile = [0] * nbins
        for mol in self.molDB.values():
            bin = int(mol.lengthKb)
            if bin >= nbins:
                continue
            sprofile[bin] += 1
        return sprofile
    
    
    def makeMolID(self, idFromFile):
        idInt = int(idFromFile)
        molID = '%04d%06d' % (self.curMolTag, idInt)
        return int(molID)

    
    def getLogSnrCutoff(self, snrRange=[3.,7.], defaultSnr=4., scaleFactor=10.):
        """Calculate SNR cutoff by log min minus one method
        Take log(snr * scaleFactor), and find minimum between bimodal
        distribution. In case the number of peaks is not equal to two,
        require the snr to be within snrRange. If not, use defaultSnr.
        """
        verbose = False #True #for get{Simple,TwoPeak}Minimum only, debugging
        ind = lambda x : int(round(math.log(x) * scaleFactor))
        inverseind = lambda x : math.e**(x/scaleFactor)
        #self.xLogSnrProfile = range(-7,45) #this seems to not be used -- deprecate
        for ii in range(self.numLabelChannels):
            #xDict = {} #not used
            freq = defaultdict(int)
            for molID in self.molDB.keys():
                mol = self.molDB[molID]
                for lab in mol.labels[ii].values():
                    try:
                        i = ind(lab.snr)
                    except:
                        continue
                    #if xDict.has_key(i):
                    freq[i] += 1 #use defaultdict so don't need to call has_key
                    #else:
                        #xDict[i] = math.e**(i/scaleFactor)
                    #freq[i] = 1
            keySort = sorted(freq.keys())
            #print "\nfreq =", freq, "\nkeySort =", keySort #debug
            #normalization isn't necessary
            #cum = float(sum(freq.values()))
            #valSort = [freq[x]/cum for x in keySort]
            valSort = map(freq.get, keySort) #[freq[x] for x in keySort]
            #print "valSort =", valSort, "\n" #debug
            #self.LogSnrProfile = []
            if keySort == []:
                return
            #indMaxOtsu = getOtsuThreshold(keySort,valSort)
            simplemin  = inverseind(keySort[getSimpleMinimum(valSort)])
            twopeakmin = inverseind(keySort[getTwoPeakMinimum(valSort, verbose=verbose)])
            #print "Simple min:", simplemin #debug
            #print "Twopeak min:", twopeakmin #debug
            if twopeakmin > snrRange[0] and twopeakmin < snrRange[1] :
                self.LogSnrSimpleThresh[ii] = twopeakmin
            elif simplemin > snrRange[0] and simplemin < snrRange[1] :
                print "Warning: getLogSnrCutoff:", twopeakmin, "outside snrRange of", snrRange, "using old algo"
                self.LogSnrSimpleThresh[ii] = simplemin
            else :
                #self.varsP.updateInfoReport #don't have varsP
                print "Warning: getLogSnrCutoff:", twopeakmin, "and", simplemin, "outside snrRange of", snrRange, "Defaulting to", defaultSnr
                self.LogSnrSimpleThresh[ii] = defaultSnr

            #what is this block for??? -- deprecate
            #for val in self.xLogSnrProfile:
            #    if freq.has_key(val):
            #        self.LogSnrProfile.append(freq[val]/cum)
            #    else:
            #        self.LogSnrProfile.append(0)

        
    def CalculateOccupancy(self, NumScans):
        if self.FOVCapacityMb == 0 : return 0 #just in case
        return self.MegabasesDetected / (self.FOVCapacityMb * 104)
    

    #simply average all molecule snr entries for all molecules in molDB
    #just take zeroth label channel
    def calculateMolAvgSnr(self, recalculate=False) :
        if not self.molDB : #need to prevent division by 0 for no entries in molDB
            return
        if not recalculate and self.molAvgSnr != -1 : #default -1: this means it has already been calculated
            return
        #need to call molecule.getAvgSnr
        avgsnr = 0
        for m in self.molDB.values() :
            m.getAvgSnr()
            avgsnr += m.avgSnr[0]
        avgsnr /= len(self.molDB)
        self.molAvgSnr = avgsnr


    #call molecules.getNLabels for each molecule in molDB
    def getTotNLabels(self) :
        tnl = [0]*self.numLabelChannels
        for mol in self.molDB.values() :
            assert mol.nLabelChannels <= self.numLabelChannels #consistency check, prevent index out of bounds on tnl
            nlab = mol.getNLabels()
            for i,ele in enumerate(nlab) :
                tnl[i] += ele
        return tnl


    def makeExperimentHeader(self):
        #output = '% 10s% 10s% 10s% 10s% 10s% 10s' %('H', 'ID', 'totalMB', 'usedMB', 'Lab/100kb', 'Occup.')
        #output = '% 10s% 10s% 10s% 10s' %('H', 'ID', 'totalMB', 'Occup.')
        #assume 2 labels
        output = "%6s  %10s  %5s  %5s  %7s  %7s  %4s\n" % ("Nmol", "Mb", "Occ", "bpp", "nLab1", "nLab2", "avgSnr")
        return output

    
    def makeExperimentReport(self):
        assert len(self.molDB) == self.nMolecules, ("Invalid n molecules: %i %i" % (len(self.molDB), self.nMolecules))
        #if self.molAvgSnr == -1 :
        self.calculateMolAvgSnr() #this is not called elsewhere--call always just to be safe
        nlab = self.getTotNLabels()
        #output =  '%i'       % len(self.molDB) #redundant with nMolecules
        output =  '%6i'      % self.nMolecules 
        output += '  %10.2f' % self.MegabasesDetected
        output += '  %5.3f' % self.CalculateOccupancy(1)
        output += '  %5.1f'  % self.bpp
        output += '  %7i'    % nlab[0]
        output += '  %7i'    % (nlab[1] if len(nlab) > 1 else 0)
        output += '  %4.2f'  % self.molAvgSnr
        #output += '  %5.1f'  % self.originalBpp
        #output += '%i'       % self.nLabelsDetected #this is a list, not sure why

        #output +=  '% 10.2f' % self.MegabasesMapInput
        # if self.MegabasesMapInput == 0.:
            # output += '% 10s' % 'N/A'
        # else:
            # output += '% 10.3f' % (0.1 *(self.nLabelsMapInput / self.MegabasesMapInput))
        output += "\n"
        return output
    

    def readMolFile(self, molFile, verbose=0):
        try :
            handle = open(molFile)
        except IOError, e :
            if verbose :
                print "ERROR in readMolFile: file not found:\n", e
            return 1 #a-la bash exit codes
        while(True):
            line = handle.readline()
            if line == '':
                handle.close()
                break
            if line[0] == 'm' or line[0] == 'M' or line[0] == 'S' or line[0] == '#':
                continue
            self.addMolecule(line)
        return 0 #a-la bash exit codes
    

    def addMolecule(self,line):
        tokens = line.split('\t')
        molID = self.makeMolID(int(tokens[0]))
        lengthPx = float(tokens[12])
        lengthKb = 0.001 * lengthPx * self.bpp
        intensity = float(tokens[15])
        isStitched = tokens[4] != tokens[5]
        mol = molecule(molID, lengthPx, lengthKb, 0, 0, 0, 0, 0, intensity, isStitched, nLabelChannels = self.numLabelChannels)
        self.addMol(mol)

        
    def addMol(self,mol):
        self.MegabasesDetected += mol.lengthKb * 0.001
        assert not self.molDB.has_key(mol.ID), "Molecule ID already in molDB: "+str(mol.ID)
        self.molDB[mol.ID] = mol   
        self.nMolecules += 1
        for i in range(self.numLabelChannels):
            self.nLabelsDetected[i] += mol.nLabels[i]
    

    def scaleTo500(self):
        k = self.bpp / 500.
        for mol in self.molDB.values():
            #mol = self.molDB[molID]
            mol.lengthPx *= k
            for i,labDict in enumerate(mol.labels):
                newLabels = {}
                for key in mol.labels[i].keys():
                    lab = mol.labels[i][key]
                    lab.relypos *= k
                    newLabels[lab.relypos] = lab
                mol.labels[i] = newLabels
        self.bpp = 500


    def annotateLabels(self, labFile, labelChannel = 0, labelSNR = -1):

        handle = open(labFile)
        newMasterDict = {}
        while(True):
            rawInput = handle.readline()
            if rawInput == '':
                handle.close()
                break
            ri = rawInput[0]
            if ri == 'l' or ri == '#' or ri == 'L' or ri == 'S':
                continue
            lab = label()
            lab.fromLabFileLine(rawInput, NanoStudio = self.NanoStudio)
            lab.molID = self.makeMolID(lab.molID)
            mol = self.molDB[lab.molID]
            if lab.snr >= labelSNR:
                #if mol.nLabelChannels != 1 : #debug
                #    print "nLabelChannels:", mol.nLabelChannels, "labelChannel:", labelChannel #debug
                mol.add_label(lab, labelChannel)
                self.nLabelsDetected[labelChannel] += 1
    

    def readBnxFile(self, bnxFile, mbOnly=False):
        # need to add color reading capabilities here
        currentBpp = float(self.bpp)
        f1 = open(bnxFile)
        while(True):
            line = f1.readline()
            if line == '':
                f1.close()
                break
            if line[0] == '#':
                ind = line.find('Bases per Pixel:')
                if not(ind == -1):
                    tokens = line.split('\t')
                    bppLocal = float(tokens[1].strip())
                    scaleFactor = self.bpp / bppLocal
                ind = line.find('Label Channels:')
                if not(ind == -1):
                    tokens = line.split('\t')
                    numLabelChannels = int(tokens[1])
                    self.numLabelChannels = numLabelChannels
            if line[0] == '0':
                tokens = line.split('\t')
                molID = int(tokens[1])
                lengthPx = float(tokens[2]) * scaleFactor / currentBpp
                lengthKb = float(tokens[2]) / 1000
                if mbOnly:
                    if lengthKb > 150:
                        self.MegabasesDetected += 0.001 * lengthKb
                    continue
                mol = molecule(molID, lengthPx, lengthKb, 0, 0, 0, 0, 0, 0, False, nLabelChannels = self.numLabelChannels)
                self.addMol(mol)
                continue
            if mbOnly:
                continue
            try:
                colorChannel = int(line[0]) - 1
            except:
                continue
            tokens = line.split('\t')
            ct = 1
            for val in tokens[1:]:
                labelPosBp = float(val)
                labelPosPx = labelPosBp * scaleFactor * (1/currentBpp)
                lab = label()
                lab.fromRawEntry(ct, labelPosPx, 10, 0, self.makeMolID(molID))
                mol.add_label(lab, colorChannel)
                
    
    def writeBnxFile(self, bnxFile, sizeLowerLimit = 0., numLabelsLowerLimit = 0, molIntensityUpperLim = 0.8, sortKeys = False, quality=False): #, sizeUpperLimit = 100000.
        self.MegabasesMapInput = 0.
        bnxData = ''
        intensityCutoff = molIntensityUpperLim * 2**14
        molWriteCount = 0
        lineCt = 0

        f1 = open(bnxFile, 'w')
        #makeBnxHeader(bpp, nLabelChannels = 1, minmollen = 0, labsnr = 0, quality=False):
        header = makeBnxHeader(self.bpp, nLabelChannels=self.numLabelChannels, minmollen = sizeLowerLimit, labsnr=self.labelSNRcutoff, quality=quality)
        f1.write(header)
        f1.close()
        
        if sortKeys:
            keys = sorted(self.molDB.keys())
        else:
            keys = self.molDB.keys()
        for molID in keys:
            mol = self.molDB[molID]
            
            bnxMolID = '%06d' % molID
            molEntry = mol.makeBnxEntry(bnxMolID, self.bpp, quality=quality)
            self.MegabasesMapInput += mol.lengthKb * 0.001
            bnxData += molEntry
            lineCt += 1
            # why so many open/close? I guess for buffer? Change 100 to 1000
            if not(lineCt%1000):
                f1 = open(bnxFile, 'a')
                f1.write(bnxData)
                f1.close()
                bnxData = ''
        f1 = open(bnxFile, 'a')
        f1.write(bnxData)
        f1.close()
            

    def annotateMapResults(self, mapFile):
        handle = open(mapFile)
        mapEntries = -1
        while(True):
            rawInput = handle.readline()
            if rawInput == '':
                break
            ri = rawInput[0]
            if ri == '#' or ri == 'M' or ri == 'S':
                continue
            mapEntries = 1
            detabString = [x for x in rawInput.strip().split('\t')]
            mapMolID = int(detabString[0])
            molID = int(detabString[1])
            molID = '%09.0d' % molID
            #vfixxID = int(detabString[1])
            #molID = self.vfixxIDLookup[vfixxID]
            #self.mappedIDLookup[mapMolID] = molID
            if not(self.molDB.has_key(molID)):
                continue
            self.mappedMolecules.add(molID)
            mol = self.molDB[molID]
            score = float(detabString[4])
            confidence = float(detabString[5])
            truepos = int(detabString[12])
            falsepos = int(detabString[13])
            falseneg = int(detabString[14])
            mol.add_alignment(score, confidence, truepos, falsepos, falseneg)
            
            mapStart = int(detabString[7])
            mapStop = int(detabString[8])
            direction = detabString[6]
            forward = direction == '+'
            mol.detail_alignment(mapStart, mapStop, direction)            
        
            origLabKeys = sorted(mol.labels.keys())
            
            nEntries = detabString.__len__()
            alignCount = 1
            for ind in range(18,nEntries,2):
                labelKey = origLabKeys[int(detabString[ind]) - 1]
                refInd = int(detabString[ind + 1])
                lab = mol.labels[labelKey]
                lab.mapResult = refInd
        
        handle.close()
        return mapEntries


#this fn is broken becuase it assumes the global maximum is the first maximum,
# and that there is a minimum which is not at the end of the distribution
def getSimpleMinimum(valSort):
    #maxVal = 0
    #i_max = 0
    #for i,val in enumerate(valSort):
    #    if val > maxVal:
    #        maxVal = val
    #        i_max = i
    #cur_val = valSort[i_max]
    i_max = valSort.index(max(valSort))
    i_cur = i_max
    while(True):
        if i_cur + 1 >= valSort.__len__():
            break
        if valSort[i_cur + 1] > valSort[i_cur]:
            break
        i_cur += 1
    return i_cur


#this fn is an improvement on above because it looks on both sides of maximum
#start to the right, same as above fn, but instead of returning on
# the first increase (positive slope), require a minimum number of consecutive increasing slopes
# start to the right because typically, there are only two peaks and the left-most
#  is the highest
#if there is no peak on the right of the maximum, look on the left
#return index of minimum
#the main limitations of this algorithm are:
# -if there is a third peak and two are to the left of the maximum, it will
#  find the first valley (immediately to left of maximum) whereas we probably want
#  the left-most one.
# -if the peak is so small it's smaller than peakthresh
# These cases should both be rare enough we should be ok.
#valSort are y-values (you don't need the x-values, you just need the order,
# which you have in valSort bc it's a list)
#peakthresh is the minimum number of points in a row with the same slope 
# which are considered a peak
def getTwoPeakMinimum(valSort, peakthresh=3, verbose=False):

    maxidx = valSort.index(max(valSort)) 
    posslope = [] #look for increasing slopes
    if verbose :
        print "maxidx:", maxidx
    #descend to right : start at peak, so slope must be negative
    for i in range(maxidx, len(valSort)-1) : #-1 bc check next index for slope
        slope = (valSort[i] < valSort[i+1]) #True is positive, False is negative
        if i == maxidx :
            assert not slope, "getTwoPeakMinimum: problem with maximum, right"
        if not slope :
            if posslope :
                posslope = [] #reset if any negative slope is found
            continue
        posslope.append(i)
        if verbose :
            print "posslope:", posslope
        #if have this many positive slopes, find minimum up to here and you're done
        if len(posslope) >= peakthresh :
            if verbose :
                print "minimum between", maxidx, "and", posslope[-1]
            minlist = valSort[maxidx:posslope[-1]] #slicing makes copy, which is more robust bc index just finds first ele
            minimum = min(minlist)
            minidx = minlist.index(minimum) + maxidx
            assert minimum == valSort[minidx], "Bad index for "+str(minimum)+" : "+str(valSort[minidx])
            if verbose :
                print "minimum:", minimum, "idx:", minidx
            return minidx

    #if above loop ends, descend to left of maximum
    negslope = [] #now, look for decreasing slopes
    for i in range(maxidx,1,-1) : #backwards from maxidx to 1, compare to previous
        slope = (valSort[i] > valSort[i-1]) #True is positive, False is negative
        if i == maxidx : #left of max means positive slope
            assert slope, "getTwoPeakMinimum: problem with maximum, left"
        if slope :
            if negslope :
                negslope = []
            continue
        negslope.append(i)
        if verbose :
            print "negslope:", negslope
        if len(negslope) >= peakthresh :
            if verbose :
                print "minimum between", negslope[-1], "and", maxidx
            minlist = valSort[negslope[-1]:maxidx]
            minimum = min(minlist)
            minidx = minlist.index(minimum) + negslope[-1]
            assert minimum == valSort[minidx], "Bad index for "+str(minimum)+" : "+str(valSort[minidx])
            if verbose :
                print "minimum:", minimum, "idx:", minidx
            return minidx
            
    #will get here for no minimum, ie, single peak--return 0, which is first ele in list, which is very low snr
    if verbose :
        print "no minimum found"
    return 0



def makeHeader(bpp, nMolecules, vfixxFile):
    header = ''
    header += '#' + tab + 'File Name:' + tab + vfixxFile + newline
    header += '#' + tab + 'N Colors:' + tab + '1' + newline
    header += '#' + tab + 'N Fragments:' + tab + str(nMolecules) + newline
    header += '#' + tab + 'Bases per Pixel:' + tab + str(bpp) + newline
    header += '#' + tab + 'Frag Len (px)	Nick Locs (px)...' + newline
    return header


def makeBnxHeader(bpp, nLabelChannels = 1, minmollen = 0, labsnr = 0, quality=False):
    header = ''
    header += '# BNX File Version:\t0.1\n'
    header += '# Label Channels:\t%d\n' % nLabelChannels
    header += '# Nickase Recognition Site 1:\tUnknown\n'
    header += '# Measurement Timecode:\tXXXX-XX-XXXXX:XX:XX\n'
    header += '# Instrument Serial:\tXXXXXXXXXXXXX\n'
    header += '# Nanometers per Pixel:\t0\n'
    header += '# Stretch Factor:\t0\n'
    header += '# Bases per Pixel:\t%d\n' % bpp
    header += '# Min Molecule Length (Kb):\t%i\n' % minmollen
    header += '# Min Label SNR:\t%.1f\n' % labsnr
    if quality :
        #header += '# Quality Score QX01:	SNR\n'
        #header += '# Quality Score QX02:	Ave Intensity\n'
        #above is apparently wrong (according to RefAligner): QX01 is nothing, QX02 is SNR, QX12 is Intensity
        #disable Intensity--it's zero anyway and confusing RefAligner
        header += '# Quality Score QX02:	SNR\n'
        #header += '# Quality Score QX12:	Ave Intensity\n'
    header += '#0h Label Channel\tMapID\tLength\n'
    header += '#0f\tint\tint\tfloat\n'
    header += '#1h Label Channel\tLabelPositions[N]\n'
    header += '#1f int\tfloat\n'
    if nLabelChannels > 1 :
        header += '#2h\tLabelChannel\tLabelPositions[N]\n'
        header += '#2h\tint\tfloat\n'
    if quality :
        header += '#Qh\tQualityScoreID\tQualityScores[N]\n'
        header += '#Qf\tstr\tfloat\n'
    return header    


class molecule():
    """Container for single molecule data, author of single bnx entry
    """
    def __init__(self, ID, lengthPx, lengthKb, x, top, bottom, row, col, intensity, isStitched, nLabelChannels = 1):
        self.ID = ID
        self.lengthPx = lengthPx
        self.lengthKb = lengthKb
        self.x = x
        self.top = top
        self.bottom = bottom
        self.row = row
        self.col = col
        self.intensity = intensity
        self.avgSnr = [0] * nLabelChannels #just one float per label channel--label objects have their own snr
        self.nLabels = [0] * nLabelChannels
        self.nLabelChannels = nLabelChannels
        self.labels = [] #size is nLabelChannels, eles are dicts of label pos : label object
        self.isStitched = isStitched
        for i in range(nLabelChannels):
            self.labels.append({})
        self.mapResult = False


    #perform some consistency checks, return self.nLabels
    def getNLabels(self) :
        assert len(self.labels) == self.nLabelChannels
        assert len(self.nLabels) == len(self.labels)
        for i,nlab in enumerate(self.nLabels) :
            assert nlab == len(self.labels[i])
        return self.nLabels


    def getAvgSnr(self) :
        '''Loop over all labels and put average snr into self.avgSnr.'''
        #first, null avgSnr, otherwise you add to the average
        for i in self.avgSnr :
            i = 0
        for i,labdict in enumerate(self.labels) :
            for lab in labdict.values() :
                self.avgSnr[i] += lab.snr
        for i in range(len(self.avgSnr)) :
            if self.nLabels[i] : #avoid divide by zero
                self.avgSnr[i] /= self.nLabels[i] #make average by dividing by n labels
        

    def makeBnxEntry(self, bnxMolID, bpp, padding = 0.2, quality=False):
        lowestPos = 100
        highestPos = -100
        for labelDict in self.labels:
            for key in labelDict.keys():
                if key < lowestPos:
                    lowestPos = key
                if key > highestPos:
                    highestPos = key
        
        extendMolLeft = min(0, lowestPos - padding)
        extendMolRight = max(0, highestPos - self.lengthPx + padding)
        extendedLength = self.lengthPx - extendMolLeft + extendMolRight
        extendedLengthBb = bpp * extendedLength
        routeLine = '0\t' + bnxMolID + '\t%3.1f\n' % extendedLengthBb
        #qualLine1 = 'QX01'
        #qualLine2 = 'QX02' #disable this--it's zero anyway and confusing RefAligner
        qualLine1 = 'QX02' #this is SNR in 0.1
        curLabelChannel = 1
        labelOutput = ''
        labCount = 0
        stitched = 0
        dataLine = ''
        for labelDict in self.labels:
            dataLine += str(curLabelChannel)
            relyposs = sorted(labelDict.keys())
            if relyposs.__len__():
                lab0 = labelDict[relyposs[0]]
                prevFOV = lab0.FOV
                for val in relyposs:
                    lab = labelDict[val]
                    relypos = bpp * (val - extendMolLeft)
                    dataLineAdd = '\t%3.1f' % relypos
                    dataLine += dataLineAdd
                    qualLine1 += '\t%3.1f' % lab.snr
                    #fovDiff = lab.FOV - prevFOV
                    #qualLine2 += '\t%d' % (fovDiff)
                    #prevFOV = lab.FOV
                
            dataLine += '\t%3.1f\n' % extendedLengthBb
            qualLine1 += '\t0.0\n'
            #qualLine2 += '\t0\n'
            curLabelChannel += 1
        labelOutput = routeLine + dataLine
        if quality:
            #labelOutput += qualLine1 + qualLine2
            labelOutput += qualLine1 #disable Intensity
        return labelOutput
    
    def add_alignment(self, score, confidence, truepos, falsepos, falseneg):
        self.mapResult = True
        self.score = score
        self.confidence = confidence
        self.truepos = truepos
        self.falsepos = falsepos
        self.falseneg = falseneg
    
    def detail_alignment(self, mapstart, mapstop, direction):
        self.mapstart = mapstart
        self.mapstop = mapstop
        self.direction = direction

    def detail_chimeric(self, hits, distant):
        self.hits = hits
        self.distant = distant

    def add_label(self, lab, colorChannel):
        if colorChannel > len(self.labels)-1 : #obviously, this shouldn't happen
            return
        if self.labels[colorChannel].has_key(lab.relypos):
            #print 'Molecule has duplicate label as Rel Y position'
            #print('Label 1 ID: ' + str(self.labels[label.relypos].ID))
            #print('Label 2 ID: ' + str(label.ID))
            labelOrig = self.labels[colorChannel][lab.relypos]
            if labelOrig.snr < lab.snr:
                self.labels[colorChannel][lab.relypos] = lab
        else:
            self.nLabels[colorChannel] += 1
            self.labels[colorChannel][lab.relypos] = lab
            
    def getDiffSet(self, colorChannel):
        labyposs = self.labels[colorChannel].keys()
        diffSet = []
        for i,val0 in enumerate(labyposs):
            for j, val1 in enumerate(labyposs):
                if i >= j:
                    continue
                diffSet.append(abs(val1-val0))
        return diffSet

def getBppFromBnx(bnxFile):
    f1 = open(bnxFile)
    bpp = 500
    found = False
    while(True):
        line = f1.readline()
        if line == '':
            f1.close()
            break
        ind = line.find('Bases per Pixel:')
        if not(ind == -1):
            tokens = line.split('\t')
            bpp = float(tokens[1])
            found = True
            f1.close()
            break
    if not(found):
        print('  BPP not found in BNX file: %s' % bnxFile)
    return bpp



class label():
    """Container of single label data, parses .lab file
    """
    def __init__(self):
        pass
#        lo.LabID	lo.MolID	lo.Length	lo.LengthKb	lo.Scan	lo.FoV	lo.Row	lo.Col	lo.Channel	lo.Left	lo.Right	lo.Top	lo.Bot	lo.Xpos	lo.Ypos	lo.RealYpos	lo.RealYposKb	lo.CMR	lo.PixelCount	lo.AvgIntensity	lo.SNR	lo.Bkgnd	lo.Noise	lo.ChannelShift	lo.GaussFitErrX	lo.GaussFitErrY	lo.GaussSigmaX	lo.GaussSigmaY	lo.GaussShiftX	lo.GaussShiftY	lo.GaussAmpX	lo.GaussAmpY	lo.Mapped
#        x = float(tokens[13])
#        y = float(tokens[14])
    def fromRawEntry(self, ID, relypos, snr, fov, molID):
        self.ID = ID
        self.relypos = relypos
        self.snr = snr
        self.mapResult = 0
        self.FOV = fov
        self.molID = molID


    def fromLabFileLine(self, line, NanoStudio = False):
        tokens = line.split('\t')
        if NanoStudio:
            ID = int(tokens[0])
            relypos = float(tokens[11])
            pixelCount = float(tokens[13])
            snr = float(tokens[15])
            fov = int(tokens[5])
            molID = int(tokens[1])
        else:
            ID = int(tokens[0])
            relypos = float(tokens[15])
            pixelCount = float(tokens[18])
            snr = float(tokens[20])
            fov = int(tokens[5])
            molID = int(tokens[1])

        self.ID = ID
        self.relypos = relypos
        self.snr = snr
        self.mapResult = 0
        self.FOV = fov
        self.molID = molID
        

if __name__ == "__main__":
    #print(sys.argv)
    
    if len(sys.argv) < 2 or sys.argv[1] == '-h':
        print('\n  Makes two color bnx file from mol, 0.lab, 1.lab\n    Usage: $python molecule.py molFile\n')
        sys.exit()
    if not(os.path.exists(sys.argv[1])):
        print('  File Not Found %s' % sys.argv[1])
        sys.exit()
    #convertDetectOutToBnxNew(sys.argv[1], 500, labelSNR = 3.0, moltag = 0, numLabelChannels = 2, sizeLowerLimit = 150., sizeUpperLimit = 100000.)
    if len(sys.argv) == 2 :
        convertDetectOutToBnxNew(sys.argv[1])
    else :
        convertDetectOutToBnxNew(sys.argv[1], sys.argv[2])
    
    
