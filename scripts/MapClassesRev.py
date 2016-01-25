import pdb
import math
import os

import utilities as util

"""@package MapClassesRev Interprate RefAligner map results for characterization

Read Cmap files and .map file and get contig statistics including mapping stats
"""


util.setVersion("$Id: MapClassesRev.py 3560 2015-02-12 18:40:59Z wandrews $")


def ContigCharacterizationNoRef(varsP, stgstr=""):
    """Report contig stats for cases without reference
    """
    qCmap = MultiCmap(varsP.latestMergedCmap)
    outStr = ''
    for i,cMap in enumerate(qCmap.cmapDB.values()):
        if i==0:
            outStr += cMap.cmapHeader() + '\n'
        outStr += cMap.reportAlignStatsUnmapped() + '\n'

    k = 1e-6
    varsP.totAssemblyLenMb = k * qCmap.totalLen
    if stgstr :
        outStr += 'Stage Summary: %s\n' % stgstr
    outStr += 'N Genome Maps: %d\n' % qCmap.n
    outStr += 'Total Genome Map Len (Mb):   %6.3f\n' % (k * qCmap.totalLen)
    outStr += 'Avg. Genome Map Len  (Mb):   %6.3f\n' % (k * qCmap.aveLen)
    outStr += 'Genome Map n50       (Mb):   %6.3f\n' % (k * qCmap.n50)
    return outStr

def TopLevelCharacterization(varsP, resultFileRoots, stgstr=""):
    """Main Analysis Function for characterization results (.cmap, .map)
    
    """
    if not varsP.curCharacterizeCmaps : #if you don't even have cmaps, no useful information can come below--bail
        print "Error in MapClassesRev.TopLevelCharacterization: varsP.curCharacterizeCmaps is empty"
        return
    mapFiles = []
    qcmapFiles = []
    rcmapFiles = []
    xmapFiles = []
    errFiles = []
    for i,curCmapFile in enumerate(varsP.curCharacterizeCmaps):
        qryCmapFile = resultFileRoots[i] + '_q.cmap'
        #it appears that curCmapFile is the input and qryCmapFile is the output
        # but the _q.cmap only has mapped contigs, while the curCmapFile should have all of them.
        if i==0:
            #qCmap = MultiCmap(curCmapFile,qCmapFile=qryCmapFile)
            qCmap = MultiCmap(curCmapFile)
            #qCmap = MultiCmap(qryCmapFile)
        else:
            qCmap.readCmapFile(curCmapFile)
            #qCmap.readCmapFile(qryCmapFile)
        if not len(qCmap.cmapDB) : #if varsP.curCharacterizeCmaps has a bad path, quit
            return
        qcmapFiles.append(qryCmapFile)
        refCmapFile = resultFileRoots[i] + '_r.cmap'
        if i==0:
            rCmap = MultiCmap(refCmapFile)
        else:
            rCmap.readCmapFile(refCmapFile)
        rcmapFiles.append(refCmapFile)
        mapFiles.append(resultFileRoots[i] + '.map')
        xmapFiles.append(resultFileRoots[i] + '.xmap')
        errFiles.append(resultFileRoots[i] + '.err')
        
    qCmap.profileMapLengths()
    rCmap.profileMapLengths()
    

    xResult = MapResults(qCmap,rCmap)
    for mapFile in mapFiles:
        xResult.readMapFile(mapFile)
    xResult.calculateStats()
    datOut = xResult.mapHeader() + '\n' 
    datOut += xResult.reportResults()
    if stgstr :
        datOut += 'Stage Summary: %s\n' % stgstr
    datOut += xResult.reportSummary()
    varsP.totAssemblyLenMb = 1e-6 * xResult.qryCmap.totalLen
    
    if not varsP.latestMergedCmap : #if no merged cmap, no need to do file manipulation below
        print "Error in MapClassesRev.TopLevelCharacterization: varsP.latestMergedCmap is empty"
        return datOut

    charOutputFileRoot = os.path.split(varsP.latestMergedCmap)[-1]
    charOutputFilePath = os.path.join(varsP.contigAlignTarget, charOutputFileRoot.rstrip('.cmap'))
    if len(varsP.curCharacterizeCmaps) > 1:
        qcmapResultFile = charOutputFilePath + '_q.cmap'
        CatFilesAfterSplitCharacterization(qcmapFiles, qcmapResultFile)
        rcmapResultFile = charOutputFilePath + '_r.cmap'
        CatFilesAfterSplitCharacterization(rcmapFiles, rcmapResultFile, uniqueIDsOnly=True)
        xmapResultFile = charOutputFilePath + '.xmap'
        CatFilesAfterSplitCharacterization(xmapFiles, xmapResultFile, uniquifyIDs=True, xmapFile=True)
        for target in varsP.curCharacterizeCmaps:
            os.remove(target)
            
        # Remove .err and .map files -- no need to remove potentially useful files
        #for target in errFiles: # + mapFiles:
        #    os.remove(target)
    
    return datOut

def CatFilesAfterSplitCharacterization(fileList, outputFile, uniqueIDsOnly = False, uniquifyIDs=False, xmapFile = False):
    '''Merge map and cmap files after distributed mapping operation
    
    Rejoin all the files that were originally split to distribute the 
    characterization job.  Remove originals when finished.  Certain 
    ID operations needed for certain files - must be Irys-CAT compatible.
    '''
    
    f2 = open(outputFile, 'w')
    usedIDs = []
    curReplaceID = 1
    for i, curFile in enumerate(fileList):
        f1 = open(curFile)
        curID = ''
        prevID = ''
        writeCurID = False
        while(True):
            line = f1.readline()
            if line == '':
                f1.close()
                break
            if line[0] == '#':
                if i == 0:
                    if xmapFile:
                        if line.startswith('# Ref'):
                            f2.write('# Reference Maps From:\t%s\n' % (outputFile.replace('.xmap', '_r.cmap')))
                            continue
                        if line.startswith('# Que'):
                            f2.write('# Query Maps From:\t%s\n' % (outputFile.replace('.xmap', '_q.cmap')))
                            continue
                    f2.write(line)
                continue
            if uniqueIDsOnly or uniquifyIDs:
                indID = line.find('\t')
                curID = line[:indID]
                if uniquifyIDs:
                    f2.write('%d' % curReplaceID)
                    f2.write(line[indID:])
                    curReplaceID += 1
                    continue
                if curID == prevID:
                    if writeCurID:
                        f2.write(line)
                else:
                    if usedIDs.__contains__(curID):
                        writeCurID = False
                    else:
                        writeCurID = True
                        usedIDs.append(curID)
                        f2.write(line)
                prevID = curID
            else:
                f2.write(line)
        os.remove(curFile)
    f2.close()
    

class MapResults():
    """Syntesizes Cmap class and map result class for statistics
    """
    def __init__(self, qryCmap, refCmap):
        self.hitDB = [] #can't use map id as dict key because a map can have multiple alignments; use list instead
        self.usedCmapIDs = set()
        self.qryCmap = qryCmap
        self.refCmap = refCmap
        self.uniqueMapLength = 0
        self.totalAlignLen = 0
    
    def mapHeader(self):
        outstr = ''
        outstr += "Genome Map Characterization:\n"
        outstr += "cID  rID  len"
        outstr += "  Cov"
        #outstr += "  rID" #ref index for either of these
        outstr += "  alignlen  alignlen/len"
        outstr += "  Conf  Conf/lenkb"
        outstr += "  M  FP  FN RMS  O  sf  sd  C%" #"  res" #--ditch res (not bpp)
        return outstr
        
    def reportResults(self):
        outStr = ''
        for curMapRes in self.hitDB :
            outStr += curMapRes.reportAlignStats() + '\n'
            self.usedCmapIDs.add(curMapRes.qryCmapID)
        for qrySingleCmap in self.qryCmap.cmapDB.values():
            if qrySingleCmap.cmapID in self.usedCmapIDs :
                continue
            outStr += qrySingleCmap.reportAlignStatsUnmapped() + '\n'
        return outStr

    def reportSummary(self):
        outStr = ''
        k = 1e-6
        outStr += 'N Genome Maps: %d\n' % self.qryCmap.n
        outStr += 'Total Genome Map Len (Mb):   %6.3f\n' % (k * self.qryCmap.totalLen)
        outStr += 'Avg. Genome Map Len  (Mb):   %6.3f\n' % (k * self.qryCmap.aveLen)
        outStr += 'Genome Map n50       (Mb):   %6.3f\n' % (k * self.qryCmap.n50)
        if not self.refCmap.totalLen : #probably none of below makes sense in this case, plus avoid raising ZeroDivisionError
            return outStr
        outStr += 'Total Ref Len    (Mb):   %6.3f\n' % (k * self.refCmap.totalLen)
        outStr += 'Total Genome Map Len / Ref Len  : %.3f\n' % (self.qryCmap.totalLen/self.refCmap.totalLen)
        if self.qryCmap.n:
            outStr += 'N Genome Maps total align       :   %d (%4.2f)\n' % (len(self.usedCmapIDs), len(self.usedCmapIDs)/float(self.qryCmap.n))
        outStr += 'Total Aligned Len             (Mb) : %6.3f\n' % (k * self.totalAlignLen)
        outStr += 'Total Aligned Len / Ref Len        : %6.3f\n' % (self.totalAlignLen/self.refCmap.totalLen)
        #outStr += 'Unique - within match group extents\n'
        outStr += 'Total Unique Aligned Len      (Mb) : %6.3f\n' % (k * self.uniqueMapLength)
        outStr += 'Total Unique Len / Ref Len         : %6.3f\n' % (self.uniqueMapLength/self.refCmap.totalLen)
        #outStr += 'Unique Aligned Len     (Mb) : Unknown\n'
        #outStr += 'Unique Aligned Len / Ref Len: Unknown\n' 
        return outStr
        
    def readMapFile(self, mapFile, verbose=0):
        commentChars = ['#','S','M']
        if not util.checkFile(mapFile, ".map") :
            print "Error in MapResults.readMapFile: missing file", mapFile
            return
        for line in open(mapFile) :
            if commentChars.__contains__(line[0]):
                continue
            curResult = SingleMapResult(line, self.qryCmap, self.refCmap)
            #if verbose > 0 and self.hitDB.has_key(curResult.qryCmapID): #if you want this back, use 'in'
            #    print "  Warning MapID %d already counted" % curResult.qryCmapID
            self.hitDB.append( curResult )

            
    def calculateStats(self):
        # Prepare a graph of overlapping matchgroups for DFS
        self.totalAlignLen = 0
        refContigDB = {}
        for curRes in self.hitDB :
            try : #this raised ZeroDivisionError -- can also potentially have ValueError due to, eg, log(0)
                curRes.fitSizingErrorModel()
            except :
                pass
            curRes.getSegmentRMS()
            if refContigDB.has_key(curRes.refCmapID):
                refContigDB[curRes.refCmapID].append((curRes.startMatch, curRes.endMatch))
            else:
                refContigDB[curRes.refCmapID] = [(curRes.startMatch, curRes.endMatch)]
            self.totalAlignLen += curRes.alignLen
        self.uniqueMapLength = 0
        for key in refContigDB.keys():
            self.uniqueMapLength += FindUniqueMatchLength(refContigDB[key])
        
        
def FindUniqueMatchLength(inList):
    """Search to find overlapping map results
    
    Make graph, nodes:contigs, edges:hits
    Depth First Search, join map regions
    """
    # Uses DFS routine to group all overlapping match groups together
    uniqueLength = 0
    matchDict = {}
    for i,pair0 in enumerate(inList):
        for j,pair1 in enumerate(inList):
            if i==0:
                matchDict[j] = set([j])
            if i <= j:
                continue
            if (pair0[1] >= pair1[0] and pair0[1] <= pair1[1]) \
            or (pair1[1] >= pair0[0] and pair1[1] <= pair0[1]):
                matchDict[j].add(i)
                matchDict[i].add(j)
    
    curGSearch = GraphSearch(matchDict)
    for connectedNodes in curGSearch.ConnectedNodeGroups:
        currentEnds = []
        for ind in connectedNodes:
            currentEnds.append(inList[ind][0])
            currentEnds.append(inList[ind][1])
        uniqueLength += max(currentEnds) - min(currentEnds)
    return uniqueLength

class GraphSearch():
    """ Finding the complete set of maps that have overlap based
    on the mapping results is acutally a search problem.
    Here we use DFS to find the connected nodes.
    """
    def __init__(self, nonDirGraph):
        # States: 0=White, 1=Gray, 2=Black
        self.Graph = nonDirGraph
        self.StateDict = {}
        self.AllWhite()
        self.ConnectedNodeGroups = []
        self.CurrentPath = set([])
        self.PerformSearch()
        
    def AllWhite(self):
        for key in self.Graph.keys():
            self.StateDict[key] = 0
    
    def PerformSearch(self):
        for key in self.Graph.keys():
            if self.StateDict[key] > 0:
                continue
            self.CurrentPath = set([])
            self.RunDFS(key)
            self.ConnectedNodeGroups.append(map(None, self.CurrentPath))
            
    def RunDFS(self, key):
        # Recursive DFS part
        self.StateDict[key] = 1
        self.CurrentPath.add(key)
        for subKey in self.Graph[key]:
            if self.StateDict[subKey] == 0:
                self.RunDFS(subKey)
        self.StateDict[key] = 2
            
    
    
    '''            
    uniqueMapLength = 0
    while(True):
        # find the total length of overlaping hits
        if matchMembers.__len__() == 0:
            break
        curMatchMembers = map(None, matchMembers)     # i
        curMatchGroup = matchDB[curMatchMembers[0]]  # [i,j]
        curEndPoints = []
        for curMatch in curMatchGroup:
            for endPoint in inList[curMatch]:
                curEndPoints.append(endPoint)
            matchMembers.remove(curMatch)
            dump = matchDB.pop(curMatch)
        uniqueMapLength += max(curEndPoints) - min(curEndPoints)
        
    for key in matchDB:
        # find the total length of non-overlapping hits
        uniqueMapLength += inList[key][1] - inList[key][0]
    return uniqueMapLength
    '''
       
        
class SingleMapResult():
    """Parsing and analysis of single map entry
    
    Generates noise characteristics of sequence comparison
    Re-implementation of noise calce from refalign2.cc in RefAligner project
    """
    def __init__(self, line, qryMap, refMap):
        self.X = []
        self.Y = []
        self.qryMap = qryMap
        self.refMap = refMap
        self.parseMapLine(line)
        self.matchPairs = []

        self.C = 1
        self.SD = 0.2
        self.SF = 0.2
        self.rms = 99.
        
    def reportAlignStats(self):
        k = 1e-6
        curQryMap = self.qryMap.cmapDB[self.qryCmapID]
        curRefMap = self.refMap.cmapDB[self.refCmapID]
        outstr = ''
        outstr += '% 3d% 3d' % (self.qryCmapID, self.refCmapID)
        outstr += '% 8.3f% 4d' % (k*curQryMap.cmapLen, curQryMap.aveCovg)
        outstr += '% 8.3f% 5.1f' % (k*self.alignLen, self.alignLen/curQryMap.cmapLen)
        confPerKb = 1000 * self.conf/curQryMap.cmapLen
        outstr += '% 6.1f% 7.3f' % (self.conf, confPerKb)
        outstr += '% 5d% 3d% 3d' % (self.trueLabels, self.falsePos, self.falseNeg)
        outstr += '% 6.3f% 3d' % (self.rms, self.out)
        outstr += '% 6.3f% 6.3f% 5.1f' % (self.SF, self.SD, (self.C-1)*100)
        return outstr

        
    def parseMapLine(self, line):
        # MappedMoleculeId	MoleculeId	MoleculeIndex	ContigId	Score	Confidence	Direction	StartLocation	
        # EndLocation	StartMatchLocation	EndMatchLocation	DetectedLabelCount	TruePositiveLableCount	FalsePositiveLableCount	FalseNegativeLabelCount	Log10Pvalue	LeftEndChim	RightEndChim	NanoLabelID1	RefLabelIDleft1	RefLabelIDright1	NanoLabelID2	RefLabelIDleft2	RefLabelIDright2
        tokens = line.split('\t')
        self.qryCmapID = int(tokens[1])
        self.refCmapID = int(tokens[3])
        self.conf = float(tokens[5])
        directionTag = tokens[6].strip()
        if directionTag == '+':
            self.forward = 1
        else:
            self.forward = 0
        self.startLoc = float(tokens[7])
        self.endLoc = float(tokens[8])
        self.startMatch = float(tokens[9])
        self.endMatch = float(tokens[10])
        self.alignLen = abs(self.endMatch - self.startMatch)
        self.detectedLabels = int(tokens[11])
        self.trueLabels = int(tokens[12])
        self.falsePos = int(tokens[13])
        self.falseNeg = int(tokens[14])
        matches = [int(x) for x in tokens[18:]]
        grpMatches = []
        while(True):
            grpMatches.append(matches[:3])
            if matches.__len__() == 3:
                break
            matches = matches[3:]
        self.getMatchSegments(grpMatches)
    
    def getSegmentRMS(self, outlier=0.5, verbose=0):
        sqrsum = 0.
        ct = 0
        self.out = 0
        for i,Y in enumerate(self.Y):
            X = self.X[i]
            err = Y - self.C * X
            if Y == 0 : continue #just in case
            deviation = abs(err/Y)
            if deviation > outlier:
                self.out += 1
                if verbose > 0 :
                    print '  Qry=%d,Ref=%d,Deviation=%6.3f' % (self.qryCmapID, self.refCmapID, deviation)
                continue
            sqrsum += err * err
            ct += 1
        self.rms = (math.sqrt((1./ct) * sqrsum) if ct > 0 else 0)
        
    
    def getMatchSegments(self, matches):
        self.cumX = 0
        self.cumY = 0
        leftMatchY = 0
        leftMatchX = 0
        for i,val in enumerate(matches):
            qrySite = val[0]
            refSite0 = val[1]
            refSite1 = val[2]
            if qrySite == -1 or refSite0 == -1:
                continue
            rightMatchX = self.qryMap.getXPos(self.qryCmapID, qrySite)
            rightMatchY = self.refMap.getYPos(self.refCmapID, refSite0, refSite1)
            if leftMatchY:
                X = 0.001 * math.fabs(rightMatchX - leftMatchX)
                Y = 0.001 * math.fabs(rightMatchY - leftMatchY)
                #print('Y=%3.1f,X=%3.1f' % (Y,X))
                #if Y > 1000:
                #    pdb.set_trace()
                self.X.append(X)
                self.Y.append(Y)
                self.cumX += X
                self.cumY += Y
            leftMatchX = rightMatchX
            leftMatchY = rightMatchY
    
    def fitSizingErrorModel(self):
        ''' Implements size model fitting from 
        refalign2.cc from RefAligner project '''
        C = self.C
        A = self.SF * self.SF
        B = self.SD * self.SD
        Lsum=0.0        # sum of log(A+By) excluding cases with perr->end */
        errsum=0.0      # sum of (Cx-y)^2/(A+By) */
        xysum = 0.0     # sum of xy/(A+By) */
        y2sum = 0.0     # sum of yy/(A+By) */
        for i,Y in enumerate(self.Y):
            X = self.X[i]
            var = (A+B*Y)
            err = Y-X
            #print('Y=%3.1f,X=%3.1f' % (Y,X))
            Ivar = 1.0/var
            errsum += err*err*Ivar
            Ivar *= Y
            y2sum += Y*Ivar
            xysum += X*Ivar
        Cscale = y2sum/xysum
        C *= Cscale
        A *= Cscale*Cscale
        B *= Cscale*Cscale
        
        Lsum = 0.0
        errsum = 0.0
        Lcnt = 0.0
        for i,Y in enumerate(self.Y):
            X = self.X[i]
            err = Y - C * X
            var = A+B*Y
            Lsum += math.log(var)
            Lcnt += 1.0
            errsum += err*err/var
        delta = errsum/Lcnt
        A *= delta
        B *= delta

        bestLL = -1e30
        bestA = A
        ii = 0
        while(ii<3):
            Lsum=0.0    # /* sum of log(A+By) */
            errsum=0.0  # /* sum of (Cx-y)^2/(A+By) */
            Isum=0.0    # /* sum of 1/(A+By) */
            errsum2=0.0 # /* sum of (Cx-y)^2/(A+By)^2 */
            Isum2=0.0   #/* sum of 1/(A+By)^2 */
            errsum3=0.0 #/* sum of (Cx-y)^2/(A+By)^3 */
            for i,Y in enumerate(self.Y):
                X = self.X[i]
                err = Y - C * X
                err *= err
                var = A+B*Y
                Ivar = 1.0/var
                Ivar2 = Ivar*Ivar
                errsum += err * Ivar;
                errsum2 += err*Ivar2;
                errsum3 += err*Ivar2*Ivar;
                
                Lsum += math.log(var)
                Isum += Ivar
                Isum2 += Ivar2
                
            LL = -(Lsum+errsum)*0.5
            LL1 = 0.5*(errsum2-Isum)
            LL2 = 0.5*Isum2-errsum3
            if LL < bestLL:
                A = (A+bestA)*0.5
                if LL < bestLL - 1e-6:
                    ii -= 1
                    #print('change iter')
                    continue
            bestA = A
            bestLL = LL
            
            maxdelta = math.fabs(A)*0.5;
            if LL2 < 0.0:
                delta = min(maxdelta,math.fabs(LL1/LL2))
            else:
                delta = maxdelta
            delta = math.copysign(delta,LL1)
            if (A + delta < 0.0):
                delta = -A
            A += delta
            ii += 1
            
        #print('SF=%3.3f->%3.3f' % (self.SF, math.sqrt(A)))
        #print('SD=%3.3f->%3.3f' % (self.SD, math.sqrt(B)))
        #print('C=%3.3f->%3.3f' % (self.C, C))
        self.SF = math.sqrt(A)
        self.SD = math.sqrt(B)
        self.C = C

        
class MultiCmap():
    """Container for collection of contigs, generates set statistics
    """
    def __init__(self, oCmapFile, qCmapFile=None):
        self.cmapDB = {}
        self.readCmapFile(oCmapFile)
        if qCmapFile:
            self.readCmapFile(qCmapFile)
        self.n = 0
        self.n50 = 0
        self.totalLen = 0
        self.aveLen = 0
        self.profileMapLengths()
        
    def profileMapLengths(self):
        cmapLens = []
        nMaps = 0
        totalLen = 0
        for curCmap in self.cmapDB.values():
            cmapLen = curCmap.cmapLen
            cmapLens.append(cmapLen)
            totalLen += cmapLen
            nMaps += 1
        
        if nMaps:
            self.n = nMaps
            self.n50 = util.getn50(cmapLens)
            self.totalLen = totalLen
            self.aveLen = totalLen / nMaps
    
    def readCmapFile(self, cmapFile):
        if util.checkFile(cmapFile) :
            f1 = open(cmapFile)
        else :
            print "Error in MapClassesRev.MultiCmap.readCmapFile: missing file", cmapFile
            return
        newCmap = True
        for line in f1 :
            if line[0] == '#':
                continue
            tokens = line.split('\t')
            if newCmap:
                cmapID = int(tokens[0])
                cmapLen = float(tokens[1])
                curCmap = Cmap(cmapID, cmapLen)
                self.cmapDB[cmapID] = curCmap
                nSites = int(tokens[2])
                newCmap = False
                #print('Adding CMAP %d' % cmapID)
            siteID = int(tokens[3])
            if siteID > nSites:
                newCmap = True
                continue
            siteLoc = float(tokens[5])
            covg = float(tokens[7])
            curCmap.addSite(siteID, siteLoc, covg)
    
    def getXPos(self, contigID, siteID):
        curCmap = self.cmapDB[contigID]
        return curCmap.siteDB[siteID]
        
    def getYPos(self, contigID, siteID0, siteID1):
        try:
            curCmap = self.cmapDB[contigID]
        except:
            pdb.set_trace()
        sitePos0 = curCmap.siteDB[siteID0]
        if siteID0 != siteID1:
            sitePos1 = curCmap.siteDB[siteID1]
            if sitePos0 > sitePos1:
                return (sitePos0 + sitePos1) / 2
            else:
                return (sitePos0 + sitePos1) / 2
        else:
            return sitePos0 
        
class Cmap():
    """Container for single contig
    """
    def __init__(self,cmapID,cmapLen):
        self.cmapID = cmapID
        self.cmapLen = cmapLen
        self.nSites = 0
        self.sumCovg = 0
        self.aveCovg = 0
        self.siteDB = {}
        self.sites = []
        self.segments = []
        
    def addSite(self, siteID, siteLoc, covg):
        self.siteDB[siteID] = siteLoc
        self.sites.append(siteLoc)
        if self.nSites:
            self.segments.append(siteLoc - self.sites[-1])
        self.nSites += 1
        self.sumCovg += covg
        self.aveCovg = self.sumCovg / self.nSites
    
    def cmapHeader(self):
        outstr = ''
        outstr += "Contig Characterization:\n"
        outstr += "cID  rID  len"
        outstr += "  Cov"
        return outstr
        
    def reportAlignStatsUnmapped(self):
        k = 1e-6
        outstr = ''
        outstr += '% 3d% 3s' % (self.cmapID, '-')
        outstr += '% 8.3f% 4d' % (k*self.cmapLen, self.aveCovg)
        return outstr
        
    
    
