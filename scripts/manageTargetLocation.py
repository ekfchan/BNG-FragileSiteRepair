import os
from xml.dom.minidom import parse, parseString
import shutil
#import pdb

"""@package manageTargetLocation Metadata for IRYS intstrument/experiment
    
Reads instrument/experiment xml file and passes relevant information to 
ImageProcessingModule.Device() 
Also finds device tiff images

"""


import utilities
utilities.setVersion("$Id: manageTargetLocation.py 2686 2014-04-15 16:51:46Z vdergachev $")


class repositoryFolder():
    """Metadata container for instrument/experiment information subset. Pass info
    to ImageProcessingModule.Device() 
    
    Also finds device tiff images
    """
    def __init__(self, remoteDirectory, localDirectory, leadTiffTag = '', stretchFactor = 0.85):
        self.remoteDirectory = remoteDirectory
        self.localDirectory = localDirectory
        self.leadTiffTag = leadTiffTag
        self.FinalCycleIsFlush = ''        
        self.nColors = 0
        self.ExpectedOverlap = 0
        self.basesPerPixel = 500
        self.isAlpha40 = False #check for ChipTypeName is "Alpha 4.0" -- use in ImageProcessing.Scan
        self.localFailedTiffs = []
        self.ScanRowCount = 4
        self.ScanColumnCount = 26

        #done with defaults--read xml, set parameters
        self.getBasesPerPixel(stretchFactor)

        self.remoteTiffList = self.getRemoteTiffPathList()
        self.localTiffList = self.setLocalTiffList()
        
        
    def getBasesPerPixel(self, stretchFactor):
        """Reads and parses instrument/experiment xml file
        """
        folderName = os.path.split(self.remoteDirectory.rstrip('/'))[1]
        xmlFile = os.path.join(self.remoteDirectory, folderName + ' Metadata.xml')
        if not os.path.exists(xmlFile):
            xmlFile = ""
            #try to find any old xml file in remoteDirectory
            if os.path.isdir(self.remoteDirectory) :
                for qfile in os.listdir(self.remoteDirectory) :
                    if qfile.endswith(".xml") :
                        xmlFile = os.path.join(self.remoteDirectory, qfile)
                        break #take first xml found
            if not xmlFile :
                print 'XML File Not Found : ' + folderName
                return -1
        dom1 = parse(xmlFile)
        a1 = dom1.getElementsByTagName('FieldOfViewSizeY')[0]
        
        b1 = a1.childNodes[0]
        fovSize = float(b1.data)
        basesPerPixel = (1./stretchFactor) * (fovSize/512.) * (1/3.4e-4) 
        self.basesPerPixel = basesPerPixel
        
        self.Pitch=0.7*512.0/fovSize

        #a2 = dom1.getElementsByTagName('OperatorUserName')[0]
        #b2 = a2.childNodes[0]
        #self.OperatorUserName = b2.data
        
        a3 = dom1.getElementsByTagName('ScanCount')[0]
        b3 = a3.childNodes[0]
        self.ScanCount = int(b3.data)
        
        a4 = dom1.getElementsByTagName('TimeStarted')[0]
        b4 = a4.childNodes[0]
        self.TimeStarted = b4.data
        
        a5 = dom1.getElementsByTagName('InstrumentSerialNumber')[0]
        b5 = a5.childNodes[0]
        self.InstrumentSerialNumber = b5.data
        
        a6 = dom1.getElementsByTagName('StageStepY')[0]
        b6 = a6.childNodes[0]
        self.StageStepY = b6.data

        a7 = dom1.getElementsByTagName('FieldOfViewSizeY')[0]
        b7 = a7.childNodes[0]
        self.FieldOfViewSizeY = b7.data
        
        a8 = dom1.getElementsByTagName('DatasetName')[0]
        b8 = a8.childNodes[0]
        self.DatasetName = b8.data
        
        a9 = dom1.getElementsByTagName('FinalCycleIsFlush')[0]
        b9 = a9.childNodes[0]
        self.FinalCycleIsFlush = b9.data
        
        a10 = dom1.getElementsByTagName('LasersCount')[0]
        b10 = a10.childNodes[0]
        self.nColors = int(b10.data)
        
        a11 = dom1.getElementsByTagName('ScanRowCount')[0]
        b11 = a11.childNodes[0]
        self.ScanRowCount = int(b11.data)
        
        a12 = dom1.getElementsByTagName('ScanColumnCount')[0]
        b12 = a12.childNodes[0]
        self.ScanColumnCount = int(b12.data)
        
        
        #print "FOVy:", self.FieldOfViewSizeY, "stagestepy:", self.StageStepY #debug
        self.ExpectedOverlap = 512 * (float(self.FieldOfViewSizeY) - float(self.StageStepY)) / float(self.FieldOfViewSizeY)
        #print "ExpectedOverlap", self.ExpectedOverlap #debug

        chiptype = dom1.getElementsByTagName('ChipTypeName')[0].childNodes[0].data
        if chiptype.find("4.") != -1 : #"4." is an attempt to keep this for whenever 4.0 is 4.1. If that's correct.
            #print "isAlpha40 -- getBasesPerPixel" #debug
            self.isAlpha40 = True


        
    def getLocalTiffList(self):
        localFiles = os.listdir(self.localDirectory)
        tiffs = []
        for localFile in localFiles:
            if localFile.endswith('.tiff'):
                tiffs.append(localFile)
        return tiffs
            
    def setLocalTiffList(self):
        localTifName = 0
        localTifPaths = []
        for tiffPath in self.remoteTiffList:
            localTifName += 1
            localTifFileName = '%02d.tiff' % localTifName
            activeTiffFile = os.path.join(self.localDirectory, self.leadTiffTag + localTifFileName)
            localTifPaths.append(activeTiffFile)
        return localTifPaths

    def getRemoteTiffPathList(self):
        remoteTiffList = os.listdir(self.remoteDirectory)
        if self.FinalCycleIsFlush == 'True':
            dump = remoteTiffList.pop()
        remoteTiffs = [os.path.join(self.remoteDirectory, x) for x in remoteTiffList if x.endswith('.tiff') and (x.find('Processed') == -1) and (x.find("_Scan") != -1)]
        return remoteTiffs
        
    def transferFilesLocally(self):
        if os.path.exists(self.localDirectory):
            shutil.rmtree(self.localDirectory)
        os.makedirs(self.localDirectory)
        nScans = self.remoteTiffList.__len__()
        for i in range(nScans):
            localTiff = self.localTiffList[i]
            remoteTiff = self.remoteTiffList[i]
            #print('Copying File from ' + tiffPath + ' to local folder')
            transferCount = 0
            while(True):
                shutil.copy(remoteTiff, localTiff)
                if os.path.getsize(remoteTiff) != os.path.getsize(localTiff):
                    transferCount += 1
                    if transferCount == 2:
                        self.localTiffList.remove(localTiff)
                        self.localFailedTiffs.append(localTiff)
                        break
                    continue
                break
                
            
            

