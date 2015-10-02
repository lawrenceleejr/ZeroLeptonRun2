import rootpy.ROOT as ROOT
from rootpy.io import root_open
#import rootpy.numpy
import root_numpy as rnp

#import ZeroLeptonNTVars

import os
from os import path
# Workaround to fix threadlock issues with GUI
ROOT.PyConfig.StartGuiThread = False
import logging
logging.basicConfig(level=logging.INFO)

logging.info("loading packages")
#import shutil                                                                                                                                       #shutil.copyfile(ROOT.gSystem.ExpandPathName('$ROOTCOREDIR/scripts/load_packages.C'), 'load_packages.C')
#lineLoadPackages = '.x ' + ROOT.gSystem.ExpandPathName('$ROOTCOREDIR/scripts/load_packages.C')
#print lineLoadPackages
ROOT.gROOT.Macro("$ROOTCOREDIR/scripts/load_packages.C")
#ROOT.gROOT.ProcessLine(lineLoadPackages)

# Initialize the xAOD infrastructure
ROOT.xAOD.Init()

files = [f for f in os.listdir(".") if os.path.isfile(f)]

print files

totalEventCounter = 0

for ifile in files:
    if( (".root" in ifile) ):
#    print ifile.isfile
        rootfile = root_open(ifile)#ROOT.TFile.Open(ifile)
        treeName = "PassThroughNT"
#        print dir(rootfile)

        if(rootfile.GetListOfKeys().Contains(treeName)):
            tree = rootfile.PassThroughNT
            nentries = tree.GetEntries()

            totalEventCounter += nentries
            print nentries
            for event in tree:

                print "jetpt size : " + str(event.jetPt.size())
#                jetpt = event.jetPt.at(0)
#                print jetpt

                triggerBits = event.triggerBits
                print "triggerBits size : " + str(triggerBits.size())
                if(triggerBits.size() > 0 ) :
                    print triggerBits.at(0)

        rootfile.Close()

print "Total number of events in the sample:" , totalEventCounter
