import rootpy.ROOT as ROOT
import rootpy.stl as stl
from rootpy.io import root_open
#import rootpy.numpy
import root_numpy as rnp

import os
from os import path
# Workaround to fix threadlock issues with GUI
import logging
logging.basicConfig(level=logging.INFO)

logging.info("loading packages")
#import shutil                                                                                                                                       #shutil.copyfile(ROOT.gSystem.ExpandPathName('$ROOTCOREDIR/scripts/load_packages.C'), 'load_packages.C')
lineLoadPackages = '.x ' + ROOT.gSystem.ExpandPathName('$ROOTCOREDIR/scripts/load_packages.C')
print lineLoadPackages
ROOT.gROOT.ProcessLine(lineLoadPackages)

# Initialize the xAOD infrastructure
ROOT.xAOD.Init()
ROOT.gStyle.SetOptStat(0)

#files = [f for f in os.listdir(".") if os.path.isfile(f)]
files = ["/afs/cern.ch/work/r/rsmith/ttbar_trigger/ttbar_alljet.root",
         "/afs/cern.ch/work/r/rsmith/ttbar_trigger/ttbar_tenjet.root"]

print files

totalEventCounter = 0

treeName = "PassThroughNT"

rootfile = {files[0] : root_open(files[0]),
            files[1] : root_open(files[1])
            }


trees = {files[0] : rootfile[files[0]].PassThroughNT,
         files[1] : rootfile[files[1]].PassThroughNT
         }

# nentries = tree.GetEntries()
# totalEventCounter += nentries
# print nentries
#vecFloatCreator = stl.vector(float)


mDeltaR_allVsTen = {
   5  :  ROOT.Hist2D( 50, 0, 1000, 50, 0, 1000, name="mDeltaR_allVsTen_GT_5jets", title="mDeltaR_allVsTen_GT_5jets"),
   7  :  ROOT.Hist2D( 50, 0, 1000, 50, 0, 1000, name="mDeltaR_allVsTen_GT_7jets", title="mDeltaR_allVsTen_GT_7jets"),
   10 :  ROOT.Hist2D( 50, 0, 1000, 50, 0, 1000, name="mDeltaR_allVsTen_GT_10jets", title="mDeltaR_allVsTen_GT_10jets"),
}

njet = ROOT.Hist( 15 , -.5 , 14.5, name = "njet" )
counter = 0

for entry in trees[files[0]] :
    trees[files[1]].GetEntry(counter)#update the other one too
    if(counter % 10000 == 0 ) : print counter
    if(counter > 300000 ) : break
    counter += 1

    # jetPtVec = trees[files[0]].jetPt# vecFloatCreator()
    nhltjet = trees[files[0]].GetLeaf("nHLTJets").GetValue(0)
    njet.Fill(nhltjet)

    mDeltaR = {}
    for ifile in files :
        mDeltaR[ifile] = trees[ifile].GetLeaf("RJVars_PP_MDeltaR").GetValue(0)/1000.

    for  njetcut in mDeltaR_allVsTen.keys() :
#        if jetPtVec.size() > njetcut :
        if nhltjet > njetcut :
            mDeltaR_allVsTen[njetcut].Fill(mDeltaR[files[0]] , mDeltaR[files[1]])

outfile = root_open("outfile_ttbar_doubleloop.root", "recreate")
for hist in mDeltaR_allVsTen.values() :
    hist.Write()
outfile.Close()
