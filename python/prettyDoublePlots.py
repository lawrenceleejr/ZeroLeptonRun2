import rootpy.ROOT as ROOT
import rootpy.stl as stl
from rootpy.io import root_open
#import rootpy.numpy
import root_numpy as rnp

import AtlasUtils
import AtlasStyle

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
ROOT.SetAtlasStyle()
ROOT.gStyle.SetOptStat(0)

totalEventCounter = 0

treeName = "PassThroughNT"

rootfile = root_open("outfile_ttbar_doubleloop.root")

mDeltaR_allVsTen = {
   5  :  rootfile.mDeltaR_allVsTen_GT_5jets,
   7  :  rootfile.mDeltaR_allVsTen_GT_7jets,
   10 :  rootfile.mDeltaR_allVsTen_GT_10jets,
}

for njetcut in mDeltaR_allVsTen.keys() :
    canvas1 = ROOT.TCanvas("c1"+str(njetcut),"c1"+str(njetcut),600, 600)
    canvas1.SetRightMargin(0.13);
    canvas1.cd()

    canvas1.SetLogz()
    leg = ROOT.TLegend(.2, .2, .5, .5)
    leg.AddEntry(mDeltaR_allVsTen[njetcut])
    leg.Draw()
    mDeltaR_allVsTen[njetcut].SetTitle("")

    mDeltaR_allVsTen[njetcut].Draw("colz")
    mDeltaR_allVsTen[njetcut].GetXaxis().SetTitle("Offline M_{#Delta}^{R}(all jets) (GeV)")
    mDeltaR_allVsTen[njetcut].GetXaxis().SetNdivisions(505)
    mDeltaR_allVsTen[njetcut].GetYaxis().SetNdivisions(505)
    mDeltaR_allVsTen[njetcut].GetYaxis().SetTitle("Offline M_{#Delta}^{R}(ten jets) (GeV)")
    AtlasUtils.myText(.3,.85,ROOT.kBlack, "N HLT Jets > " + str(njetcut))
    canvas1.SetGrid(1,1)

    canvas1.SaveAs("plots/2d_alljet_tenjet_njetcut"+str(njetcut)+".eps")

    canvas2 = ROOT.TCanvas("c2"+str(njetcut),"c2"+str(njetcut),600,600)
    canvas2.SetRightMargin(0.13);
    canvas2.cd()

    canvas2.SetLogy()
    alljet = mDeltaR_allVsTen[njetcut].ProjectionX()
    tenjet = mDeltaR_allVsTen[njetcut].ProjectionY()
    alljet.Draw()
    alljet.GetYaxis().SetTitle("Entries")
    tenjet.SetMarkerColor(ROOT.kRed)
    tenjet.Draw("same")
    AtlasUtils.myText(.3,.85,ROOT.kBlack, "N HLT Jets > " + str(njetcut))
    canvas2.SetGrid(1,1)
    canvas2.SaveAs("plots/oneAxis_"+str(njetcut)+".eps")

import time
time.sleep(120)

print "Total number of events in the sample:" , totalEventCounter
