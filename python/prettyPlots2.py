import rootpy.ROOT as ROOT
from rootpy.io import root_open
#import rootpy.numpy
import root_numpy as rnp

import AtlasUtils
import AtlasStyle

import os
from os import path
# Workaround to fix threadlock issues with GUI
ROOT.PyConfig.StartGuiThread = False
import logging
logging.basicConfig(level=logging.INFO)

logging.info("loading packages")
#import shutil                                                                                                                                       \
#shutil.copyfile(ROOT.gSystem.ExpandPathName('$ROOTCOREDIR/scripts/load_packages.C'), 'load_packages.C')
#lineLoadPackages = '.x ' + ROOT.gSystem.ExpandPathName('$ROOTCOREDIR/scripts/load_packages.C')
#print lineLoadPackages
ROOT.gROOT.Macro("$ROOTCOREDIR/scripts/load_packages.C")
#ROOT.gROOT.ProcessLine(lineLoadPackages)

# Initialize the xAOD infrastructure
ROOT.xAOD.Init()

myfiles = {
    "" : root_open("outfile.ttbar.user.rsmith.6764481._000055.o.root")
}

jetPtCuts = { 30 : myfiles[""].pt30_vs_pt40.ProjectionX(),
              40 : myfiles[""].pt30_vs_pt40.ProjectionY(),
              50 : myfiles[""].pt30_vs_pt50.ProjectionY(),

    }

ROOT.SetAtlasStyle()
ROOT.gStyle.SetOptStat(0)

canvas1 = ROOT.TCanvas("MDR_jetPtCuts", "MDR_jetPtCuts", 600, 600)
canvas1.SetRightMargin(0.13);
canvas1.SetLogy()

jetPtCuts[30].SetTitle("")
jetPtCuts[30].GetYaxis().SetTitleOffset(1.4)
jetPtCuts[30].GetXaxis().SetTitle("HLT M_{#Delta}^{R} (GeV)")
jetPtCuts[30].GetXaxis().SetRangeUser(0,1500)
jetPtCuts[30].GetYaxis().SetTitle("Entries")


for color,value in enumerate(jetPtCuts.values()) :
#    value.Scale(1./value.Integral()) #a.u.
    value.SetMarkerColor(color+2 )
    print "jetPtCut " + str(color) + "  "  +  str(value.Integral())

jetPtCuts[30].Draw()
jetPtCuts[40].Draw("same")
jetPtCuts[50].Draw("same")

leg3 = ROOT.TLegend(.5, .5, 0.7 , 0.7)
leg3.AddEntry(jetPtCuts[30], "jet pt > 30 GeV" )
leg3.AddEntry(jetPtCuts[40], "jet pt > 40 GeV" )
leg3.AddEntry(jetPtCuts[50], "jet pt > 50 GeV" )
leg3.Draw("same")

#AtlasUtils.myText(.3,.85,ROOT.kBlack, "Ten Jets")
AtlasUtils.ATLAS_LABELInternal(.3,.75, ROOT.kBlack)
canvas1.SetGrid(1,1)

canvas2 = ROOT.TCanvas("nHLTJets", "nHLTJets", 600, 600)
canvas2.SetRightMargin(0.13);
canvas2.SetGrid(1,1)
canvas2.SetLogy()

hltjets = {
    ""     : myfiles[""].nHLTjets,
    "xe70" : myfiles[""].nHLTjets_xe70
}

hltjets[""].SetTitle("")
hltjets[""].GetYaxis().SetTitleOffset(1.4)
hltjets[""].GetXaxis().SetTitle("Number of HLT Jets")
hltjets[""].GetYaxis().SetTitle("a.u.")


for color,value in enumerate(hltjets.values()) :
    value.Scale(1./value.Integral()) #a.u.
    value.SetMarkerColor(color+2 )

hltjets[""]    .Draw()
hltjets["xe70"].Draw("same")

leg2 = ROOT.TLegend(.7, .5, 0.9 , 0.7)
leg2.AddEntry(hltjets[""], "No trigger selection" )
leg2.AddEntry(hltjets["xe70"], "HLT_xe70" )
leg2.Draw("same")

AtlasUtils.ATLAS_LABELInternal(.6,.75, ROOT.kBlack)

print "% with more than 10 jets " + str(hltjets["xe70"].Integral(hltjets["xe70"].FindBin(10.5) , 100000) / hltjets["xe70"].Integral())


canvas1.Print('plots/'+canvas1.GetName()+".eps")
canvas2.Print('plots/'+canvas2.GetName()+".eps")

import time
time.sleep(150)
