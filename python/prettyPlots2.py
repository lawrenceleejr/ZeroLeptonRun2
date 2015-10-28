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

leg3 = ROOT.TLegend(.4, .4, 0.7 , 0.7)
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


from collections import OrderedDict
hltjets=OrderedDict()

#   hltjets["30"]    = myfiles[""].nHLTjets_30,
#hltjets["30_xe70"]= myfiles[""].nHLTjets_xe70_30
hltjets["30_L1"]  = myfiles[""].nHLTjets_L1_30

#   hltjets["40"]    = myfiles[""].nHLTjets_40,
#hltjets["40_xe70"]= myfiles[""].nHLTjets_xe70_40
hltjets["40_L1"]  = myfiles[""].nHLTjets_L1_40

#   hltjets["50"]    = myfiles[""].nHLTjets_50,
#hltjets["50_xe70"]= myfiles[""].nHLTjets_xe70_50
hltjets["50_L1"]  = myfiles[""].nHLTjets_L1_50


firstkey = hltjets.keys()[0]

hltjets[firstkey].SetTitle("")
hltjets[firstkey].GetYaxis().SetTitleOffset(1.4)
hltjets[firstkey].GetXaxis().SetTitle("Number of HLT Jets")
hltjets[firstkey].GetYaxis().SetTitle("a.u.")


for color,value in enumerate(hltjets.values()) :
    value.Scale(1./value.Integral()) #a.u.
    value.SetMarkerColor(color+2 )

hltjets[firstkey]    .Draw()
leg2 = ROOT.TLegend(.5, .4, 0.8 , 0.7)

for key , value in hltjets.items() :
    value.Draw("same")
    if ( "xe70" in key ) :
        leg2.AddEntry(value , "HLT_xe70, jet pt > " + key.split("_")[0] + " GeV ")
    elif ( "L1" in key ) :
        leg2.AddEntry(value , "L1 seed , jet pt > " + key.split("_")[0]+ " GeV ")
    else :
        leg2.AddEntry(value , "No trigger selection " + key)

leg2.Draw("same")

AtlasUtils.ATLAS_LABELInternal(.5,.75, ROOT.kBlack)


for key , value in hltjets.items() :
    print "% with more than 10 jets for " + key + " : " + str(value.Integral(value.FindBin(10.5) , 100000) / value.Integral())


canvas1.Print('plots/'+canvas1.GetName()+".eps")
canvas2.Print('plots/'+canvas2.GetName()+".eps")

import time
time.sleep(150)
