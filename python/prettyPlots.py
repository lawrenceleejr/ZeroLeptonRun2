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

files = {
    "tenjet" : root_open("outfile.ttbar.ttbar_tenjet.root"),
    "alljet" : root_open("outfile.ttbar.ttbar_alljet.root"),
}

onlineMDR_vs_offlineMDR = {
    "tenjet" : files['tenjet'].onlineMDR_vs_offlineMDR,
    "alljet" : files['alljet'].onlineMDR_vs_offlineMDR,
}
eff_xe10_razor170_off = {
    "tenjet" : files['tenjet'].eff_xe10_razor170_off,
    "alljet" : files['alljet'].eff_xe10_razor170_off,
}
eff_xe10_razor170_off_metcut = {
    "tenjet" : files['tenjet'].eff_xe10_razor170_off_metcut,
    "alljet" : files['alljet'].eff_xe10_razor170_off_metcut,
}

eff_xe10_razor170_0L_off = {
    "tenjet" : files['tenjet'].eff_xe10_razor170_0L_off,
    "alljet" : files['alljet'].eff_xe10_razor170_0L_off,
}

ROOT.SetAtlasStyle()
ROOT.gStyle.SetOptStat(0)

canvas1 = ROOT.TCanvas("tenjet2d", "tenjet2d", 600, 600)
canvas1.SetRightMargin(0.13);

onlineMDR_vs_offlineMDR["tenjet"].Draw("colz")
onlineMDR_vs_offlineMDR["tenjet"].SetTitle("")
onlineMDR_vs_offlineMDR["tenjet"].GetYaxis().SetTitleOffset(1.4)
onlineMDR_vs_offlineMDR["tenjet"].GetXaxis().SetTitle("Offline M_{#Delta}^{R} (GeV)")
onlineMDR_vs_offlineMDR["tenjet"].GetYaxis().SetTitle("HLT M_{#Delta}^{R} (GeV) ")
onlineMDR_vs_offlineMDR["tenjet"].Draw("colz")
AtlasUtils.myText(.3,.85,ROOT.kBlack, "Ten Jets")
AtlasUtils.ATLAS_LABELInternal(.3,.75, ROOT.kBlack)
canvas1.SetGrid(1,1)

canvas2 = ROOT.TCanvas("alljet2d", "alljet2d", 600, 600)
canvas2.SetRightMargin(0.13);
onlineMDR_vs_offlineMDR["alljet"].Draw("colz")
onlineMDR_vs_offlineMDR["alljet"].SetTitle("")
onlineMDR_vs_offlineMDR["alljet"].GetXaxis().SetTitle("Offline M_{#Delta}^{R} (GeV)")
onlineMDR_vs_offlineMDR["alljet"].GetYaxis().SetTitleOffset(1.4)
onlineMDR_vs_offlineMDR["alljet"].GetYaxis().SetTitle("HLT M_{#Delta}^{R} (GeV) ")
onlineMDR_vs_offlineMDR["alljet"].Draw("colz")
AtlasUtils.myText(.3,.85,ROOT.kBlack, "All Jets")
AtlasUtils.ATLAS_LABELInternal(.3,.75, ROOT.kBlack)
canvas2.SetGrid(1,1)

canvas3 = ROOT.TCanvas("eff_nometcut", "eff_nometcut" , 600 , 600)
canvas3.SetRightMargin(0.13);
canvas3.SetGrid(1,1)

eff_xe10_razor170_off["alljet"].Draw()
eff_xe10_razor170_off["tenjet"].Draw()
eff_xe10_razor170_off["alljet"].SetTitle("")

eff_xe10_razor170_off["alljet"].GetTotalHistogram().GetXaxis().SetTitle("Offline M_{#Delta}^{R} (GeV)")
eff_xe10_razor170_off["alljet"].GetTotalHistogram().GetXaxis().SetLabelSize(.02)

eff_xe10_razor170_off["alljet"].GetTotalHistogram().GetYaxis().SetTitle("Efficiency wrt L1 seed")
eff_xe10_razor170_off["alljet"].GetTotalHistogram().GetYaxis().SetTitleOffset(1.4)

eff_xe10_razor170_off["alljet"].Draw()
eff_xe10_razor170_off["tenjet"].SetMarkerColor(ROOT.kRed)
eff_xe10_razor170_off["tenjet"].Draw("same")

#leg3 = ROOT.TLegend(.6, 0.2, 0.8 , 0.4)
leg3 = ROOT.TLegend(.2, .3, 0.4 , 0.5)
leg3.AddEntry(eff_xe10_razor170_off["alljet"] , "all jets")
leg3.AddEntry(eff_xe10_razor170_off["tenjet"] , "ten jets")
leg3.Draw("same")
AtlasUtils.ATLAS_LABELInternal(.2,.2,ROOT.kBlack)

canvas4 = ROOT.TCanvas("mDeltaR", "mDeltaR" , 600 , 600)
canvas4.Draw()

pad1 = ROOT.TPad("pad1","pad1",0. ,0.3, 1.,1.0  )
pad1.SetBottomMargin(0);
pad1.SetGrid(1,1)
pad1.Draw()
pad1.cd()
hlt_MDR_tenjet = onlineMDR_vs_offlineMDR["tenjet"].ProjectionY()
hlt_MDR_alljet = onlineMDR_vs_offlineMDR["alljet"].ProjectionY()
hlt_MDR_alljet.SetMarkerColor(ROOT.kRed)

hlt_MDR_tenjet.Draw()
hlt_MDR_alljet.Draw("same")
leg4 = ROOT.TLegend(.5, 0.2, 0.7 , 0.4)
leg4.AddEntry(eff_xe10_razor170_off["alljet"] , "all jets")
leg4.AddEntry(eff_xe10_razor170_off["tenjet"] , "ten jets")
leg4.Draw("same")
AtlasUtils.ATLAS_LABELInternal(.6,.8,ROOT.kBlack)

canvas4.cd()
pad2 = ROOT.TPad("pad2","pad2",0. ,0.0, 1.,0.3  )
pad2.SetTopMargin(0);
pad2.SetGrid(1,1)
pad2.Draw()
pad2.cd()

ratio = hlt_MDR_tenjet.Clone()
ratio.GetYaxis().SetTitle("ten jet / all jet");
ratio.GetYaxis().SetNdivisions(505);
ratio.GetYaxis().SetTitleSize(20);
ratio.GetYaxis().SetTitleFont(43);
ratio.GetYaxis().SetTitleOffset(1.55);
ratio.GetYaxis().SetLabelFont(43);
ratio.GetYaxis().SetLabelSize(15);
ratio.GetXaxis().SetTitleSize(20);
ratio.GetXaxis().SetTitleFont(43);
ratio.GetXaxis().SetTitleOffset(4.);
ratio.GetXaxis().SetLabelFont(43);
ratio.GetXaxis().SetLabelSize(15);

ratio.SetMinimum(.9999)
ratio.SetMaximum(1.0001)
ratio.Sumw2()

ratio.Divide(hlt_MDR_alljet)
ratio.Draw()

canvas5 = ROOT.TCanvas("eff_0L", "eff_0L" , 600 , 600)
canvas5.SetRightMargin(0.13);
canvas5.SetGrid(1,1)

eff_xe10_razor170_0L_off["alljet"].Draw()
eff_xe10_razor170_0L_off["tenjet"].Draw()
eff_xe10_razor170_0L_off["alljet"].SetTitle("")

eff_xe10_razor170_0L_off["alljet"].GetTotalHistogram().GetXaxis().SetTitle("Offline M_{#Delta}^{R} (GeV)")
eff_xe10_razor170_0L_off["alljet"].GetTotalHistogram().GetXaxis().SetLabelSize(.02)

eff_xe10_razor170_0L_off["alljet"].GetTotalHistogram().GetYaxis().SetTitle("Efficiency wrt L1 seed")
eff_xe10_razor170_0L_off["alljet"].GetTotalHistogram().GetYaxis().SetTitleOffset(1.4)

eff_xe10_razor170_0L_off["alljet"].Draw()
eff_xe10_razor170_0L_off["tenjet"].SetMarkerColor(ROOT.kRed)
eff_xe10_razor170_0L_off["tenjet"].Draw("same")

#leg3 = ROOT.TLegend(.6, 0.2, 0.8 , 0.4)
leg5 = ROOT.TLegend(.2, .3, 0.4 , 0.5)
leg5.AddEntry(eff_xe10_razor170_0L_off["alljet"] , "all jets")
leg5.AddEntry(eff_xe10_razor170_0L_off["tenjet"] , "ten jets")
leg5.Draw("same")
AtlasUtils.ATLAS_LABELInternal(.2,.2,ROOT.kBlack)
AtlasUtils.myText(.2,.6,ROOT.kBlack, "0L selection")


canvas1.Print('plots/'+canvas1.GetName()+".eps")
canvas2.Print('plots/'+canvas2.GetName()+".eps")
canvas3.Print('plots/'+canvas3.GetName()+".eps")
canvas4.Print('plots/'+canvas4.GetName()+".eps")
canvas5.Print('plots/'+canvas5.GetName()+".eps")

import time
time.sleep(60)
