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

datadirs =[
#"/data/users/rsmith/razor_trigger/user.rsmith.trig.v7.alljetRJ.mc15_13TeV.370911.MadGraphPythia8EvtGen_GG_direct_800_600.SUSY1.e3962_a766_a777_r6282_p2419_o.root",
#"/data/users/rsmith/razor_trigger/user.rsmith.trig.v7.tenjetRJ.mc15_13TeV.370911.MadGraphPythia8EvtGen_GG_direct_800_600.SUSY1.e3962_a766_a777_r6282_p2419_o.root",
#"/data/users/rsmith/razor_trigger/user.rsmith.trig.v7.alljetRJ.mc15_13TeV.370938.MadGraphPythia8EvtGen_GG_direct_1200_200.SUSY1.e3962_a766_a777_r6282_p2419_o.root",
     #"/data/users/rsmith/razor_trigger/user.rsmith.trig.v7.tenjetRJ.mc15_13TeV.370938.MadGraphPythia8EvtGen_GG_direct_1200_200.SUSY1.e3962_a766_a777_r6282_p2419_o.root",
#     "/data/users/rsmith/razor_trigger/user.rsmith.trig.v7.tenjetRJ.mc15_13TeV.410000.PowHPEvG_ttbar_nonallhad.SUSY1.e3698_s2608_s2183_r6765_r6282_p2419_o.root/",
#     "/data/users/rsmith/razor_trigger/user.rsmith.trig.v7.alljetRJ.mc15_13TeV.410000.PowHPEvG_ttbar_nonallhad.SUSY1.e3698_s2608_s2183_r6765_r6282_p2419_o.root",
#"/data/users/rsmith/razor_trigger_data/user.rsmith.trig.v13.tenjetRJ.data15_13TeV.00279984.physics_Main.merge.DAOD_SUSY1.f629_m1504_p2425_o.root",
#"/data/users/rsmith/razor_trigger_data/user.rsmith.trig.v13.alljetRJ.data15_13TeV.00279984.physics_Main.merge.DAOD_SUSY1.f629_m1504_p2425_o.root"
#"/afs/cern.ch/work/r/rsmith/ttbar_trigger/ttbar_alljet.root",
#"/afs/cern.ch/work/r/rsmith/data_trigger_test/"
"/afs/cern.ch/work/r/rsmith/data_trigger/"
]

# def findPt( mydict,  ptCut ) :
#       return (v for (k,v) in mDeltaR.iteritems() if ptCut in k).next()

for datadir in datadirs :
     files = []
     #files.append(datadir)
     print "searching " +datadir
     for root, _, filenames in os.walk(datadir):
          for filename in filenames:
              files.append(os.path.join(root, filename))
     print files

     mDeltaR_hlt_jetpt               = { "pt30_vs_pt40" : ROOT.TH2F("pt30_vs_pt40", "pt30_vs_pt40", 100 , 0 , 3000 , 100 , 0 , 3000),
                                         "pt30_vs_pt50" : ROOT.TH2F("pt30_vs_pt50", "pt30_vs_pt50", 100 , 0 , 3000 , 100 , 0 , 3000),
                                         "pt40_vs_pt50" : ROOT.TH2F("pt40_vs_pt50", "pt40_vs_pt50", 100 , 0 , 3000 , 100 , 0 , 3000),
                                         "pt40_vs_pt100" : ROOT.TH2F("pt40_vs_pt100", "pt40_vs_pt100", 100 , 0 , 3000 , 100 , 0 , 3000),

          }

     eff_xe10_razor170_off        = ROOT.TEfficiency("eff_xe10_razor170_off" , "xe10_razor170_eff_off; mDeltaR (GeV) ; wrt L1 seed", 100 , 0 , 3000)
     eff_xe10_razor170_off_metcut = ROOT.TEfficiency("eff_xe10_razor170_off_metcut" , "xe10_razor170_eff_off_metcut; mDeltaR (GeV) ; wrt L1 seed", 100 , 0 , 3000)
     eff_xe10_razor170_off_0L     = ROOT.TEfficiency("eff_xe10_razor170_0L_off" , "xe10_razor170_0L_off; mDeltaR (GeV) ; wrt L1 seed", 100 , 0 , 3000)
     njets          = ROOT.TH1F("njet", "njet", 13 , -.5 , 12.5 )

     nHLTJets       ={ 30: ROOT.TH1F("nHLTjets_30"     , "nHLTjets_30"     , 25 , -.5 , 24.5 ),
                       40: ROOT.TH1F("nHLTjets_40"     , "nHLTjets_40"     , 25 , -.5 , 24.5 ),
                       50: ROOT.TH1F("nHLTjets_50"     , "nHLTjets_50"     , 25 , -.5 , 24.5 ),
                       }

     nHLTJets_L1       ={ 30: ROOT.TH1F("nHLTjets_L1_30"     , "nHLTjets_L1_30"     , 25 , -.5 , 24.5 ),
                          40: ROOT.TH1F("nHLTjets_L1_40"     , "nHLTjets_L1_40"     , 25 , -.5 , 24.5 ),
                          50: ROOT.TH1F("nHLTjets_L1_50"     , "nHLTjets_L1_50"     , 25 , -.5 , 24.5 ),
                       }
     nHLTJets_xe70       ={ 30: ROOT.TH1F("nHLTjets_xe70_30"     , "nHLTjets_xe70_30"     , 25 , -.5 , 24.5 ),
                       40: ROOT.TH1F("nHLTjets_xe70_40"     , "nHLTjets_xe70_40"     , 25 , -.5 , 24.5 ),
                       50: ROOT.TH1F("nHLTjets_xe70_50"     , "nHLTjets_xe70_50"     , 25 , -.5 , 24.5 ),
                       }

     for filecount, ifile in enumerate(files):
     #    print ifile.isfile
          totalEventCounter = 0
          rootfile = root_open(ifile)
          treeName = "PassThroughNT"

          if(rootfile.GetListOfKeys().Contains(treeName)) :
             tree = rootfile.PassThroughNT
             nentries = tree.GetEntries()

             totalEventCounter += nentries
             print nentries

             counter = 0 #todo fix
             for event in tree:
                 counter = counter + 1

                 if ( counter % 1000 == 0 ) : print counter
     #             RJvars = event.NTRJigsawVars #ROOT.NTRJigsawVars()

                 triggerBits = int( event.triggerBits.at(0))
     #            print triggerBits
     #             print "triggerBits size : " + str(triggerBits.size())
                 triggers = [#this is ordered!
                 "L1_XE50",
                 "L1_XE70",
                 "HLT_xe70",
                 #    "HLT_xe70_pueta",
                 "HLT_xe100",
                 #    "HLT_xe100_pueta",
                 "HLT_e28_tight_iloose",
                 "HLT_e60_medium",
                 "HLT_mu26_imedium",
                 #    "HLT_j30_xe10_razor170",
                 "HLT_xe70_tc_em",
                 "HLT_xe70_tc_lcw",
                 #    "HLT_xe70_mht",
                 #    "HLT_xe70_pufit",
                 "HLT_xe100_tc_em",
                 "HLT_xe100_tc_lcw",
                 #    "HLT_xe100_mht",
                 #    "HLT_xe100_pufit",
                 "HLT_3j175",
                 "HLT_4j85",
                 "HLT_5j85",
                 "HLT_6j25",
                 "HLT_6j45_0eta240",
                 #    "HLT_6j55_0eta240_L14J20",
                 "HLT_7j45",
                 "L1_2J15",
                 "HLT_2j55_bloose",
                 "HLT_j80_xe80",
                 "HLT_e24_lhmedium_iloose_L1EM18VH",
                 "HLT_e60_lhmedium",
                 "HLT_mu20_iloose_L1MU15",
                 "HLT_mu40",
                 "HLT_mu50",
                 "HLT_g120_loose",
                 "HLT_g120_lhloose",
                 "HLT_mu18",
                 "HLT_e17_lhloose_L1EM15",
                 "HLT_e17_loose_L1EM15",
                 "HLT_mu14_iloose",
                 "HLT_j30_xe10_razor100",
                 "HLT_j30_xe10_razor170",
                 "HLT_j30_xe10_razor185",
                 "HLT_j30_xe10_razor195",
                 "HLT_j30_xe60_razor100",
                 "HLT_j30_xe60_razor170",
                 "HLT_j30_xe60_razor185",
                 "HLT_j30_xe60_razor195",
                 "L1_2J15_XE55",
                 ];

                 triggerDict = {}

                 for count, trig in enumerate(triggers) :
                     triggerDict[trig] = ( (triggerBits&(1<<count) ) > 0.)

                 listOfBranches = tree.GetListOfBranches()

                 missingEt = tree.GetLeaf("met").GetValue(0);
                 njets.Fill( event.jetPt.size())




                 nlepton = event.elPt.size() + event.muPt.size()

                 mDeltaR = {}
                 njets_hlt = {}# tree.GetLeaf("nHLTJets").GetValue(0)
                 for branch in listOfBranches :
                      if( "NTRJigsaw" in branch.GetName() ) :
                           mDeltaR[branch.GetName()]   = branch.GetLeaf("RJVars_PP_MDeltaR").GetValue(0)/1000.
                           njets_hlt[branch.GetName()] = branch.GetLeaf("RJVars_V1_N").GetValue(0) + branch.GetLeaf("RJVars_V2_N").GetValue(0)

                 # if( tree.GetLeaf("nHLTJets").GetValue(0) > .5 ) :
                 #      print "from HLT jet counting : " +  str(tree.GetLeaf("nHLTJets").GetValue(0))
                 #      print "from RJ vars          : " +  str(njets_hlt["NTRJigsawVars_hlt_jetPt30"])

                 nHLTJets[30].Fill(njets_hlt["NTRJigsawVars_hlt_jetPt30"])
                 nHLTJets[40].Fill(njets_hlt["NTRJigsawVars_hlt_jetPt40"])
                 nHLTJets[50].Fill(njets_hlt["NTRJigsawVars_hlt_jetPt50"])

                 if(triggerDict["L1_2J15_XE55"]) :
                      nHLTJets_L1[30].Fill(njets_hlt["NTRJigsawVars_hlt_jetPt30"])
                      nHLTJets_L1[40].Fill(njets_hlt["NTRJigsawVars_hlt_jetPt40"])
                      nHLTJets_L1[50].Fill(njets_hlt["NTRJigsawVars_hlt_jetPt50"])
                 if(triggerDict["HLT_xe70"]) :
                      nHLTJets_xe70[30].Fill(njets_hlt["NTRJigsawVars_hlt_jetPt30"])
                      nHLTJets_xe70[40].Fill(njets_hlt["NTRJigsawVars_hlt_jetPt40"])
                      nHLTJets_xe70[50].Fill(njets_hlt["NTRJigsawVars_hlt_jetPt50"])



                 mDeltaR_hlt_jetpt["pt30_vs_pt50"].Fill(mDeltaR["NTRJigsawVars_hlt_jetPt30"], mDeltaR["NTRJigsawVars_hlt_jetPt50"])
                 mDeltaR_hlt_jetpt["pt30_vs_pt40"].Fill(mDeltaR["NTRJigsawVars_hlt_jetPt30"], mDeltaR["NTRJigsawVars_hlt_jetPt40"])
                 mDeltaR_hlt_jetpt["pt40_vs_pt50"].Fill(mDeltaR["NTRJigsawVars_hlt_jetPt40"], mDeltaR["NTRJigsawVars_hlt_jetPt50"])
                 mDeltaR_hlt_jetpt["pt40_vs_pt100"].Fill(mDeltaR["NTRJigsawVars_hlt_jetPt40"], mDeltaR["NTRJigsawVars_hlt_jetPt100"])

     print "nentries :"  + str(eff_xe10_razor170_off.GetTotalHistogram().GetEntries())
     rootfile.Close()
     outfile  = root_open('outfile.ttbar.'+ifile.split('/')[-1].replace('.root','')+'.root', 'recreate')
     outfile.cd()
     for value in mDeltaR_hlt_jetpt.values() :
          value.Write()
     njets.Write()

     for key in nHLTJets.keys() :
          nHLTJets     [key].Write()
          nHLTJets_L1  [key].Write()
          nHLTJets_xe70[key].Write()

     eff_xe10_razor170_off.Write()
     eff_xe10_razor170_off_metcut.Write()
     eff_xe10_razor170_off_0L.Write()
#     onlineMDR_vs_offlineMDR_min8jet.Write()
#     onlineMDR_vs_offlineMDR.Write()
     outfile .Close()

print "Total number of events in the sample:" , totalEventCounter
