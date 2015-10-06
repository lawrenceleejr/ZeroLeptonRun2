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
#"/data/users/rsmith/razor_trigger/user.rsmith.trig.v7.tenjetRJ.mc15_13TeV.370911.MadGraphPythia8EvtGen_GG_direct_800_600.SUSY1.e3962_a766_a777_r6282_p2419_o.root"
#"/data/users/rsmith/razor_trigger/user.rsmith.trig.v7.alljetRJ.mc15_13TeV.370938.MadGraphPythia8EvtGen_GG_direct_1200_200.SUSY1.e3962_a766_a777_r6282_p2419_o.root",
#"/data/users/rsmith/razor_trigger/user.rsmith.trig.v7.tenjetRJ.mc15_13TeV.370938.MadGraphPythia8EvtGen_GG_direct_1200_200.SUSY1.e3962_a766_a777_r6282_p2419_o.root",
 "/data/users/rsmith/razor_trigger/user.rsmith.trig.v7.tenjetRJ.mc15_13TeV.410000.PowHPEvG_ttbar_nonallhad.SUSY1.e3698_s2608_s2183_r6765_r6282_p2419_o.root/",
 "/data/users/rsmith/razor_trigger/user.rsmith.trig.v7.alljetRJ.mc15_13TeV.410000.PowHPEvG_ttbar_nonallhad.SUSY1.e3698_s2608_s2183_r6765_r6282_p2419_o.root/"
]
files = []

for datadir in datadirs :
     for root, _, filenames in os.walk(datadir):
          for filename in filenames:
              files.append(os.path.join(root, filename))
     print files

     totalEventCounter = 0
     eff_xe10_razor170 = ROOT.TEfficiency("eff_xe10_razor170" , "xe10_razor170_eff", 100 , 0 , 3000)

     for ifile in files:
     #    print ifile.isfile
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
     # print count
                     triggerDict[trig] = ( (triggerBits&(1<<count) ) > 0.)

                 branch = tree.GetBranch( "NTRJigsawVars" );
                 mDeltaR = branch.GetLeaf("RJVars_PP_MDeltaR").GetValue(0);
                 if tree.jetPt.size() > 8 : print "jet size : " + str(tree.jetPt.size())

                 if(triggerDict["L1_XE50"]) :
                     eff_xe10_razor170.Fill(triggerDict["HLT_j30_xe10_razor170"], mDeltaR/1000.)

             #        outfile .Write()
         rootfile.Close()

     print "exiting early"
     exit()

     outfile  = root_open('outfile'+ifile.split('v7')[1].split("mc15")[0]+'root', 'recreate')
     outfile.cd()
     eff_xe10_razor170.Write()
     outfile .Close()


print "Total number of events in the sample:" , totalEventCounter
