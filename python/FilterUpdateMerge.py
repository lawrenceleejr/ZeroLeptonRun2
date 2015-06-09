#!/bin/env python

import sys, os, copy

import ROOT


def parseCmdLine():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--input_sig", dest="inputfile_signal", 
                      help="List of input filenames for signal stored in a text file (ls <myfolder>/*output.root* > mylist.txt)", 
                      default="tree_signal.list")
    parser.add_option("--input_data", dest="inputfile_data", 
                      help="List of input filenames for data", 
                      default="tree_data.list")
    parser.add_option("--input_mc", dest="inputfile_mc", 
                      help="List of input filenames for mc", 
                      default="tree_mc.list")
    parser.add_option("--input_qcd", dest="inputfile_qcd", 
                      help="List of input filenames for qcd", 
                      default="tree_qcd.list")
    parser.add_option("--doSignal", dest="doSignal", 
                      help="do Signal (default=%default)", 
                      action='store_true', default=False)
    parser.add_option("--doQCD", dest="doQCD", 
                      help="do QCD (default=%default)", 
                      action='store_true', default=False)
    parser.add_option("--doData", dest="doData", 
                      help="do Data (default=%default)", 
                      action='store_true', default=False)
    parser.add_option("--doBackground", dest="doBackground", 
                      help="do Background (default=%default)", 
                      action='store_true', default=False)
    parser.add_option("--doNormWeight", dest="doNormWeight", 
                      help="(Re)Compute normWeight variable (default=%default)", 
                      action='store_true', default=False)
    parser.add_option("--dbName", dest="dbName", 
                      help="Weights database (default=%default)", 
                      default="ZeroLeptonRun2/data/MCBackgroundDB.dat")
    parser.add_option("--filter", dest="filter",
                      help="Select events as done in UpdateTree::SelectEvent ? (default=%default)", 
                      action='store_true', default=False)
    parser.add_option("--mergeExtraVars", dest="mergeExtraVars",
                      help="Also merge the ExtraVars block ?", 
                      action='store_true', default=True)
    parser.add_option("--mergeCRWTVars", dest="mergeCRWTVars",
                      help="Also merge the CRWTVars block ?", 
                      action='store_true', default=True)
    parser.add_option("--mergeCRZVars", dest="mergeCRZVars",
                      help="Also merge the CRZVars block ?", 
                      action='store_true', default=True)
    parser.add_option("--mergeCRYVars", dest="mergeCRYVars",
                      help="Also merge the CRYVars block ?", 
                      action='store_true', default=True)
    parser.add_option("--verbose", dest="verbose", type='int', 
                      help="Verbose level (0=minimum, default=%default)", default=0)
    parser.add_option("--prefix", dest="prefix", default="mc15_13TeV",
                      help="Prefix to identify mc production",) 

    (config, args) = parser.parse_args()
    return config

#------------------------------------------------------------------------

class Sample:
    def __init__(self, name, myarg, inputfile, config, 
                 Extraname="", treename="", nokfactor=False):
        self.name = name
        self.config = config
        self.files=[]
        if treename=="":
            self.treename=self.name
        else:            
            self.treename=treename
        self.extraName=Extraname
        self.signaldb=None
        self.nokfactor=nokfactor
        try:
            for file in open(inputfile,'read'):
                file.strip()
                file = file[:-1] # get rid of the \n
                if len(file) == 0: continue
                for id in myarg:
                    if not str(id) in file: continue
                    fname = file
                    if fname.startswith('/eos'): fname = 'root://eosatlas/'+fname
                    if self.emptyFile(fname):
                        print 'Empty file, ignored',fname
                    else:
                        self.files.append(fname)
        except:
            print "can't read file",treename,inputfile

    def emptyFile(self,fname):
        """ file must have at least 1D/2D counter histograms """
        ftemp = ROOT.TFile.Open(fname)
        if not ftemp or ftemp.IsZombie():
            print 'Could not open input file ',fname
            return True
        keys=ftemp.GetListOfKeys()
        for k in keys:
            kname=k.GetName()
            if not kname.startswith('Counter_'): continue
            obj=ftemp.Get(kname)
            if obj.IsA().InheritsFrom(ROOT.TH1D.Class()) or obj.IsA().InheritsFrom(ROOT.TH2D.Class()): 
                ftemp.Close()
                return False
        ftemp.Close()
        return True

    def getList(self):
        return self.files

    def getNewTreeName(self,kname,ds=None):
        """Rename tree removing NT and also SRAll for control regions."""
        if type(self.treename) == dict and ds in self.treename:
            newname=self.treename[ds]+"_"+kname
        else:
            newname=self.treename+"_"+kname
        newname=newname.replace("NT","")
        if newname.find("TopMcAtNlo")>=0: newname=newname.replace("McAtNlo","")+"_McAtNlo"
        if newname.find("TopPowheg")>=0: newname=newname.replace("Powheg","")+"_Powheg"
        if newname.find("WAlpgen")>=0: newname=newname.replace("Alpgen","")+"_Alpgen"
        if newname.find("WSherpa")>=0: newname=newname.replace("Sherpa","")+"_Sherpa"
        #if newname.find("ZAlpgen")>=0: newname=newname.replace("Alpgen","")+"_Alpgen"        
        if newname.find("CR")>=0 or newname.find("VR")>=0:
            newname=newname.replace("_SRAll","_")
        else:
            newname=newname.replace("SRAll","SRAll_")
        if newname[-1]=="_": newname=newname[:-1]
        newname=newname.replace("__","_")
        return newname

    def getListOfTrees(self,fname):
        """Find list of trees."""
        try:
            ftemp = ROOT.TFile.Open(fname)
        except:
            "ERROR can't open",fname
            return
        keys=ftemp.GetListOfKeys()
        mykeys=[]
        for k in keys:
            kname=k.GetName()            
            obj=ftemp.Get(kname)
            if not obj.IsA().InheritsFrom( ROOT.TTree.Class() ): continue
            mykeys.append(kname)
        ftemp.Close()
        return copy.deepcopy(mykeys)

    def extractChannelFromFilename(self,filename):
        for prefix in ('mc11_7TeV.', 'mc12_8TeV.', 'mc14_8TeV.','mc14_13TeV.', 'mc15_13TeV.'):
            if prefix in filename:
                return int(filename.split(prefix)[1].split('.')[0])
        return None

    def normWeight(self,filename,xsecDB,isSignal):
        if isSignal: return 1.
        if self.config.doData: return 1.

        ds = self.extractChannelFromFilename(filename)
        if not ds:
            print 'Could not identify channel number for ',filename
            sys.exit(1)
            
        xsec = xsecDB.xsectTimesEff(ds)
        sumW = 0.
        if ds >= 147910 and ds <= 147917:
            sumW = xsecDB.process(ds).stat()
        else:
            sumW = xsecDB.sumweight(ds)

        if sumW != 0: return xsec/sumW
        return 0.


    def MergeCounters(self,outFile,inFiles,xsecDB,isSignal):
        """Merge histrograms corresponding to counters"""
        if len(inFiles) == 0: return

        # get the list of histograms from the first file
        scale = 1.
        if self.config.doNormWeight: scale = self.normWeight(inFiles[0],xsecDB,isSignal)
        ftemp = ROOT.TFile.Open(inFiles[0])
        if not ftemp or ftemp.IsZombie():
            print 'Could not open input file ',inFiles[0]
            sys.exit(1)
        keys=ftemp.GetListOfKeys()
        hkeys=[]
        outhists = {}
        for k in keys:
            kname=k.GetName()
            if not kname.startswith('Counter_'): continue
            obj=ftemp.Get(kname)
            if not obj.IsA().InheritsFrom(ROOT.TH1D.Class()) and  not obj.IsA().InheritsFrom(ROOT.TH2D.Class()): continue
            hkeys.append(kname)
            outobj = obj.Clone()
            outobj.SetDirectory(outFile)
            outobj.Scale(scale)
            outhists[kname] = outobj
        ftemp.Close()
        del ftemp
        inFiles = inFiles[1:]

        # merge with histograms from other files
        for fname in inFiles:
            ftemp = ROOT.TFile.Open(fname)
            if not ftemp or ftemp.IsZombie():
                print 'Could not open input file ',fname
                sys.exit(1)
            scale = 1.
            if self.config.doNormWeight: scale = self.normWeight(fname,xsecDB,isSignal)
            for k in hkeys:
                h = ftemp.Get(k)
                if not h:
                    print 'Could not extract histogram',k,'in file',fname
                    sys.exit(1)
                outhists[k].Add(h,scale)
                pass
            ftemp.Close()
            del ftemp
        outFile.cd()
        for h in outhists.itervalues():
            if not isSignal:
                # number of events and sum weight no longer have meanings for merged and normalized counters
                h.SetBinContent(1,0.)
                h.SetBinContent(2,0.)
            h.Write()
        pass

    def ModifyAndMerge(self,mergename="",dbName="ZeroLeptonRun2/data/MCBackgroundDB.dat",isSignal=False,update=False):
        """Merge all input files and rename trees with modifications in NTVars."""
        if len(self.files) == 0: return
        # list of trees in the input files
        treenames = self.getListOfTrees(self.files[0])

        # create C++ merger
        xsecDB = ROOT.SUSY.CrossSectionDB(dbName, True)
        merger = ROOT.FilterUpdateMerge(xsecDB)

        # open output file
        if mergename=="": 
            mergename=self.name+self.extraName+".root"
        if update:
            newfile = ROOT.TFile.Open(mergename,"UPDATE")
        else:
            newfile = ROOT.TFile.Open(mergename,"RECREATE")

        self.MergeCounters(newfile, self.files,xsecDB,isSignal)

        # Build a list of output trees and files that should be merged into it
        # Needed for signal grid as we have one tree name per mc channel
        # while for background MC and signal we merge all files

        # key = output tree name  value = (TTree, input-tree-name, list-of-files)
        outTreeDict = {} 
        for inTree in treenames:
            for filename in self.files:
                ds = None
                if self.config.doQCD or self.config.doData: 
                    ds = 1
                else:
                    ds = self.extractChannelFromFilename(filename)
                if not ds:
                    print 'Could not identify channel number for ',filename
                    sys.exit(1)
                newtname=self.getNewTreeName(inTree, ds)
                if outTreeDict.has_key(newtname):
                    outTreeDict[newtname][2].append(filename)
                else:
                    outTree = ROOT.TTree(newtname,'0-lepton small ntuple')
                    outTreeDict[newtname] = (outTree, inTree, [filename])

        # now do the real merging
        for outTreeName, (outTree, inTreeName, filelist) in outTreeDict.iteritems():
            inList = ROOT.FilterUpdateMergeFileList()
            for name in filelist:
                inList.add(name)

            doExtraVars = self.config.mergeExtraVars
            doCRWTVars = self.config.mergeCRWTVars
            doCRZVars = self.config.mergeCRZVars
            doCRYVars = self.config.mergeCRYVars
            # check that the first file in the list has the requested block
            for name in filelist:
                testf = ROOT.TFile.Open(name)
                if not testf or testf.IsZombie(): continue
                testt = testf.Get(inTreeName)
                if not testt: continue
                if doExtraVars and not testt.GetBranch('NTExtraVars'):
                    doExtraVars = False
                if doCRWTVars and not testt.GetBranch('NTCRWTVars'):
                    doCRWTVars = False
                if doCRZVars and not testt.GetBranch('NTCRZVars'):
                    doCRZVars = False
                if doCRYVars and not testt.GetBranch('NTCRYVars'):
                    doCRYVars = False
                break
            merger.process(outTree, inTreeName, inList, isSignal, self.config.doNormWeight, self.config.filter, doExtraVars, doCRWTVars, doCRZVars, doCRYVars)
            newfile.cd()
            outTree.Write()

        newfile.Close()

        return

if __name__ == '__main__':
    config = parseCmdLine()
    #config = parseCmdLine(sys.argv[1:])

    # Load RootCore libs
    ROOT.gROOT.ProcessLine(".x $ROOTCOREDIR/scripts/load_packages.C")
    print 'Packages loaded'

    # List of dataset used in analysis
    if config.prefix == 'mc12_8TeV' or  config.prefix == 'mc14_8TeV':
        from mc12_8TeV_MCSampleList import *
    elif config.prefix == 'mc14_13TeV':
        from mc14_13TeV_MCSampleList import *
    elif config.prefix == 'mc15_13TeV':
        from mc15_13TeV_MCSampleList import *
    elif config.prefix == 'mc15_week1':
        from mc15_13TeV_week1_MCSampleList import *
    else:
        print 'Unsupported mc production type',config.prefix
        sys.exit(1)


    AllSamples=[]
    if config.doBackground==True:
        Top = Sample('Top',lTop,config.inputfile_mc,config,treename="Top")
        WMassiveCB = Sample('WMassiveCB',lWjets,config.inputfile_mc,config,treename="W")
        ZMassiveCB = Sample('ZMassiveCB',lZjets,config.inputfile_mc,config,treename="Z")
        DB = Sample('DibosonMassiveCB',lDiBoson,config.inputfile_mc,config,treename="Diboson")
        QCD = Sample('QCD',lQCDMC,config.inputfile_mc,config)
        GAMMAMassiveCB = Sample('GAMMAMassiveCB',lYjets,config.inputfile_mc,config,treename="GAMMA")

        AllSamples.append(Top)
        AllSamples.append(WMassiveCB)
        AllSamples.append(ZMassiveCB)
        AllSamples.append(QCD)
        AllSamples.append(GAMMAMassiveCB)


    if config.doData==True:
        # Run 1 streams
        DataJetTauEtmiss = Sample('Data',['JetTauEtmiss'],config.inputfile_data,config,Extraname="JetTauEtmiss")
        DataMuon = Sample('Data',['Muon'],config.inputfile_data,config,Extraname="Muon")
        DataEgamma = Sample('Data',['Egamma'],config.inputfile_data,config,Extraname="Egamma")
        # Run 2 stream
        DataMain = Sample('Data',['Main'],config.inputfile_data,config,Extraname="Main")
        AllSamples += [DataJetTauEtmiss,DataMuon,DataEgamma,DataMain]
        
    if config.doQCD==True:
        QCDdd1 = Sample('QCDdd',['JetTauEtmiss'],config.inputfile_qcd,config)
        AllSamples += [QCDdd1]

    # do merging
    for s in AllSamples:
        s.ModifyAndMerge(dbName=config.dbName)
        

    
    if config.doSignal==True:
        import pickle
        if not os.path.isfile('signalPointPickle.pkl'):
            print 'pickled signal grids not found, exiting...'
            sys.exit()

        picklefile = open('signalPointPickle.pkl','rb')
        pointdict = pickle.load(picklefile)
        picklefile.close()

        for grid,mapgrid in pointdict.items():
            myarg = []
            mytprefix = {}
            for key,info in mapgrid.items():
                name=grid+"_"+str(info[0])
                if(len(info)>=2):name+="_"+str(info[1])
                if(len(info)>=3):name+="_"+str(info[2])
                print name," ",key
                if type(key) == int:
                    myarg.append("."+str(key)+".")
                else:
                    myarg.append(str(key))
                mytprefix[key] = name
            signal = Sample(grid,myarg,config.inputfile_signal,config,treename=mytprefix)
            signal.ModifyAndMerge(mergename=grid+".root",dbName=config.dbName,isSignal=True,update=False)
