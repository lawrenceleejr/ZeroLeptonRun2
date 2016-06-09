#!/bin/env python

import sys, os, copy, re, time

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
    parser.add_option("--skipNormWeight", dest="skipNormWeight",
                      help="Don't (re)Compute normWeight variable (default=%default)",
                      action='store_true', default=False)
    parser.add_option("--dbName", dest="dbName",
                      help="Weights database (default=%default)",
                      default="ZeroLeptonRun2/data/MCBackgroundDB.dat")
    parser.add_option("--filter", dest="filter",
                      help="Select events as done in UpdateTree::SelectEvent ? (default=%default)",
                      action='store_true', default=True)
    parser.add_option("--skipExtraVars", dest="skipExtraVars",
                      help="Don't merge the ExtraVars block ?",
                      action='store_true', default=False)
    parser.add_option("--skipCRWTVars", dest="skipCRWTVars",
                      help="Don't merge the CRWTVars block ?",
                      action='store_true', default=False)
    parser.add_option("--skipCRZVars", dest="skipCRZVars",
                      help="Don't merge the CRZVars block ?",
                      action='store_true', default=False)
    parser.add_option("--skipCRYVars", dest="skipCRYVars",
                      help="Don't merge the CRYVars block ?",
                      action='store_true', default=False)
    parser.add_option("--mergeRJigsawVars", dest="mergeRJigsawVars",
                      help="Also merge the RJigsawVars block ?",
                      action='store_true', default=False)
    parser.add_option("--allBackgroundsNT", dest="allBackgroundsNT",
                      help="merge all BGs",
                      action='store_true', default=False)
    parser.add_option("--DBNT", dest="DBNT",
                      help="Merge for Dibosons ?",
                      action='store_true', default=False)
    parser.add_option("--QCDNT", dest="QCDNT",
                      help="Merge for QCD MC ?",
                      action='store_true', default=False)
    parser.add_option("--SplitQCDNT", dest="SplitQCDNT",
                      help="Merge for QCD MC process per process ?",
                      action='store_true', default=False)
    parser.add_option("--TopNT", dest="TopNT",
                      help="Merge for Top ?",
                      action='store_true', default=False)
    parser.add_option("--WNT", dest="WNT",
                      help="Merge for W+jets ?",
                      action='store_true', default=False)
    parser.add_option("--ZNT", dest="ZNT",
                      help="Merge for Z+jets ?",
                      action='store_true', default=False)
    parser.add_option("--GAMNT", dest="GAMNT",
                      help="Merge for gamma+jets ?",
                      action='store_true', default=False)
    parser.add_option("--WMadP8NT", dest="WMadP8NT",
                      help="Merge for W+jets Madgraph+P8 ?",
                      action='store_true', default=False)
    parser.add_option("--ZMadP8NT", dest="ZMadP8NT",
                      help="Merge for Z+jets Madgraph+P8 ?",
                      action='store_true', default=False)
    parser.add_option("--TopPowP8NT", dest="TopPowP8NT",
                      help="Merge for Top Powheg+P8?",
                      action='store_true', default=False)
    parser.add_option("--TopPowHNT", dest="TopPowHNT",
                      help="Merge for Top Powheg+Herwig ?",
                      action='store_true', default=False)
    parser.add_option("--TopMCNLONT", dest="TopMCNLONT",
                      help="Merge for Top MC@NLO?",
                      action='store_true', default=False)
    parser.add_option("--whichTrees", dest="treePatterns", default="",
                      help="comma separated python pattern to match tree names selected for merging",)
    parser.add_option("--excludeTrees", dest="treeExcludePatterns", default="",
                      help="comma separated python pattern to match tree names to exclude from merging",)
    parser.add_option("--verbose", dest="verbose", type='int',
                      help="Verbose level (0=minimum, default=%default)", default=0)
    parser.add_option("--prefix", dest="prefix", default="mc15_13TeV",
                      help="Prefix to identify mc production",)
    parser.add_option("--outputDir", dest="outputDir", default="./",
                      help="Output directory location",)

    (config, args) = parser.parse_args()
    return config

#------------------------------------------------------------------------

class Sample:
    def __init__(self, name, myarg, inputfile, config,
                 Extraname="", treename="", nokfactor=False):
        self.name = name
        self.config = config
        if treename=="":
            self.treename=self.name
        else:
            self.treename=treename
        self.extraName=Extraname
        self.signaldb=None
        self.nokfactor=nokfactor

        self.treePatterns = None
        if len(config.treePatterns)>0:
            self.treePatterns = config.treePatterns.split(',')

        self.treeExcludePatterns = None
        if len(config.treeExcludePatterns)>0:
            self.treeExcludePatterns = config.treeExcludePatterns.split(',')

        # if the input is a pickle file, it is supposed to contain a dictionary
        # of file lists indexed by the dataset name. In this case the channel
        # number should be extracted from the dataset name and not the files
        # name since in general file URL in rucio do not contain the channel nb
        if inputfile.endswith('.pkl'):
            self.files={}
            import pickle
            if not os.path.isfile(inputfile):
                print "can't read file",treename,inputfile
                sys.exit()
            picklefile = open(inputfile,'rb')
            myfiles = pickle.load(picklefile)
            picklefile.close()

            for (did,pfns) in myfiles.iteritems():
                for id in myarg:
                    if not str(id) in did: continue
                    self.files[did] = pfns
            pass
        else:
            self.files=[]
            try:
                for myfile in open(inputfile,'read'):
                    myfile.strip()
                    myfile = myfile[:-1] # get rid of the \n
                    print myfile
                    if len(myfile) == 0: continue
                    for id in myarg:
                        if not str(id) in myfile: continue
                        fname = myfile
                        if fname.startswith('/eos'): fname = 'root://eosatlas/'+fname
                        if self.emptyFile(fname):
                            print 'Empty file, ignored',fname
                        else:
                            self.files.append(fname)
            except:
                print sys.exc_info()[0]
                print "can't read file",inputfile,'for',treename
            pass

    def OpenWithRetries(self,fname):
        # only retry remote files
        if ( ':' in fname ):
            for delay in [10,30,120,300,1200]:
                f = ROOT.TFile.Open(fname,'READ')
                if f and not f.IsZombie(): return f
                print 'Problem opening file, wait for',delay,'seconds and retry'
                time.sleep(delay)
            return None
        else:
            f = ROOT.TFile.Open(fname,'READ')
            return f

    def emptyFile(self,fname):
        """ file must have at least 1D/2D counter histograms """
        ftemp = self.OpenWithRetries(fname)
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
            ftemp = self.OpenWithRetries(fname)
        except:
            "ERROR can't open",fname
            return
        keys=ftemp.GetListOfKeys()
        mykeys=[]
        for k in keys:
            kname=k.GetName()
            # could be have been already found with different cycle number
            if kname in mykeys: continue
            if self.treePatterns:
                accept = False
                for rexpr in self.treePatterns:
                    if re.match(rexpr,kname):
                        accept = True
                        break
                if not accept : continue
            if self.treeExcludePatterns:
                accept = True
                for rexpr in self.treeExcludePatterns:
                    if rexpr in kname:
                        accept = False
                        break
                if not accept : continue
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


    def MergeCounters(self,outFile,InFiles,xsecDB,isSignal):
        """Merge histrograms corresponding to counters"""
        if len(InFiles) == 0: return
        inFiles = copy.deepcopy(InFiles)

        inFileIsDict = False
        if type(inFiles) == dict: inFileIsDict = True

        if inFileIsDict:
            refname = inFiles.keys()[0]
            firstfile = inFiles[refname][0]
        else:
            firstfile = inFiles[0]
            refname = firstfile

        # get the list of histograms from the first file
        scale = 1.
        if not self.config.skipNormWeight: scale = self.normWeight(refname,xsecDB,isSignal)
        ftemp = self.OpenWithRetries(firstfile)
        if not ftemp or ftemp.IsZombie():
            print 'Could not open input file ',firstfile
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
        # also copy the metadata (TObjArrays)
        filesmeta = ftemp.Get('FileList_JobBookeeping_JobBookeeping')
        if filesmeta:
            filsemeta = filesmeta.Clone()
        btag_syst_names = ftemp.Get('btag_weights_names')
        if btag_syst_names :
            btag_syst_names = btag_syst_names.Clone()
        event_syst_names = ftemp.Get('event_weights_names')
        if event_syst_names :
            event_syst_names = event_syst_names.Clone()
        ftemp.Close()
        del ftemp
        if inFileIsDict:
            inFiles[refname] = inFiles[refname][1:]
        else:
            inFiles = inFiles[1:]

        # merge with histograms from other files
        if inFileIsDict:
            for (did, flist ) in inFiles.iteritems():
                scale = 1.
                if not self.config.skipNormWeight: scale = self.normWeight(did,xsecDB,isSignal)
                for fname in flist:
                    ftemp = self.OpenWithRetries(fname)
                    if not ftemp or ftemp.IsZombie():
                        print 'Could not open input file ',fname
                        sys.exit(1)
                        pass
                    for k in hkeys:
                        h = ftemp.Get(k)
                        if not h:
                            print 'Could not extract histogram',k,'in file',fname
                            # sys.exit(1)
                            continue
                        outhists[k].Add(h,scale)
                        pass
                    # append list of files metadata
                    if filesmeta:
                        fm = ftemp.Get('FileList_JobBookeeping_JobBookeeping')
                        if fm:
                            for i in range(fm.GetEntries()):
                                filesmeta.Add(fm.At(i).Clone())
                    ftemp.Close()
                pass
            pass
        else:
            for fname in inFiles:
                ftemp = self.OpenWithRetries(fname)
                if not ftemp or ftemp.IsZombie():
                    print 'Could not open input file ',fname
                    sys.exit(1)
                    pass
                scale = 1.
                if not self.config.skipNormWeight: scale = self.normWeight(fname,xsecDB,isSignal)
                for k in hkeys:
                    h = ftemp.Get(k)
                    if not h:
                        print 'Could not extract histogram',k,'in file',fname
                        # sys.exit(1)
                        # pass
                        continue
                    outhists[k].Add(h,scale)
                    pass
                # append list of files metadata
                if filesmeta:
                    fm = ftemp.Get('FileList_JobBookeeping_JobBookeeping')
                    if fm:
                        for i in range(fm.GetEntries()):
                            filesmeta.Add(fm.At(i).Clone())
                ftemp.Close()
                del ftemp
            pass

        outFile.cd()
        for h in outhists.itervalues():
            if not isSignal:
                # number of events and sum weight no longer have meanings for merged and normalized counters
                h.SetBinContent(1,0.)
                h.SetBinContent(2,0.)
            h.Write()
        if filesmeta: filesmeta.Write('FileList_JobBookeeping_JobBookeeping',ROOT.TObject.kSingleKey)
        if event_syst_names: event_syst_names.Write('event_weights_names',ROOT.TObject.kSingleKey)
        if btag_syst_names: btag_syst_names.Write('btag_weights_names',ROOT.TObject.kSingleKey)

        pass

    def ModifyAndMerge(self,mergename="",outputDir="",dbName="ZeroLeptonRun2/data/MCBackgroundDB.dat",isSignal=False,update=False):
        """Merge all input files and rename trees with modifications in NTVars."""
        if len(self.files) == 0: return
        filesIsDict = False
        if type(self.files) == dict: filesIsDict = True

        # list of all files
        if filesIsDict:
            allfiles = []
            for flist in self.files.itervalues(): allfiles += flist
        else:
            allfiles = self.files

        # list of trees in the input files
        treenames = self.getListOfTrees(allfiles[0])

        print treenames

        # create C++ merger
        xsecDB = None
        try :
            xsecDB = ROOT.SUSY.CrossSectionDB(dbName, False, True)#no path resolver, but use extendedDB
        except :
            print "You are using a SUSYTools version without the extended database set to on by default."
            print "This default was on in AnalysisBase versions <= 2.4.6, and will be fixed for releases >= 2.4.9."
            print "You can switch to these releases to use FilterUpdateMerge.py"
            print "Exiting"
            exit()


        merger = ROOT.FilterUpdateMerge(xsecDB)

        # open output file
        if mergename=="":
            mergename=self.name+self.extraName+".root"
        if update:
            newfile = ROOT.TFile.Open(outputDir+"/"+mergename,"UPDATE")
        else:
            newfile = ROOT.TFile.Open(outputDir+"/"+mergename,"RECREATE")

        # merge counters
        self.MergeCounters(newfile, self.files, xsecDB, isSignal)

        # Build a list of output trees and files that should be merged into it
        # Needed for signal grid as we have one tree name per mc channel
        # while for background MC and signal we merge all files

        # key = output tree name  value = (TTree, input-tree-name, list-of-files)
        outTreeDict = {}
        for inTree in treenames:
            print inTree
            if filesIsDict:
                # we read a dictionary of (did, list-of-files), extract
                # mc channel from dataset name
                for (did,filelist) in self.files.iteritems():
                    ds = None
                    if self.config.doQCD or self.config.doData:
                        ds = 1
                    else:
                        ds = self.extractChannelFromFilename(did)
                        pass

                    if not ds:
                        print 'Could not identify channel number for ',did
                        sys.exit(1)
                        pass

                    if not self.config.doData and not isSignal and xsecDB.rawxsect(ds) < 0.:
                        print 'No cross-section for sample',ds,'skip file',did
                        continue

                    newtname=self.getNewTreeName(inTree, ds)
                    if outTreeDict.has_key(newtname):
                        for f in filelist: outTreeDict[newtname][2].append(f)
                    else:
                        outTree = ROOT.TTree(newtname,'0-lepton small ntuple')
                        outTreeDict[newtname] = (outTree, inTree, copy.deepcopy(filelist))

                    pass
            else:
                # we read files from a list, extract channel number from
                # filename
                for filename in self.files:
                    ds = None
                    if self.config.doQCD or self.config.doData:
                        ds = 1
                    else:
                        ds = self.extractChannelFromFilename(filename)
                        pass

                    if not ds:
                        print 'Could not identify channel number for ',filename
                        sys.exit(1)
                        pass

                    if not self.config.doData and not isSignal and xsecDB.rawxsect(ds) < 0.:
                        print 'No cross-section for sample',ds,'skip file',filename
                        continue

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

            doExtraVars = not self.config.skipExtraVars
            doRJigsawVars = self.config.mergeRJigsawVars
            doCRWTVars = not self.config.skipCRWTVars
            doCRZVars = not self.config.skipCRZVars
            doCRYVars = not self.config.skipCRYVars
            # check that the first file in the list has the requested block
            for name in filelist:
                testf = self.OpenWithRetries(name)
                if not testf or testf.IsZombie(): continue
                testt = testf.Get(inTreeName)
                if not testt: continue
                if doExtraVars and not testt.GetBranch('NTExtraVars'):
                    doExtraVars = False
                if doRJigsawVars and not testt.GetBranch('NTRJigsawVars'):
                    doRJigsawVars = False
                if doCRWTVars and not testt.GetBranch('NTCRWTVars'):
                    doCRWTVars = False
                if doCRZVars and not testt.GetBranch('NTCRZVars'):
                    doCRZVars = False
                if doCRYVars and not testt.GetBranch('NTCRYVars'):
                    doCRYVars = False
                break
            merger.process(outTree, inTreeName, inList, isSignal, not self.config.skipNormWeight, self.config.filter, doExtraVars, doRJigsawVars, doCRWTVars, doCRZVars, doCRYVars)
            newfile.cd()
            outTree.Write("",ROOT.TObject.kOverwrite)

        newfile.Close()

        return

if __name__ == '__main__':
    config = parseCmdLine()
    #config = parseCmdLine(sys.argv[1:])

    # Load RootCore libs
    ROOT.gROOT.SetBatch()
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

        if config.allBackgroundsNT or config.TopNT: AllSamples.append(Top)
        if config.allBackgroundsNT or config.WNT: AllSamples.append(WMassiveCB)
        if config.allBackgroundsNT or config.ZNT: AllSamples.append(ZMassiveCB)
        if config.allBackgroundsNT or config.DBNT: AllSamples.append(DB)
        if config.allBackgroundsNT or config.GAMNT: AllSamples.append(GAMMAMassiveCB)
        if config.allBackgroundsNT or config.QCDNT: AllSamples.append(QCD)
        if config.allBackgroundsNT or config.SplitQCDNT:
            for id in lQCDMC:
                AllSamples.append(Sample(name='QCD'+str(id),myarg=[id],inputfile=config.inputfile_mc,config=config,treename='QCD'))

        # define alternative samples
        WMadgraph = Sample('WMadgraphPythia8',lWjetsMadgraph,config.inputfile_mc,config,treename="WMadgraphPythia8")
        ZMadgraph = Sample('ZMadgraphPythia8',lZjetsMadgraph,config.inputfile_mc,config,treename="ZMadgraphPythia8")
        TopPowhegP8 = Sample('TopPowhegPythia8',lttbarPowhegP8[:]+lsinglet+lttbarX,config.inputfile_mc,config,treename="TopPowhegPythia8")
        TopPowhegHerwig = Sample('TopPowhegHerwig',lttbarPowhegHerwing[:]+lsinglet+lttbarX,config.inputfile_mc,config,treename="TopPowhegHerwig")
        TopMCatNLO = Sample('TopMCatNLO',lttbarMCatNLO[:]+lsinglet+lttbarX,config.inputfile_mc,config,treename="TopMCatNLO")
        if config.allBackgroundsNT or config.TopPowP8NT: AllSamples.append(TopPowhegP8)
        if config.allBackgroundsNT or config.TopPowHNT: AllSamples.append(TopPowhegHerwig)
        if config.allBackgroundsNT or config.TopMCNLONT: AllSamples.append(TopMCatNLO)
        if config.allBackgroundsNT or config.WMadP8NT: AllSamples.append(WMadgraph)
        if config.allBackgroundsNT or config.ZMadP8NT: AllSamples.append(ZMadgraph)

        print AllSamples
    #if config.doZLO==True:
    #    ZLO = Sample('ZLO',lZjetsLO,config.inputfile_mc,config,treename="Z")
    #    AllSamples.append(ZLO)

    if config.doData==True:
        # Run 1 streams
        DataJetTauEtmiss = Sample('Data',['JetTauEtmiss'],config.inputfile_data,config,Extraname="JetTauEtmiss")
        DataMuon = Sample('Data',['Muon'],config.inputfile_data,config,Extraname="Muon")
        DataEgamma = Sample('Data',['Egamma'],config.inputfile_data,config,Extraname="Egamma")
        # Run 2 stream
        DataMain = Sample('Data',['Main'],config.inputfile_data,config,Extraname="Main")
        DataDebug = Sample('Data',['debugrec_hlt'],config.inputfile_data,config,Extraname="Debug")
        AllSamples += [DataJetTauEtmiss,DataMuon,DataEgamma,DataMain,DataDebug]

    if config.doQCD==True:
        QCDdd1 = Sample('QCDdd',['JetTauEtmiss'],config.inputfile_qcd,config)
        AllSamples += [QCDdd1]

    # do merging
    for s in AllSamples:
        s.ModifyAndMerge(dbName=config.dbName, outputDir = config.outputDir)



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
            signal.ModifyAndMerge(mergename=grid+".root",outputDir = config.outputDir,dbName=config.dbName,isSignal=True,update=False)
