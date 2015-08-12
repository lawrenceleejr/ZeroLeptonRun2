#!/usr/bin/env python
"""
############################################################
#
# Create an extended version of susy_crosssections.txt
# including information on MC statistics and sum of weights
#
#FIXME: Should define a way to use files on the grid storage
#       for which the file path does not allow to infer the
#       dataset ID
############################################################
"""

import ROOT
import sys, os, commands, re, shutil, subprocess, copy
from math import *
ROOT.gROOT.SetBatch(True)

#####################################################


#####################################################
#
#####################################################

def getListFiles(indir):
    """Get list of folders and files for an input directory."""
    if not os.path.isdir(indir): return
    fout = open('files.list','w')
    for subdir in os.listdir(indir):
        if os.path.isdir(indir+'/'+subdir):
            for fname in os.listdir(indir+'/'+subdir):
                if os.path.isfile(indir+'/'+subdir+'/'+fname) and '.root' in fname:
                    fout.write(indir+'/'+subdir+'/'+fname+'\n')
        elif os.path.isfile(indir+'/'+subdir) and '.root' in subdir:
            fout.write(indir+'/'+subdir+'\n')
    fout.close()
    return

def getMCSampleList(strlds):
    """Get list of specific dataset numbers to check via --lds field."""
    lds = []
    if strlds == "": return lds
    for strds in strlds.split(','):
        if hasattr(MCSampleList,strds):
            print "getMCSampleList:",strds,'found in MCSampleList'
            lds += copy.deepcopy(getattr(MCSampleList,strds))
        elif strds.isdigit():
            lds.append(int(strds))
    print "getMCSampleList: list of dataset numbers found",lds
    return lds

def updateMap(processMap,xsecDB,key,nb_of_events,weight):
    if key in processMap.keys():
        processMap[key][0]+=nb_of_events
        processMap[key][1]+=weight
    else:
        xsec = -1.
        if isinstance(key,tuple):
            xsec = xsecDB.xsectTimesEff(key[0],key[1])
        else:
            xsec = xsecDB.xsectTimesEff(key)

        if xsec < 0:
            print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
            print 'No cross-section found for sample ID '+str(key),"!!!!!!!!!!!!!!!"
            print key,xsec
            print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
            xsec=0
        # get number of events that the job saw
        processMap[key]=[nb_of_events,weight,xsec]
        pass
    pass
    return processMap


def processFile(config,fname,channel,xsecDB,myMap):
    if fname.startswith('/eos'): fname = 'root://eosatlas/'+fname
    tfile=ROOT.TFile.Open(fname)
    hist=tfile.Get('Counter_JobBookeeping_JobBookeeping')
    if hist!=None and hist.GetNbinsX()>=2:
        # if not channel in hsum:
        #     hsum[channel] = copy.deepcopy(ROOT.TH1D(hist))
        # else:
        #     hsum[channel].Add(hist)
        nb_of_events=0
        weight=0
        if config.isSignal==False:
            nb_of_events=hist.GetBinContent(1)
            weight=hist.GetBinContent(2)
            updateMap(myMap,xsecDB,channel,nb_of_events,weight)
        else:
            """
            gets complicated with derivations: from JobBookeeping we
            can extract the sum of weights for all events before skimming
            but we don't know for which hardproc. All we can do is calculate
            the average efficiency of the derivation over hardprocs and
            rescale for that. Would be OK for simplified models but can
            be really wrong if many very different processes are present.
            """
            hist2=tfile.Get('Counter_for_ZeroLeptonCounterSRAll')
            sumwini = hist.GetBinContent(2)
            sumwkept = 0.
            for ybin in range(1,hist2.GetNbinsY()+1):
                sumwkept+=hist2.GetBinContent(2,ybin)
                pass
            scale = sumwini/sumwkept
            for ybin in range(1,hist2.GetNbinsY()+1):
                nb_of_events=hist2.GetBinContent(1,ybin)
                weight=hist2.GetBinContent(2,ybin)
                if nb_of_events>0:
                    updateMap(myMap,xsecDB,(channel,ybin-1),nb_of_events*scale,weight*scale)
                    pass
                pass
            pass
        pass
    else:
        print "WARNING: file ignored ", fname 
        pass
    tfile.Close()

def extractChannel(name):
    for prefix in ('mc11_7TeV.', 'mc12_8TeV.', 'mc14_8TeV.','mc14_13TeV.', 'mc15_13TeV.'):
        if prefix in name:
            return int(name.split(prefix)[1].split('.')[0])
    return None


def getMCInfo(config):
    myMap={}
    xsecDB = ROOT.SUSY.CrossSectionDB(config.xsecfile)

    lds = getMCSampleList(config.lds)
    hsum = {}

    if config.inputfilename.endswith('.pkl'):
        import pickle
        if not os.path.isfile(config.inputfilename):
            print "can't read file",config.inputfilename
            sys.exit()
            pass
        picklefile = open(config.inputfilename,'rb')
        myfiles = pickle.load(picklefile)
        picklefile.close()

        for (did,flist) in myfiles.iteritems():
            channel = extractChannel(did)
            for fname in flist:
                processFile(config,fname,channel,xsecDB,myMap)
    else:
        for fname in open(config.inputfilename,'read'):    
            fname.strip()
            if len(fname) == 0: continue
            if fname.startswith('#'): continue
            fname = fname.strip()

            # get dataset number from input file
            channel = int(fname.split(".")[int(config.NB)])
            if len(lds) and not channel in lds: continue

            if ( not config.isSignal and xsecDB.rawxsect(channel) < 0. ): 
                print 'No cross section for sample',channel,'skip file',fname
                continue

            # get weight from histogram
            processFile(config,fname,channel,xsecDB,myMap)


    f=open(config.outputfilename,'w')
    ftex=open(config.outputfilename.split('.')[0]+'.tex','w')
    if config.isSignal==False:
        f.write("id/I:name/C:xsec/F:kfac/F:eff/F:relunc/F:sumw/F:stat/F\n")
        ftex.write(r"""\begin{table}
\scriptsize
\begin{center}
\begin{tabular}{|l|l|r|r|r|r|r|}
\hline
Dataset ID & Dataset name & $\sigma \times \epsilon$ [pb] & k-factor & $N_{gen}$ & $\Sigma w$ & $\mathcal{L}_{int}\ [\mathrm{fb}^{-1}]$ \\
\hline
""")
        for channel,info in sorted(myMap.items()):
            line=str(channel)+" "+str(xsecDB.name(channel))+" "+str(xsecDB.rawxsect(channel))+" "+str(xsecDB.kfactor(channel))+" "+str(xsecDB.efficiency(channel))+" "+str(xsecDB.rel_uncertainty(channel))+" "+str(info[1])+" "+str(info[0])#not very optimal
            #print line
            f.write(line+"\n")
            if  channel == 110070 or channel == 110071 :
                # Powheg weighted events, sum(weight) = xsec independent of
                # number of generated events so use number of generated events
                # as an approximation of the statistics
                lint = info[0] / (xsecDB.rawxsect(channel) * xsecDB.efficiency(channel) * xsecDB.kfactor(channel) * 1000.)
            else:
                lint = info[1] / (xsecDB.rawxsect(channel) * xsecDB.efficiency(channel) * xsecDB.kfactor(channel) * 1000.)
            ftex.write("%d & %s & %.3g & %3.2f & %d & %.4g & %.3g \\\\\n" % (channel,xsecDB.name(channel).replace('_','\_'),xsecDB.rawxsect(channel) * xsecDB.efficiency(channel),xsecDB.kfactor(channel),info[0],info[1],lint))
            #print channel," 0 ",info[1]
            pass
        ftex.write(r"""\hline
\end{tabular}
\end{center}
\caption{Summary of samples used in the analysis. $\sigma \times \epsilon$ corresponds to the cross-section provided by the Monte Carlo generator times the possible truth level filter efficiency. k-factor is used to normalize the generator cross-section to the best known cross-section for a process. $\Sigma w$ includes the generator level weights and the pileup reweighting effect. $\mathcal{L}_{int}$ is the generated equivalent integrated luminosity.\label{tab:crosssection_%s}}
\end{table}
""" % config.outputfilename.split('.')[0])
    else:
        f.write("id/I:hardproc/I:xsec/F:kfac/F:eff/F:relunc/F:sumw/F:stat/F\n")
        ftex.write(r"""\begin{table}
\scriptsize
\begin{center}
\begin{tabular}{|l|l|l|r|r|r|r|}
\hline
Dataset ID & Dataset name & $\sigma \times \epsilon$ [pb] & $N_{gen}$ & $\mathcal{L}_{int}\ [\mathrm{fb}^{-1}]$ \\
\hline
""")
        for key,info in sorted(myMap.items()):
            channel = key[0]
            hardproc = key[1]
            line=str(channel)+" "+str(hardproc)+" "+str(xsecDB.rawxsect(channel,hardproc))+" "+str(xsecDB.kfactor(channel,hardproc))+" "+str(xsecDB.efficiency(channel,hardproc))+" "+str(xsecDB.rel_uncertainty(channel,hardproc))+" "+str(info[1])+" "+str(int(round(info[0])))#not very optimal
            #print line
            f.write(line+"\n")

            lint = info[1] / (xsecDB.rawxsect(channel,hardproc) * xsecDB.efficiency(channel,hardproc) * xsecDB.kfactor(channel,hardproc) * 1000.)
            ftex.write("%d & %s & %.3g & %d & %.5g \\\\\n" % (channel,xsecDB.name(channel).replace('_','\_'),xsecDB.rawxsect(channel,hardproc) * xsecDB.efficiency(channel,hardproc),info[0],lint))
        ftex.write(r"""\hline
\end{tabular}
\end{center}
\caption{Summary of signal samples used in the analysis. $\sigma \times \epsilon$ corresponds to the NLO cross-section times the possible truth level filter efficiency. $\mathcal{L}_{int}$ is the generated equivalent integrated luminosity.\label{tab:crosssection_%s}}
\end{table}
""" % config.outputfilename.split('.')[0])

    f.close()
    ftex.close()

    return myMap

#------------------------------------------------------------------------
def parseCmdLine(args):
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--xsecfile", dest="xsecfile", help="cross section file", default="SUSYTools/data/susy_crosssections_13TeV.txt")
    parser.add_option("--input", dest="inputfilename", help="List of input filenames", default="")
    parser.add_option("--output", dest="outputfilename", help="output filename", default="MCBackgroundDB.dat")
    parser.add_option("--no", dest="NB", help="position of mc_channel_id in the string (splitted by the character '.')", default=3)
    parser.add_option("--lds", dest="lds", help="List of dataset numbers like 105200, lW, lTop", default="")
    parser.add_option("--indir", dest="indir", help="Input directory (can be used instead of --input)", default="")
    parser.add_option("--isSignal", dest="isSignal", help="isSignal",action='store_true', default=False)
    parser.add_option("--prefix", dest="prefix", help="Prefix to identify mc production",default="mc15_13TeV") 
    (config, args) = parser.parse_args(args)
    return config


if __name__ == '__main__':
    config = parseCmdLine(sys.argv[1:])

    ## Load RootCore libs
    ROOT.gROOT.SetBatch()
    ROOT.gROOT.ProcessLine(".x $ROOTCOREDIR/scripts/load_packages.C")
    print 'Packages loaded'

    if config.inputfilename == "" and config.indir != "":
        getListFiles(config.indir)
        config.inputfilename = 'files.list'

    # List of dataset used in analysis
    if config.prefix == 'mc12_8TeV':
        import mc12_8TeV_MCSampleList as MCSampleList
    elif config.prefix == 'mc14_13TeV':
        import mc14_13TeV_MCSampleList as MCSampleList
    elif config.prefix == 'mc15_13TeV':
        import mc15_13TeV_MCSampleList as MCSampleList
    elif config.prefix == 'mc15_week1':
        import mc15_13TeV_week1_MCSampleList as MCSampleList
    else:
        print 'Unsupported mc production type',config.prefix
        sys.exit(1)


    myMap=getMCInfo(config)
