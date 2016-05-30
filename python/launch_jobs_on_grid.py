#!/usr/bin/env python

__doc__ = """
This script can be used to launch jobs on the grid
Usage:
setup RootCore and pathena 
chmod u+x ZeroLeptonRun2/python/launch_jobs_on_grid.py;
ZeroLeptonRun2/python/launch_jobs_on_grid.py --help
ZeroLeptonRun2/python/launch_jobs_on_grid.py
"""

import os,sys,subprocess,datetime,re

def parseCmdLine(args):
    """ Parse input command line to optdict.
    To get the whole list of options type : launch_jobs_on_grid.py -h"""
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--inDS", dest="inDS", help="Name of single input dataset",default="") 
    parser.add_option("--inDSfile", dest="inDSfile", help="Names of input datasets stored in a file",default="") 
    parser.add_option("--prefix", dest="prefix", help="Prefix appended to the output DS name, should be something like user.yourname",default="user.%s." % os.environ['USER']) 
    parser.add_option("--suffix", dest="suffix", help="Suffix appended to the output DS name in case you need to add a version for example",default="") 
    parser.add_option("--debug", dest="debug", help="Debug mode", action='store_true', default=False)
    parser.add_option("--config", dest="configfile", help="cafe configuration file", default="")
    parser.add_option("--runopts", dest="runopts", help="Command line arguments to cafe executable", default="")
    parser.add_option("--container", dest="container", help="If the input datasets are in a container with multiple run numbers, configuure prun to keep them distinct", action='store_true', default=False)
    parser.add_option("--prunopts", dest="prunopts", help="Command line arguments to prun", default="")
    parser.add_option("--signal", dest="signal", help="Input is signal MC", action='store_true', default=False)
    parser.add_option("--noGRL", dest="noGRL", help="If true no GRL is applied", action='store_true', default=False)

    parser.add_option("--tmpDir", dest="tmpDir", help="Tmp dir", default="")
    (config, args) = parser.parse_args(args)
    return config


def findFiles(topdir,extensions):
    matches = []
    for root, dirnames, filenames in os.walk(topdir):
        if '.svn' in root: continue
        for filename in filenames:
            for ext in extensions:
                if filename.endswith(ext):
                    matches.append(os.path.join(root, filename))
    return matches

def main():
    config = parseCmdLine(sys.argv[1:])

    ## Get list of input datasets    
    linDS = []
    if config.inDS != "": 
        linDS.append(config.inDS)
    elif config.inDSfile != "":
        finDS = open(config.inDSfile,'r')
        for inDS in finDS:
            inDS = inDS.strip()
            if len(inDS)>0 :
                inDS = inDS.split()[0]
                if not inDS.startswith('#'): linDS.append(inDS)
        finDS.close()
    if config.debug: print 'list of datasets',linDS
    if len(linDS) == 0:
        print "You did not specify any input dataset or file containing a list of input datasets.\nPlease type launch_jobs_on_grid.py -h"
        sys.exit(0)

    # find local file to explicitely include them on the prun command line
    rootCoreDir = os.environ['ROOTCOREBIN']
    localfiles = []
    localfiles += findFiles('CDIFiles',('.root'))
    flpkgs = open(rootCoreDir+'/packages','r')
    for pkgname in flpkgs:
        pkgname = pkgname.split(':')[0]
        if not os.path.isdir(pkgname): continue
        localfiles += findFiles(pkgname,('.root','.txt','.xml','.env','tarball'))

    ## Loop over input datasets
    first = True
    theproc = ''
    lastproc = ''
    lastpdf = ''
    systlist = []
    for inDS in linDS:
        if inDS.startswith("#"): 
            continue
        
        print 'launching grid job for',inDS

        isSignal = config.signal
        # check for known signal MC patterns
        if not isSignal:
            for pattern in ['SM_SS', 'SS_direct', 'SM_GG', 'GG_direct', 'GG_onestep', 'SM_SG', 'Gtt', 'Gluino_Stop_charm', 'Tt', 'compsusy', 'msugra', 'MSUGRA', '_0_10_P', 'NUHMG', 'NUHM2', 'bRPV']:
                if pattern in inDS:
                    print 'override option to signal = True'
                    isSignal = True
                    break

        outDS = re.sub('^(user\.[a-zA-Z]*\.)?',config.prefix,inDS)
        if ',' in outDS: outDS = outDS.split(',')[0]
        if config.suffix != "":
            suffix = config.suffix
            if config.debug: suffix += '.debug'
            if outDS.endswith('/'):
                outDS = outDS.split('/')[0]+suffix+'/'
            else:
                outDS += suffix
        if len(outDS) > 100: # panda will complain dataset name is too long
            outDS = outDS.replace('.merge.AOD','')
            outDS = outDS.replace('.merge.DAOD_SUSY1','.SUSY1')
            outDS = outDS.replace('.merge','')
            outDS = outDS.replace('.NTUP_SUSY','')
            outDS = outDS.replace('Sherpa_CT10_','Sherpa_')
            outDS = outDS.replace('MassiveCB','MCB')
            outDS = outDS.replace('CJetFilter','CF')
            outDS = outDS.replace('BFilter','BF')
            outDS = outDS.replace('CVeto','CV')
            outDS = outDS.replace('BVeto','BV')
            outDS = outDS.replace('PowhegPythia8EvtGen_','PowHP8EvG_')
            outDS = outDS.replace('PowhegPythiaEvtGen_','PowHPEvG_')
            outDS = outDS.replace('PowhegHerwigppEvtGen_','PowHHppEvG_')
            outDS = outDS.replace('MadGraphPythiaEvtGen_', 'MG_PEvG_')
            outDS = outDS.replace('MadGraphPythia8EvtGen_','MG_PEvG_')
            outDS = outDS.replace('SinglePhoton','1Gam_')
            if len(outDS) > 100: outDS = outDS.replace('_AUET2CTEQ6L1MPI','')
            if len(outDS) > 100: outDS = outDS.replace('_AUET2BCTEQ6L1','')
            if len(outDS) > 100: outDS = outDS.replace('_AUET2_CTEQ6L1','')
            if len(outDS) > 100: outDS = outDS.replace('_AUET2','')
            if len(outDS) > 100: outDS = outDS.replace('_CTEQ6L1','')
            if len(outDS) > 100: outDS = outDS.replace('_MPI','')
            if len(outDS) > 100: outDS = outDS.replace('_AU2CT10','')
            if len(outDS) > 100: outDS = outDS.replace('_CT10','')
            if len(outDS) > 100: outDS = outDS.replace('_UEEE3_CTEQ6L','')
            if len(outDS) > 100: outDS = outDS.replace('_AZNLOCTEQ6L1','')
            if len(outDS) > 100: outDS = outDS.replace('_A14NNPDF23LO','')
            if len(outDS) > 100: outDS = outDS.replace('_P2012','')
            if len(outDS) > 100: outDS = outDS.replace('_hdamp172p5','')
            if len(outDS) > 100: outDS = outDS.replace('_NNPDF30NNLO','')
            if '_tid' in outDS and len(outDS) > 100:
                rstr = outDS.split(config.suffix)[0].split('_tid')[1]
                outDS = outDS.replace('_tid%s' % rstr,'')

        print 'output DS: ' + outDS, len(outDS+"_o.root/") 

        # which cafe config file should be used ?
        if config.configfile == "":
            if ("mc15" in inDS or "data15" in inDS) and not("TRUTH1" in inDS):
                cafeconfig = "ZeroLeptonRun2/config/zerolepton.config"
            elif "mc14" in inDS or "data12" in inDS:
                cafeconfig = "ZeroLeptonRun2/config/zerolepton_DC14.config"
            elif "TRUTH1" in inDS:
                cafeconfig = "ZeroLeptonRun2/config/zeroleptontruth.config"
            else:
                print "Unexpected dataset name, could not figure out which cafe config file to use"
                sys.exit(1)
        else:
            cafeconfig = config.configfile

        outputs = "o.root" 
        scriptcmd = r"""cp ../PoolFileCatalog.xml . ; ZeroLeptonRun2/python/pfc2txt.py; cat pfc.txt; echo %IN| sed 's/\,/\n/g'>inputfiles; unset ROOT_TTREECACHE_SIZE; ln -sf \${TestArea}/CDIFiles/13TeV . ; cafe """ + cafeconfig + """ Events: -1  Input: filelist:inputfiles  Output: o.root """ 

        # Real data ?
        if "data11" in inDS or "data12" in inDS or "data15" in inDS: 
            scriptcmd += " Global.IsData: TRUE Global.IsSignal: FALSE "
            if 'physics_Egamma' in inDS:
                 scriptcmd += "crwt.IsElectronChannel: TRUE vrwt.IsElectronChannel: TRUE "
            elif 'physics_Muons' in inDS:
                 scriptcmd += "crwt.IsMuonChannel: TRUE vrwt.IsMuonChannel: TRUE "

        # Fast or full simulation ?
        if "mc11" in inDS or "mc12" in inDS or "mc14" in inDS  or "mc15" in inDS :
            tag = inDS.split(".")[-1]
            tags = tag.split("_")
            if len(tags)>2 and tags[2].startswith("a"):
                scriptcmd += " Global.IsAtlfast: TRUE"
                scriptcmd += " Global.doFake: FALSE"


        # Run period
        if "_7TeV." in inDS:
            print "7TeV data/MC not supported in ZeroLeptonRun2"
            sys.exit(1)
        elif "_8TeV." in inDS:
            scriptcmd += " Global.Period: p8tev"
        elif "_13TeV." in inDS or "_14TeV." in inDS :
            scriptcmd += " Global.Period: p13tev"
        else:
            print "Could not identify the run period (7/8/13/14 TeV) from the input dataset ",inDS
            sys.exit(1)

        # test special derivation tags
        if not 'Global.DerivationTag' in config.runopts:
            tag = inDS.split(".")[-1]
            knowntags = ['p2353', 'p2363', 'p2372', 'p2375', 'p2377', 'p2384', 'p2419', 'p2425' , 'p2436', 'p2452', 'p2470']
            found = False
            if not('TRUTH1' in inDS):
                for t in knowntags:
                    if t in tag:
                        scriptcmd += " Global.DerivationTag: "+t
                        found = True
                        break
                    pass
                if not found:
                    scriptcmd += " Global.DerivationTag: NA "

        # signal events
        if isSignal:
            scriptcmd += " Global.IsSignal: TRUE "

        # GRL
        if config.noGRL:
            scriptcmd += " grl.passAll: TRUE "

        extfiles=' --extFile '
        for f in localfiles:
            extfiles += f+','
        extfiles = extfiles[:-1]+' '

        scriptcmd += " " + config.runopts
        cmd = r"""prun --excludeFile=\*/.svn/\*,\*oxbridgekinetics-0.6/\*,\*OxbridgeKinetics/lib\* --exec "%(scriptcmd)s" --useRootCore --inDS %(inDS)s --outputs %(outputs)s  --outDS %(outDS)s  %(prunopts)s""" % {'inDS':inDS,'outDS':outDS,'scriptcmd':scriptcmd,'prunopts':config.prunopts,'outputs':outputs}
        cmd += extfiles
        if config.container and not 'data11' in inDS:
            cmd += ' --useContElementBoundary'
        if not ( "--nGBPerJob" in config.prunopts or "--nFilesPerJob"  in config.prunopts):
            cmd += ' --nGBPerJob=MAX '

        if len(linDS) > 1:
            if config.tmpDir != "":
                tarfile =  config.tmpDir+'/TarBall_prun.tar.gz'
            else:
                tarfile = '/tmp/'+os.environ['USER']+'/TarBall_prun.tar.gz'
                pass
            # remove tarfile the first time around, unless it might be the one
            # specified in the user command line
            if not 'inTarBall' in cmd:
                if first and os.path.exists(tarfile):
                    ret = subprocess.call('rm -f %s' % tarfile, shell=True)
                if first:
                    cmd += ' --outTarBall=%s' % tarfile
                else:
                    cmd += ' --inTarBall=%s' % tarfile
        if not "--excludedSite" in config.prunopts:
            cmd += ' --excludedSite=TEST'
        if config.tmpDir != "":
            cmd += ' --tmpDir %s' % config.tmpDir
        if config.debug:
            cmd += ' --nFiles 1'
            print cmd
        first = False
        ret = 0
        #print cmd
        #sys.exit(0)
        ret = subprocess.call(cmd, shell=True)
        if ret != 0:
            print 'command failed => stop script'
            #sys.exit(1)
    return

if __name__ == "__main__":
    main()
