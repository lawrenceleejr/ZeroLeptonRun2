#!/bin/env python

#
# localSetupPyAMI
#

import os,sys,subprocess,datetime,re,string
try:
    import pyAMI.client
    import pyAMI.atlas.api
except:
    print 'Could not import AMI modules, did you "localSetupPyAMI" ?'


# List of susy signals
lsignals = ['Herwigpp_UEEE3_CTEQ6L1_MSUGRA','Herwigpp_UEEE4_CTEQ6L1_Gtt','SM_SS','SM_GG','SM_SG','pMSSM_qL_to_h']

# List of keywords used to remove non-used datasets
filters = ['Pythia_GGM','PG11','DGemt','DGnoL','_Wprime','_Zprime','PythiaB','_3Leptons','LeptonPhotonFilter','ADDGraviton','Sherpa_CT10_SingleTop','BlackMaxPythia8','Charybdis2Pythia8','QBHPythia8','_WimpPair_','_ZH125_','_ggH125_','_VBFH125_','_WH125_','ParticleGenerator_']

# bad tags
# https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/DC14DataMCSampleInfo#13_TeV_Monte_Carlo
# https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/AtlasProductionGroupMC14a#DC14_Datasets
# what about r5869, r5899, r5894, r5910, r5898, r5812
# preproduction samples +  no pileup 
badtags = ['_r5720_','_r5721_'] + ['_r5803_','_r5804_','_r5862_']

# mc15 reco tag lists
mc15_rtags = {}
mc15_rtags['week1']  = ['r6633', 'r6853' ]
mc15_rtags['50ns']  = [ 'r6630', 'r6655', 'r6647', 'r6793', 'r6828', 'r6767', 'r6802' ]
mc15_rtags['25ns']  = [ 'r6725', 'r6765' ]


def genParamsFromParents(client,datasetName,datasetNumber):
    from pyAMI.atlas.api import get_dataset_info
    approx_GenFiltEff = None
    xsec = None
    prov =  pyAMI.atlas.api.get_dataset_prov(client, datasetName)
    for parent in prov['node']:
        # minbias overlays are also parents so need to 
        # check the channel number
        if int( parent[u'logicalDatasetName'].split(".")[1]) == datasetNumber:
            parentinfo = get_dataset_info(client,parent[u'logicalDatasetName'])[0]
            if parentinfo.has_key(u'approx_GenFiltEff'):
                approx_GenFiltEff = parentinfo[u'approx_GenFiltEff']
                pass
            if parentinfo.has_key(u'crossSection') and parentinfo[u'crossSection'] != u'NULL':
                xsec = float(parentinfo[u'crossSection'])
                pass
            if approx_GenFiltEff and xsec:
                break
        pass
    return (xsec,approx_GenFiltEff)

    

def badDataset(datasetName,generatorString,version):
    dsname = datasetName
    if '/' in dsname: dsname = dsname[:-1]
    # any of the veto'ed tags
    for t in badtags:
        if t in version: return True
    # for now veto AtlFast II
    tags = string.split(version,'_')
    for t in tags:
        if t.startswith('a'): return True
    return False

def parseCmdLine(args):
    """ Parse input command line to optdict.
    To get the whole list of options type : get_list_dataset.py -h"""
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--prefix", dest="prefix", help="Prefix to search datasets in AMI",default="mc15_13TeV") 
    parser.add_option("--whichMC15", dest="whichMC15", help="week1, 50ns or 25ns",default="") 
    parser.add_option("--datatype", dest="datatype", help="datatype",default="%.merge.DAOD_SUSY1%") 
    parser.add_option("--suffix", dest="suffix", help="Suffix appended to the output file name",default="") 
    parser.add_option("--tag", dest="tag", help="Production tag",default="_r6630_r6264_p2353") 
    parser.add_option("--server", dest="server", help="AMI server (main or replica)",default="main") 
    parser.add_option("--onlyComplete", dest="onlyComplete", help="skip datasets not marked with all events available", action='store_true', default=False)
    parser.add_option("--official", dest="official", help="Select only baseline or alternative samples from MCSampleList", action='store_true', default=False)
    parser.add_option("--baseline", dest="baseline", help="Select only baseline samples from MCSampleList", action='store_true', default=False)
    parser.add_option("--sample", dest="sample", help="Selected any sample list defined in MCSampleList, e.g. lTop", default=None)
    parser.add_option("--signal", dest="signal", help="For MC select only signal samples", action='store_true', default=False)
    (config, args) = parser.parse_args(args)
    return config

def main():
    # configurable options
    config = parseCmdLine(sys.argv[1:])

    if (config.baseline or config.official ) and config.sample:
        print "--baseline, --official and --sample are mutually exclusive"
        sys.exit(1)

    # AMI client connection
    client = pyAMI.client.Client('atlas')
    pyAMI.client.endpoint = config.server
    pyAMI.atlas.api.init()

    # consistency checks
    if config.whichMC15 != '':
        if config.whichMC15 == 'week1' and config.prefix != 'mc15_week1':
            print 'prefix changed to mc15_week1 in agrement with whichMC15'
            config.prefix = 'mc15_week1'
        elif config.whichMC15 == '50ns'  and config.prefix != 'mc15_13TeV':
            print 'prefix changed to mc15_13TeV in agrement with whichMC15'
            config.prefix = 'mc15_13TeV'
        elif config.whichMC15 == '25ns'  and config.prefix != 'mc15_13TeV':
            print 'prefix changed to mc15_13TeV in agrement with whichMC15'
            config.prefix = 'mc15_13TeV'

    # data type is NTUP_SUSY for 2011/2012 and AOD for 2014 on
    datatype = config.datatype
    if 'mc11_' in config.prefix or 'mc12_' in config.prefix or 'data11_' in config.prefix or 'data12_' in config.prefix : datatype = '%.merge.NTUP_SUSY%'

    # make list of official datasets (baseline+alt)
    officialids = []
    if config.official or config.baseline or config.sample:
        if 'mc12_8TeV' in config.prefix or 'mc14_8TeV' in config.prefix:
            import mc12_8TeV_MCSampleList as mcsl
        elif 'mc14_13TeV' in config.prefix:
            import mc14_13TeV_MCSampleList as mcsl
        elif 'mc15_13TeV' in config.prefix:
            import mc15_13TeV_MCSampleList as mcsl
        elif 'mc15_week1' in config.prefix:
            import mc15_13TeV_week1_MCSampleList as mcsl
        else:
            print '--official is only supported for mc12_8TeV, mc14_8TeV, mc14_13TeV, mc15_13TeV and mc15_week1'
            sys.exit(1)
        if config.sample:
            officialids = mcsl.__dict__[str(config.sample)]
        else:
            officialids = mcsl.__dict__["lbaseline"]
            if config.official:
                officialids += mcsl.__dict__["lalt"]

    # get all datasets matching prefix & tag and then filter them
    from pyAMI.atlas.api import get_dataset_info, list_datasets

    alldatasets = []
    if config.whichMC15 != '':
        prefix = config.prefix
        if prefix == 'mc15_week1': prefix = 'mc15_13TeV'
        for tag in mc15_rtags[config.whichMC15]:
            dskey = prefix+datatype+tag+config.tag
            print 'Querying AMI for datasets matching pattern',dskey
            alldatasets += list_datasets(client,dskey)
    else:
        prefix = config.prefix
        if prefix == 'mc15_week1': prefix = 'mc15_13TeV'
        dskey = config.prefix+datatype+config.tag
        print 'Querying AMI for datasets matching pattern',dskey
        alldatasets = list_datasets(client,dskey)

    acceptedDS = []
    for DSlist in alldatasets:
        dsname = DSlist['ldn']
        cut = False
        for filter in filters:
            if filter in dsname.split('.')[2]: cut = True
        if (config.official or config.baseline or config.sample) and not int(dsname.split('.')[1]) in officialids: cut = True
        if config.signal :
            cut = True
            for pattern in lsignals:
                if pattern in dsname: cut = False
        if cut: continue
        acceptedDS.append(dsname)
        pass
    acceptedDS.sort()

    # get informations for all accepted datasets
    dsinfos = []
    for dsname in acceptedDS:
        dsinfos.append(get_dataset_info(client,dsname)[0])
        pass

    # write file
    coveredids = set()
    fout = open('datasets.txt','w')
    for info in dsinfos:
        try:
            dsname = info['logicalDatasetName']
            generatorString  = info['generatorName']
            version  = info['version']
            if badDataset(dsname,generatorString,version): continue
            availability = info['prodsysStatus']
            if config.onlyComplete and availability != u'ALL EVENTS AVAILABLE':
                print 'Skip incomplete dataset',dsname,availability
                continue
            nFiles = int(info['nFiles'])
            if nFiles>0:
                period = 'MC'
                xsec = 0.
                effic = 1.
                if info.has_key('period'):
                    period = info['period']
                else:
                    datasetNumber = int(info[u'datasetNumber'])
                    coveredids.add(datasetNumber)
                    # confirmed with AMI team that this should be enought, no need
                    # to re-implement get_dataset_xsec_effic for PyAMI5


                    # there are sometime problems in the propagation of these
                    # properties to the xAOD/derived datasets so go back in
                    # parentage to find the information
                    xsec = info[u'crossSection']
                    if info.has_key(u'approx_GenFiltEff'):
                        effic =  info[u'approx_GenFiltEff']
                    if xsec == u'NULL' or not info.has_key(u'approx_GenFiltEff'):
                        xsec,effic = genParamsFromParents(client,dsname,datasetNumber)

                    if not xsec: xsec = 0
                    if not effic:
                        print 'No approx_GenFiltEff found for',dsname,'set to 0 !!!!'
                        effic = 0
                    pass
                nevts = info['totalEvents']
                nfiles = info['nFiles']
                if not dsname.endswith('/'): dsname += '/'
                fout.write("%s %s %s %s %s %s\n" % (dsname,nevts,nfiles,period,xsec,effic))
        except KeyError as prop:
            print 'Missing property',prop,'for dataset ',dsname,'in AMI, skip'
    fout.close()

    if len(coveredids) == 0:
        if not config.prefix.startswith('data'): print 'Could not extract any channel IDs from datasets found, this is OK for data but suspicious for MC'
    else:
        for id in officialids:
            if not id in coveredids:
                print 'No dataset found for channel ',id

    pass


if __name__ == "__main__":
    main()
