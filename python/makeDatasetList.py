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
lsignals = ['Herwigpp_UEEE3_CTEQ6L1_MSUGRA','Herwigpp_UEEE4_CTEQ6L1_Gtt','SM_SS','SM_GG','SM_SG']

# List of keywords used to remove non-used datasets
filters = ['Pythia_GGM','PG11','DGemt','DGnoL','_Wprime','_Zprime','PythiaB','_3Leptons','LeptonPhotonFilter','ADDGraviton','Sherpa_CT10_SingleTop','BlackMaxPythia8','Charybdis2Pythia8','QBHPythia8','_WimpPair_','_ZH125_','_ggH125_','_VBFH125_','_WH125_','ParticleGenerator_']

# bad tags
# https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/DC14DataMCSampleInfo#13_TeV_Monte_Carlo
# https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/AtlasProductionGroupMC14a#DC14_Datasets
# what about r5869, r5899, r5894, r5910, r5898, r5812
# preproduction samples +  no pileup 
badtags = ['_r5720_','_r5721_'] + ['_r5803_','_r5804_','_r5862_']


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
    parser.add_option("--prefix", dest="prefix", help="Prefix to search datasets in AMI",default="mc14_13TeV") 
    parser.add_option("--datatype", dest="datatype", help="data type",default="%.merge.AOD%") 
    parser.add_option("--suffix", dest="suffix", help="Suffix appended to the output file name",default="") 
    parser.add_option("--tag", dest="tag", help="Production tag",default="r5853") 
    parser.add_option("--server", dest="server", help="AMI server (main or replica)",default="main") 
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
        else:
            print '--official is only supported for mc12_8TeV, mc14_8TeV and mc14_13TeV'
            sys.exit(1)
        if config.sample:
            officialids = mcsl.__dict__[str(config.sample)]
        else:
            officialids = mcsl.__dict__["lbaseline"]
            if config.official:
                officialids += mcsl.__dict__["lalt"]

    # get all datasets matching prefix & tag and then filter them
    from pyAMI.atlas.api import get_dataset_info, list_datasets

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
    fout = open('datasets.txt','w')
    for info in dsinfos:
        try:
            dsname = info['logicalDatasetName']
            generatorString  = info['generatorName']
            version  = info['version']
            if badDataset(dsname,generatorString,version): continue
            availability = info['prodsysStatus']
            nFiles = int(info['nFiles'])
            if nFiles>0:
                period = 'MC'
                xsec = 0.
                effic = 1.
                if info.has_key('period'):
                    period = info['period']
                else:
                    #(xsec, effic) = get_dataset_xsec_effic(client,info.info['logicalDatasetName'])
                    # confirmed with AMI team that this should be enought, no need
                    # to re-implement get_dataset_xsec_effic for PyAMI5
                    xsec = info[u'crossSection']
                    effic =  info[u'approx_GenFiltEff']
                nevts = info['totalEvents']
                nfiles = info['nFiles']
                if not dsname.endswith('/'): dsname += '/'
                fout.write("%s %s %s %s %s %s\n" % (dsname,nevts,nfiles,period,xsec,effic))
        except KeyError as prop:
            print 'Missing property',prop,'for dataset ',dsname,'in AMI, skip'
    fout.close()
    pass


if __name__ == "__main__":
    main()
