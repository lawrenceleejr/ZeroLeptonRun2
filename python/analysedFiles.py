#!/bin/env python
# --------------------------------------------------------------------------
#
# Extract the list of files that were processed to produce the mini-tuples
#
# need to setup root and possibly rucio (depending on options) like
# lsetup root rucio
#
# If the argument is a pickle file (ends with .pkl):
#   -  Writes a pickle file with a dictionary of analysed files { dataset: {pfn: name}}
# Else, the argument is supposed to be a merged root file and the input
# file list is taken directly from the file metadata.
#
# if the --inputDS option is used, the file is supposed to contain the parent 
# dataset (used for mini-tuple production) and the script will compare 
# the analysed files to the files in the parent dataset (needs rucio)
#
# TODO: allow as input a list of files on disk 
# --------------------------------------------------------------------------

import os,sys,pickle,copy
from optparse import OptionParser
import ROOT

from rucio.client import Client
from rucio.common.exception import (DataIdentifierAlreadyExists, Duplicate, FileAlreadyExists, AccessDenied, ResourceTemporaryUnavailable,
                                    DataIdentifierNotFound, InvalidObject, RSENotFound, InvalidRSEExpression, DuplicateContent, RSEProtocolNotSupported,
                                    RuleNotFound, CannotAuthenticate, MissingDependency, RSEOperationNotSupported, UnsupportedOperation)


def printHelp():
    print 'analysedFiles.py [--inputDS file-with-parent-datasets] xxx.pkl'
    sys.exit(0)


def get_client(args):
    """
    Returns a new client object.
    """
    if args.auth_strategy == 'userpass':
        creds = {'username': args.username, 'password': args.password}
    else:
        creds = None

    try:
        client = Client(rucio_host=args.host, auth_host=args.auth_host,
                        account=args.account,
                        auth_type=args.auth_strategy, creds=creds,
                        ca_cert=args.ca_certificate, timeout=args.timeout)
    except CannotAuthenticate, error:
        logger.error(error)
        if not args.auth_strategy:
            if 'RUCIO_AUTH_TYPE' in os.environ:
                auth_type = os.environ['RUCIO_AUTH_TYPE']
            else:
                try:
                    auth_type = config_get('client', 'auth_type')
                except (NoOptionError, NoSectionError):
                    logger.error('Cannot get AUTH_TYPE')
                    sys.exit(FAILURE)
        if auth_type == 'x509_proxy':
            logger.error('Please verify that your proxy is still valid and renew it if needed.')
        sys.exit(FAILURE)
    return client

def extract_scope(did):
    # Try to extract the scope from the DSN
    if did.find(':') > -1:
        scope, name = did.split(':')[0], did.split(':')[1]
        if name.endswith('/'):
            name = name[:-1]
        return scope, name
    else:
        scope = did.split('.')[0]
        if did.startswith('user') or did.startswith('group'):
            scope = ".".join(did.split('.')[0:2])
        if did.endswith('/'):
            did = did[:-1]
        return scope, did

def extractAnalysedFiles(fileurl):
    analysedFiles = set()
    tfile = ROOT.TFile.Open(fileurl)
    if not tfile or tfile.IsZombie():
        print 'could not open',fileurl
        sys.exit(2)
        pass
    infilelist = tfile.Get('FileList_JobBookeeping_JobBookeeping')
    for f in  infilelist: analysedFiles.add(str(f).split()[0])
    tfile.Close()
    return copy.deepcopy(analysedFiles)

def listAnalysedFiles(filedict):
    analysedFiles = {}
    for (did,pfns) in filedict.iteritems():
        if not analysedFiles.has_key(did):
            analysedFiles[did] = {}
            pass
        for pfn in pfns:
            analysedFiles[did][pfn] = []
            tfile = ROOT.TFile.Open(pfn)
            if not tfile or tfile.IsZombie():
                print 'could not open',pfn
                sys.exit(2)
                pass
            infilelist = tfile.Get('FileList_JobBookeeping_JobBookeeping')
            for f in  infilelist: analysedFiles[did][pfn].append(str(f).split()[0])
            tfile.Close()
            pass
        pass
    return copy.deepcopy(analysedFiles)

def  filesInDS(inDS):
    class Struct:
        def __init__(self,**fields): self.__dict__=fields

    template = Struct(account=None, auth_host=None, auth_strategy=None, ca_certificate=None, certificate=None, host=None, password=None, timeout=None, username=None, verbose=False, which='')

    fargs = template
    fargs.all_states=False
    client = get_client(fargs)

    files = {}
    fileSet = set()
    for dsname in inDS:
        files[dsname]=[]
        scope, name = extract_scope(dsname)
        for f in client.list_files(scope=scope, name=name):
            files[dsname].append(f['name'])
            fileSet.add(f['name'])
            pass
        pass
    return (copy.deepcopy(files),copy.deepcopy(fileSet))

def matchFiles(analysedDict, rucioDict, config):
    for userDS,analysedFilesDict in analysedDict.iteritems():
        userDSname = userDS
        if ':' in userDS: userDSname = userDS.split(':')[1]
        analysedFilesList = []
        for flist in analysedFilesDict.values():
            for f in flist:
                analysedFilesList.append(f)
                pass
            pass
        analysedFilesList.sort()
        #print userDSname,analysedFilesList

        # now try to match to a dataset in the rucio list
        fields = userDSname.split('.')[2:]
        matchedDS = ''
        for parentDS,rucioFiles in rucioDict.iteritems():
            parentDSname = parentDS
            if ':' in parentDS: parentDSname = parentDS.split(':')[1]
            parentFields = parentDSname.split('.')
            #print fields,parentFields
            if config.strictMatch:
                # strict match: try to match exactly the first five field and loosely the sixth
                if fields[:5] == parentFields[:5] and parentFields[5] in fields[5]:
                    if matchedDS != '':
                        print 'error, this dataset was already matched to ',matchedDS
                        pass
                    matchedDS = parentDS
                    pass
                pass
            else:
                # for MC the description field is in general altered to limit the length of the output dataset do don't require the third field to match
                if fields[:2] == parentFields[:2] and fields[3:5] == parentFields[3:5] and parentFields[5] in fields[5]:
                    #print 'matched ',userDSname,parentDSname
                    if matchedDS != '':
                        print 'error, this dataset was already matched to ',matchedDS
                        pass
                    matchedDS = parentDS
                    pass
                pass
            pass
        
        if matchedDS =='':
            print 'could not find a matching dataset for ',userDSname
            sys.exit(2)
        rucioFileList = rucioDict[matchedDS]
        rucioFileList.sort()
        if rucioFileList != analysedFilesList:
            print 'analysed files do not match the file list from rucio !',userDSname
            print analysedFilesList
            print rucioFileList
        else:
            print userDSname,'validated'
        pass
    pass

# extract arguments    
parser = OptionParser()
parser.add_option("--inputDS", dest="inputDS", 
                  help="file containing the list of input datasets, triggers comparison with rucio",default="") 
parser.add_option("--strictMatch", dest="strictMatch", 
                      help="stricter matching criteria between parent and user dataset (default=%default)", 
                  action='store_true', default=False)
 
(config, args) = parser.parse_args(sys.argv[1:])
if len(args)<1:
    printHelp()

# read dictionary from input file
inputfile = args[0]

if inputfile.endswith('.pkl'):
    if not os.path.isfile(inputfile):
        print "can't read file",inputfile
        sys.exit(1)
    picklefile = open(inputfile,'rb')
    filedict = pickle.load(picklefile)
    picklefile.close()

    # extract list of files processed to create the mini-tuples
    analysedFiles = listAnalysedFiles(filedict)
    outfile = open('analysedFiles_'+inputfile,'wb')
    pickle.dump(analysedFiles,outfile)
    outfile.close()
else:
    analysedFiles = extractAnalysedFiles(inputfile)

if config.inputDS != "":
    inDS = []
    for line in open(config.inputDS,'r'):
        line = line.strip()
        if line.startswith('#'): continue
        if len(line)>0: inDS.append(line)
        pass
    rucioFiles, rucioFileSet = filesInDS(inDS)
    if type(analysedFiles) is dict:
        matchFiles(analysedFiles,rucioFiles,config)
    else:
        print 'len',len(analysedFiles),len(rucioFileSet),analysedFiles==rucioFileSet
        print 'analysed but not in datasets',analysedFiles-rucioFileSet
        print 'in datasets but not analysed',rucioFileSet-analysedFiles
        pass



