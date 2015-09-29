#!/usr/bin/env python
"""
A script to validate a mini-tuple production by comparing the files that
were run on to files known to AMI
"""


import os,sys
import ROOT
try:
    import pyAMI.client
    import pyAMI.atlas.api
except:
    print 'Could not import AMI modules, did you "localSetupPyAMI" ?'

# check arguments
if len(sys.argv) < 3:
    print 'zlProdXcheck.py directory-with-mini-tuple  input-sample-dataset'
    sys.exit(0)

dir = sys.argv[1]
ds  = sys.argv[2]

if not os.path.isdir(dir):
    print dir,'is not a directory !'
    sys.exit(1)

#-------------------------- extract information from mini-tuples root files
# build list of root files
fnames = []
for root, dirnames, filenames in os.walk(dir):
    for fname in filenames:
        fnames.append(fname)

# extract the list of files read from the bookeeping information
readfiles = []
readGUIDs = []
for fname in fnames:
    tfile = ROOT.TFile.Open(dir+'/'+fname)
    if not tfile or tfile.IsZombie():
        print 'could not open root file',fname
        sys.exit(2)
    flist = tfile.Get('FileList_JobBookeeping_JobBookeeping')
    if not flist:
        print 'Could not extract file list from ',fname
        continue
    for x in flist:
        (name, guid, evts) =str(x).split()
        readfiles.append(name)
        readGUIDs.append(guid)
    tfile.Close()
#    break

readfiles.sort()
readGUIDs.sort()
# check for duplicates
for i in range(len(readfiles)-1):
    if readfiles[i+1] == readfiles[i]:
        print 'duplicate input file ',readfiles[i]

#print readfiles
#print readGUIDs

#-------------------------- extract information from AMI
# AMI client connection
client = pyAMI.client.Client('atlas')
pyAMI.client.endpoint = 'main'
pyAMI.atlas.api.init()
infos = pyAMI.atlas.api.get_dataset_info(client,ds)
amifileinfos = pyAMI.atlas.api.list_files(client,ds)
amifiles = []
amiGUIDs = []
for fileinfos in amifileinfos:
    amifiles.append(fileinfos[u'LFN'])
    amiGUIDs.append(fileinfos[u'fileGUID'])
amifiles.sort()
amiGUIDs.sort()

#-------------------------- compare information for root files and from AMI
print 'file sizes',len(readfiles),len(amifiles),readfiles==amifiles
for f in readfiles:
    if not f in amifiles:
        print f,'read but not found in AMI'
for f in amifiles:
    if not f in readfiles:
        print f,'in AMI but not processed '

