#!/usr/bin/env python
"""
Find duplicate event based on RunNumber and EventNumber

Note that for MC it might give false positives since a prodsys bug lead to 
production of events with different truth kinematics with the same run and
event numbers.
"""
import sys,copy,collections
import ROOT

def help():
    """Print usage and exit"""
    print 'usage:  findDuplicateEvent.py input-xAOD-files'
    print ''
    print 'Note that it is known that mc15 has duplicate (RunNumber,EventNumber) events with different kinematics.'
    sys.exit(0)

inputfiles = []
if len(sys.argv)>1:
    arg = sys.argv[1]
    if arg == '-h' or arg =='--help':
        help()
    else:
        inputfiles = sys.argv[1:]
else:
    help()

rundict = {}
for fname in inputfiles:
    tfile = ROOT.TFile.Open(fname)
    if not tfile or tfile.IsZombie():
        print 'could not open root file',fname
        sys.exit(2)
    tree = tfile.Get('CollectionTree')
    if not tree:
        print 'Could not find tree CollectionTree in',fname
        sys.exit(2)

    b_EventInfoAux = tree.GetBranch('EventInfoAux.')
    l_run = tree.GetLeaf('EventInfoAux.runNumber')
    l_evt = tree.GetLeaf('EventInfoAux.eventNumber')

    for i in range(tree.GetEntries()):
        b_EventInfoAux.GetEntry(i)
        run =  int(l_run.GetValue())
        evt = int(l_evt.GetValue())
        if not rundict.has_key(run):
            rundict[run] = []
        rundict[run].append(evt)
    tfile.Close()

for (run,evts) in rundict.iteritems():
    evts.sort()
    for i in range(len(evts)-1):
        if evts[i+1] == evts[i]:
            print 'duplicate ',run,evts[i]



