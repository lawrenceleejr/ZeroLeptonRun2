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
    print 'usage:  findDuplicateEvent.py input-mini-tuple'
    print ''
    print 'Note that it is known that mc15 has duplicate (RunNumber,EventNumber) events with different kinematics.'
    sys.exit(0)

def getListOfTrees(tfile):
    """Find list of trees."""
    keys=tfile.GetListOfKeys()
    mykeys=[]
    for k in keys:
        kname=k.GetName()            
        obj=tfile.Get(kname)
        if not obj.IsA().InheritsFrom( ROOT.TTree.Class() ): continue
        mykeys.append(kname)
    # return list made unique
    return list(collections.OrderedDict.fromkeys(mykeys))

inputfile = ''
if len(sys.argv)>1:
    arg = sys.argv[1]
    if arg == '-h' or arg =='--help':
        help()
    else:
        inputfile = arg
else:
    help()


tfile = ROOT.TFile.Open(inputfile)
if not tfile or tfile.IsZombie():
    print 'could not open file ',inputfile
    sys.exit(0)

for treename in getListOfTrees(tfile):
    print 'Checking tree ',treename
    rundict = {}
    tree = tfile.Get(treename)
    if not tree:
        print 'Could not file a tree named',treename,'in',inputfile
        continue
    b_NTVars = tree.GetBranch('NTVars')
    l_run = tree.GetLeaf('RunNumber')
    l_evt = tree.GetLeaf('EventNumber')
    for i in range(tree.GetEntries()):
        b_NTVars.GetEntry(i)
        run =  int(l_run.GetValue())
        evt = int(l_evt.GetValue())
        if not rundict.has_key(run):
            rundict[run] = []
        rundict[run].append(evt)

    for (run,evts) in rundict.iteritems():
        evts.sort()
        for i in range(len(evts)-1):
            if evts[i+1] == evts[i]:
                print 'duplicate ',run,evts[i]



