#!/usr/bin/env python

########### Initialization ######################################
##
##

import ROOT
import logging
import shutil
import os, sys
import itertools
import array

import mc15_13TeV_MCSampleList as samplelist
import discoverInput

logging.basicConfig(level=logging.INFO)
from optparse import OptionParser

import atexit
@atexit.register
def quite_exit():
	ROOT.gSystem.Exit(0)

ROOT.gROOT.SetBatch()
# from multiprocessing import Pool
import multiprocessing as mp

def main():

	logging.info("loading packages...")
	ROOT.gROOT.Macro("$ROOTCOREDIR/scripts/load_packages.C")

	parser = OptionParser()
	parser.add_option("--inDir", help   = "dir with output", default="/afs/cern.ch/user/l/leejr/work/public/fromGrid/")
#	parser.add_option("--nproc", help     = "number of parallel processes", default="4"  )
	parser.add_option("--selection", help     = "selection string for skimming", default="1"  )
	parser.add_option("--outputDir", help ="where to write the output merged ntuples" , default = "output")
	parser.add_option("--dryRun",  help = 'do a dry run, without actually launching the merge job', action='store_true', default=False)

	(options, args) = parser.parse_args()


	print options
	search_directories = [options.inDir]
	print search_directories

#these names should be as they are called in your samplelist
	outputSampleNames = [
		"QCD",
		"Top",
		"Wjets",
		"ZMassiveCB",
		"DibosonMassiveCB",
		"GammaMassiveCB",

	        "SS_direct",
	        "GG_direct",
		"Data2015",
		"Data2016",
	]

	try:
		os.stat(options.outputDir)
	except:
		os.mkdir(options.outputDir)
	try:
		os.stat(options.outputDir+"/signal")
	except:
		os.mkdir(options.outputDir+"/signal")


	logging.info("creating new sample handler")
	sh_all = ROOT.SH.SampleHandler()

	discoverInput.discover(sh_all, search_directories, "*"  )#todo add options
	print 'length of sh_all : ', len(sh_all)

	logging.info("adding my tags defined in discoverInput.py")
	discoverInput.addTags(sh_all)

	print sh_all

	ROOT.SH.readSusyMetaDir(sh_all,"$ROOTCOREBIN/data/SUSYTools")
	ROOT.SH.readSusyMetaDir(sh_all,"$ROOTCOREBIN/data/SUSYTools/mc15_13TeV/")



	processes = []
	for outputSampleName in outputSampleNames:
		sh = sh_all.find(outputSampleName)
		p = mp.Process(target=processTheSH,
			       args=(sh,
			     	  options.outputDir.strip('/'),
			     	  outputSampleName,
			     	  options.selection,
				  )
			       )
		processes.append(p)

	print "will run processes : " , processes
	print "over sample handlers :", outputSampleNames

	if not options.dryRun :
		for p in processes :
			p.start()

		for p in processes :
			p.join()


	return


def processTheSH( sh,
		  outputDirectory = "output",
		  sampleName = "OTHER.root",
		  selection  = "1.",
		  treePrefix = "",
		  ) :
#	print len(sh)

	## Split up samplehandler into per-BG SH's based on tag metadata

	treesToProcess = []
	filesToEventuallyHadd = []

	for sample in sh :
		sample_name = sample.getMetaString("sample_name")
		# print sample_name
		dsid = int(sample.getMetaDouble("dsid"))

		if len(treesToProcess) == 0:
			treesToProcess = getListOfTreeNames(sample, treePrefix)

		attachCounters(sample)

		tmpOutputDirectory = os.path.join(outputDirectory, "tmpOutput")
#		print tmpOutputDirectory

		try:
			os.stat(tmpOutputDirectory)
		except:
			os.mkdir(tmpOutputDirectory)

		outputSampleFileName = "%s/%s.root"%(tmpOutputDirectory, dsid)
		filesToEventuallyHadd.append(outputSampleFileName)
#		print 'in loop, files to hadd',  filesToEventuallyHadd
		outputSampleFile = ROOT.TFile(outputSampleFileName,"RECREATE")

		print "Starting"
#		print os.stat(outputSampleFileName).st_size
#		print treesToProcess

		for itree in treesToProcess:
#			print itree
			if ("SRAllNT" not in itree) : continue
			sh.setMetaString("nc_tree", itree)
			outputSampleFile.cd()
			mytree = sample.makeTChain().Clone(itree)

#			print mytree, mytree.GetEntries()
			if mytree.GetEntries() and getNormFactor(sample, bool(sh.find(ROOT.SH.TagList(':'.join(['Data2015','Data2016' ]))))):#second arg to check if data
				try:
					outputTree = ROOT.addBranch( mytree, getNormFactor(sample) , selection)
				except:
					print 'failed to add branch'
					continue
#				print outputTree.GetEntries()
				outputTree.Write()
#				print "Saved tree %s with %s events . . ." % ( outputTree.GetName(), outputTree.GetEntries() )
#			print os.stat(outputSampleFileName).st_size

#		print os.stat(outputSampleFileName).st_size
		print "WRITING FILE...."
		outputSampleFile.Write()
#		print os.stat(outputSampleFileName).st_size
		outputSampleFile.Close()
#		print os.stat(outputSampleFileName).st_size

#	print 'filesToEventuallyHadd' , filesToEventuallyHadd

	if sampleName == "GG_direct" or sampleName == "SS_direct":
		for myfile in filesToEventuallyHadd:
			os.system('cp %s %s/signal/.'% (myfile,outputDirectory)  )
 	else :
		os.system('hadd -O -f %s/%s.root %s'%
			  (outputDirectory, sampleName, " ".join(filesToEventuallyHadd) )
			  )

	return

#To scale the histograms in the files after the event loop is done...
def getNormFactor(sample, isData = False):
	if isData : return 1.

	tempxs = sample.getMetaDouble("nc_xs") * sample.getMetaDouble("kfactor") * sample.getMetaDouble("filter_efficiency")

#	print "Norm weight for %s is %f/(%f or %f)"%(sample.getMetaString("short_name"), tempxs, sample.getMetaDouble("nc_nevt"), sample.getMetaDouble("nc_sumw"))
	m_eventscaling = tempxs
	if sample.getMetaDouble("nc_nevt"):
		m_eventscaling /= sample.getMetaDouble("nc_nevt") if "jetjet" in sample.getMetaString("short_name") else sample.getMetaDouble("nc_sumw")
	else:
		print "nevt not SET! normweight = 0 !!!"
		m_eventscaling = 0.
	return m_eventscaling


addBranchCode = """
TTree * addBranch(TTree* tree, float normalization, TString selection="1"){

		TTree * newtree = tree->CopyTree(selection);
		float normweight = normalization;
		TBranch * bnormweight = newtree->Branch("normweight",&normweight,"normweight/F");
		int nevents = newtree->GetEntries();

		for (Long64_t i=0;i<nevents;i++) {
			newtree->GetEntry(i);
			bnormweight->Fill();
			if(i%10000==0) cout<< i << " of " << nevents << endl;
		}

		return newtree;
}
"""

ROOT.gInterpreter.Declare(addBranchCode)

# This python function is replaced by the c++ function above for speed
# def addLeaves(tree,normalization,selection="1"):
# 	return 0
# 	events = tree.GetEntries()
# 	leaves = "normweight/D"
# 	leafValues = array.array("d", [0.])
# 	# newtree = tree.CloneTree(0)
# 	newtree = tree.CopyTree(selection)
# 	newBranch = newtree.Branch( "normweight" , leafValues, leaves )
# 	for i in xrange(events):
# 		tree.GetEntry(i)
# 		leafValues[0] = ROOT.Double(normalization)
# 		newtree.Fill()
# 		if i % 10000 == 0:
# 			print "%s of %s: %s" % (i,events,leafValues)
# 	return newtree


def attachCounters(sample):
#	print sample
#	print sample.makeFileList()

	nevt = 0
	sumw = 0

	for fname in  sample.makeFileList() :
#		print fname
		f = ROOT.TFile(fname )
#		print f
		nevt += f.Get("Counter_for_ZeroLeptonCounterSRAll").GetBinContent(1)
		sumw += f.Get("Counter_for_ZeroLeptonCounterSRAll").GetBinContent(2)

#	print nevt
#	print sumw


	sample.setMetaDouble("nc_nevt",nevt)
	sample.setMetaDouble("nc_sumw",sumw)

def getListOfTreeNames(sample, treePrefix = "" ):
	f = ROOT.TFile(sample.fileName(0) )
	listOfTrees = [key.GetName() for key in f.GetListOfKeys() if treePrefix in key.GetName()]
	return listOfTrees




if __name__ == "__main__":
	main()

