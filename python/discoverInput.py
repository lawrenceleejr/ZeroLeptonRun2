import ROOT
import logging
import shutil
import os
import mc15_13TeV_MCSampleList as samplelist

def discover(sh, search_directories, pattern="*"):

	# scan for datasets in the given directories
	for directory in search_directories:
	    ROOT.SH.ScanDir().samplePattern(pattern).scan(sh, directory)


	logging.info("%d different datasets found scanning all directories", len(sh))

	return sh



def addTags(sh_all):
	for sample in sh_all:
		sample_name = sample.getMetaString("sample_name")
#		print sample_name

		if 'data15_13TeV' in sample_name :
			sample.addTag('Data2015')
			return
		if 'data16_13TeV' in sample_name :
			sample.addTag('Data2016')
			return

                dsid        = int( sample_name.split(".")[3])
		print dsid
		sample.setMetaDouble("dsid", dsid)

#		print samplelist

		for name, dsidrange in (vars(samplelist)).iteritems() :
			if isinstance(dsidrange, list ) :
				if not name.startswith('__') :
					if dsid in dsidrange :
						print 'for ', dsid ,' , adding tag ' , name
						sample.addTag(name)

		print sample
