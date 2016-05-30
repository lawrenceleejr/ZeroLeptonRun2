#!/usr/bin/env python

import os, sys
import os.path
from optparse import OptionParser
import subprocess
import math
import time

def parseCommonOptions() :
    parser = OptionParser()
    parser.add_option("--submitDir", help   = "dir to store the output", default="submit_dir")
    parser.add_option("--inputDS", help     = "must be the name of an existing .ds file in ZeroLeptonRun2/data/", default="")
    parser.add_option("--nChunks", help  = "number of chunks to split this .ds file", default="1")
    parser.add_option("--iChunk", help   = "which chunk of this .ds file do you want to run?", default="0")

    parser.add_option("--gridUser", help    = "gridUser"  , default= os.environ.get("USER"))
    parser.add_option("--gridTag", help     = "gridTag"   , default= time.strftime("%y%m%d"))

    parser.add_option("--tag", help     = "factory tag"   , default= '')
    parser.add_option("--trunk", help     = "factory trunk"  , action="store_true" , default=False)

    parser.add_option('--doOverwrite', help = "Overwrite submit dir if it already exists",action="store_true", default=False)
    return parser


def chopUpTheDataSetFile(inputFile, nChunks, iChunk):

    allDatasets = []
    with open(inputFile) as f:
        allDatasets = f.read().splitlines()

    nDatasets = len(allDatasets)
    nDatasetsPerChunk = int(math.ceil(nDatasets/float(nChunks) ) )

    subsetDatasets = allDatasets[iChunk*nDatasetsPerChunk:(iChunk+1)*nDatasetsPerChunk]

    f = file("tmpInputDataSets.ds", "w")
    f.write("\n".join(subsetDatasets)+"\n")


def main():

    myCWD = os.getcwd()

    #########################################################################
    # Parse the command line options

    parser = parseCommonOptions()

    (options, args) = parser.parse_args()


    #########################################################################
    # Set up sandboxed directory structure

    if not os.path.exists(options.submitDir):
        os.makedirs(options.submitDir)
    else:
        raise Exception('There already exists a run directory named %s. Either get rid of it or change the name with --submitDir.'%options.submitDir)

    os.makedirs(options.submitDir+"/run/")
    os.makedirs(options.submitDir+"/tmp/")


    #########################################################################
    # Configure what factory to checkout into sandbox

    factoryPackageToCheckOut = ""
    compileCommand = ""
    if options.trunk:
        factoryPackageToCheckOut = "trunk/"
        compileCommand = "rc compile"
    else:
        factoryPackageToCheckOut = "tags/%s"%options.tag

    if factoryPackageToCheckOut == "":
        raise Exception('You either need to set --trunk or --tag to run!')

    # Get the desired factory code and install
    os.chdir(options.submitDir + '/run/')
    os.system('svn co svn+ssh://svn.cern.ch/reps/atlasoff/PhysicsAnalysis/SUSYPhys/Factory/ZeroLeptonRun2/%s ZeroLeptonRun2'%factoryPackageToCheckOut)



    #########################################################################
    # Let's split up the input dataset files

    if os.path.isfile("ZeroLeptonRun2/data/%s"%options.inputDS):
        chopUpTheDataSetFile("ZeroLeptonRun2/data/%s"%options.inputDS,
            int(options.nChunks),
            int(options.iChunk)
            )
    elif os.path.isfile("../../%s"%options.inputDS):
        chopUpTheDataSetFile("../../%s"%options.inputDS,
            int(options.nChunks),
            int(options.iChunk)
            )
    else:
        os.system("echo %s > tmpInputDataSets.ds"%option.inputDS)

    ##############################################################
    # The heavy lifting...

    os.system(
        """

        #########################################################################
        # Setup the environment from scratch or else things get angry
        # Make sure you do this on a system with cvmfs access!


        alias setupATLAS='source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh'


        #########################################################################
        # Get a new grid certificate

        voms-proxy-init -voms atlas


        #########################################################################
        # Get and compile all the dependencies for this tag/trunk
        # It will try to re-checkout ZeroLeptonRun2 but will fail because it's done above

        source ZeroLeptonRun2/INSTALL

        # Depending on whether you're using a tag or not - it'll compile before sending.
        {0}

        echo "Going to run on..."
        cat tmpInputDataSets.ds

        read -p "Press [Enter] to start submission if everything looks good. Else break (ctrl-C)."


        #########################################################################
        # Setup Panda Client

        localSetupPandaClient


        ##########################################################################
        # Submit!

        echo "ZeroLeptonRun2/python/launch_jobs_on_grid.py --prunopts='--nGBPerJob=10 --excludedSite=ANALY_UIO,SIGNET,LUNARC ' --prefix='user.{1}.'  --suffix='_{2}' --tmpDir=../tmp/ --inDSfile tmpInputDataSets.ds"
        ZeroLeptonRun2/python/launch_jobs_on_grid.py --prunopts='--nGBPerJob=10 --excludedSite=ANALY_UIO,SIGNET,LUNARC ' --prefix='user.{1}.'  --suffix='_{2}' --tmpDir=../tmp/ --inDSfile tmpInputDataSets.ds


        """.format(
            compileCommand,
            options.gridUser,
            options.gridTag,
            options.inputDS
            )
        )






if __name__ == "__main__":
    main()

