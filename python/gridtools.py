#from Brad Axen

import os
import re
import types
import argparse
import datetime

from rucio_down import download

# Import pbook
try:
    tmp, __name__ = __name__, 'pbook'
    execfile(os.path.join(os.environ["ATLAS_LOCAL_ROOT_BASE"], "x86_64/PandaClient/current/bin/pbook") )
    __name__ = tmp
    pbook = PBookCore(False, False, False) #enforceEnter, verbose, restoreDB
except ImportError:
    print("Failed to load PandaClient, please set up locally")
    sys.exit(1)

class bcolors:
    BLUE = '\033[95m'
    LIGHT = '\033[94m'
    ORANGE = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    ENDC = '\033[0m'

def is_listy(x):
    return any(isinstance(x, t) for t in [types.TupleType,types.ListType])


def apply_each(func, args, status='finished'):
    pbook.sync()

    start_time = (datetime.datetime.utcnow() - datetime.timedelta(days=args.history)).strftime('%Y-%m-%d %H:%M:%S')

    jobs = Client.getJobIDsJediTasksInTimeRange(start_time)[2]
    for jobid, info in jobs.iteritems():
        if is_listy(status):
            if not any(info['status'] == s for s in status):
                continue
        elif info['status'] != status:
            continue
        job = pbook.getJobInfo(jobid)
        if not job: continue
        if args.match in job.jobName and re.match(args.regex, job.jobName):
            func(job)

def retry_job(job):
    if 'failed' in job.jobStatus:
        print bcolors.LIGHT + bcolors.BOLD + "Found finished job with failed subjob(s): " + bcolors.ENDC + bcolors.ORANGE + job.jediTaskID + bcolors.ENDC + " " + job.jobName
        pbook.retry(job.groupID)

def kill_job(job):
    print bcolors.LIGHT + bcolors.BOLD + "Found running job: " + bcolors.ENDC + bcolors.ORANGE + job.jediTaskID + bcolors.ENDC + " " + job.jobName
    pbook.kill(job.groupID)

def download_job(job):
    print bcolors.LIGHT + bcolors.BOLD + "Found completed job: " + bcolors.ENDC + bcolors.ORANGE + job.jediTaskID + bcolors.ENDC + " " + job.jobName
    # Look up output datasets and download
    for ds in job.outDS.split(','):
        download(ds)

def print_job(job):
    print bcolors.LIGHT + bcolors.BOLD + job.jediTaskID + bcolors.ENDC + job.jobName
    print "    ", job.inDS
    print "    ", job.outDS


def main(args):
    parser = argparse.ArgumentParser('Grid tool to kill, retry, or download jobs in batches.')
    subparsers = parser.add_subparsers(help='')

    subs = []
    subs.append(subparsers.add_parser('retry', help='Retry failed jobs.'))
    subs[-1].set_defaults(func=lambda args: apply_each(retry_job, args, status='finished'))

    subs.append(subparsers.add_parser('kill', help='Kill running jobs.'))
    subs[-1].set_defaults(func=lambda args: apply_each(kill_job, args, status=['running','scouting','ready','defined']))

    subs.append(subparsers.add_parser('download', help='Download jobs marked "done" using rucio.'))
    subs[-1].set_defaults(func=lambda args: apply_each(download_job, args, status=['finished','running','done']))

    subs.append(subparsers.add_parser('print', help='Print job information.'))
    subs[-1].set_defaults(func=lambda args: apply_each(print_job, args, status=['finished','running','done']))

    # Set common arguments on each so that they print help correctly.
    for p in subs:
        p.add_argument('--history', default=7, type=int, help='How far back in the history to check for jobs.')
        p.add_argument('--match', default='', help='Apply subcommand to jobs which contain this string.')
        p.add_argument('--regex', default='.*', help='Apply subcommand to jobs which match this regex.')

    args = parser.parse_args(args)
    args.func(args)

if __name__ == "__main__":
    main(sys.argv[1:])
