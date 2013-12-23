#!/usr/bin/env python
""" Run htseq-count on a sam file
Support for parallel threads

Usage: -i input_dir -o outdir -g gtf_file [-s STRANDED] [-o OPTIONS] [-N num_processors (default 1)]

works as of HTSeq 0.5.4p5

"""

from __future__ import with_statement
import os
import sys
import glob
import re
import subprocess
from optparse import OptionParser
from multiprocessing import Process, Queue, current_process

from demultiplex import *

## Function to run htseq on one file
#def run_htseq(infile, gtf, outfile, strandedness=None,options=None):


def call_htseq(infile, gtf, outfile, strandedness=None,options=''):
    cmd = ["htseq-count",infile, gtf]

    strandedness = {'yes': '-s yes', 'no': '-s no', 'reverse': '-s reverse'  }.get(strandedness, '')
    if strandedness != '':
        cmd += strandedness.split(' ')
    
    if options != '' and options != None:
        cmd += options.split(' ')
    
    return [cmd, outfile]


## Functions to try to make this parallelized ...
def worker(work_queue, done_queue):
    try:
        for cmd in iter(work_queue.get, 'STOP'):
            print cmd
            f = open(cmd[1],'w') # stdout will be mapped to file
            p = subprocess.Popen(cmd[0], stdout=f)
            print "pid: %s" % p.pid
            retval = p.wait()
            f.close()
            done_queue.put("%s - %s got %s, %s." % (current_process().name, cmd, p.pid, retval))

    except Exception, e:
        done_queue.put("%s failed on %s with: %s" % (current_process().name, cmd, e.message))
    return True


def main(options):

    outdir = options.outdir
    indir = options.input_dir

    try:
        num_processors = int(options.num_processors)
    except ValueError, TypeError:
        num_processors = 1
    if num_processors < 1:
        num_processors = 1

    if options.option_string is None:
        options.option_string = '' 

    if not os.path.exists(outdir):
        try:
            os.makedirs(outdir)
        except OSError:
            assert os.path.isdir(outdir)

    assert os.path.isdir(indir)
    
    # Look through the input directory for .sam files
    sams = glob.glob("%s/*.sam" % indir)
    
    # Get sample prefix
    #prefix = [re.sub('\.sam$','',os.path.basename(sam)) for sam in sams]

    ## parallel stuff:
    work_queue = Queue()
    done_queue = Queue()
    processes = []

    for sam in sams:
        prefix = re.sub('\.sam$','',os.path.basename(sam))
        cmd = call_htseq(sam, options.gtf, "%s/%s.counts.tab" % (outdir, prefix), options.stranded,options.option_string)
        work_queue.put(cmd)

    for w in xrange(num_processors):
        p = Process(target=worker, args=(work_queue, done_queue))
        p.start()
        processes.append(p)
        work_queue.put('STOP')

    for p in processes:
        p.join()

    done_queue.put('STOP')

    for status in iter(done_queue.get, 'STOP'):
        print status

    print "Completed counting"



# if __name__ == "__main__":
#     parser = OptionParser()
#     parser.add_option("-i", "--input_dir", dest="input_dir", action="store",
#                       default=False)
#     parser.add_option("-o", "--outdir", dest="outdir", action="store",
#                       default=None)
#     parser.add_option("-g", "--gtf_file", dest="gtf", action="store",
#                       default=None)
#     parser.add_option("-s", "--stranded", dest="stranded", action="store",
#                       default=None)
#     parser.add_option("-O", "--Options", dest="option_string", action="store",
#                       default=None)
#     parser.add_option("-N", "--num_processors", dest="num_processors", action="store",
#                       default=None)
#     parser.add_option("-H", "--Help", dest="Help", action="store",
#                       default=None)
#     (options, args) = parser.parse_args()
#     if options.Help or options.outdir is None or options.input_dir is None or options.gtf is None:
#         print __doc__
#         sys.exit()
#     main(options)

