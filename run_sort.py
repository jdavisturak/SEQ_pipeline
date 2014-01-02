#!/usr/bin/env python
""" Run sam_to_sam.py on a (directory of) .sam files
Support for parallel threads

Usage: -i input_dir -o outdir [-O OPTIONS] [-N num_processors (default 1)]


"""

from __future__ import with_statement
import os
import sys
import glob
import re
import subprocess
from optparse import OptionParser
from multiprocessing import Process, Queue, current_process

## Function to generate cmd to run on one file

def call_sort(infile,  outdir,ref,options=''):
    cmd = ["sam_to_bam.py",'--convert',infile, outdir,'--ref=%s' % ref]

    if options != '' and options != None:
        cmd += options.split(' ')
    
    return [cmd]


## Worker function for parallelization
def worker(work_queue, done_queue):
    try:
        for cmd in iter(work_queue.get, 'STOP'):
            print cmd
            #f = open(cmd[1],'w') # stdout will be mapped to file
            p = subprocess.Popen(cmd[0])
            print "pid: %s" % p.pid
            retval = p.wait()
            #f.close()
            done_queue.put("%s - %s got %s, %s." % (current_process().name, cmd, p.pid, retval))

    except Exception, e:
        done_queue.put("%s failed on %s with: %s" % (current_process().name, cmd, e.message))
    return True


def main(indir, outdir, ref, option_string, num_processors):


    try:
        num_processors = int(num_processors)
    except ValueError, TypeError:
        num_processors = 1
    if num_processors < 1:
        num_processors = 1

    if option_string is None:
        option_string = '' 

    if not os.path.exists(outdir):
        try:
            os.makedirs(outdir)
        except OSError:
            assert os.path.isdir(outdir)

    assert os.path.isdir(indir)
    
    # Look through the input directory for .sam files
    sams = glob.glob("%s/*.sam" % indir)
    
    ## parallel stuff:
    work_queue = Queue()
    done_queue = Queue()
    processes = []

    for sam in sams:
        prefix = re.sub('\.sam$','',os.path.basename(sam))
        cmd = call_sort(sam, outdir, ref, option_string)
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
                        
    print "Completed splice junction counting"


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-i", "--input_dir", dest="input_dir", action="store",
                      default=False)
    parser.add_option("-o", "--outdir", dest="outdir", action="store",
                      default=None)
    parser.add_option("-O", "--Options", dest="option_string", action="store",
                      default=None)
    parser.add_option("-N", "--num_processors", dest="num_processors", action="store",
                      default=None)
    parser.add_option("-r", "--ref", dest="ref", action="store",
                      default=None)
    parser.add_option("-H", "--Help", dest="Help", action="store",
                      default=None)
    (options, args) = parser.parse_args()
    if options.Help or options.outdir is None or options.input_dir is None:
        print __doc__
        sys.exit()
    main(options.input_dir, options.outdir,options.ref, options.option_string, options.num_processors)

