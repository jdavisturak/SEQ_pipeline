#!/usr/bin/env python
""" Run bam to bigwig script on a (directory of) .sorted.bam files
Support for parallel threads

Usage: -i input_dir -o outdir -r ChromInfo [-e 'chr'] [-O OPTIONS] [-N num_processors (default 1)]
-e 'chr' converts from Ensembl to UCSC notation  (of somatic and sex chromosomes in human and mouse, anyway)

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

def call_bw(infile,  outbase,ChromInfo, extraChrString,options=''):
    
    myAwk = "{chr=\"%s\"$1; printf(\"%s\\t%d\\t%d\\t0\\t0\\t%s\\n\",chr,$2,$3,$6) | \" genomeCoverageBed -i stdin -bg -g %s | wigToBigWig stdin %s %s%s.bw\"s;}" % (extraChrString,ChromInfo, ChromInfo,outbase, extraChrString)
    cmd=["bamToBed","-i",infile, "-splitD"]

    print myAwk

    if options != '' and options != None:
        cmd += options.split(' ')

    cmd = [cmd, ['awk', myAwk]]
    return cmd


## Worker function for parallelization
def worker_pipe1(work_queue, done_queue):
    try:
        for cmd in iter(work_queue.get, 'STOP'):
            print cmd            
            p1 = subprocess.Popen(cmd[0], stdout=subprocess.PIPE)
            print "pid: %s" % p1.pid

            done_queue.put("%s - %s got %s." % (current_process().name, cmd, p1.pid))

            p2 = subprocess.Popen(cmd[1], stdin=p1.stdout)
            p1.stdout.close()
            print "pid2: %s" % p2.pid
            
            retval = p2.wait()

            done_queue.put("%s - %s got %s, %s." % (current_process().name, cmd, p2.pid, retval))

    except Exception, e:
        done_queue.put("%s failed on %s with: %s" % (current_process().name, cmd, e.message))
    return True


def main(indir, outdir, ref, extraChrString, option_string, num_processors):


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
    
    # Look through the input directory for .bam files
    bams = glob.glob("%s/*.sorted.bam" % indir)
    
    ## parallel stuff:
    work_queue = Queue()
    done_queue = Queue()
    processes = []

    for bam in bams:
        prefix = re.sub('\.sorted.bam$','',os.path.basename(bam))
        cmd = call_bw(bam, outdir, ref, extraChrString, option_string)
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
    parser.add_option("-e", "--extraChrString", dest="extraChrString", action="store",
                      default='')
    parser.add_option("-H", "--Help", dest="Help", action="store",
                      default=None)
    (options, args) = parser.parse_args()
    if options.Help or options.outdir is None or options.input_dir is None or options.ref is None:
        print __doc__
        sys.exit()
    main(options.input_dir, options.outdir,options.ref, options.extraChrString,options.option_string, options.num_processors)

