#!/usr/bin/env python
''' run_STAR
Usage:
python run_STAR.py -i <indir> -o <outdir> -r <ref> [-N num_processors, --paired=[(yes)\|no], --skip_regex='^ambig']

'''

from __future__ import with_statement
import os
import sys
import glob
import re
from optparse import OptionParser
import subprocess

# cd /home/RNAseq/Kim/RelA_mutants/fastq
# mkdir ../aligned
# mkdir ../aligned/sams
# mkdir ../aligned/SJs
# mkdir ../aligned/logs

def make_or_open_dir(Dir):
    if not os.path.exists(Dir):
        try:
            os.makedirs(Dir)
        except OSError:
            assert os.path.isdir(Dir)


def main(indir, outdir, ref, paired, option_string, num_processors, skip):

    # Check inputs:
    try:
        num_processors = int(num_processors)
    except ValueError, TypeError:
        num_processors = 1
    if num_processors < 1:
        num_processors = 1

    if option_string is None:
        option_string = '' 
    if skip == ''
        skip = None

    make_or_open_dir(outdir)
    sam_dir = os.path.abspath(outdir+'/sams')
    SJs_dir = os.path.abspath(outdir+'/SJs')
    log_dir = os.path.abspath(outdir+'/logs')

    make_or_open_dir(sam_dir)
    make_or_open_dir(SJs_dir)
    make_or_open_dir(log_dir)
    
    indir = os.path.abspath(indir)

    assert os.path.isdir(indir)

    # Make and change to a new directory:
    make_or_open_dir(indir + 'tempOut')
    try:
        os.chdir(indir + 'tempOut')
    

    # Check if it's paired:
    isPaired=True
    if paired is None or paired == False or paired=='no' paired=='unpaired':
        isPaired=False

    # Get list of files
    if isPaired:
        myGlob = "*_1.fastq"
    else:
        myGlob = "*.fastq"

    fastqIn = glob.glob("%s/%s" % (indir, myGlob))

    # Loop through all files and compose a command:
    cmd = ['STAR','--genomeDir',ref, 'runThreadN', num_processors] + ' '.split(option_string) + ['--readFilesIn']

    for fq in fastqIn:

        if  continue

        if isPaired:
            prefix = re.sub('_1\.fastq$','', fq)
            cmd += ["%s_1.fastq" % prefix, "%s_2.fastq" % prefix ]
        else:
            prefix = re.sub('\.fastq$','', fq)
            cmd += [fq]
        print cmd

        try:
            p = subprocess.Popen(cmd)
            print "pid: %s" % p.pid
            retval = p.wait()
            done_queue.put("%s - %s got %s, %s." % (current_process().name, cmd, p.pid, retval))

        except Exception, e:
            done_queue.put("%s failed on %s with: %s" % (current_process().name, cmd, e.message))

        # Once the file is finished processing, move outputfiles to new location:
        os.rename('Aligned.out',"%s/%s.sam" % (sam_dir, prefix) )
        os.rename('SJ.out.tab',"%s/%s.SJ.tab" % (SJs_dir, prefix) )
        logs = glob.glob("*log")
        [[os.rename(LOGFILE,"%s/%s.%s" % (log_dir, prefix, LOGFILE) )] for LOGFILE in logs]
            
    # now erase temp directory
    os.chdir("..")
    os.rmdir('tempOut')        



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
    parser.add_option("-p", "--paired", dest="paired", action="store",
                      default=False)
    parser.add_option("-s", "--skip_regex, dest="skip", action="store",
                      default='ambig*')

    (options, args) = parser.parse_args()
    if options.Help or options.outdir is None or options.input_dir or options.ref is None:
        print __doc__
        sys.exit()
    main(options.input_dir, options.outdir,options.ref, options.paired,options.option_string, options.num_processors,options.skip)












