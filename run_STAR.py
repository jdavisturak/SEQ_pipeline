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
    if skip == '':
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

    # Summary:
    print "input directory:\t%s\nsam directory:\t%s\nSJs directory:\t%s\nLogs directory:\t%s\n" % (indir, sam_dir, SJs_dir, log_dir)

    # Make and change to a new directory:
    make_or_open_dir('tempOut')
    os.chdir('tempOut')

    

    # Check if it's paired:
    isPaired = True
    if paired is None or paired == False or paired=='no' or paired=='unpaired':
        isPaired=False

    # Get list of files
    if isPaired:
        myGlob = "*_1.fastq"
    else:
        myGlob = "*.fastq"

    fastqIn = glob.glob("%s/%s" % (indir, myGlob))
    
    ## Check for regular expression match 'skip', and then skip those matches
    if not (skip is None or skip==''):
        skip = re.sub("^\^","^%s/" % indir ,skip); # if skip starts with '^' make sure to add the input directory 
        fastqIn = [f for f in fastqIn if not re.match(skip,f)]

    print("%s: processing %d samples" % ("Paired" if isPaired else "Unpaired", len(fastqIn)))


    # Loop through all files and compose a command:
    
    for fq in fastqIn:

        cmd = ['STAR','--genomeLoad', 'LoadAndKeep', '--genomeDir',ref, '--runThreadN', str(num_processors)]  
        if not option_string is None:
            cmd +=  option_string.split(' ')

        cmd += ['--readFilesIn']


        if isPaired:
            prefix = re.sub('_1\.fastq$','', fq)
            cmd += ["%s_1.fastq" % prefix, "%s_2.fastq" % prefix ]
        else:
            prefix = re.sub('\.fastq$','', fq)
            cmd += [fq]
        print ' '.join(cmd)

        try:
            p = subprocess.Popen(cmd)
            print "pid: %s" % p.pid
            retval = p.wait()
            print("%s: got %s, %s." % (cmd, p.pid, retval))

        except Exception, e:
            print("failed on %s with: %s" % (cmd, e.message))
            exit()

        # Once the file is finished processing, move outputfiles to new location:
        prefix_short = os.path.basename(prefix)
        target = "%s/%s.sam" % (sam_dir, prefix_short)     
        try:
            os.rename('Aligned.out.sam',target)
        except:
            print "could not move to %s" % target
            exit()

        os.rename('SJ.out.tab',"%s/%s.SJ.tab" % (SJs_dir, prefix_short) )
        logs = glob.glob("Log*")
        [[os.rename(LOGFILE,"%s/%s.%s" % (log_dir, prefix_short, LOGFILE) )] for LOGFILE in logs]
            
    # now erase temp directory
    try:
        os.chdir("..")
        os.rmdir('tempOut')  
    except Exception, e:
        print "Warning: could not erase temp. directory. No biggie."
    


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-i", "--input_dir", dest="input_dir", action="store",
                      default=".")
    parser.add_option("-o", "--outdir", dest="outdir", action="store",
                      default=None)
    parser.add_option("-O", "--Options", dest="option_string", action="store",
                      default=None)
    parser.add_option("-N", "--num_processors", dest="num_processors", action="store",
                      default=None)
    parser.add_option("-r", "--ref", dest="ref", action="store",
                      default=None)
    parser.add_option("-p", "--paired", dest="paired", action="store_true",
                      default=False)
    parser.add_option("-s", "--skip_regex", dest="skip", action="store",
                      default='ambig*')

    (options, args) = parser.parse_args()
    if options.outdir is None or options.input_dir is None or options.ref is None:
        print __doc__
        sys.exit()
    main(options.input_dir, options.outdir, options.ref, options.paired, options.option_string, options.num_processors, options.skip)












