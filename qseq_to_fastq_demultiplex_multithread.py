#!/usr/bin/env python
"""Convert output solexa qseq files into fastq format, handling multiplexing.

Works with qseq output from Illumina's on-machine base caller

Usage:
    qseq_to_fastq_demultiplex_multithread.py <run name>  <targets file name> -o <outdir> [-N num_processors][-r ] [-R] [-f]

Output files will be in the <outdir> directory as <run_name><sample>_l<lane>[_1/2].fastq

Illumina barcoded samples contain barcodes in a separate qseq lane, which are
identified by being much shorter than the primary read. 

Optional arguments:
    --failed (-f): Write out reads failing the Illumina quality checks also, to outdir/failed.
    --reverse (-r): Reads map to anti-sense strand of RNA.  This flag has different consequences for paired or unpaired reads:
            Unpaired (single-end) reads get reverse-complemented, whereas Paired reads swap their identities 
    
"""

from __future__ import with_statement
import os
import sys
import glob
from optparse import OptionParser
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from demultiplex import *
import multiprocessing as mp

def main(run_name, target_name, do_fail=False, inputDir='.', outdir=None, reverse=False,RevBarcodes=False, num_processors=1):
    
    try:
        num_processors = int(options.num_processors)
    except ValueError, TypeError:
        num_processors = 1
    if num_processors < 1:
        num_processors = 1

    print target_name
    assert os.path.exists(target_name)

    if outdir is None:
        outdir = os.path.join(os.getcwd(), "fastq")
    if not os.path.exists(outdir):
        try:
            os.makedirs(outdir)
        except OSError:
            assert os.path.isdir(outdir)
    if do_fail:
        fail_dir = os.path.join(outdir, "failed")
        if not os.path.exists(fail_dir):
            try:
                os.makedirs(fail_dir)
            except OSError:
                assert os.path.isdir(fail_dir)
    else:
        fail_dir = None
    targets = read_targets(filename=target_name)

    barcodes2 = []
 
    for lane_num in list(set(targets[:,1])):
        print "Demultiplexing lane %s" % lane_num
        lane_prefix = os.path.join(inputDir,"s_%s" % lane_num)
        out_prefix = run_name
        samples = [target[0] for target in targets[targets[:,1]==lane_num]]
        barcodes = [target[2] for target in targets[targets[:,1]==lane_num]]
        
        if RevBarcodes:
            barcodes = [str(Seq(b, generic_dna).reverse_complement()) for b in barcodes]

        if (targets.shape[1] > 3):
            barcodes2 = [target[3] for target in targets[targets[:,1]==lane_num]]
            if RevBarcodes:
                barcodes2 = [str(Seq(b, generic_dna).reverse_complement()) for b in barcodes2]

        write_lane(lane_prefix, out_prefix, outdir, fail_dir,target_name,samples,barcodes,barcodes2, lane_num, reverse, num_processors)

def write_lane(lane_prefix, out_prefix, outdir, fail_dir,target_name,samples,barcodes,barcodes2, lane_num, reverse=False, num_processors=1):

    ## Determine barcodes, etc
    bc_len = len(barcodes[0])
    bc2_len = len(barcodes2[0]) if barcodes2 else None
    qseq_files = glob.glob("%s_*qseq.txt" % lane_prefix)
    one_files, two_files, bc_files, bc2_files = _split_paired(qseq_files)
    is_paired = len(two_files) > 0
    is_double_bc = len(bc2_files) > 0
    print "Detecting %s end reads and %d barcode%s for lane %s\n" % \
        ('paired' if is_paired else 'single', 2 if is_double_bc else 1,  's' if is_double_bc else '',lane_num)

    # Get paths to all files
    out_files = (_get_outfiles(out_prefix, outdir, is_paired, samples,lane_num)
                 if not fail_dir else None)
    fail_files = (_get_outfiles(out_prefix, fail_dir, is_paired, samples, lane_num)
                  if fail_dir else None)
    bc_map = create_barcode_map(barcodes,barcodes2,samples)
    

    for (num, files) in [("1", one_files), ("2", two_files)]:

        if num=='2' and not is_paired:
            break
 
        # Set up multi-threaded processing

        worker_queue = mp.Queue()
        writer_queue = [mp.Queue() for s in xrange(len(out_files[num])) ]
        
        Worker_proc = [mp.Process(target=Work, args= (worker_queue, writer_queue) ) for n in xrange(num_processors)]
        if fail_files:
            Writer_proc = [mp.Process(target=Write, args=(writer_queue[s], fail_files[num][s]) ) for s in xrange(len(fail_files[num])) ]
        else:
            Writer_proc = [mp.Process(target=Write, args=(writer_queue[s], out_files[num][s]) ) for s in xrange(len(out_files[num]))]

        for i, fname in enumerate(files):
            bc_file = _get_associated_barcode(i, fname, bc_files)
            bc2_file = _get_associated_barcode(i ,fname, bc2_files)
            paramsList = [fname, num, bc_file, bc2_file, out_files[num], bc_map,  bc_len, bc2_len,  is_paired,fail_files, reverse]
            
            # Send the file name and other info to the worker queue
            worker_queue.put(paramsList)

        ## Start up the parallel processes
        for p in Worker_proc:
            p.start()
       
        for p in Writer_proc:
            p.start()
        
        for p in Worker_proc:
            p.join()
        
        for p in Writer_proc:
            p.join()
        
        print ("Ok now time to quit??")
        # I'm not sure when to add the quitting signal:
        for p in worker_queue:
            p.put('kill')
        for p in writer_queue:
            p.put('kill')
    
        #convert_qseq_to_fastq(fname, num, bc_file, bc2_file, out_files, bc_map, ambig_seqs, ambig_seqs2, ambig_bcs, bc_len, bc2_len,  is_paired,fail_files, reverse)


def Work(worker_queue, writer_queue):
    while True:
        params = worker_queue.get() 
        if (params=='kill'):
            print 'quitting processing after  file %s' % params[0]
            break

        print "Processing input file ", params[0]
        res = convert_qseq_to_fastq(params, writer_queue)
            

def Write(write_queue, fhandle):
    while True:
    
        res = write_queue.get()
        
        if(res == 'kill'):
            print "Closing file handle %s" % fhandle
            break
        fhandle.write(res)

    fhandle.close()
    



def _get_associated_barcode(file_num, fname, bc_files):
    """Get barcodes for the first read if present.
    """
    if len(bc_files) > 0:
        bc_file = bc_files[file_num]
        bc_parts = bc_file.split("_")
        read_parts = fname.split("_")
        assert (bc_parts[1] == read_parts[1] and
                bc_parts[3] == read_parts[3]), (bc_parts, read_parts)
        return bc_file
    return None

def convert_qseq_to_fastq(params, writer_queue):

    """Convert a qseq file into the appropriate fastq output.
    """

    (fname, num, bc_file, bc2_file, out_files, bc_map, bc_len, bc2_len, is_paired, fail_files, reverse) \
        = params;
    if fail_files:
        fail_files = fail_files[num]

    bc_iterator = _qseq_iterator(bc_file, fail_files is None) if bc_file else None
    bc2_iterator = _qseq_iterator(bc2_file, fail_files is None) if bc2_file else None

    # If paired-end and reverse==True, swap the identity of the read (1 v 2) 
    # "num" is used both to add the flag at the end of the read NAME and to determine the target file name
    if reverse and is_paired:
        num = str(3-int(num))
    
    ambig_queue     = writer_queue[ len(out_files) - 2]
    ambig_bcs_queue = writer_queue[ len(out_files) - 1]

    if fail_files:
        ambig_seqs = fail_files[len(out_files) - 2]
        ambig_bcs  = fail_files[len(out_files) - 2]
    else:
        ambig_seqs = out_files[len(out_files) - 2]
        ambig_bcs  = out_files[len(out_files) - 2]
      

    # if fail_files:
    #     fambig_queue     = fwriter_queue[ len(fail_files) - 2]
    #     fambig_bcs_queue = fwriter_queue[ len(fail_files)  1]
    count = 0

    for basename, seq, qual, passed in _qseq_iterator(fname, fail_files is None):
        count+=1;
        #On read %d for file $%s" % (count, fname)


        # If single-end and reverse==True, rev comp the read
        if reverse and not is_paired:
            # reverse quality string
            qual = qual[::-1] 
            # reverse COMPLEMENT DNA string
            seq = str(Seq(seq, generic_dna).reverse_complement())


        if bc_iterator:
          (_, bc_seq, _, _) = bc_iterator.next()

        bc2_seq = []
        if bc2_iterator:
            (_, bc2_seq, _, _) = bc2_iterator.next()
            
        BC_SEQ =  bc_seq[0:bc_len]
        BC2_SEQ = bc2_seq[0:bc2_len] if bc2_file else []
            
        full_seq = "%s\t%s" % (BC_SEQ,BC2_SEQ)
        name = "%s/%s" % (basename, num)
        out = "@%s\n%s\n+\n%s\n" % (name, seq, qual)

        # This figures out which file to write to
        my_index = get_sample_index(BC_SEQ, BC2_SEQ, bc_map)       

        if passed and not fail_files:
            if not my_index ==[]:
                writer_queue[my_index].put(out)
            else:
                ambig_queue.put(out)
                if not ambig_bcs is None:
                    ambig_bcs_queue.put(full_seq+"\n")
        
        elif fail_files and not passed:                
            if not my_index == []:
                writer_queue[my_index].put(out)
            else:
                ambig_queue.put(out)
                if not ambig_bcs is None:
                    ambig_bcs_queue.put(full_seq+"\n") 
    print "Finished processing file %s (%d reads)" % (fname, count)  


def _qseq_iterator(fname, pass_wanted):
    """Return the name, sequence, quality, and pass info of qseq reads.

    Names look like:

    HWI-EAS264:4:1:1111:3114#0/1
    """
    with open(fname) as qseq_handle:
        for line in qseq_handle:
            parts = line.strip().split("\t")
            passed = int(parts[-1]) == 1
            if passed is pass_wanted:
                name = ":".join([parts[0]] +  parts[2:6]) + "#" + parts[6]
                seq = parts[8].replace(".", "N")
                qual = parts[9]
                assert len(seq) == len(qual)
                yield name, seq, qual, passed

def _get_outfiles(out_prefix, outdir, has_paired_files, samples,lane_num):

    out_files = {}
    if has_paired_files:
        for num in ("1", "2"):
            out_files[num] = []
            for sample in samples:
                out_files[num].append(os.path.join(outdir, "%s%s_l%s_%s.fastq" % (
                    out_prefix, sample, lane_num,num)))
            #add ambigugous files:
            out_files[num].append(os.path.join(outdir, "%sambig_l%s_%s.fastq" % (
                    out_prefix, lane_num,num)))
        out_files['1'].append(os.path.join(outdir, "%sambig_bcs_l%s.fastq" % (
                    out_prefix, lane_num)))
        out_files['2'].append('')


    else:
        out_files["1"] = []
        for sample in samples:
            out_files["1"].append(os.path.join(outdir, "%s%s_l%s.fastq" % (
                out_prefix, sample,lane_num)))
        out_files["1"].append(os.path.join(outdir, "%sambig_l%s.fastq" % (
                out_prefix,lane_num)))
        out_files["1"].append(os.path.join(outdir, "%sambig_bcs_l%s.fastq" % (
                out_prefix,lane_num)))
                
    for index, flist in out_files.items():
        for fname in out_files[index]:
          if os.path.isfile(fname):
            raise ValueError("File exists: %s" % fname) 
        out_files[index] = [open(fname, "w") if not fname=='' else None for fname in out_files[index]]
    return out_files

 

def _split_paired(files):
    """Identify first read, second read and barcode sequences in qseqs.

    Barcoded sequences are identified by being much shorter than reads
    in the first lane.
    """
    files.sort()
    one = []
    two = []
    bcs = []
    bcs2 = []
    ref_size = None
    
    for f in files:
        parts = f.split("_")
        
        if parts[2] == "1":
            one.append(f)
            if ref_size is None:
                ref_size = _get_qseq_seq_size(f) // 2
    
        elif parts[2] == "2":
            cur_size = _get_qseq_seq_size(f)
            assert ref_size is not None
            if cur_size < ref_size:
                bcs.append(f)
            else:
                two.append(f)
    
        elif parts[2] == "3":
            cur_size = _get_qseq_seq_size(f)
            assert ref_size is not None
            if cur_size < ref_size:
                bcs2.append(f)
            else:
                two.append(f)
        elif parts[2] == "4":
            two.append(f)
        else:
            raise ValueError("Unexpected part: %s" % f)
    one.sort()
    two.sort()
    bcs.sort()
    bcs2.sort()
    if len(two) > 0: assert len(two) == len(one)
    if len(bcs) > 0: assert len(bcs) == len(one)
    if len(bcs2) > 0: assert len(bcs2) == len(one)
    return one, two, bcs, bcs2

def _get_qseq_seq_size(fname):
    with open(fname) as in_handle:
        parts = in_handle.readline().split("\t")
        return len(parts[8])

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-f", "--failed", dest="do_fail", action="store_true",
                      default=False)
    parser.add_option("-o", "--outdir", dest="outdir", action="store",
                      default=None)
    parser.add_option("-r", "--reverse", dest="reverse", action="store_true",
                      default=False)
    parser.add_option("-R", "--RevBarcodes", dest="RevBarcodes", action="store_true",
                      default=False)
    parser.add_option("-N", "--num_processors", dest="num_processors", action="store",
                      default=1)
    parser.add_option("-i", "--inputDir", dest="inputDir", action="store",
                      default=".")
    (options, args) = parser.parse_args()
    if len(args) < 2:
        print __doc__
        sys.exit()
    main(args[0], args[1], options.do_fail, options.inputDir, options.outdir, options.reverse, options.RevBarcodes, options.num_processors)

