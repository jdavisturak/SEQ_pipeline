## Utility to reverse complement a .sam file


import sys
import os
import pysam

def reverse_reads(infile,outfile):
  # do stuff!
  
  # Read in sam data
  samfile = pysam.Samfile(infile, "r" )

  # start creating output file
  outreads = pysam.Samfile(outfile, "wh", template=samfile) # copies the header over too

  #
  for read in samfile.fetch():
     
     read.is_reverse = not read.is_reverse 
     
     if read.is_paired:
        read.mate_is_reverse = not read.mate_is_reverse    
        #read.is_read1 = not read.is_read1
        #read.is_read2 = not read.is_read2
              
     outreads.write(read)
  
  return(0)



def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-i", dest="IN", nargs=1, default=None,
                      help="Input file")
    parser.add_option("-o", dest="OUT", nargs=1, default="RC.sam",
                      help="Output file")
    (options, args) = parser.parse_args()

    if options.IN == None:
        print >> sys.stderr, "Error: need -i"
        sys.exit(1)
        
    inputFile = os.path.abspath(os.path.expanduser(options.IN));
    outputFile = os.path.abspath(os.path.expanduser(options.OUT));
   
    if os.path.exists(outputFile):
      print >> sys.stderr, "Error: output file exists already"
      sys.exit(1)


    reverse_reads(inputFile,outputFile)

if __name__ == '__main__':
    main()
 