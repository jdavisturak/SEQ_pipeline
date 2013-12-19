from __future__ import division
from itertools import chain
import numpy as np
import csv
import warnings

# Hamming distance between two strings
def hamming_distance(s1, s2):
  assert len(s1) == len(s2)
  return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

# Compute all pairwise hamming distances between all strings in a list, and return the minimum distance
def min_hamming(list1):
  min_dist = len(list1[0])
  for i in range(0,len(list1)-1):
    for j in range(i+1,len(list1)):
      min_dist = min(hamming_distance(list1[i],list1[j]),min_dist)
  return (min_dist)


## Generate all 6-mers: recursively!
def generate_Kmers(num=6, alphabet=['A','C','G','T','N']):
  if (num <= 1):
    return(alphabet)
  # Otherwise:
  # Get all k-1 mers
  Km1mers = generate_Kmers(num-1,alphabet)
  
  # Loop through all characters
  # apppend EACH to all (k-1)-mers of   
  
  a=[[x+alphabet[i] for x in Km1mers] for i in range(0,len(alphabet))]
  return list(chain(*a))

    
### Read in targets file: looks for column names 'Sample' and 'Barcode'
def read_targets(filename='targets.txt',samples_name='Sample', barcodes_name='Barcode', lanes_name='Lane', barcodes2_name = 'Barcode2'):
  
  infile=open(filename,'r')
  r = csv.DictReader(infile,dialect=csv.Sniffer().sniff(infile.read(1000)))
  infile.seek(0)
  
  samples, barcodes, lanes, barcodes2 = [],[],[],[]
  for row in r:
    samples.append(row[samples_name])
    barcodes.append(row[barcodes_name])
    lanes.append(row[lanes_name])
    if barcodes2_name in row.keys():
      barcodes2.append(row[barcodes2_name])     
  
  assert all(len(x) == len(barcodes[0]) for x in barcodes)
  assert all(len(x) == len(barcodes2[0]) for x in barcodes2)
  
  if barcodes2:
    targets = np.array([samples,lanes,barcodes,barcodes2]).transpose()
  else:
    targets = np.array([samples,lanes,barcodes]).transpose()
  return(targets)


def find_barcode(kmer, barcodes, max_distance, Map):
  # Distance from this kmer to all barcodes:
  dists = [hamming_distance(kmer, code) for code in barcodes]  
  min_dist = min(dists)
  
  # Find those that equal the minimum distance:
  which_min = [x==min_dist for x in dists]  
  myMins = [int(x) for x in list(chain(*np.where(which_min)))] # this is like the 'where' function in R
  
  if((len(myMins)) == 1 & (min_dist <= max_distance)):
    Map[kmer] = myMins[0]
  else:
    Map[kmer] = []

  # The Map object is passed back by reference to the scope that declared it


 
    
### Create dictionary:
def create_barcode_dictionary(barcodes, max_distance=2):

  print barcodes
 
  if not barcodes:
      return None
  print "calculating Hamming distances"

  # Generate kmers                                    
  num = len(barcodes[0]);
  kmers = generate_Kmers(num=num)
  
  # Get minumum barcode distance
  min_dist = min_hamming(barcodes)
  
  # If the user specified a max_distance that is greater than min_dist, then give a warning and reset max_distance
  if (min_dist < max_distance):
    warnings.warn("max_distance (%s) is bigger than minimum Hamming distance(%s).  Resetting the max_distance" % (max_distance,min_dist))
    max_distance = min_dist
  
  Map = {}
  out=[find_barcode(K,barcodes, max_distance,Map) for K in kmers]
  print "%d %d-mers match" % (len(np.where([type(Map[key]) == int  for key in Map.keys()])[0]), num)
  
  return Map
  

### Create a dict of dicts so that we can decipher the samples...
def samples_map(barcodes1,barcodes2,samples):
  
  myDict = {}
  b1 = set(barcodes1)
  b2 = set(barcodes2)
  
  # Set up the initial list:
  for one in b1:
      myDict[one] = {}
      for two in b2:
        myDict[one][two]=[]
        
  for s in range(0,len(samples)):
    myDict[barcodes1[s]][barcodes2[s]] = samples[s]
      
  return myDict
 


### Function to decide what type of map to return:
def create_barcode_map(barcodes1, barcodes2, samples, max_distance=2):
  if barcodes2 == []:
    return create_barcode_dictionary(barcodes1,max_distance)
  
  ## if we have two barcodes:  
  # Get unique barcodes for each list; ordered 
  i1 = list(set(barcodes1))
  i2 = list(set(barcodes2))
  i1.sort()
  i2.sort()  
  
  # Now create dictionaries
  dict1 = create_barcode_dictionary(i1,max_distance)
  dict2 = create_barcode_dictionary(i2,max_distance)    
 
  # so these are hash tables into the closest INDEX (into unique list) to the key of the dict  
  ## For each list of barcodes, figure out where they are in their respective lists
 
  index1 = []
  for b in barcodes1:
    for a in range(0,len(i1)):
      if i1[a] == b:
        index1.append(a)
        continue
  
 
  index2 = []
  for b in barcodes2:
    for a in range(0,len(i2)):
      if i2[a] == b:
        index2.append(a)
        continue
  
  # Return a list containing the two dictionaries and the samples list
  return [samples_map(index1, index2, range(0,len(samples))), dict1, dict2]
  
 
      

## Function to return the INDEX into a samples given two barcodes:
def get_sample_index(bc1, bc2, bc_map):
  if bc2 == []:
    return bc_map[bc1]
  
  # Otherwise:
  (samples_map,dict1, dict2) = bc_map
 
  # My barcodes indices:
  i1 = dict1[bc1]
  i2 = dict2[bc2] 
   
  if i1 == [] or i2 == []:
    return []
  return samples_map[i1][i2]
   
  

     
  
## Testing: check out how many matches there are to various error tolerances in the first qseq file
#infile=open('s_2_2_1101_qseq.txt','r')
#reader=csv.reader(infile,delimiter='\t')
#codes = []
#for row in reader:
#  codes.append(row[8])

# 6-mers
#codes2 = [code[0:6] for code in codes]

# Read targets
#(samples, barcodes) = read_targets(filename='Pilot_targets.txt')

# For various max_distance values, check how many 6-mers in my sample match one of the allowed 6-mers


#Map=create_barcode_dictionary(barcodes,0)
#len(np.where([(type(Map[code])==int) for code in codes2])[0])/len(codes2)

#Map=create_barcode_dictionary(barcodes,1)
#len(np.where([(type(Map[code])==int) for code in codes2])[0])/len(codes2)

#Map=create_barcode_dictionary(barcodes,2)
#len(np.where([(type(Map[code])==int) for code in codes2])[0])/len(codes2)




  
