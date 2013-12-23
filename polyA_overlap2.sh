#!/bin/bash
# This is a script to combine lots of existing tools to compute the overlap of a .sam or .bam file 

# Inputs: 1) a partial file path 
#         2) path to polyA bed file
#         3) Output prefix: appended to file path
#         4) Boolean: whether to first convert FROM .sam
#         5) Boolean: whether to convert to .bw format
#         6) Path to ChromInfo.txt, which has the lenghts of the chromosomes, 
#         7) Boolean: whether to ALSO convert to .bw format using the 'chr' UCSC convention  
#         8) 'paired' 

F=$1
A_bed=$2
prefix=$3;
Is_sam=$4 # convert from .sam format first
convert=$5 # convert to bigWig format first
ChromInfo=$6
convert2=$7 # Also convert to bigWig format where I add 'chr' to all the chromosomes.
paired=$8     # 'paired'

echo "file :$F"
echo "Poly(A) file :$A_bed"
echo "prefix :$prefix"
echo "Sam? :$Is_sam"
echo "Convert to bigwig? :$convert"
echo "Convert to bigwig2? :$convert2"


# Convert to sam.   Doesn't execute if the .sorted.bam already exists.
if [[ ! -e "${Hits_prefix}.sorted.bam" && "$Is_sam" -eq "1" ]]
  then 
    echo 'Converting to .bam'    #If also converting to bigWig, don't bother outputting the .bam b/c we are gonna need to sort it anyway.
    if [ "$convert" -eq "1" ]
      then 
        samtools view -b -S ${F}.sam | samtools sort - $F.sorted
    else  
      samtools view -b -S ${F}.sam -o ${F}.bam; 
    fi
fi


# Convert .bam to .bw
if [ "$convert" -eq "1" ]
  then
    echo 'coverting to bigWig format using '$ChromInfo
    sh /home/home/jeremy/Code/make_bigWig_generic.sh $F $ChromInfo $paired
fi 

if [ "$convert2" -eq "1" ]
  then
    echo 'coverting to bigWig format with "chr" using '$ChromInfo.UCSC
    sh /home/home/jeremy/Code/make_bigWig_generic_addChr.sh $F $ChromInfo.UCSC $paired
fi 


# Do poly-A overlap

  rm ${F}.${prefix}.A.plus.Down
  rm ${F}.${prefix}.A.plus.Up
  rm ${F}.${prefix}.A.minus.Down
  rm ${F}.${prefix}.A.minus.Up
      
	bigWigAverageOverBed ${F}.bw+ $A_bed ${F}.plus.${prefix}.A
  awk '$1 ~ /Up/ {a="Up"} $1 ~ /Down/ {a="Down"} {printf($0 "\n") >> "'${F}'.'${prefix}'.A.plus."a}' ${F}.plus.${prefix}.A
  rm ${F}.plus.${prefix}.A

  bigWigAverageOverBed ${F}.bw- $A_bed ${F}.minus.${prefix}.A 
  awk '$1 ~ /Up/ {a="Up"} $1 ~ /Down/ {a="Down"} {printf($0 "\n") >> "'${F}'.'${prefix}'.A.minus."a}' ${F}.minus.${prefix}.A
  rm ${F}.minus.${prefix}.A
