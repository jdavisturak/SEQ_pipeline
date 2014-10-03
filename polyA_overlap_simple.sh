#!/bin/bash
# This is a script to combine lots of existing tools to compute the overlap of a .sam or .bam file 

# Inputs: 1) a  file path 
#         2) path to polyA bed file
#         3) Output prefix: appended to file path

f=$1
A_bed=$2
prefix=$3;

F=$(basename $f .bw+)

echo "file :$F"
echo "Poly(A) file :$A_bed"
echo "prefix :$prefix"



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
