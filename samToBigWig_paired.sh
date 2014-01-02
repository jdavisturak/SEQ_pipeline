#!/bin/bash
#
# Jeremy Davis-Turak
#
# Convert .sam to bigWig
# Input: sam file: X.sam
# Output: X.bw+, X.bw-
#
# Restrictions: .sam is paired end, with all pairs on adjacent lines

Sam=$1
ChromInfo=$2
extraChrString=$3

## Note to self: check for '1' or '2' at end of the .tag name!
F=$(basename $Sam .sam)
cmd="samtools view -b -S $F.sam | bamToBed -splitD  | awk  '{s=\$6; str=substr(\$4,length(\$4),1); if(str==2){if(s==\"+\"_) s=\"-\"; else s=\"+\"} chr=\"$extraChrString\"\$1; printf(\"%s\t%d\t%d\t0\t0\t%s\n\",chr,\$2,\$3,s) | \"sort -k1,1 -k2,2n | genomeCoverageBed -i stdin -bg -g $ChromInfo | wigToBigWig stdin $ChromInfo $F$extraChrString.bw\"s;}'"
echo $cmd
eval $cmd   
