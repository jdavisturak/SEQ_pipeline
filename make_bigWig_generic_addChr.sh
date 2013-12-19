# Argument 1: Hits_prefix='accepted_hits'
# Argument 2: ChromInfo='/home/jeremy/RNAseq/indices/Mus_musculus/UCSC/mm10/Sequence/Chromosomes/chromInfo.txt'

Hits_prefix=$1
ChromInfo=$2
paired=$3

echo "ChromInfo: "$ChromInfo

# First sort and convert to .bed (and get rid of the name of the read)
                                                                                                                               
if [ ! -e "${Hits_prefix}.sorted.bam" ]
  then
    samtools sort $Hits_prefix.bam  $Hits_prefix.sorted
  fi

if [ ! -e "${Hits_prefix}.bed2+" ]
  then  
    bamToBed -split -i  $Hits_prefix.sorted.bam | awk  '{printf("chr%s\t%d\t%d\t0\t0\t%s\n",$1,$2,$3,$6) >> "'${Hits_prefix}'.bed2"$6;}'  
  fi
  
# Convert to bedGraph then .wig 
bedtools genomecov -i $Hits_prefix.bed2+ -bg -g $ChromInfo | wigToBigWig stdin $ChromInfo $Hits_prefix.plus2.bw

bedtools genomecov -i $Hits_prefix.bed2- -bg -g $ChromInfo | wigToBigWig stdin $ChromInfo $Hits_prefix.minus2.bw
#####
