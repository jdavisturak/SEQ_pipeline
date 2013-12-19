IN=$1
OUT=$2 

awk '{if(NR%4==1){a=$0} if(NR%4==2) {b=$0} if(NR%4==3) {c=$0}   if(NR%4==0 && b ~ /^TTTTT/) {printf("%s\n%s\n%s\n%s\n",a,b,c,$0 )  >> "'${OUT}'T.fastq";} if(NR%4==0 && b ~ /AAAAA$/) {printf("%s\n%s\n%s\n%s\n",a,b,c,$0 )  >> "'${OUT}'A.fastq";} }' $IN



