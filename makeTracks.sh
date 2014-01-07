#!/bin/sh

inputDir=$1
outputDir=$2
CWD=`pwd`
cd $inputDir
fileDest=$CWD/$3


if [ "$#" -lt 3 ]; then
    echo 'no third argument'
    fileDest='trackDb.txt'
fi

# echo "$fileDest"

#scp *bw? signal@signalingsystems.ucsd.edu:$outputDir



touch $fileDest

for bw in *bw?; do

    Name=$(basename $bw .bw?)
    Name=${Name/plus/"+"}
    Name=${Name/minus/"-"}
    Location=${outputDir/\/Library\/WebServer\/Documents\/signal/http://signalingsystems.ucsd.edu}

    #echo "track name=$Name description=$Name type=bigWig bigDataUrl=$Location/$bw visibility=full maxHeightPixels=80:30:15"
    echo "track $Name" >> $fileDest
    echo "bigDataUrl $Location/$bw" >> $fileDest
    echo "shortLabel $Name" >> $fileDest
    echo "longLabel $Name" >> $fileDest
    echo "type bigWig" >> $fileDest
    echo " 
    " >> $fileDest      

done