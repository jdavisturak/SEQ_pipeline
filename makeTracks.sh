#!/bin/sh

inputDir=$1
outputDir=$2
cd $inputDir

scp *bw? signal@signalingsystems.ucsd.edu:$outputDir

for bw in *bw?; do

	Name=$(basename $bw .bw?)
	Name=${Name/plus/"+"}
	Name=${Name/minus/"-"}
	Location=${outputDir/\/Library\/WebServer\/Documents\/signal/http://signalingsystems.ucsd.edu}

	echo "track name=$Name description=$Name type=bigWig bigDataUrl=$Location/$bw visibility=full maxHeightPixels=80:30:15"

done


