# Script to take a list of codes from Suhua (UCLA seq core) and make a password file and commands for downloading:
#!bin/bash

outputDir=$1 # this will be the prefix to "."

while read line
do
	#echo $line
	arr=(${line//:/ }) # each line should look like SxaQSEQsWA080L5:WYawMlDt with the semicolon separating the Dataset and Password
    D=${arr[0]}
    P=${arr[1]}
	#echo $P
	# Create a password file:
	echo $P > pwd_$D
	# set permissions:
	chmod 660 pwd_$D

	# create command of interest:
	cmd="rsync --password-file=pwd_$D --recursive --times --verbose --stats --progress --itemize-changes rsync://$D@pan.pellegrini.mcdb.ucla.edu/$D/  $outputDir."
	echo $cmd

done

