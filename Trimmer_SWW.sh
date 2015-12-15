#!/bin/bash

input_path_R1="/home/shahrokh/bio720/FinalProject/data/*f.fastq.gz"
out_path="/home/shahrokh/bio720/FinalProject/results"

for x in $input_path_R1
do  
	java -jar /usr/local/trimmomatic/trimmomatic-0.33.jar PE -threads 20 -basein $x -baseout $out_path/`basename $x`\
	ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done
