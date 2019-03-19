#!/bin/bash

#launchs eqtl calculation

# regions to analyze
regions=("chr4:103052643-104182151")
# groups to analyze
groups=("CEU" "GBR" "FIN" "TSI" "CEU,GBR,FIN,TSI")

# path to HlaQtl jar
eqtl_jar=/mnt/ngsengine/HlaQtl-1.0-SNAPSHOT-bin/dist/HlaQtl-1.0-SNAPSHOT.jar
# configuration file (.properties)
properties=/home/ubuntu/ngsengine.properties
# output directory
out_dir=/mnt/ngsengine/output/
# cufflinks data to use: cufflinks_ensOnly_ or cufflinks_refOnly_
cuff_pre=cufflinks_ensOnly_ 

for i in "${regions[@]}"
do
	for j in "${groups[@]}"
	do
		echo "Performing eqtl for locus" $i "and group" $j
		echo mkdir $out_dir
		echo java -Xmx5g -jar $eqtl_jar eqtl --props=$properties --output=$out_dir --locus=$i --group=$j --cuffPre=$cuff_pre > $i"_"$j.out
		echo rm -rf $out_dir
	done
done
