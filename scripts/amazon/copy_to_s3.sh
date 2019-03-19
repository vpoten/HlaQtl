#!/bin/bash

#usage: copy_to_s3.sh <bucket> <directory>

files=$2/*

for f in $files
do
if [ -f "$f" ]; then
	filename=$(basename $f)
	echo "copying $filename file to $1..."
	echo aws put $1$filename $f
fi
done
