#!/bin/bash
# Purpose: For each FASTA file in the directory, this script will merge them into a single FASTA file

for file in *.fa
do
	# Merge FASTA records *within* a file
	cat $file | sed -e '1!{/^>.*/d;}' | grep -v ">" | tr -d "\n\r" | fold -w 50 -s > temp
	# Rename header to match file name (i.e. sample name)
	echo ">$file" >> combined_samples.fa
	# Add it to the combined samples FASTA
   	tail -n +2 temp >> combined_samples.fa
   	echo >> combined_samples.fa
done
