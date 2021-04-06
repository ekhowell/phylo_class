for file in *.fa
do
	cat $file | sed -e '1!{/^>.*/d;}' | grep -v ">" | tr -d "\n\r" | fold -w 50 -s > temp
	echo ">$file" >> combined_samples.fa
   	tail -n +2 temp >> combined_samples.fa
   	echo >> combined_samples.fa
done
