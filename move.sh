# Remove the annoying EB prefix.
#!/bin/sh 


#for file in *.bam;
#do
#    mv "$file" "${file#EB}"
#done

ALL_SAMPLES=$1

while read line ; do

SAMPLE_ID=$line

if [ ! -d /home/ubuntu/media/data/adipose_samples/"$SAMPLE_ID" ]; then

	mkdir /home/ubuntu/media/data/adipose_samples/"$SAMPLE_ID"
	mkdir /home/ubuntu/media/data/adipose_samples/"$SAMPLE_ID"/reference

	mv /home/ubuntu/media/data/FAT/"$SAMPLE_ID"_sorted.bam /home/ubuntu/media/data/adipose_samples/"$SAMPLE_ID"/"$SAMPLE_ID"_sorted.bam

else

mv /home/ubuntu/media/data/FAT/"$SAMPLE_ID"_sorted.bam /home/ubuntu/media/data/adipose_samples/"$SAMPLE_ID"/"$SAMPLE_ID"_sorted.bam

fi

done < ${ALL_SAMPLES}
