# Remove the annoying EB prefix.
#!/bin/sh 


#for file in *.bam;
#do
#    mv "$file" "${file#EB}"
#done

ALL_SAMPLES=$1

while read line ; do

SAMPLE_ID=$line

if [ ! -d /media/data/blood_samples/"$SAMPLE_ID" ]; then

	mkdir /media/data/blood_samples/"$SAMPLE_ID"
	mkdir /media/data/blood_samples/"$SAMPLE_ID"/reference

	mv /media/data/BLOOD/"$SAMPLE_ID"_sorted.bam /media/data/blood_samples/"$SAMPLE_ID"/"$SAMPLE_ID"_sorted.bam

else

mv /media/data/BLOOD/"$SAMPLE_ID"_sorted.bam /media/data/blood_samples/"$SAMPLE_ID"/"$SAMPLE_ID"_sorted.bam

fi

done < ${ALL_SAMPLES}
