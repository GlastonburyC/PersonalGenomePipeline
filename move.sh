# Remove the annoying EB prefix.
#!/bin/sh 


#for file in *.bam;
#do
#    mv "$file" "${file#EB}"
#done

ALL_SAMPLES=$1

while read line ; do

SAMPLE_ID=$line

if [ ! -d /media/shared_data/data/blood_samples/"$SAMPLE_ID" ]; then

	mkdir /media/shared_data/data/blood_samples/"$SAMPLE_ID"
	mkdir /media/shared_data/data/blood_samples/"$SAMPLE_ID"/reference

	mv /media/shared_data/data/BLOOD/"$SAMPLE_ID"_sorted.bam /media/shared_data/data/blood_samples/"$SAMPLE_ID"/"$SAMPLE_ID"_sorted.bam

else

mv /media/shared_data/data/BLOOD/"$SAMPLE_ID"_sorted.bam /media/shared_data/data/blood_samples/"$SAMPLE_ID"/"$SAMPLE_ID"_sorted.bam

fi

done < ${ALL_SAMPLES}
