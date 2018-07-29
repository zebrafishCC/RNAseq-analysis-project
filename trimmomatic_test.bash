
#trim adapter for the first time
for SAMPLE in $(cat samples_remained.txt)

do
	#echo ${SAMPLE} 
	trimmomatic PE -phred33 -summary ./trimmomatic_results/${SAMPLE}.summary.txt ./zebrafish/${SAMPLE}_L003_R1_001.fastq.gz ./zebrafish/${SAMPLE}_L003_R2_001.fastq.gz ./trimmomatic_results/${SAMPLE}_forward_paired.fq.gz ./trimmomatic_results/${SAMPLE}_forward_unpaired.fq.gz ./trimmomatic_results/${SAMPLE}_reverse_paired.fq.gz ./trimmomatic_results/${SAMPLE}_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

done



#trim the overrepresented sequence for the second time on hoffman2 cluster
for SAMPLE in $(cat sample_names.txt)

do
        #echo ${SAMPLE} 
        trimmomatic PE -phred33 -summary ./trimmomatic_OR_results/${SAMPLE}.summary.txt ./paired/${SAMPLE}_forward_paired.fq.gz ./paired/${SAMPLE}_reverse_paired.fq.gz ./trimmomatic_OR_results/${SAMPLE}_F_paired.fq.gz ./trimmomatic_OR_results/${SAMPLE}_F_unpaired.fq.gz ./trimmomatic_OR_results/${SAMPLE}_R_paired.fq.gz ./trimmomatic_OR_results/${SAMPLE}_R_unpaired.fq.gz ILLUMINACLIP:uniq_overrepresented_sequence.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

done
