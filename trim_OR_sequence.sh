for SAMPLE in $(cat sample_names.txt)

do
        #echo ${SAMPLE} 
        trimmomatic PE -phred33 -summary ./trimmomatic_OR_results/${SAMPLE}.summary.txt ./paired/${SAMPLE}_forward_paired.fq.gz ./paired/${SAMPLE}_reverse_paired.fq.gz ./trimmomatic_OR_results/${SAMPLE}_F_paired.fq.gz ./trimmomatic_OR_results/${SAMPLE}_F_unpaired.fq.gz ./trimmomatic_OR_results/${SAMPLE}_R_paired.fq.gz ./trimmomatic_OR_results/${SAMPLE}_R_unpaired.fq.gz ILLUMINACLIP:uniq_overrepresented_sequence.fa:0:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

done
