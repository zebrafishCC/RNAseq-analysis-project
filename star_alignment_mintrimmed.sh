set -ueo pipefail

#--alignEndsType EndToEnd is only used for mim_trimmed data to facilitate rMATS alternative expression analysis. Be careful!!!

#shared memory
STAR --genomeLoad LoadAndExit --genomeDir ./genome/

for SAMPLE in $(cat samples_remained.txt)
do 
	STAR --runThreadN 16 --alignEndsType EndToEnd --genomeDir ./genome/ --genomeLoad LoadAndKeep --limitBAMsortRAM 20000000000 --sjdbOverhang 100 --readFilesIn ./min_trimmed/${SAMPLE}_F_paired.fq.gz ./min_trimmed/${SAMPLE}_R_paired.fq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ./min_trimmed_bamfiles/${SAMPLE}

done
STAR --genomeLoad Remove --genomeDir ./genome/

