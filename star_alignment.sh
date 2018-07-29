for SAMPLE in $(cat samples_remained.txt)
do 
	STAR --runThreadN 16 --alignEndsType EndToEnd --genomeDir ./genome/ --sjdbGTFfile ./genome/ensembl_zv11_rRNArm.gtf --sjdbOverhang 100 --readFilesIn ./paired/${SAMPLE}_forward_paired.fq.gz ./paired/${SAMPLE}_reverse_paired.fq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ./bamfiles/${SAMPLE}
Tdone
