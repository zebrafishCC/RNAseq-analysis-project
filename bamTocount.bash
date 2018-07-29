for SAMPLE in $(cat sample_names.txt)

do
	htseq-count --quiet --format=bam --stranded=no --order=pos ./bamfiles/${SAMPLE}Aligned.sortedByCoord.out.bam /home/chengchen/data/RNAseq/zv11_ensembl.gtf > ./countfile/${SAMPLE}_count.txt

done
