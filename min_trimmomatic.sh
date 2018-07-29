set -ueo pipefail
cat sample_names.txt| xargs -n 1 -P 4 -I{} trimmomatic PE -phred33 ./trimmomatic_results/paired/{}_forward_paired.fq.gz ./trimmomatic_results/paired/{}_reverse_paired.fq.gz  ./min_trimmed/{}_F_paired.fq.gz ./min_trimmed/{}_F_unpaired.fq.gz ./min_trimmed/{}_R_paired.fq.gz ./min_trimmed/{}_R_unpaired.fq.gz CROP:120 MINLEN:120

