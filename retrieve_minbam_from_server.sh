#!/bin/bash
set -ueo pipefail
mkdir min_trimmed_bamfiles
for SAMPLE in $(cat samples_remained.txt)
do
	scp woaiyaoy@hoffman2.idre.ucla.edu:/u/flashscratch/w/woaiyaoy/min_trimmed_bamfiles/${SAMPLE}Aligned.sortedByCoord.out.bam ./min_trimmed_bamfiles
done
