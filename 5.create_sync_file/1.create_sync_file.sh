PILEUP_INPUT=/media/raid/home/schmielewski/sex_ratio_EE/results/4.filter_indels/EE_SR_IndelsRm.mpileup # with filtered indels
SYNC_OUTPUT=/sex_ratio_EE/results/5.create_sync_file/EE_SR_IndelsRm.sync

java -ea -Xmx60g -jar ~/software/popoolation2_1201/mpileup2sync.jar \
	--input ${PILEUP_INPUT} \
	--output ${SYNC_OUTPUT} \
	--fastq-type sanger \
	--min-qual 20 \
	--threads 50
