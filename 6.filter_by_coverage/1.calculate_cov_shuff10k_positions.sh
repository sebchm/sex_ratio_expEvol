#! /bin/bash
# calculate coverage for subsampled sync file (1/10k lines sampled)

SYNC_FILE=/media/raid/home/schmielewski/sex_ratio_EE/results/5.create_sync_file/EE_SR_IndelsRm.sync

awk 'NR % 10000  == 0' ${SYNC_FILE} | python /media/raid/home/schmielewski/sex_ratio_EE/bin/6.filter_by_coverage/GetSumCovFromSync.py > /media/raid/home/schmielewski/sex_ratio_EE/results/6.filter_by_coverage/EE_SR_IndelsRm_SumCovShuf.OneEverytenthous.txt
