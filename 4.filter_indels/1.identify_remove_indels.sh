      # identify indels
perl ~/software/popoolation_1.2.2/basic-pipeline/identify-genomic-indel-regions.pl --indel-window 5 --min-count 5 --input /media/raid/home/schmielewski/sex_ratio_EE/results/3.create_mpileup_files/EE_SR.mpileup --output /media/raid/home/schmielewski/sex_ratio_EE/results/4.filter_indels/EE_SR_indels.gtf
        # not sure that --min-count to use. Barghi et al. 2017 used --min-count 5 or --min-count 2% (?) and a value between 2 and 5 is most usually used. However, Franssen et al. 2015 used –min-count 90
       
        # remove indels
perl ~/software/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input /media/raid/home/schmielewski/sex_ratio_EE/results/3.create_mpileup_files/EE_SR.mpileup --gtf /media/raid/home/schmielewski/sex_ratio_EE/results/4.filter_indels/EE_SR_indels.gtf --output /media/raid/home/schmielewski/sex_ratio_EE/results/4.filter_indels/EE_SR_IndelsRm.mpileup
