#!/bin/bash

# pipeline for trimming reads with Trimmomatic, mapping the reads with bwa mem, masking and removing duplicates with samtools
# modified by SCh based on JP scripts

##prepare refernce genome before running (e.g. index) - to avoid doing it each time for each sample

FQPATH=/media/nas/EE_SR                             #<--- PATH TO FASTQ FILES
REFSEQ=~/sex_ratio_EE/data/Rhrob_anchored.fa     #<--- INDEXED (SAMTOOLS), GATK DICTIONARY, BWA INDEX

ADAPTERS=TruSeq3-PE.fa
ADAPTERS_PATH=~/software/Trimmomatic-0.39/adapters

OUTDIR=~/sex_ratio_EE/results/2.trim_map_markDup_reads

###INFORMATION PROVIDED FROM COMMAND LINE:

IDNAME=$1     #INDIVIDUAL ID    $1
FASTQ_1=$2    #FASTQ R1 FILE    $2
FASTQ_2=$3    #FASTQ R2 FILE    $3


##COMPUTATIONAL OPTIONS (MEM, CPUs etc.)

NTHREADS=50

##LINKS TO SOFTWARE NOT IN $PATH
TRIMMOMATIC=~/software/Trimmomatic-0.39/trimmomatic-0.39.jar
PICARD=~/software/picard/build/libs/picard.jar


##!!!!!!!!!!!!!
#########################
#######END OF EDITING####
#########################


echo Starting at `date`

#=======DIR FOR SAMPLE===============
##create tmp working directory and copy fastq files
if [ ! -d "${OUTDIR}/${IDNAME}" ]
then
    mkdir ${OUTDIR}/${IDNAME}
    mkdir ${OUTDIR}/${IDNAME}/tmp_${IDNAME}
fi

cd ${OUTDIR}/${IDNAME}


#====================================
#====TRIM READS WITH TRIMMOMATIC=====
#====================================
echo Starting trimming reads at `date`

cp ${ADAPTERS_PATH}/${ADAPTERS} .
java -jar ${TRIMMOMATIC} PE -threads ${NTHREADS} ${FQPATH}/${FASTQ_1} ${FQPATH}/${FASTQ_2} ${FASTQ_1}_trimmed.fastq_P.gz ${FASTQ_1}_trimmed.fastq_U.gz \
   ${FASTQ_2}_trimmed.fastq_P.gz ${FASTQ_2}_trimmed.fastq_U.gz ILLUMINACLIP:${ADAPTERS}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 


echo Finished trimming reads at `date`

#====================================
#===========MAPPING==================
#====================================

echo Starting mapping at `date`

R1=${FASTQ_1}_trimmed.fastq_P.gz
R2=${FASTQ_2}_trimmed.fastq_P.gz

bwa mem -t ${NTHREADS} -R "@RG\tID:${IDNAME}\tSM:${IDNAME}\tLB:library1" ${REFSEQ} ${R1} ${R2} | samtools sort -@${NTHREADS} -n -O BAM -o ${IDNAME}.out_sorted.bam -
# rm -f ./*.gz  #remove it or not


#====================================
#=======ADD TAGS USING FIXMATE ======
#====================================
echo Starting adding tags for duplicate marking at `date`

samtools fixmate --threads ${NTHREADS} -m ${IDNAME}.out_sorted.bam ${IDNAME}.out_sorted.fixmate.bam
samtools sort ${IDNAME}.out_sorted.fixmate.bam --threads ${NTHREADS} -O BAM -o ${IDNAME}.out_sorted.fixmate_sorted.bam
#rm ${IDNAME}.out_sorted.fixmate.bam

#====================================
#=======MARK DUPLICATES =============
#====================================
echo Starting marking duplicates at `date`

samtools markdup --threads ${NTHREADS} ${IDNAME}.out_sorted.fixmate_sorted.bam ${IDNAME}.dupmarked.bam

#==========================================
#=======REMOVE AMBIGOUS READS =============
#==========================================
echo Starting removing ambigous reads at `date`

samtools view --threads ${NTHREADS} -q 20 -bS ${IDNAME}.dupmarked.bam | samtools sort --threads ${NTHREADS} -o ${IDNAME}.dupmarked.AmbigRm.bam

#======================
#=======INDEX BAM======
#======================

samtools index ${IDNAME}.dupmarked.AmbigRm.bam

echo Finished mapping at `date`

cd ~


