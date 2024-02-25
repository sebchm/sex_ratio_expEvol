
# activate conda
source ~/../anaconda/anaconda3/bin/activate
conda activate btk

cd ~/sex_ratio_EE/results/1_A.blobtools

# get taxdump
mkdir -p taxdump
cd taxdump
curl -L ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz | tar xzf -

# make blast db
makeblastdb -in (...)/Rhrob_anchored.fa \
  -input_type fasta -dbtype nucl \
  -title Rhrob_anchored -parse_seqids \
  -taxid 223528 \
  -out Rhrob_anchored_blastdb

# blast the ref. genome
blastn -db Rhrob_anchored_blastdb \
       -query (...)Rhrob_anchored.fa \
       -outfmt "6 qseqid staxids bitscore std" \
       -max_target_seqs 10 \
       -max_hsps 1 \
       -evalue 1e-25 \
       -num_threads 16 \
       -out blast.out


  
