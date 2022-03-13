
BASEDIR=/scratch/gguo094/batch_v3
OUTPUT=$BASEDIR/output_anvio
CONTIGS=$BASEDIR/megahit_coassembly_cutadapt/co-assembly_contigs.fa
READS=/scratch/gguo094/raw_data
MAPPING=$OUTPUT/mapping




#### IPRS_interest
for annotation in `awk '{print $1}' /scratch/gguo094/batch_v3/output_anvio/SUMMARY-refined_bins_2_FINAL/gene_call_db/iprs_interest.txt`; 
do grep $annotation $OUTPUT/SUMMARY-refined_bins_2_FINAL/bin_by_bin/Bin_*/*-gene_calls.txt | awk -v x=${annotation} -F "/|\t|:" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28}' > $OUTPUT/SUMMARY-refined_bins_2_FINAL/iprs_interest/${annotation}_gene_calls.tmp & done;

cat $OUTPUT/SUMMARY-refined_bins_2_FINAL/iprs_interest/$*gene_calls.tmp > $OUTPUT/SUMMARY-refined_bins_2_FINAL/iprs_interest/iprs_interest-calls.txt    


for annotation in `awk '{print $1}' /scratch/gguo094/batch_v3/output_anvio/SUMMARY-refined_bins_2_FINAL/gene_call_db/iprs_interest.txt`; 
do grep $annotation $OUTPUT/SUMMARY-refined_bins_2_FINAL/bin_by_bin/Bin_*/*-gene_calls.txt | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28}' > $OUTPUT/SUMMARY-refined_bins_2_FINAL/iprs_interest2/${annotation}_gene_calls.tmp & done;

cat $OUTPUT/SUMMARY-refined_bins_2_FINAL/iprs_interest2/*gene_calls.tmp > $OUTPUT/SUMMARY-refined_bins_2_FINAL/iprs_interest2/iprs_interest2-calls.txt    
 
#### GO_TERMS


for annotation in `awk '{print $1}' /scratch/gguo094/batch_v3/output_anvio/SUMMARY-refined_bins_2_FINAL/gene_call_db/go_terms_call.txt`; 
do grep $annotation $OUTPUT/SUMMARY-refined_bins_2_FINAL/bin_by_bin/Bin_*/*-gene_calls.txt | awk -v x=${annotation} -F "/|\t|:" '{print  $0}' > $OUTPUT/SUMMARY-refined_bins_2_FINAL/go_terms/${annotation}_go_terms_calls.tmp & done;

cat $OUTPUT/SUMMARY-refined_bins_2_FINAL/go_terms/*go_terms_calls.tmp > $OUTPUT/SUMMARY-refined_bins_2_FINAL/go_terms/go_terms_calls.txt     



###########################################################################################################################################################################################################################

## antimicrobial resistant gene using fargene


module load seqtk
module load emboss

for bins in $OUTPUT/SUMMARY-refined_bins_2_FINAL/bin_by_bin/Bins/*.fa
do
echo "$(basename "$bins")"
fargene -i $OUTPUT/SUMMARY-refined_bins_2_FINAL/bin_by_bin/Bins/"$(basename "$bins")" --hmm-model class_a -o $OUTPUT/SUMMARY-refined_bins_2_FINAL/bin_by_bin/Bins/class_a/"$(basename "$bins")"-class_a --force
fargene -i $OUTPUT/SUMMARY-refined_bins_2_FINAL/bin_by_bin/Bins/"$(basename "$bins")" --hmm-model class_b_1_2 -o $OUTPUT/SUMMARY-refined_bins_2_FINAL/bin_by_bin/Bins/class_b_1_2/"$(basename "$bins")"-class_b_1_2 --force
fargene -i $OUTPUT/SUMMARY-refined_bins_2_FINAL/bin_by_bin/Bins/"$(basename "$bins")" --hmm-model class_b_3 -o $OUTPUT/SUMMARY-refined_bins_2_FINAL/bin_by_bin/Bins/class_3/"$(basename "$bins")"-class_b_3 --force
fargene -i $OUTPUT/SUMMARY-refined_bins_2_FINAL/bin_by_bin/Bins/"$(basename "$bins")" --hmm-model class_c -o $OUTPUT/SUMMARY-refined_bins_2_FINAL/bin_by_bin/Bins/class_c/"$(basename "$bins")"-class_c --force
fargene -i $OUTPUT/SUMMARY-refined_bins_2_FINAL/bin_by_bin/Bins/"$(basename "$bins")" --hmm-model class_d_1 -o $OUTPUT/SUMMARY-refined_bins_2_FINAL/bin_by_bin/Bins/class_d_1/"$(basename "$bins")"-class_d_1 --force
fargene -i $OUTPUT/SUMMARY-refined_bins_2_FINAL/bin_by_bin/Bins/"$(basename "$bins")" --hmm-model class_d_2 -o $OUTPUT/SUMMARY-refined_bins_2_FINAL/bin_by_bin/Bins/class_d_2/"$(basename "$bins")"-class_d_2 --force
fargene -i $OUTPUT/SUMMARY-refined_bins_2_FINAL/bin_by_bin/Bins/"$(basename "$bins")" --hmm-model qnr -o $OUTPUT/SUMMARY-refined_bins_2_FINAL/bin_by_bin/Bins/qnr/"$(basename "$bins")"-qnr --force
done

#put all aa file of amr hits in the same folder (has bin name and amr class ID)
cp /scratch/gguo094/batch_v3/output_anvio/SUMMARY-refined_bins_2_FINAL/bin_by_bin/Bins/*/*/predictedGenes/*-filtered-peptides.fasta fargene_amr_all/



##### TRIAL / NOT USED 



for annotation in `awk '{print $1}' /scratch/gguo094/batch_v3/output_anvio/SUMMARY-refined_bins_2_FINAL/gene_call_db/gene_database_all.txt`; 
do grep $annotation $OUTPUT/SUMMARY-refined_bins_2_FINAL/bin_by_bin/Bin_*/*-gene_calls.txt | cut -f 1,2,6,7| awk -v x=${annotation} -F "/|\t|:" '{print x,$9,$10, $11,$15,$16}' > $OUTPUT/SUMMARY-refined_bins_2_FINAL/tigrfram_marker_genes/${annotation}_gene_calls.tmp & done;
 
 
cat $OUTPUT/SUMMARY-refined_bins_2_FINAL/cogs_marker_genes/*_gene_calls.tmp > $OUTPUT/SUMMARY-refined_bins_2_FINAL/marker_gene_compiled/cogs_gene_calls.txt                                         


for annotation in `awk '{print $1}' /scratch/gguo094/batch_v3/output_anvio/SUMMARY-refined_bins_2_FINAL/gene_call_db/gene_database_all_cogs.txt`; do grep $annotation $OUTPUT/SUMMARY-refined_bins_2_FINAL/bin_by_bin/Bin_*/*-gene_calls.txt | cut -f 1,2,6,7| awk -v x=${annotation} -F "/|\t|:" '{print x,$9,$10, $11,$15,$16}' > $OUTPUT/SUMMARY-refined_bins_2_FINAL/cogs_marker_genes/${annotation}_gene_calls.tmp & done;



for annotation in `awk '{print $1}' /scratch/gguo094/batch_v3/output_anvio/SUMMARY-refined_bins_2_FINAL/gene_call_db/gene_database_all_cogs.txt`; do grep $annotation $OUTPUT/SUMMARY-refined_bins_2_FINAL/bin_by_bin/Bin_*/*-gene_calls.txt | cut -f 1,2,3,4,5,6,7 | awk -v x=${annotation} -F "/|\t|:" '{print x,$9,$10, $11,12, $13}' > $OUTPUT/SUMMARY-refined_bins_2_FINAL/trial_marker_genes/${annotation}_gene_calls.tmp & done;


for annotation in `awk '{print $1}' /scratch/gguo094/batch_v3/output_anvio/SUMMARY-refined_bins_2_FINAL/gene_call_db/gene_database_all_cogs.txt`; do grep $annotation $OUTPUT/SUMMARY-refined_bins_2_FINAL/bin_by_bin/Bin_*/*-gene_calls.txt | cut -f1,2,3,4,5,6,7 | awk -v x=${annotation} -F "/|\t|:" '{print x, $10, $11, $13, $14, $15, $16, $17,$23}' > $OUTPUT/SUMMARY-refined_bins_2_FINAL/trial_marker_genes/${annotation}_gene_calls.tmp & done;


for annotation in `awk '{print $1}' /scratch/gguo094/batch_v3/output_anvio/SUMMARY-refined_bins_2_FINAL/gene_call_db/gene_database_draft2.txt`; do grep $annotation $OUTPUT/SUMMARY-refined_bins_2_FINAL/bin_by_bin/Bin_*/*-gene_calls.txt | cut -f 1,2,3,4,5,6,7| awk -v x=${annotation} -F "/|\t|:" '{print x,$1,$2,$3,$4,$5,$6,$7,$8, $9,$10, $11,$15,$16}' > $OUTPUT/SUMMARY-refined_bins_2_FINAL/gene_database_draft2/${annotation}_gene_calls.tmp & done;

 
cat $OUTPUT/SUMMARY-refined_bins_2_FINAL/gene_database_draft2/*_gene_calls.tmp > $OUTPUT/SUMMARY-refined_bins_2_FINAL/gene_database_draft2/gene_database_draft2-calls.txt    