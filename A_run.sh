#!/bin/bash
#$ -S /bin/bash
#$ -o /isilon/ottawa-rdc/users/shared/chenw_lab/galen/temp_files
#$ -e /isilon/ottawa-rdc/users/shared/chenw_lab/galen/temp_files


source ~/miniconda3/bin/activate
conda activate anvio-6

#directories

GALEN=/isilon/ottawa-rdc/users/shared/chenw_lab/galen
WORK=$GALEN/batch
RAW=$WORK/raw
output=$WORK/output
MAPPING=$output/mapping
TRIM=$WORK/output/trimmomatic
MEGAHIT=$output/megahit_coassembly
ANVIO_A=$output/anvio/indA
ANVIO_B=$output/anvio/indB


bowtie2-build $ANVIO_A/IndA_contigs.fa $MAPPING/IndA/IndA_contigs --threads 25

cd $TRIM

for sample in `awk '{print $1}' $MAPPING/sample_name`;
do
R1="${sample}.pair1.fq.gz"
R2="${sample}.pair2.fq.gz"
sam="${sample}.sam" 
bowtie2 -x $MAPPING/IndA/IndA_contigs -1 $R1 -2 $R2 -S $MAPPING/$sam --threads 25

done
conda deactivate


echo convert sam to bam

for sample in `awk '{print $1}' $MAPPING/sample_name`;
do
sam=".sam"
bam=".bam"
samtools view -S -b $MAPPING/${sample}.sam > $MAPPING/${sample}.bam
done


conda activate anvio-7


## bam to anvio_bam profiling

for sample in `awk '{print $1}' $MAPPING/sample_name`;
do
anvi-init-bam $MAPPING/${sample}.bam -o $MAPPING/${sample}_anvi.bam
done


echo mapping_done

anvi-gen-contigs-database -f $ANVIO_A/IndA_contigs.fa -o $ANVIO_A/IndA_contigs.db -n IndA --num-threads 25
anvi-run-hmms -c $ANVIO_A/IndA_contigs.db -T 50
anvi-run-ncbi-cogs -c $ANVIO_A/IndA_contigs.db --cog-data-dir $GALEN/COG --num-threads 20 --search-with diamond 
anvi-get-sequences-for-gene-calls -c $ANVIO_A/IndA_contigs.db \
                                    --get-aa-sequences \
                                    -o $ANVIO_A/IndA_amino-acid-sequences.fa
anvi-display-contigs-stats $ANVIO_A/IndA_contigs.db  --report-as-text  --output-file $ANVIO_A/IndA_contigs_post-hmm-cogs.txt
anvi-import-functions -c $ANVIO_A/IndA_contigs.db -i $ANVIO_A/interpro_output/indA_iprs2anvio.tsv
file_ext="_anvi.bam"
for sample in `awk '{print $1}' $ANVIO_A/profile/sample_name`;
do
echo "$sample$file_ext"
anvi-profile -i $MAPPING/"$sample$file_ext" -c $ANVIO_A/IndA_contigs.db --output-dir $ANVIO_A/profile/$sample --sample-name $sample -T 100 --min-contig-length 2500
done
