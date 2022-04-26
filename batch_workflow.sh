#!/bin/bash
#$ -S /bin/bash


source ~/miniconda3/bin/activate
conda activate anvio7

#directories

GALEN=/isilon/ottawa-rdc/users/shared/chenw_lab/galen
WORK=$GALEN/batch
RAW=$WORK/raw
output=$WORK/output
MAPPING=$output/mapping
TRIM=$WORK/output/trimmomatic
MEGAHIT=$output/megahit_coassembly
ANVIO=$output/anvio
ANVIO_A=$ANVIO/indA
ANVIO_B=$ANVIO/indB

####################################################################
### QC of fastqc pre-trim
####################################################################

cd $RAW
mkdir fastqc
fastqc $RAW/*.gz -t 10 -o $RAW/fastqc/ -f fastq
multiqc $RAW/fastqc/ -o $RAW/fastqc/

###################################################################
### removal of adapter using trimmomatic
###################################################################

echo Trimming
mkdir $OUTPUT
mkdir $OUTPUT/trimmomatic

# ls $RAW/* > $RAW/sample_name

## manually remove extension


for i in `awk '{print $1}' $RAW/sample_name`;
do
echo $i
    file1=$RAW/${i}_R1.fastq.gz
    file2=$RAW/${i}_R2.fastq.gz
echo $file1
echo $file2
	java -jar $GALEN/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
	-threads 20 \
	$file1 $file2 \
	$TRIM/$i.pair1.fq.gz $TRIM/$i.unpair1.fq.gz \
	$TRIM/$i.pair2.fq.gz $TRIM/$i.unpair2.fq.gz \
	ILLUMINACLIP:$GALEN/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10 \
	LEADING:2 TRAILING:2 SLIDINGWINDOW:3:15 MINLEN:36 HEADCROP:6
done

echo Trimming_done


####################################################################
### QC of fastqc post-trim
####################################################################
echo QC_of_Trimming

mkdir $TRIM/fastqc
fastqc -t 20 $TRIM/*fq.gz -o $TRIM/fastqc/
multiqc $TRIM/fastqc/*.pair* -o $TRIM/fastqc/

echo QC_of_Trimming_done



####################################################################
### co-assembly megahit
###
####################################################################

echo co-assembling
mkdir $MEGAHIT/
indA_F=$(echo $(ls $TRIM/IndA_*.pair1.fq.gz)  | sed "s/ /,/g")
indA_R=$(echo $(ls $TRIM/IndA_*.pair2.fq.gz)  | sed "s/ /,/g")


megahit -1 $indA_F -2 $indA_R -o $MEGAHIT/indA -t 20 --min-contig-len 1000 --continue


cp $MEGAHIT/indA/final.contigs.fa $MEGAHIT/indA/IndA_contigs.fa

indB_F=$(echo $(ls $TRIM/IndB_*.pair1.fq.gz)  | sed "s/ /,/g")
indB_R=$(echo $(ls $TRIM/IndB_*.pair2.fq.gz)  | sed "s/ /,/g")


megahit -1 $indB_F -2 $indB_R -o $MEGAHIT/indB -t 20 --min-contig-len 1000 --continue

cp $MEGAHIT/indB/final.contigs.fa $MEGAHIT/indB/IndB_contigs.fa


echo co-assembling_done



####################################################################
### loading into anvio
###
####################################################################

mkdir $ANVIO
mkdir $ANVIO_A
mkdir $ANVIO_B
cd $ANVIO
cp $MEGAHIT/IndA_contigs.fa $ANVIO_A
cp $MEGAHIT/IndB_contigs.fa $ANVIO_B

anvi-script-reformat-fasta $ANVIO_A/IndA_contigs.fa -o $ANVIO_A/IndA_contigs.fa -l 0 --simplify-names
anvi-script-reformat-fasta $ANVIO_B/IndB_contigs.fa -o $ANVIO_B/IndB_contigs.fa -l 0 --simplify-names



###create anvio database
GRDI

anvi-gen-contigs-database -f $ANVIO_A/IndA_contigs.fa -o $ANVIO_A/IndA_contigs.db -n IndA --num-threads 25
anvi-gen-contigs-database -f $ANVIO_B/IndB_contigs.fa -o $ANVIO_B/IndB_contigs.db -n IndB --num-threads 25


anvi-run-hmms -c $ANVIO_A/IndA_contigs.db -T 24
anvi-run-hmms -c $ANVIO_B/IndB_contigs.db -T 24

anvi-run-kegg-kofams -c $ANVIO_A/IndA_contigs.db -T 24 --hmmer-program hmmsearch  --just-do-it
anvi-run-kegg-kofams -c $ANVIO_B/IndB_contigs.db -T 24 --hmmer-program hmmsearch  --just-do-it

anvi-run-ncbi-cogs -c $ANVIO_A/IndA_contigs.db --cog-data-dir $GALEN/COG --num-threads 24 --search-with diamond
anvi-run-ncbi-cogs -c $ANVIO_B/IndB_contigs.db --cog-data-dir $GALEN/COG --num-threads 24 --search-with diamond

anvi-get-sequences-for-gene-calls -c $ANVIO_A/IndA_contigs.db \
                                    --get-aa-sequences \
                                    -o $ANVIO_A/IndA_amino-acid-sequences.fa
anvi-get-sequences-for-gene-calls -c $ANVIO_B/IndB_contigs.db \
                                    --get-aa-sequences \
                                    -o $ANVIO_B/IndB_amino-acid-sequences.fa


anvi-run-pfams -c $ANVIO_A/IndA_contigs.db --pfam-data-dir $GALEN/database/pfam --hmmer-program hmmscan --num-threads 24
anvi-run-pfams -c $ANVIO_B/IndB_contigs.db --pfam-data-dir $GALEN/database/pfam --hmmer-program hmmscan --num-threads 24
### estimate taxonomy of bins (using GTDB)
anvi-run-scg-taxonomy -c $ANVIO_A/IndA_contigs.db -T 24
anvi-run-scg-taxonomy -c $ANVIO_B/IndB_contigs.db -T 24

anvi-display-contigs-stats $ANVIO_A/IndA_contigs.db  --report-as-text  --output-file $ANVIO_A/IndA_stats.txt
anvi-display-contigs-stats $ANVIO_B/IndB_contigs.db  --report-as-text  --output-file $ANVIO_B/IndB_stats.txt

#
# ####################################################################
# ### GENE FUNCTION CALLS
# ### interproscan on amino-acid seq from anvio
# ####################################################################
#
# mkdir $ANVIO_A/iprs/
#
# perl $GALEN/split_fasta.pl -i $ANVIO_A/IndA_amino-acid-sequences.fa -o $ANVIO_A/iprs/iprs -n 5000
#
# mkdir $ANVIO_B/iprs/
#
# perl $GALEN/split_fasta.pl -i $ANVIO_B/IndB_amino-acid-sequences.fa -o $ANVIO_B/iprs/iprs -n 5000
#
# ## Create sample list of all file created
#
# find $ANVIO_A/iprs -printf "%f\n" > $ANVIO_A/iprs/iprs_list
# find $ANVIO_B/iprs -printf "%f\n" > $ANVIO_B/iprs/iprs_list
#
#
# ## repeat below for each list.
#
# mkdir $ANVIO_A/interpro_output
#
# for sample in `awk '{print $1}' $ANVIO_A/iprs/iprs_list`;
# do
# /isilon/ottawa-rdc/users/shared/chenw_lab/galen/interproscan-5.39-77.0/interproscan.sh -i $ANVIO_A/iprs/${sample} -f tsv \
# 	-d $ANVIO_A/interpro_output/ \
# 	--tempdir $WORK/temp/ \
# 	--disable-precalc \
# 	-appl Pfam,PIRSF,SUPERFAMILY,TIGRFAM \
# 	--iprlookup \
# 	--goterms \
# 	--pathways --cpu 10
# done
#
# mkdir $ANVIO_B/interpro_output
#
# for sample in `awk '{print $1}' $ANVIO_A/iprs/iprs_list`;
# do
# /isilon/ottawa-rdc/users/shared/chenw_lab/galen/interproscan-5.39-77.0/interproscan.sh -i $ANVIO_B/iprs/${sample} -f tsv \
# 	-d $ANVIO_B/interpro_output/ \
# 	--tempdir $WORK/temp/ \
# 	--disable-precalc \
# 	-appl Pfam,PIRSF,SUPERFAMILY,TIGRFAM \
# 	--iprlookup \
# 	--goterms \
# 	--pathways --cpu 10
# done
#
# ### create new folder to house concatenate of all smaller fxnal annotation output (the script dont like a folder with too much clutter, im guessing)
#
# cat $ANVIO_A/interpro_output/*.tsv > $ANVIO_A/interpro_output/indA_all_iprs.tsv
# cat $ANVIO_B/interpro_output/*.tsv > $ANVIO_B/interpro_output/indB_all_iprs.tsv
#
# ### script to clean up and allow import to anvio
#
# ## create iprs2anvio.sh file found here: https://github.com/xvazquezc/stuff/blob/master/iprs2anvio.sh
#
# /isilon/ottawa-rdc/users/shared/chenw_lab/galen/interproscan-5.48-83.0/iprs2anvio.sh -i $ANVIO_A/interpro_output/indA_all_iprs.tsv -o $ANVIO_A/interpro_output/indA -g -p -ipr
#
# /isilon/ottawa-rdc/users/shared/chenw_lab/galen/interproscan-5.48-83.0/iprs2anvio.sh -i $ANVIO_B/interpro_output/indB_all_iprs.tsv -o $ANVIO_B/interpro_output/indB -g -p -r
#
#
# ### importing functional annotation to anvio
# cp $ANVIO_B/IndB_contigs_backup.db $ANVIO_B/IndB_contigs.db
# cp $ANVIO_A/IndA_contigs_backup.db $ANVIO_A/IndA_contigs.db
#
# anvi-import-functions -c $ANVIO_A/IndA_contigs.db -i $ANVIO_A/interpro_output/indA_iprs2anvio.tsv
# anvi-import-functions -c $ANVIO_B/IndB_contigs.db -i $ANVIO_B/interpro_output/indB_iprs2anvio.tsv
#
#
#
#
# ####################################################################
# ### TAXONOMY CALLS
# ### centrifuge  on amino-acid seq from anvio
# ####################################################################
# mkdir $ANVIO/centrifuge
#
#
# centrifuge -f -x $CENTRIFUGE_BASE/p+h+v/p+h+v $ANVIO/amino-acid-sequences.fa -S $ANVIO/centrifuge/centrifuge_hits.tsv
#
# ## Make sure there is two files in the work directory ($ANVIO/centrifuge/)
#
# anvi-import-taxonomy-for-genes -c $ANVIO/contigs.db -i $ANVIO/centrifuge/centrifuge_report.tsv $ANVIO/centrifuge/centrifuge_hits.tsv -p centrifuge

### adding Single copy gene into database

# run this first, only once.

####################################################################
### mapping with bowtie2
###
####################################################################

echo mapping

## create mapping directory

## building mapping files
mkdir $MAPPING
mkdir $MAPPING/IndA
mkdir $MAPPING/IndB
cd $MAPPING

bowtie2-build $ANVIO/indA/IndA_contigs.fa $MAPPING/IndA/IndA_contigs --threads 25
bowtie2-build $ANVIO/indB/IndB_contigs.fa $MAPPING/IndB/IndB_contigs --threads 25

## transfering sample name file to mapping folder.

cp $RAW/sample_name $MAPPING

#CHANGE SAMPLE NAME FROM A1-A6 TO IndA.... etc

### mapping read to contigs

cd $TRIM

for sample in `awk '{print $1}' $MAPPING/sample_name`;
do
R1="${sample}.pair1.fq.gz"
R2="${sample}.pair2.fq.gz"
sam="${sample}.sam"
bowtie2 -x $MAPPING/IndA/IndA_contigs -1 $R1 -2 $R2 -S $MAPPING/$sam --threads 25
bowtie2 -x $MAPPING/IndB/IndB_contigs -1 $R1 -2 $R2 -S $MAPPING/$sam --threads 25

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

####################################################################
### Profiling cont'd very long. bam --> anvio
### PROFILE DONT LIKE SAMPLE NAME WITH "-"
####################################################################

### change bam file name from "-" to "_"

cd $MAPPING
find . -depth -name '*-*' -exec rename '-' '_' {} +


## ### do everything in one file, using more cpu. (5 days)

# create folder to store single anvio profile before merging and transfering list of all samples to the appropriate folders.
mkdir $ANVIO_A/profile
cp $MAPPING/sample_name $ANVIO_A/profile

mkdir $ANVIO_B/profile
cp $MAPPING/sample_name $ANVIO_B/profile

# manually delete the sample for indB in the in indA folder and vice versa.
file_ext="_anvi.bam"
for sample in `awk '{print $1}' $ANVIO_A/profile/sample_name`;
do
echo "$sample$file_ext"
anvi-profile -i $MAPPING/"$sample$file_ext" -c $ANVIO_A/IndA_contigs.db --output-dir $ANVIO_A/profile/$sample --sample-name $sample -T 100 --min-contig-length 2500
done
# manually delete the sample for indB in the in indA folder and vice versa.
file_ext="_anvi.bam"
for sample in `awk '{print $1}' $ANVIO_B/profile/sample_name`;
do
echo "$sample$file_ext"
anvi-profile -i $MAPPING/"$sample$file_ext" -c $ANVIO_B/IndB_contigs.db --output-dir $ANVIO_B/profile/$sample --sample-name $sample -T 100 --min-contig-length 2500
done


## merge profile

anvi-merge $ANVIO_A/profile/*/PROFILE.db -o $ANVIO_A/profile_merged -c $ANVIO_A/IndA_contigs.db -S indA --skip-hierarchical-clustering -W
anvi-merge $ANVIO_B/profile/*/PROFILE.db -o $ANVIO_B/profile_merged -c $ANVIO_B/IndB_contigs.db -S indB --skip-hierarchical-clustering -W



####################################################################
### BINNING! Finally
### Concoct
####################################################################


anvi-cluster-contigs -p $ANVIO_A/profile_merged/PROFILE.db -c $ANVIO_A/IndA_contigs.db -C indA_concoct --driver concoct -T 100 --just-do-it
anvi-cluster-contigs -p $ANVIO_B/profile_merged/PROFILE.db -c $ANVIO_B/IndB_contigs.db -C indB_concoct --driver concoct -T 100 --just-do-it

####################################################################
### metabat2
####################################################################


anvi-cluster-contigs -p $ANVIO_A/profile_merged/PROFILE.db -c $ANVIO_A/IndA_contigs.db -C indA_metabat2 --driver metabat2 -T 100 --just-do-it
anvi-cluster-contigs -p $ANVIO_B/profile_merged/PROFILE.db -c $ANVIO_B/IndB_contigs.db -C indB_metabat2 --driver metabat2 -T 100 --just-do-it


####################################################################
### maxbin2
####################################################################


anvi-cluster-contigs -p $ANVIO_A/profile_merged/PROFILE.db -c $ANVIO_A/IndA_contigs.db -C indA_metabat2 --driver maxbin2 -T 100 --just-do-it
anvi-cluster-contigs -p $ANVIO_B/profile_merged/PROFILE.db -c $ANVIO_B/IndB_contigs.db -C indB_metabat2 --driver maxbin2 -T 100 --just-do-it

####################################################################
### dastool
####################################################################

anvi-cluster-contigs -p $ANVIO_A/profile_merged/PROFILE.db -c $ANVIO_A/IndA_contigs.db -S indA_concoct,indA_metabat2 --search-engine usearch --driver dastool -C indA_dastool -T 100 --just-do-it
anvi-cluster-contigs -p $ANVIO_B/profile_merged/PROFILE.db -c $ANVIO_B/IndB_contigs.db -S indB_concoct,indB_metabat2 --search-engine usearch --driver dastool -C indB_dastool -T 100 --just-do-it

####################################################################
### SUMMARIZE!
####################################################################

anvi-summarize -p $ANVIO_A/profile_merged/PROFILE.db -c $ANVIO_A/IndA_contigs.db -o $ANVIO_A/sample_summary_indA_dastool -C indA_dastool --init-gene-coverages
anvi-summarize -p $ANVIO_B/profile_merged/PROFILE.db -c $ANVIO_B/IndB_contigs.db -o $ANVIO_B/sample_summary_indB_dastool -C indB_dastool --init-gene-coverages

####################################################################
### Transfer contig.db and profile.db (and auxiliary) file to local machine to manually refine!
####################################################################

anvi-refine ....

####################################################################
### change bin name to make final collection
####################################################################


## export and change bin name manually
cd $ANVIO_A/

anvi-export-collection -C indA_dastool -p $ANVIO_A/profile_merged/PROFILE.db
anvi-import-collection $ANVIO_A/collection-indA_dastool.txt --bins-info $ANVIO_A/collection-indA_dastool-info.txt -C FINAL -p $ANVIO_A/profile_merged/PROFILE.db  -c $ANVIO_A/IndA_contigs.db


cd $ANVIO_B/

anvi-export-collection -C indB_dastool -p $ANVIO_B/profile_merged/PROFILE.db
anvi-import-collection $ANVIO_B/collection-indB_dastool.txt --bins-info $ANVIO_B/collection-indB_dastool-info.txt -C FINAL -p $ANVIO_B/profile_merged/PROFILE.db  -c $ANVIO_A/IndA_contigs.db



####################################################################
### add estimate metabolism (KEGG)
####################################################################
# estimate metabolism (KEGG) of collection, must run anvi-setup-kegg-kofams first
anvi-run-kegg-kofams -p $ANVIO_A/profile_merged/PROFILE.db  -c $ANVIO_A/IndA_contigs.db   -C FINAL -T 75 --kegg-data-dir $GALEN/KOFAM/
anvi-run-kegg-kofams -p $ANVIO_B/profile_merged/PROFILE.db  -c $ANVIO_B/IndB_contigs.db   -C FINAL -T 75 --kegg-data-dir $GALEN/KOFAM/

cd $ANVIO_A
anvi-estimate-metabolism -p $ANVIO_A/profile_merged/PROFILE.db -c $ANVIO_A/IndA_contigs.db -C FINAL --kegg-data-dir $GALEN/KOFAM/ --output-file-prefix $ANVIO_A/IndA_KEGG_FINAL --add-coverage
cd $ANVIO_B
anvi-estimate-metabolism -p $ANVIO_B/profile_merged/PROFILE.db -c $ANVIO_B/IndB_contigs.db -C FINAL --kegg-data-dir $GALEN/KOFAM/ --output-file-prefix $ANVIO_B/IndB_KEGG_FINAL --add-coverage


####################################################################
### estimate taxonomy of bins (using GTDB)
####################################################################

anvi-run-scg-taxonomy -c $ANVIO_A/IndA_contigs.db -T 100 --all-hits-output-file $ANVIO_A/SCG_out -P 5
anvi-run-scg-taxonomy -c $ANVIO_B/IndB_contigs.db -T 100 --all-hits-output-file $ANVIO_B/SCG_out -P 5

cd $ANVIO_A/
anvi-estimate-scg-taxonomy -p $ANVIO_A/profile_merged/PROFILE.db -c $ANVIO_A/IndA_contigs.db -C FINAL -T 100


cd $ANVIO_B/
anvi-estimate-scg-taxonomy -p $ANVIO_B/profile_merged/PROFILE.db -c $ANVIO_B/IndB_contigs.db -C FINAL -T 100


####################################################################
### summarize! again..
####################################################################
anvi-summarize -p $ANVIO_A/profile_merged/PROFILE.db -c $ANVIO_A/IndA_contigs.db -o $ANVIO_A/sample_summary_indA_FINAL -C FINAL --init-gene-coverages
anvi-summarize -p $ANVIO_B/profile_merged/PROFILE.db -c $ANVIO_B/IndB_contigs.db -o $ANVIO_B/sample_summary_indB_FINAL -C FINAL --init-gene-coverages



####################################################################
### checkm
####################################################################
mkdir $ANVIO_A/sample_summary_indA_FINAL/bin_by_bin/all_fas
cp $ANVIO_A/sample_summary_indA_FINAL/bin_by_bin/*/*.fa $ANVIO_A/sample_summary_indA_FINAL/bin_by_bin/all_fas

mkdir $ANVIO_B/sample_summary_indB_FINAL/bin_by_bin/all_fas
cp $ANVIO_B/sample_summary_indB_FINAL/bin_by_bin/*/*.fa $ANVIO_B/sample_summary_indB_FINAL/bin_by_bin/all_fas

checkm lineage_wf $ANVIO_A/sample_summary_indA_FINAL/bin_by_bin/all_fas $ANVIO_A/sample_summary_indA_FINAL/checkm -x fa -t 100
checkm lineage_wf $ANVIO_B/sample_summary_indB_FINAL/bin_by_bin/all_fas $ANVIO_B/sample_summary_indB_FINAL/checkm -x fa -t 100

####################################################################
### DRAM
####################################################################


conda activate DRAM

cd $ANVIO_A/sample_summary_indA_FINAL/bin_by_bin/
DRAM.py annotate --input_fasta 'all_fas/*.fa' -o $ANVIO_A/sample_summary_indA_FINAL/DRAM --threads 60

### Distilling (summarizing all data into a very nice graph)

DRAM.py distill \
	-i $ANVIO_A/sample_summary_indA_FINAL/DRAM/annotations.tsv \
	-o $ANVIO_A/sample_summary_indA_FINAL/DRAM/genome_summaries \
	--trna_path $ANVIO_A/sample_summary_indA_FINAL/DRAM/trnas.tsv \
	--rrna_path $ANVIO_A/sample_summary_indA_FINAL/DRAM/rrnas.tsv


cd $ANVIO_B/sample_summary_indB_FINAL/bin_by_bin/all_fas
DRAM.py annotate \
	-i '*.fa' \
	-o $ANVIO_B/sample_summary_indB_FINAL/DRAM \
	--threads 100

### Distilling (summarizing all data into a very nice graph)

DRAM.py distill \
	-i $ANVIO_B/sample_summary_indB_FINAL/DRAM/annotations.tsv \
	-o $ANVIO_B/sample_summary_indB_FINAL/DRAM/genome_summaries \
	--trna_path $ANVIO_B/sample_summary_indB_FINAL/DRAM/trnas.tsv \
	--rrna_path $ANVIO_B/sample_summary_indB_FINAL/DRAM/rrnas.tsv

####################################################################
### prepare fasta file for dereplication (dereplication allows the comparison between different genome from different sample sets)
####################################################################

#add prefix for each year (so no same filename occurs)
mkdir $output/anvio_derep/
mkdir $output/anvio_derep/all_fa

cp $ANVIO_A/sample_summary_indA_FINAL/bin_by_bin/all_fas/*.fa $output/anvio_derep/all_fa

for file in $output/anvio_derep/all_fa/*;
do mv $output/anvio_derep/all_fa/$file $output/anvio_derep/all_fa/indA_$file ;
done

cp $ANVIO_B/sample_summary_indB_FINAL/bin_by_bin/all_fas/*.fa $output/anvio_derep/all_fa

for file in $ANVIO_B/sample_summary_indB_FINAL/bin_by_bin/*;
do mv $ANVIO_B/sample_summary_indB_FINAL/bin_by_bin/$file $ANVIO_B/sample_summary_indB_FINAL/bin_by_bin/indB_$file ;
done

#derep (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5702732/)

dRep dereplicate $output/anvio_derep/derep -g $output/anvio_derep/all_fa/*.fa -p 50

#run DRAM on derep'ed fasta

conda activate DRAM

cd $output/anvio_derep/derep/dereplicated_genomes

DRAM.py annotate -i '*.fa' -o $output/anvio_derep/derep/DRAM --threads 100

### Distilling (summarizing all data into a very nice graph)

DRAM.py distill -i $output/anvio_derep/derep/DRAM/annotations.tsv -o $output/anvio_derep/derep/DRAM/genome_summaries --trna_path $output/anvio_derep/derep/DRAM/trnas.tsv --rrna_path $ANVI$output/anvio_derep/derep/DRAM/rrnas.tsv





####################################################################
### gtdb-tk
####################################################################
conda activate gtdb-tk

GTDBTK_DATA_PATH="/isilon/ottawa-rdc/users/shared/chenw_lab/galen/GTDB_release95/"

gtdbtk de_novo_wf --genome_dir $ANVIO_A/sample_summary_indA_FINAL/bin_by_bin/all_fas/ --bacteria -x fa --outgroup_taxon p__Chloroflexi --out_dir $ANVIO_A/sample_summary_indA_FINAL/gtdb_output --cpus 20
gtdbtk ani_rep --genome_dir $ANVIO_A/sample_summary_indA_FINAL/bin_by_bin/all_fas --out_dir $ANVIO_A/sample_summary_indA_FINAL/gtdb_output -x fa --cpus 20
gtdbtk classify_wf --genome_dir $ANVIO_A/sample_summary_indA_FINAL/bin_by_bin/all_fas --out_dir $ANVIO_A/sample_summary_indA_FINAL/gtdb_output --outgroup_taxon p__Chloroflexi -x fa --cpus 20

gtdbtk de_novo_wf --genome_dir $ANVIO_B/sample_summary_indB_FINAL/bin_by_bin/all_fas --bacteria -x fa --outgroup_taxon p__Chlamydiae --out_dir $ANVIO_B/sample_summary_indB_FINAL/gtdb_output --cpus 20
gtdbtk ani_rep --genome_dir $ANVIO_B/sample_summary_indB_FINAL/bin_by_bin/all_fas --out_dir $ANVIO_B/sample_summary_indB_FINAL/gtdb_output -x fa --cpus 20
gtdbtk classify_wf --genome_dir $ANVIO_B/sample_summary_indB_FINAL/bin_by_bin/all_fas --out_dir $ANVIO_B/sample_summary_indB_FINAL/gtdb_output --outgroup_taxon p__Chlamydiae -x fa --cpus 20
