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
ANVIO=$output/anvio
ANVIO_A=$ANVIO/indA
ANVIO_B=$ANVIO/indB


mkdir $ANVIO_B/interpro_output

for sample in `awk '{print $1}' $ANVIO_B/iprs/iprs_list`;
do
/isilon/ottawa-rdc/users/shared/chenw_lab/galen/interproscan-5.39-77.0/interproscan.sh -i $ANVIO_B/iprs/${sample} -f tsv \
	-d $ANVIO_B/interpro_output/ \
	--tempdir $WORK/temp/ \
	--disable-precalc \
	-appl Pfam,PIRSF,SUPERFAMILY,TIGRFAM \
	--iprlookup \
	--goterms \
	--pathways --cpu 10
done

