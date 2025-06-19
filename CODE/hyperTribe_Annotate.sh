#! /usr/bin/bash

#----- Activate conda environment (if any)
source /optnfs/common/miniconda3/etc/profile.d/conda.sh
conda activate /dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/owen/sharedconda/miniconda/envs/xls2csv

#----- SET ARGUMENTS
HT_RESULTS=$1 # Results from hyperTribe in the "results/" folder (usually has a .xls extension)
HT_CSV=$2 # An output name you want for the long format data (this could be whatever ending in .txt)
PREFIX=$3 # Prefix for output sample names
GTFFILE=$4 # Path to organism GTF file for annotations

#----- OUTPUT DIR WILL ALWAYS BE THE SAME
OUTPUT_DIR="Annotated/"
UTIL_DIR="Utils/"
mkdir -p "$OUTPUT_DIR"
mkdir -p "$UTIL_DIR"

#----- CONVERT RESULTS FROM XLS TO CSV
#in2csv "$HT_RESULTS" > "$HT_CSV"

#----- Might want to add some check here " if file is binary, convert to csv first..."

#--------------------------------------------------------#
# PIVOT THE DATA AND CONVERT TO BED FORMAT
#--------------------------------------------------------#

#----- Tidy the data
awk -F'\t' '{ 
    n = split($4, a, ";")
    split($5, b, ";")
    split($6, c, ";")
    for (i = 1; i <= n; i++) {
        print $1 "," $2 "," $3 "," a[i] "," b[i] "," c[i]
    }
}' "$HT_RESULTS" > "$HT_CSV"

#----- Remove intron lines
awk -F',' 'NR==1 || $4 != "INTRON"' "$HT_CSV" > "${OUTPUT_DIR}""${PREFIX}".noIntrons.csv

#----- Get exon edits in bed file format
awk -F',' '{ 
    split($6, a, "_");
    chrom = a[1];
    pos = a[2];
    if (pos ~ /^[0-9]+$/) {
        print chrom, pos-1, pos, $1, $5;
    }
}' OFS='\t' "${OUTPUT_DIR}""${PREFIX}".noIntrons.csv > "${OUTPUT_DIR}""${PREFIX}".exonEdits.bed

#--------------------------------------------------------#
# ANNOTATE
#--------------------------------------------------------#

#----- How many edits are we looking for
echo -e "Number of expected edits: "
wc -l "${OUTPUT_DIR}""${PREFIX}".exonEdits.bed
echo "\n"

#----- Filter the GTF for CDS and UTRs
awk '$3 == "CDS" || $3 == "three_prime_utr" || $3 == "five_prime_utr"' "$GTFFILE" > ${UTIL_DIR}CDS_5_3.gtf

#----- Convert filtered GTF to BED, extracting gene name
awk 'BEGIN{OFS="\t"}
{
    match($0, /gene_name "([^"]+)"/, arr);
    gene = (arr[1] != "") ? arr[1] : ".";
    print $1, $4-1, $5, gene, ".", $7, ".", $3
}' ${UTIL_DIR}CDS_5_3.gtf > ${UTIL_DIR}filt1.bed

#----- Add "chr" prefix to chromosome if needed
awk 'BEGIN{OFS="\t"} {$1 = "chr"$1; print}' ${UTIL_DIR}filt1.bed > ${UTIL_DIR}CDS_5_3.bed

#----- Clean up
rm ${UTIL_DIR}filt1.bed

awk 'BEGIN{OFS="\t"}
{
    # Extract gene_name from the attributes field
    match($0, /gene_name "([^"]+)"/, arr);
    gene = (arr[1] != "") ? arr[1] : ".";

    # Adjust chromosome format
    chrom = ($1 == "MT") ? "chrM" : "chr"$1;

    # Print fields with gene_name included in 4th column
    print chrom, $4-1, $5, gene, ".", $7, ".", $3
}' "$GTFFILE" > ${UTIL_DIR}fullGTF.bed


bedtools intersect -a "${OUTPUT_DIR}""${PREFIX}".exonEdits.bed -b CDS_5_3.bed -wa -wb | \
awk 'BEGIN{OFS="\t"} {
    if ($4 == $9) {
        query = $1 FS $2 FS $3 FS $4;
        feat = $NF;
        if (!(query in seen)) {
            print $0;
            seen[query] = feat;
        }
    }
}' > ${OUTPUT_DIR}primary_unique_edits.bed



#----- Debugging
echo -e "Number of unique edits belonging to 5'UTR, 3'UTR or CDS: "
wc -l ${OUTPUT_DIR}primary_unique_edits.bed


#----- Figure out what edits are missing from the primary unique edits
awk 'FNR==NR {seen[$4]; next} !($4 in seen)' ${OUTPUT_DIR}primary_unique_edits.bed "${OUTPUT_DIR}""${PREFIX}".exonEdits.bed > ${OUTPUT_DIR}missing_from_primary.bed
echo -e "Number of edits missing from the primary unique edits: "
wc -l ${OUTPUT_DIR}missing_from_primary.bed
echo "\n"

bedtools intersect -a ${OUTPUT_DIR}missing_from_primary.bed -b ${UTIL_DIR}fullGTF.bed -wa -wb | \
awk 'BEGIN{OFS="\t"} {
    if ($4 == $9) {
        query = $1 FS $2 FS $3 FS $4;
        feat = $NF;
        if (!(query in seen)) {
            print $0;
            seen[query] = feat;
        }
    }
}' > ${OUTPUT_DIR}secondary_unique_edits.bed

#----- Debugging
echo -e "Number of edits found by intersecting against full GTF: "
wc -l ${OUTPUT_DIR}secondary_unique_edits.bed


#----- Find what is still missing
awk 'FNR==NR {seen[$4]; next} !($4 in seen)' ${OUTPUT_DIR}secondary_unique_edits.bed "${OUTPUT_DIR}""${PREFIX}".exonEdits.bed > ${OUTPUT_DIR}still_missing.bed
awk 'FNR==NR {seen[$4]; next} !($4 in seen)' ${OUTPUT_DIR}secondary_unique_edits.bed ${OUTPUT_DIR}missing_from_primary.bed 

#----- Concatenate secondary_unique_edits.bed to primary_unique_edits.bed
echo -e "Appending secondary to primary..."
cat  ${OUTPUT_DIR}secondary_unique_edits.bed >>  ${OUTPUT_DIR}primary_unique_edits.bed
echo -e "Number of resulting lines: "
wc -l  ${OUTPUT_DIR}primary_unique_edits.bed

awk 'FNR==NR {seen[$4]; next} !($4 in seen)'  ${OUTPUT_DIR}primary_unique_edits.bed "$INPUT" >  ${OUTPUT_DIR}still_missing.bed


#----- Add header to edits
awk 'BEGIN{OFS="\t"} {print $1, $3, $4, $5, $6, $7, $8, $9, $11, $13}'  ${OUTPUT_DIR}primary_unique_edits.bed >  ${OUTPUT_DIR}hyperTRIBE_annotated_results.txt
(echo -e "Chromosome\tEdit_Site\tGene\tScore\tGTF_Chrom\tGTF_Start\tGTF_End\tStrand\tAnnotation" && cat  ${OUTPUT_DIR}hyperTRIBE_annotated_results.txt) >  ${OUTPUT_DIR}hyperTribe_Results_annotated.tsv

