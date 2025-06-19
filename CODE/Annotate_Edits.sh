#!/usr/bin/bash

#!/bin/bash

#SBATCH --job-name=annoEdits
#SBATCH --nodes=1
#SBATCH --partition=preempt1
#SBATCH --account=dac
#SBATCH --time=60:00:00
#SBATCH --mail-user=f007qps@dartmouth.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=annoEdits%j.out

#----- START
echo -e "#-----------------------------------------------------------#\n"
echo "Starting job: $SLURM_JOB_NAME (Job ID: $SLURM_JOB_ID)"
echo "Running on node: $(hostname)"
echo "Start time: $(date)"
echo -e "#-----------------------------------------------------------#\n"

#----- Activate conda environment (if any)
source /optnfs/common/miniconda3/etc/profile.d/conda.sh
conda activate /dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/owen/sharedconda/miniconda/envs/xls2csv


#----- Input = the exon edit sites bed file
GTFFILE = $1
INPUT = $2

#----- How many edits are we looking for
echo -e "Number of expected edits: "
wc -l "$INPUT"
echo "\n"


#----- Filter the GTF for CDS and UTRs
awk '$3 == "CDS" || $3 == "three_prime_utr" || $3 == "five_prime_utr"' "$GTFFILE" > CDS_5_3.gtf

#----- Convert filtered GTF to BED, extracting gene name
awk 'BEGIN{OFS="\t"}
{
    match($0, /gene_name "([^"]+)"/, arr);
    gene = (arr[1] != "") ? arr[1] : ".";
    print $1, $4-1, $5, gene, ".", $7, ".", $3
}' CDS_5_3.gtf > filt1.bed

#----- Add "chr" prefix to chromosome if needed
awk 'BEGIN{OFS="\t"} {$1 = "chr"$1; print}' filt1.bed > CDS_5_3.bed

#----- Clean up
rm filt1.bed

awk 'BEGIN{OFS="\t"}
{
    # Extract gene_name from the attributes field
    match($0, /gene_name "([^"]+)"/, arr);
    gene = (arr[1] != "") ? arr[1] : ".";

    # Adjust chromosome format
    chrom = ($1 == "MT") ? "chrM" : "chr"$1;

    # Print fields with gene_name included in 4th column
    print chrom, $4-1, $5, gene, ".", $7, ".", $3
}' "$GTFFILE" > fullGTF.bed


bedtools intersect -a "$INPUT" -b CDS_5_3.bed -wa -wb | \
awk 'BEGIN{OFS="\t"} {
    if ($4 == $9) {
        query = $1 FS $2 FS $3 FS $4;
        feat = $NF;
        if (!(query in seen)) {
            print $0;
            seen[query] = feat;
        }
    }
}' > primary_unique_edits.bed



#----- Debugging
echo -e "Number of unique edits belonging to 5'UTR, 3'UTR or CDS: "
wc -l primary_unique_edits.bed


#----- Figure out what edits are missing from the primary unique edits
awk 'FNR==NR {seen[$4]; next} !($4 in seen)' primary_unique_edits.bed "$INPUT" > missing_from_primary.bed
echo -e "Number of edits missing from the primary unique edits: "
wc -l missing_from_primary.bed
echo "\n"

bedtools intersect -a missing_from_primary.bed -b fullGTF.bed -wa -wb | \
awk 'BEGIN{OFS="\t"} {
    if ($4 == $9) {
        query = $1 FS $2 FS $3 FS $4;
        feat = $NF;
        if (!(query in seen)) {
            print $0;
            seen[query] = feat;
        }
    }
}' > secondary_unique_edits.bed

#----- Debugging
echo -e "Number of edits found by intersecting against full GTF: "
wc -l secondary_unique_edits.bed


#----- Find what is still missing
awk 'FNR==NR {seen[$4]; next} !($4 in seen)' secondary_unique_edits.bed "$INPUT" > still_missing.bed
awk 'FNR==NR {seen[$4]; next} !($4 in seen)' secondary_unique_edits.bed missing_from_primary.bed 

#----- Concatenate secondary_unique_edits.bed to primary_unique_edits.bed
echo -e "Appending secondary to primary..."
cat secondary_unique_edits.bed >> primary_unique_edits.bed
echo -e "Number of resulting lines: "
wc -l primary_unique_edits.bed

awk 'FNR==NR {seen[$4]; next} !($4 in seen)' primary_unique_edits.bed "$INPUT" > still_missing.bed


#----- Add header to edits
awk 'BEGIN{OFS="\t"} {print $1, $3, $4, $5, $6, $7, $8, $9, $11, $13}' primary_unique_edits.bed > hyperTRIBE_annotated_results.txt
(echo -e "Chromosome\tEdit_Site\tGene\tScore\tGTF_Chrom\tGTF_Start\tGTF_End\tStrand\tAnnotation" && cat hyperTRIBE_annotated_results.txt) > hyperTribe_Results_annotated.tsv









