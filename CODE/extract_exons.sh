#!/usr/bin/env bash

# Exit on error
set -euo pipefail

# Check arguments
if [[ "$#" -ne 2 ]]; then
    echo "Usage: $0 input.xls output_prefix"
    exit 1
fi

input_xls="$1"
out_prefix="$2"

# Define intermediate and final outputs
csv_file="${out_prefix}.csv"
tidy_csv="${out_prefix}.tidy.csv"
no_introns="${out_prefix}.noIntrons.csv"
exon_edits_bed="${out_prefix}.HyperTRIBE_exonEdits.bed"

# Convert XLS to CSV
echo "Converting $input_xls to CSV..."
in2csv "$input_xls" > "$csv_file"

# Tidy data
echo "Tidying data..."
awk -F',' '{
    n = split($4, a, ";")
    split($5, b, ";")
    split($6, c, ";")
    for (i = 1; i <= n; i++) {
        print $1 "," $2 "," $3 "," a[i] "," b[i] "," c[i]
    }
}' "$csv_file" > "$tidy_csv"

# Filter out intron lines
echo "Filtering out INTRON lines..."
awk -F',' 'NR==1 || $4 != "INTRON"' "$tidy_csv" > "$no_introns"

# Convert to BED format
echo "Generating BED file..."
awk -F',' '{
    split($6, a, "_");
    chrom = a[1];
    pos = a[2];
    if (pos ~ /^[0-9]+$/) {
        print chrom, pos-1, pos, $1, $5;
    }
}' OFS='\t' "$no_introns" > "$exon_edits_bed"

echo "Done. Output files:"
echo "  CSV:             $csv_file"
echo "  Tidy CSV:        $tidy_csv"
echo "  No Introns CSV:  $no_introns"
echo "  BED File:        $exon_edits_bed"