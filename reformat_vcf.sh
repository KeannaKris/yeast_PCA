#!/bin/bash

# Create new header
cat > vcf_header.txt << 'EOL'
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=SGDref#1#chrI,length=230218>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	DBVPG6044	DBVPG6765	S288C	SK1	UWOPS034614	Y12	YPS128
EOL

# Process the VCF file
awk -F'\t' 'BEGIN {OFS="\t"}
!/^#/ {
    # Replace "." with "PASS" in FILTER column
    $7 = "PASS"
    # Simplify FORMAT field to just GT
    $9 = "GT"
    # Remove QI from genotype fields
    for(i=10; i<=NF; i++) {
        split($i,a,":")
        $i = a[1]
    }
    print
}' ../reformatted_output.fixed.vcf > temp_body.vcf

# Combine header and body
cat vcf_header.txt temp_body.vcf > winpca_input.vcf

# Clean up temporary files
rm vcf_header.txt temp_body.vcf
