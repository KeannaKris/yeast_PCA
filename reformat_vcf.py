import os

def reformat_vcf(input_vcf, output_vcf):
    """Reformat VCF file to standard format"""
    with open(input_vcf, 'r') as infile, open(output_vcf, 'w') as outfile:
        # Copy header lines
        for line in infile:
            if line.startswith('##'):
                outfile.write(line)
                continue
            
            if line.startswith('#CHROM'):
                # Write the header line with sample names
                header_parts = line.strip().split('\t')
                # Add sample names
                samples = ['DBVPG6044', 'DBVPG6765', 'S288C', 'SK1', 
                         'UWOPS034614', 'Y12', 'YPS128']
                new_header = header_parts[:-1] + samples
                outfile.write('\t'.join(new_header) + '\n')
                continue
            
            # Process variant lines
            parts = line.strip().split('\t')
            if len(parts) < 10:  # Skip malformed lines
                continue
                
            # Extract sample info
            sample_info = parts[9]  # This is your sample field
            genotype, sample_data = sample_info.split(':')
            sample_name = sample_data.split('#')[0]
            
            # Create genotype fields for all samples
            sample_fields = []
            for sample in samples:
                if sample == sample_name:
                    sample_fields.append(genotype)
                else:
                    sample_fields.append('0|0')  # Reference genotype for other samples
            
            # Write reformatted line
            new_line = parts[:9] + sample_fields
            outfile.write('\t'.join(new_line) + '\n')

# Use the function
input_vcf = "/home/kjohnwill/yeast_PCA/output02.vcf"
output_vcf = "/home/kjohnwill/yeast_PCA/reformatted_output.vcf"

print("Reformatting VCF file...")
reformat_vcf(input_vcf, output_vcf)
print("Done! Reformatted VCF saved to:", output_vcf)
