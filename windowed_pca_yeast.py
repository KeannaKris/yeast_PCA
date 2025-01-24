# Dependencies
import os
import numpy as np
import pandas as pd
import allel
import matplotlib.pyplot as plt

# Set output directory
output_dir = "/home/kjohnwill/yeast_PCA"
os.makedirs(output_dir, exist_ok=True)

# Set file paths
vcf = "/home/kjohnwill/yeast_PCA/output02.vcf"

# Add validation function
def validate_vcf(vcf_path):
    try:
        # Try reading with allel
        callset = allel.read_vcf(vcf_path, fields=['variants/CHROM', 'variants/POS', 'calldata/GT'])
        print("VCF successfully loaded")
        
        # Check if required fields are present
        if 'calldata/GT' not in callset:
            print("Warning: No genotype (GT) data found in VCF")
        
        # Print first few variants
        if 'variants/CHROM' in callset and 'variants/POS' in callset:
            for i in range(min(5, len(callset['variants/CHROM']))):
                print(f"Variant {i+1}: {callset['variants/CHROM'][i]} - {callset['variants/POS'][i]}")
                
        # Print available fields
        print("Available keys:", list(callset.keys()))
        if 'variants/CHROM' in callset:
            print("Number of variants:", len(callset['variants/CHROM']))
            print("First few chromosomes:", callset['variants/CHROM'][:5])
        if 'calldata/GT' in callset:
            print("Genotype array shape:", callset['calldata/GT'].shape)
    except Exception as e:
        print(f"Error reading VCF: {str(e)}")

# Windowed PCA function and load VCF file
def windowed_PCA(vcf, window_size=1000000, window_step=10000, min_variants=10):
    # Load VCF file
    print("Loading VCF file...")
    callset = allel.read_vcf(vcf, fields=['variants/CHROM', 'variants/POS', 
                                         'calldata/GT', 'variants/SVTYPE'])
    
    print("VCF Contents:")
    print("Available keys:", list(callset.keys()))
    if 'variants/CHROM' in callset:
        print("Number of variants:", len(callset['variants/CHROM']))
        print("First few chromosomes:", callset['variants/CHROM'][:5])
    if 'calldata/GT' in callset:
        print("Genotype array shape:", callset['calldata/GT'].shape)
    else:
        print("No GT field found in VCF")
    
    # Extract chromosome and position data
    chrom = np.array([c.split('#')[-1] for c in callset['variants/CHROM']])  # Modified to handle chromosome names
    pos = np.array(callset['variants/POS'])
    
    # Extract genotype data
    print("Extracting genotype data...")
    genotypes = allel.GenotypeArray(callset['calldata/GT'])
    genotype_matrix = genotypes.to_n_alt().T  # Convert genotype data into allele counts
    
    # Initialize PC1 results list and window midpoints
    pc1_values = []
    window_midpoints = []
    chromosomes = []
    
    # Perform windowed PCA
    print("Performing windowed PCA...")
    for start in range(0, np.max(pos), window_step):
        end = start + window_size
        # Get indices of variants within window
        mask = (pos >= start) & (pos < end)
        if np.sum(mask) >= min_variants:
            windowed_genotypes = genotype_matrix[:, mask]
            # Replace NaN or Inf values with zeros
            windowed_genotypes = np.nan_to_num(windowed_genotypes, nan=0.0, posinf=0.0, neginf=0.0)
            # Check for invalid values
            if np.any(np.isnan(windowed_genotypes)) or np.any(np.isinf(windowed_genotypes)):
                print(f"Skipping window {start}-{end} due to invalid data.")
                continue
            # Perform PCA
            try:
                coords, model = allel.pca(windowed_genotypes, n_components=1)
                pc1_values.extend(coords[:, 0])
                window_mid = (start + end) // 2
                window_midpoints.extend([window_mid] * len(coords[:, 0]))
                chromosomes.extend([chrom[mask][0]] * len(coords[:, 0]))
            except ValueError as e:
                print(f"Error performing PCA on window {start}-{end}: {e}")
                continue
        else:
            print(f"Skipping window {start}-{end}: Insufficient variants.")
            pc1_values.extend([np.nan] * np.sum(mask))
            window_mid = (start + end) // 2
            window_midpoints.extend([window_mid] * np.sum(mask))
            chromosomes.extend(['NA'] * np.sum(mask))
    
    return np.array(window_midpoints), np.array(pc1_values), np.array(chromosomes)

# Visualize results across genome
def plot_windowed_pca(window_midpoints, pc1_values, chromosomes):
    print("Plotting PC1 values across genome...")
    # Plot PC1 values for each chromosome
    unique_chromosomes = np.unique(chromosomes)
    for chrom in unique_chromosomes:
        if chrom == 'NA':
            continue
        mask = chromosomes == chrom
        plt.plot(window_midpoints[mask], pc1_values[mask], label=f'Chr {chrom}', marker='o', markersize=3)
    plt.xlabel('Genomic position (Mb)')
    plt.ylabel('PC1')
    plt.title('Windowed PCA: PC1 across genomic position')
    plt.legend(loc='best')
    plt.grid(True)
    
    # Save plot
    plt.savefig(os.path.join(output_dir, 'windowed_pca_plot.png'))
    plt.close()

# Main workflow
if __name__ == "__main__":
    # Run validation first
    print("Validating VCF file...")
    validate_vcf(vcf)
    
    print("\nStarting PCA analysis...")
    window_midpoints, pc1_values, chromosomes = windowed_PCA(vcf, window_size=1000000, window_step=10000, min_variants=10)
    plot_windowed_pca(window_midpoints, pc1_values, chromosomes)
