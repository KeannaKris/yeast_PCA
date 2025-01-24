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

def windowed_PCA(vcf, window_size=1000000, window_step=10000, min_variants=10):
    # Load VCF file
    print("Loading VCF file...")
    callset = allel.read_vcf(vcf)
    
    # Extract chromosome and position data
    chrom = np.array(callset['variants/CHROM'])
    pos = np.array(callset['variants/POS'])
    
    # Extract genotype data
    print("Extracting genotype data...")
    genotypes = allel.GenotypeArray(callset['calldata/GT'])
    
    # Convert to allele counts and handle missing data
    genotype_matrix = genotypes.to_n_alt()
    
    # If we only have one sample, we need to reshape the data
    if len(genotype_matrix.shape) == 1:
        print("Only one sample detected, reshaping data...")
        genotype_matrix = genotype_matrix.reshape(-1, 1)
    
    # Initialize results
    pc1_values = []
    window_midpoints = []
    chromosomes = []
    
    # Get unique chromosomes
    unique_chroms = np.unique(chrom)
    
    print("Performing windowed PCA...")
    # Process each chromosome separately
    for current_chrom in unique_chroms:
        print(f"Processing chromosome {current_chrom}")
        # Get indices for current chromosome
        chrom_mask = chrom == current_chrom
        chrom_pos = pos[chrom_mask]
        chrom_genotypes = genotype_matrix[chrom_mask]
        
        if len(chrom_pos) == 0:
            continue
            
        # Process windows for this chromosome
        for start in range(min(chrom_pos), max(chrom_pos), window_step):
            end = start + window_size
            # Get indices of variants within window
            window_mask = (chrom_pos >= start) & (chrom_pos < end)
            
            if np.sum(window_mask) >= min_variants:
                windowed_genotypes = chrom_genotypes[window_mask]
                
                # Skip windows with no variation
                if np.all(windowed_genotypes == windowed_genotypes[0]):
                    continue
                    
                try:
                    # Normalize the data manually
                    geno_mean = np.nanmean(windowed_genotypes, axis=0)
                    geno_std = np.nanstd(windowed_genotypes, axis=0)
                    geno_std[geno_std == 0] = 1  # Avoid division by zero
                    normalized_genotypes = (windowed_genotypes - geno_mean) / geno_std
                    
                    # Replace any remaining NaN values with 0
                    normalized_genotypes = np.nan_to_num(normalized_genotypes)
                    
                    # Perform PCA
                    try:
                        coords, model = allel.pca(normalized_genotypes, n_components=1)
                        pc1_values.extend(coords[:, 0])
                        window_mid = (start + end) // 2
                        window_midpoints.extend([window_mid] * len(coords[:, 0]))
                        chromosomes.extend([current_chrom] * len(coords[:, 0]))
                    except Exception as e:
                        print(f"Skipping window {start}-{end} due to PCA error: {str(e)}")
                except Exception as e:
                    print(f"Error processing window {start}-{end}: {str(e)}")
                    continue
            else:
                print(f"Skipping window {start}-{end}: Insufficient variants ({np.sum(window_mask)} found, {min_variants} required)")
    
    return np.array(window_midpoints), np.array(pc1_values), np.array(chromosomes)

def plot_windowed_pca(window_midpoints, pc1_values, chromosomes):
    print("Plotting PC1 values across genome...")
    plt.figure(figsize=(15, 8))
    
    # Plot PC1 values for each chromosome
    unique_chromosomes = np.unique(chromosomes)
    colors = plt.cm.tab10(np.linspace(0, 1, len(unique_chromosomes)))
    
    for i, chrom in enumerate(unique_chromosomes):
        mask = chromosomes == chrom
        if np.sum(mask) > 0:  # Only plot if we have data for this chromosome
            plt.scatter(window_midpoints[mask], pc1_values[mask], 
                       label=f'Chr {chrom}', color=colors[i], alpha=0.6, s=20)
    
    plt.xlabel('Genomic position')
    plt.ylabel('PC1')
    plt.title('Windowed PCA: PC1 across genomic position')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    
    # Save plot with tight layout
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'windowed_pca_plot.png'), bbox_inches='tight', dpi=300)
    plt.close()

# Main workflow
if __name__ == "__main__":
    print("\nStarting PCA analysis...")
    # Adjust these parameters if needed
    window_midpoints, pc1_values, chromosomes = windowed_PCA(
        vcf, 
        window_size=100000,  # Reduced window size
        window_step=10000,   # Smaller step size
        min_variants=5       # Reduced minimum variants
    )
    
    if len(pc1_values) > 0:
        plot_windowed_pca(window_midpoints, pc1_values, chromosomes)
        print("Analysis complete! Check the output directory for the plot.")
    else:
        print("No PCA results were generated. Try adjusting the window size or minimum variants parameter.")
