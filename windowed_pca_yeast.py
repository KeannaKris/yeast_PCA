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
vcf = "/home/kjohnwill/yeast_PCA/reformatted_output.vcf"

def windowed_PCA(vcf, window_size=10000, window_step=2000, min_variants=2):
    # Load VCF file
    print("Loading VCF file...")
    callset = allel.read_vcf(vcf, fields=['variants/CHROM', 'variants/POS', 
                                         'calldata/GT'])
    
    if callset is None:
        print("Error: Could not load VCF file")
        return None
        
    # Extract chromosome and position data
    chrom = np.array(callset['variants/CHROM'])
    pos = np.array(callset['variants/POS'])
    
    # Convert genotypes to allele counts matrix
    genotypes = allel.GenotypeArray(callset['calldata/GT'])
    # Sum the allele counts across both haplotypes
    allele_counts = genotypes.to_n_alt()
    
    # Print diagnostic information
    print(f"\nData Summary:")
    print(f"Total variants: {len(pos)}")
    print(f"Genotype shape: {genotypes.shape}")
    print(f"Number of samples: {genotypes.shape[1]}")
    print(f"Unique chromosomes: {np.unique(chrom)}")
    
    # Initialize results
    pc1_values = []
    window_midpoints = []
    chromosomes = []
    
    # Process each chromosome separately
    for current_chrom in np.unique(chrom):
        print(f"\nProcessing {current_chrom}")
        
        # Get variants for this chromosome
        chrom_mask = chrom == current_chrom
        chrom_pos = pos[chrom_mask]
        chrom_genotypes = allele_counts[chrom_mask]
        
        if len(chrom_pos) == 0:
            continue
        
        # Process windows
        for start in range(min(chrom_pos), max(chrom_pos), window_step):
            end = start + window_size
            window_mask = (chrom_pos >= start) & (chrom_pos < end)
            
            if np.sum(window_mask) >= min_variants:
                windowed_genotypes = chrom_genotypes[window_mask]
                
                # Check for sufficient variation
                if np.all(windowed_genotypes == windowed_genotypes[0]):
                    continue
                    
                # Replace any missing values with 0
                windowed_genotypes = np.nan_to_num(windowed_genotypes, nan=0)
                
                try:
                    # Calculate variance for each sample
                    sample_variances = np.var(windowed_genotypes, axis=0)
                    
                    # Skip if no variation
                    if np.all(sample_variances == 0):
                        continue
                    
                    # Normalize only if there's variation
                    geno_mean = np.mean(windowed_genotypes, axis=0)
                    geno_std = np.std(windowed_genotypes, axis=0)
                    geno_std[geno_std == 0] = 1  # Avoid division by zero
                    
                    normalized_genotypes = (windowed_genotypes - geno_mean) / geno_std
                    
                    # Verify no NaN or inf values
                    if not np.any(np.isnan(normalized_genotypes)) and not np.any(np.isinf(normalized_genotypes)):
                        # Perform PCA
                        coords, model = allel.pca(normalized_genotypes, n_components=1)
                        pc1_values.append(coords[0, 0])
                        window_midpoints.append(start + window_size//2)
                        chromosomes.append(current_chrom)
                        print(f"Successfully processed window {start}-{end}")
                        
                except Exception as e:
                    print(f"Skipping window {start}-{end}: {str(e)}")
                    continue
    
    return np.array(window_midpoints), np.array(pc1_values), np.array(chromosomes)

def plot_windowed_pca(window_midpoints, pc1_values, chromosomes):
    if len(pc1_values) == 0:
        print("No data to plot")
        return
        
    print("Plotting PC1 values across genome...")
    plt.figure(figsize=(15, 8))
    
    # Plot for each chromosome
    unique_chromosomes = np.unique(chromosomes)
    colors = plt.cm.tab10(np.linspace(0, 1, len(unique_chromosomes)))
    
    for i, chrom in enumerate(unique_chromosomes):
        mask = chromosomes == chrom
        if np.sum(mask) > 0:
            chrom_label = chrom.split('#')[-1]
            plt.scatter(window_midpoints[mask]/1000000, pc1_values[mask], 
                       label=chrom_label, color=colors[i], alpha=0.6, s=20)
    
    plt.xlabel('Genomic position (Mb)')
    plt.ylabel('PC1')
    plt.title('Windowed PCA: PC1 across genomic position')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    
    # Save plot
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'windowed_pca_plot.png'), 
                bbox_inches='tight', dpi=300)
    plt.close()

# Main workflow
if __name__ == "__main__":
    print("\nStarting PCA analysis...")
    try:
        window_midpoints, pc1_values, chromosomes = windowed_PCA(
            vcf, 
            window_size=50000,    # Increased window size
            window_step=10000,    # Increased step size
            min_variants=5        # Increased minimum variants
        )
        
        if window_midpoints is not None and len(pc1_values) > 0:
            plot_windowed_pca(window_midpoints, pc1_values, chromosomes)
            print("\nAnalysis complete! Check the output directory for the plot.")
        else:
            print("\nNo PCA results were generated. Try adjusting parameters.")
            
    except Exception as e:
        print(f"Error during analysis: {str(e)}")
