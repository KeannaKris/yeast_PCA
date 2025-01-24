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

def process_genotypes(callset, mask):
    """Process genotypes comparing reference to samples"""
    # Get the relevant variants
    gt_data = callset['calldata/GT'][mask]
    n_variants = np.sum(mask)
    
    # Create matrix: variants x 7 samples
    genotype_matrix = np.zeros((n_variants, 7))
    
    # If GT is 1|1, it means the sample has the ALT allele (different from reference)
    genotype_matrix = (gt_data[:, 0, 0] == 1).astype(int)[:, np.newaxis]
    # Repeat the column 7 times for each sample
    genotype_matrix = np.repeat(genotype_matrix, 7, axis=1)
    
    print(f"Processed genotype matrix shape: {genotype_matrix.shape}")
    print(f"Sample of genotype values:\n{genotype_matrix[:5, :]}")
    
    return genotype_matrix

def windowed_PCA(vcf, window_size=10000, window_step=2000, min_variants=2):
    # Load VCF file
    print("Loading VCF file...")
    callset = allel.read_vcf(vcf, fields=['variants/CHROM', 'variants/POS', 
                                         'variants/REF', 'variants/ALT',
                                         'calldata/GT'])
    
    if callset is None:
        print("Error: Could not load VCF file")
        return None
        
    # Extract chromosome and position data
    chrom = np.array(callset['variants/CHROM'])
    pos = np.array(callset['variants/POS'])
    
    # Print diagnostic information
    print(f"\nData Summary:")
    print(f"Total variants: {len(pos)}")
    print(f"Genotype shape: {callset['calldata/GT'].shape}")
    print(f"Unique chromosomes: {np.unique(chrom)}")
    
    # Initialize results
    pc1_values = []
    window_midpoints = []
    chromosomes = []
    
    # Perform windowed PCA
    print("\nPerforming windowed PCA...")
    
    # Process each chromosome separately
    for current_chrom in np.unique(chrom):
        print(f"\nProcessing chromosome {current_chrom}")
        
        # Get variants for this chromosome
        chrom_mask = chrom == current_chrom
        chrom_pos = pos[chrom_mask]
        
        if len(chrom_pos) == 0:
            continue
            
        # Get genotype matrix for this chromosome
        genotype_matrix = process_genotypes(callset, chrom_mask)
        
        # Process windows
        for start in range(min(chrom_pos), max(chrom_pos), window_step):
            end = start + window_size
            # Get indices of variants within window
            window_mask = (chrom_pos >= start) & (chrom_pos < end)
            
            if np.sum(window_mask) >= min_variants:
                windowed_genotypes = genotype_matrix[window_mask]
                
                # Skip windows with no variation
                if np.all(windowed_genotypes == windowed_genotypes[0]):
                    print(f"Skipping window {start}-{end}: No variation")
                    continue
                
                try:
                    # Normalize the data
                    geno_mean = np.nanmean(windowed_genotypes, axis=0)
                    geno_std = np.nanstd(windowed_genotypes, axis=0)
                    geno_std[geno_std == 0] = 1
                    normalized_genotypes = (windowed_genotypes - geno_mean) / geno_std
                    normalized_genotypes = np.nan_to_num(normalized_genotypes)
                    
                    # Only perform PCA if we have variation
                    if np.any(normalized_genotypes != normalized_genotypes[0]):
                        coords, model = allel.pca(normalized_genotypes, n_components=1)
                        pc1_values.extend([coords[0, 0]])  # Take first PC only
                        window_midpoints.extend([start + window_size//2])
                        chromosomes.extend([current_chrom])
                        print(f"Successfully processed window {start}-{end}")
                        
                except Exception as e:
                    print(f"Error processing window {start}-{end}: {str(e)}")
                    continue
    
    return np.array(window_midpoints), np.array(pc1_values), np.array(chromosomes)

def plot_windowed_pca(window_midpoints, pc1_values, chromosomes):
    if len(pc1_values) == 0:
        print("No data to plot")
        return
        
    print("Plotting PC1 values across genome...")
    plt.figure(figsize=(15, 8))
    
    # Plot PC1 values for each chromosome
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
            window_size=10000,    # 10kb windows
            window_step=2000,     # 2kb steps
            min_variants=2        # Minimum 2 variants per window
        )
        
        if window_midpoints is not None and len(pc1_values) > 0:
            plot_windowed_pca(window_midpoints, pc1_values, chromosomes)
            print("\nAnalysis complete! Check the output directory for the plot.")
        else:
            print("\nNo PCA results were generated. Try adjusting the window size or minimum variants parameter.")
            
    except Exception as e:
        print(f"Error during analysis: {str(e)}")
