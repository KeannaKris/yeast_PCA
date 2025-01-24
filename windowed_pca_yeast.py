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
    """Process genotypes and extract information from the complex format"""
    # Get the FORMAT field which contains the variant information
    genotypes = np.zeros((np.sum(mask), 7))  # Matrix for 7 samples
    gt_data = callset['calldata/GT'][mask]
    
    # Convert presence of variant (1|1) to 1, absence to 0
    genotypes = (gt_data[:, 0, 0] == 1).astype(int)
    # Expand to 7 samples
    genotypes = np.tile(genotypes.reshape(-1, 1), (1, 7))
    
    return genotypes

# Windowed PCA function and load VCF file
def windowed_PCA(vcf, window_size=10000, window_step=2000, min_variants=2):
    # Load VCF file
    print("Loading VCF file...")
    callset = allel.read_vcf(vcf, fields=['variants/CHROM', 'variants/POS', 
                                         'variants/SVTYPE', 'calldata/GT'])
    
    if callset is None:
        print("Error: Could not load VCF file")
        return None
        
    # Extract chromosome and position data
    chrom = np.array(callset['variants/CHROM'])
    pos = np.array(callset['variants/POS'])
    svtype = np.array(callset.get('variants/SVTYPE', [None] * len(chrom)))
    
    # Print diagnostic information
    print(f"\nTotal variants: {len(pos)}")
    print(f"Genotype shape: {callset['calldata/GT'].shape}")
    print(f"Unique chromosomes: {np.unique(chrom)}")
    
    # Initialize PC1 results list and window midpoints
    pc1_values = []
    window_midpoints = []
    chromosomes = []
    
    # Perform windowed PCA
    print("\nPerforming windowed PCA...")
    
    # Process each chromosome separately
    for current_chrom in np.unique(chrom):
        print(f"Processing chromosome {current_chrom}")
        
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
                        pc1_values.extend(coords[:, 0])
                        window_midpoints.extend([start + window_size//2] * len(coords[:, 0]))
                        chromosomes.extend([current_chrom] * len(coords[:, 0]))
                        
                except Exception as e:
                    print(f"Error processing window {start}-{end}: {str(e)}")
                    continue
            else:
                print(f"Skipping window {start}-{end}: Insufficient variants ({np.sum(window_mask)} found)")
    
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
            # Extract just the chromosome number for cleaner labels
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
        
        if len(pc1_values) > 0:
            plot_windowed_pca(window_midpoints, pc1_values, chromosomes)
            print("\nAnalysis complete! Check the output directory for the plot.")
        else:
            print("\nNo PCA results were generated. Try adjusting the window size or minimum variants parameter.")
            
    except Exception as e:
        print(f"Error during analysis: {str(e)}")
