# Dependencies
import os
import numpy as np
import pandas as pd
import allel
import matplotlib.pyplot as plt
from collections import defaultdict

# Set output directory
output_dir = "/home/kjohnwill/yeast_PCA"
os.makedirs(output_dir, exist_ok=True)

# Set file paths
vcf = "/home/kjohnwill/yeast_PCA/output02.vcf"

def extract_sample_info(format_str):
    """Extract sample name from the complex format string"""
    # Example format: 1|1:DBVPG6044#1#chrI@15833@15834@P
    try:
        # Split on ':' and get the second part
        sample_part = format_str.split(':')[1]
        # Get the sample name (everything before the first '#')
        sample_name = sample_part.split('#')[0]
        return sample_name
    except:
        return None

def process_variants_by_position(callset, mask):
    """Process variants grouping by position and sample"""
    pos = callset['variants/POS'][mask]
    format_strings = callset['FORMAT'][mask]
    
    # Create a dictionary to store variants by position
    pos_dict = defaultdict(lambda: np.zeros(7))
    sample_to_idx = {name: idx for idx, name in enumerate([
        'DBVPG6044', 'DBVPG6765', 'S288C', 'SK1', 
        'UWOPS034614', 'Y12', 'YPS128'
    ])}
    
    # Process each variant
    for i, (p, fmt) in enumerate(zip(pos, format_strings)):
        sample_name = extract_sample_info(fmt)
        if sample_name in sample_to_idx:
            pos_dict[p][sample_to_idx[sample_name]] = 1
    
    # Convert to matrix
    unique_pos = sorted(pos_dict.keys())
    genotype_matrix = np.array([pos_dict[p] for p in unique_pos])
    
    print(f"Generated matrix shape: {genotype_matrix.shape}")
    if len(genotype_matrix) > 0:
        print(f"Sample of first few variants:\n{genotype_matrix[:5]}")
    
    return genotype_matrix, unique_pos

def windowed_PCA(vcf, window_size=10000, window_step=2000, min_variants=2):
    # Load VCF file
    print("Loading VCF file...")
    callset = allel.read_vcf(vcf, fields=['variants/CHROM', 'variants/POS',
                                         'variants/FORMAT'])
    
    if callset is None:
        print("Error: Could not load VCF file")
        return None
        
    # Extract basic data
    chrom = np.array(callset['variants/CHROM'])
    pos = np.array(callset['variants/POS'])
    
    # Print diagnostic information
    print(f"\nData Summary:")
    print(f"Total variants: {len(pos)}")
    
    # Initialize results
    pc1_values = []
    window_midpoints = []
    chromosomes = []
    
    # Process each chromosome separately
    for current_chrom in np.unique(chrom):
        print(f"\nProcessing chromosome {current_chrom}")
        
        # Get variants for this chromosome
        chrom_mask = chrom == current_chrom
        
        # Process variants for this chromosome
        genotype_matrix, unique_pos = process_variants_by_position(callset, chrom_mask)
        
        if len(unique_pos) == 0:
            continue
            
        # Process windows
        for start in range(min(unique_pos), max(unique_pos), window_step):
            end = start + window_size
            window_mask = (unique_pos >= start) & (unique_pos < end)
            
            if np.sum(window_mask) >= min_variants:
                windowed_genotypes = genotype_matrix[window_mask]
                
                # Skip windows with no variation between samples
                if np.all(windowed_genotypes == windowed_genotypes[0]):
                    continue
                
                try:
                    # Normalize the data
                    geno_mean = np.nanmean(windowed_genotypes, axis=0)
                    geno_std = np.nanstd(windowed_genotypes, axis=0)
                    geno_std[geno_std == 0] = 1
                    normalized_genotypes = (windowed_genotypes - geno_mean) / geno_std
                    normalized_genotypes = np.nan_to_num(normalized_genotypes)
                    
                    # Perform PCA if we have variation
                    coords, model = allel.pca(normalized_genotypes, n_components=1)
                    pc1_values.extend([coords[0, 0]])
                    window_midpoints.extend([start + window_size//2])
                    chromosomes.extend([current_chrom])
                    
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
