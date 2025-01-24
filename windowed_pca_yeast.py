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

def read_vcf(vcf_path):
    """Read VCF file while preserving header information"""
    # Store header lines
    header_lines = []
    data_lines = []
    
    with open(vcf_path, 'r') as file:
        for line in file:
            if line.startswith('##'):
                header_lines.append(line.strip())
            else:
                data_lines.append(line.strip())
    
    # Get column names from the #CHROM line
    columns = data_lines[0].strip('#').split('\t')
    
    # Create DataFrame from remaining lines
    data = [line.split('\t') for line in data_lines[1:]]
    df = pd.DataFrame(data, columns=columns)
    
    # Convert POS to integer
    df['POS'] = pd.to_numeric(df['POS'])
    
    return df, header_lines

def extract_sample_info(sample_str):
    """Extract sample name from the sample field"""
    try:
        if '|' not in sample_str:
            return None
            
        parts = sample_str.split(':')
        if len(parts) < 2:
            return None
            
        sample_part = parts[1]
        sample_name = sample_part.split('#')[0]
        return sample_name
    except:
        return None

def process_variants_by_position(df, chromosome):
    """Process variants grouping by position and sample"""
    # Filter for chromosome
    chrom_df = df[df['#CHROM'] == chromosome].copy()
    
    # Define sample names and create mapping
    sample_names = ['DBVPG6044', 'DBVPG6765', 'S288C', 'SK1', 
                   'UWOPS034614', 'Y12', 'YPS128']
    sample_to_idx = {name: idx for idx, name in enumerate(sample_names)}
    
    # Group by position
    pos_dict = defaultdict(lambda: np.zeros(len(sample_names)))
    variant_count = 0
    
    # Process each variant
    for _, row in chrom_df.iterrows():
        pos = row['POS']
        sample_info = row['sample']
        
        if '1|1' in sample_info:  # Check if variant is present
            sample_name = extract_sample_info(sample_info)
            if sample_name in sample_to_idx:
                pos_dict[pos][sample_to_idx[sample_name]] = 1
                variant_count += 1
    
    # Convert to matrix
    unique_pos = sorted(pos_dict.keys())
    genotype_matrix = np.array([pos_dict[p] for p in unique_pos])
    
    print(f"\nChromosome {chromosome}:")
    print(f"Total variants processed: {variant_count}")
    print(f"Unique positions: {len(unique_pos)}")
    print(f"Matrix shape: {genotype_matrix.shape}")
    
    return genotype_matrix, unique_pos

def windowed_PCA(vcf_path, window_size=10000, window_step=2000, min_variants=2):
    """Perform windowed PCA analysis on VCF data"""
    # Read VCF file
    print("Loading VCF file...")
    vcf_df, header_lines = read_vcf(vcf_path)
    
    print(f"\nData Summary:")
    print(f"Total variants: {len(vcf_df)}")
    print(f"Chromosomes: {vcf_df['#CHROM'].unique()}")
    
    # Initialize results
    pc1_values = []
    window_midpoints = []
    chromosomes = []
    
    # Process each chromosome
    for current_chrom in sorted(vcf_df['#CHROM'].unique()):
        print(f"\nProcessing {current_chrom}")
        
        # Get variant matrix for this chromosome
        genotype_matrix, unique_pos = process_variants_by_position(vcf_df, current_chrom)
        
        if len(unique_pos) == 0:
            continue
        
        # Process windows
        for start in range(min(unique_pos), max(unique_pos), window_step):
            end = start + window_size
            window_mask = (np.array(unique_pos) >= start) & (np.array(unique_pos) < end)
            
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
                    
                    # Perform PCA
                    coords, model = allel.pca(normalized_genotypes, n_components=1)
                    pc1_values.append(coords[0, 0])
                    window_midpoints.append(start + window_size//2)
                    chromosomes.append(current_chrom)
                    
                except Exception as e:
                    print(f"Error in window {start}-{end}: {str(e)}")
                    continue
    
    return np.array(window_midpoints), np.array(pc1_values), np.array(chromosomes)

def plot_windowed_pca(window_midpoints, pc1_values, chromosomes):
    """Plot PCA results"""
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
            window_size=10000,    # 10kb windows
            window_step=2000,     # 2kb steps
            min_variants=2        # Minimum 2 variants per window
        )
        
        if window_midpoints is not None and len(pc1_values) > 0:
            plot_windowed_pca(window_midpoints, pc1_values, chromosomes)
            print("\nAnalysis complete! Check the output directory for the plot.")
        else:
            print("\nNo PCA results were generated. Try adjusting parameters.")
            
    except Exception as e:
        print(f"Error during analysis: {str(e)}")
