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

def parse_vcf_line(line):
    """Parse a single VCF data line"""
    if line.startswith('#'):
        return None
    parts = line.strip().split('\t')
    if len(parts) < 10:  # Need at least 10 columns
        return None
    
    return {
        'chrom': parts[0],
        'pos': int(parts[1]),
        'ref': parts[3],
        'alt': parts[4],
        'info': parts[7],
        'format': parts[8],
        'sample': parts[9]
    }

def read_vcf(vcf_path):
    """Read VCF file line by line"""
    variants = []
    with open(vcf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            variant = parse_vcf_line(line)
            if variant:
                variants.append(variant)
    return variants

def extract_sample_name(sample_field):
    """Extract sample name from the sample field"""
    try:
        # Format is like "1|1:DBVPG6044#1#chrI@15833@15834@P"
        genotype, info = sample_field.split(':')
        sample_name = info.split('#')[0]
        return sample_name if genotype == "1|1" else None
    except:
        return None

def process_chromosome(variants, chrom):
    """Process variants for a single chromosome"""
    # Filter variants for this chromosome
    chrom_variants = [v for v in variants if v['chrom'] == chrom]
    
    # Define our samples
    samples = ['DBVPG6044', 'DBVPG6765', 'S288C', 'SK1', 
              'UWOPS034614', 'Y12', 'YPS128']
    sample_dict = {s: i for i, s in enumerate(samples)}
    
    # Create position-based variant matrix
    pos_dict = defaultdict(lambda: np.zeros(len(samples)))
    
    # Process each variant
    for variant in chrom_variants:
        pos = variant['pos']
        sample_name = extract_sample_name(variant['sample'])
        if sample_name in sample_dict:
            pos_dict[pos][sample_dict[sample_name]] = 1
    
    # Convert to matrix
    positions = sorted(pos_dict.keys())
    if not positions:
        return None, None
        
    matrix = np.array([pos_dict[p] for p in positions])
    return matrix, positions

def windowed_pca(vcf_path, window_size=10000, window_step=2000, min_variants=2):
    """Perform windowed PCA analysis"""
    print("Reading VCF file...")
    variants = read_vcf(vcf_path)
    print(f"Total variants read: {len(variants)}")
    
    # Get unique chromosomes
    chroms = sorted(set(v['chrom'] for v in variants))
    print(f"Found chromosomes: {chroms}")
    
    # Results storage
    all_pc1 = []
    all_positions = []
    all_chroms = []
    
    # Process each chromosome
    for chrom in chroms:
        print(f"\nProcessing {chrom}")
        
        # Get variant matrix for this chromosome
        geno_matrix, positions = process_chromosome(variants, chrom)
        if geno_matrix is None:
            continue
            
        print(f"Generated matrix of shape {geno_matrix.shape}")
        
        # Process windows
        for start in range(min(positions), max(positions), window_step):
            end = start + window_size
            window_mask = (np.array(positions) >= start) & (np.array(positions) < end)
            
            if sum(window_mask) >= min_variants:
                window_geno = geno_matrix[window_mask]
                
                # Skip if no variation
                if np.all(window_geno == window_geno[0]):
                    continue
                    
                try:
                    # Normalize
                    geno_mean = np.mean(window_geno, axis=0)
                    geno_std = np.std(window_geno, axis=0)
                    geno_std[geno_std == 0] = 1
                    norm_geno = (window_geno - geno_mean) / geno_std
                    
                    # PCA
                    coords, model = allel.pca(norm_geno, n_components=1)
                    
                    # Store results
                    all_pc1.append(coords[0, 0])
                    all_positions.append((start + end) // 2)
                    all_chroms.append(chrom)
                    
                except Exception as e:
                    print(f"Error in window {start}-{end}: {str(e)}")
                    continue
    
    return np.array(all_positions), np.array(all_pc1), np.array(all_chroms)

def plot_results(positions, pc1_values, chroms):
    """Plot PCA results"""
    if len(pc1_values) == 0:
        print("No results to plot")
        return
        
    plt.figure(figsize=(15, 8))
    
    # Plot each chromosome
    unique_chroms = sorted(set(chroms))
    colors = plt.cm.tab10(np.linspace(0, 1, len(unique_chroms)))
    
    for i, chrom in enumerate(unique_chroms):
        mask = chroms == chrom
        if sum(mask) > 0:
            label = chrom.split('#')[-1]  # Just show chrI instead of SGDref#1#chrI
            plt.scatter(positions[mask]/1000000, pc1_values[mask],
                       label=label, color=colors[i], alpha=0.6, s=20)
    
    plt.xlabel('Genomic Position (Mb)')
    plt.ylabel('PC1')
    plt.title('Windowed PCA across Genome')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'pca_results.png'),
                bbox_inches='tight', dpi=300)
    plt.close()

if __name__ == "__main__":
    print("Starting analysis...")
    try:
        positions, pc1_values, chroms = windowed_pca(
            vcf,
            window_size=10000,
            window_step=2000,
            min_variants=2
        )
        
        if len(pc1_values) > 0:
            plot_results(positions, pc1_values, chroms)
            print("Analysis complete! Check output directory for plot.")
        else:
            print("No PCA results generated. Try adjusting parameters.")
    except Exception as e:
        print(f"Error during analysis: {e}")
