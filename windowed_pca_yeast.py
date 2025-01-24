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

def analyze_variant_distribution(callset):
    """Analyze the distribution of variants across chromosomes and types"""
    chrom = callset['variants/CHROM']
    pos = callset['variants/POS']
    variant_types = defaultdict(int)
    
    # Count variants by type and chromosome
    for i, v in enumerate(callset.get('variants/SVTYPE', [None] * len(chrom))):
        if v is None:
            variant_types[f"{chrom[i]}_SNP"] += 1
        else:
            variant_types[f"{chrom[i]}_{v}"] += 1
    
    # Print summary
    print("\nVariant Distribution:")
    for key, count in sorted(variant_types.items()):
        print(f"{key}: {count}")
    
    return variant_types

def process_genotypes(genotypes, variant_type=None):
    """Process genotypes based on variant type"""
    # Convert to presence/absence for structural variants
    if variant_type in ['INS', 'DEL']:
        return (genotypes.to_n_alt() > 0).astype(int)
    # Use actual allele counts for SNPs
    else:
        return genotypes.to_n_alt()

def windowed_PCA(vcf, window_size=50000, window_step=5000, min_variants=3):
    # Load VCF file with specific fields
    print("Loading VCF file...")
    callset = allel.read_vcf(vcf, fields=['variants/CHROM', 'variants/POS', 
                                         'variants/SVTYPE', 'calldata/GT',
                                         'variants/ID', 'variants/ALT', 'variants/REF'])
    
    # Debug print
    print("\nFirst few variants:")
    for i in range(5):
        print(f"\nVariant {i+1}:")
        print(f"CHROM: {callset['variants/CHROM'][i]}")
        print(f"POS: {callset['variants/POS'][i]}")
        print(f"GT: {callset['calldata/GT'][i]}")
        print(f"SVTYPE: {callset.get('variants/SVTYPE', [None]*len(callset['variants/CHROM']))[i]}")
        print(f"REF: {callset['variants/REF'][i]}")
        print(f"ALT: {callset['variants/ALT'][i]}")

    # Check if we have any data
    if not callset or 'calldata/GT' not in callset:
        print("No genotype data found in VCF")
        return None
        
    # Print raw genotype shape
    print(f"\nRaw genotype shape: {callset['calldata/GT'].shape}")
    
    # Print first few genotypes
    print("\nFirst few genotypes:")
    print(callset['calldata/GT'][:5])
    
    if callset is None:
        raise ValueError("Failed to load VCF file")
    
    # Analyze variant distribution
    variant_dist = analyze_variant_distribution(callset)
    
    # Extract basic data
    chrom = np.array(callset['variants/CHROM'])
    pos = np.array(callset['variants/POS'])
    svtype = np.array(callset.get('variants/SVTYPE', [None] * len(chrom)))
    
    # Print diagnostic information
    print(f"\nData Summary:")
    print(f"Total variants: {len(pos)}")
    print(f"Genotype shape: {callset['calldata/GT'].shape}")
    print(f"Unique chromosomes: {np.unique(chrom)}")
    
    # Initialize results dictionaries for different variant types
    results = {
        'SNP': {'pc1_values': [], 'window_midpoints': [], 'chromosomes': []},
        'SV': {'pc1_values': [], 'window_midpoints': [], 'chromosomes': []}
    }
    
    # Process each chromosome separately
    for current_chrom in np.unique(chrom):
        print(f"\nProcessing chromosome {current_chrom}")
        
        # Process SNPs and SVs separately
        for variant_class in ['SNP', 'SV']:
            if variant_class == 'SNP':
                mask = (chrom == current_chrom) & (svtype == None)
            else:
                mask = (chrom == current_chrom) & (svtype != None)
            
            chrom_pos = pos[mask]
            if len(chrom_pos) == 0:
                continue
                
            # Get genotypes for this subset
            genotypes = allel.GenotypeArray(callset['calldata/GT'][mask])
            genotype_matrix = process_genotypes(genotypes, 
                                              'SV' if variant_class == 'SV' else None)
            
            print(f"{variant_class} variants in {current_chrom}: {len(chrom_pos)}")
            
            # Process windows
            for start in range(min(chrom_pos), max(chrom_pos), window_step):
                end = start + window_size
                window_mask = (chrom_pos >= start) & (chrom_pos < end)
                
                if np.sum(window_mask) >= min_variants:
                    windowed_genotypes = genotype_matrix[window_mask]
                    
                    try:
                        # Normalize the data
                        geno_mean = np.nanmean(windowed_genotypes, axis=0)
                        geno_std = np.nanstd(windowed_genotypes, axis=0)
                        geno_std[geno_std == 0] = 1
                        normalized_genotypes = (windowed_genotypes - geno_mean) / geno_std
                        normalized_genotypes = np.nan_to_num(normalized_genotypes)
                        
                        # Perform PCA
                        coords, model = allel.pca(normalized_genotypes, n_components=1)
                        
                        # Store results
                        results[variant_class]['pc1_values'].append(coords[0, 0])
                        results[variant_class]['window_midpoints'].append(start + window_size//2)
                        results[variant_class]['chromosomes'].append(current_chrom)
                        
                    except Exception as e:
                        print(f"Error processing {variant_class} window {start}-{end}: {str(e)}")
                        continue
    
    return results

def plot_pca_results(results, output_dir):
    """Plot PCA results for SNPs and SVs separately"""
    for variant_class in results:
        if len(results[variant_class]['pc1_values']) == 0:
            print(f"No results to plot for {variant_class}")
            continue
            
        plt.figure(figsize=(15, 8))
        
        # Convert to numpy arrays
        window_midpoints = np.array(results[variant_class]['window_midpoints'])
        pc1_values = np.array(results[variant_class]['pc1_values'])
        chromosomes = np.array(results[variant_class]['chromosomes'])
        
        # Plot for each chromosome
        unique_chromosomes = np.unique(chromosomes)
        colors = plt.cm.tab10(np.linspace(0, 1, len(unique_chromosomes)))
        
        for i, chrom in enumerate(unique_chromosomes):
            mask = chromosomes == chrom
            if np.sum(mask) > 0:
                plt.scatter(window_midpoints[mask]/1000000, pc1_values[mask], 
                           label=chrom.split('#')[-1], color=colors[i], alpha=0.6, s=20)
        
        plt.xlabel('Genomic position (Mb)')
        plt.ylabel('PC1')
        plt.title(f'Windowed PCA: PC1 across genomic position\n({variant_class} variants)')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.grid(True, alpha=0.3)
        
        # Save plot
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'windowed_pca_{variant_class.lower()}.png'), 
                    bbox_inches='tight', dpi=300)
        plt.close()

# Main workflow
if __name__ == "__main__":
    print("\nStarting PCA analysis...")
    try:
        # Run PCA analysis
        results = windowed_PCA(
            vcf, 
            window_size=50000,    # 50kb windows
            window_step=5000,     # 5kb steps
            min_variants=3        # Minimum 3 variants per window
        )
        
        # Plot results
        plot_pca_results(results, output_dir)
        print("\nAnalysis complete! Check the output directory for plots.")
        
    except Exception as e:
        print(f"Error during analysis: {str(e)}")
