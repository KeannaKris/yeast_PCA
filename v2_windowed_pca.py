# Dependencies
import os
import numpy as np
import pandas as pd
import allel
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from collections import defaultdict

# Configuration
MIN_VARIANTS_PER_WINDOW = 10
FLOAT_PRECISION = 4
MEAN_THRESHOLD = 0.6

def process_window_pca(genotypes, scaler='patterson'):
    """Perform PCA on a window of genotypes"""
    try:
        # Convert to alternate allele counts
        ac = genotypes.to_n_alt()
        
        # Skip windows with no variation
        if np.all(ac == ac[0]):
            return None, None, 0
            
        # Perform PCA
        coords, model = allel.pca(ac, n_components=2, scaler=scaler)
        
        return (coords[:, 0], # PC1
                model.explained_variance_ratio_[0] * 100,
                ac.shape[0]) # number of variants
    except Exception as e:
        print(f"Error in PCA: {str(e)}")
        return None, None, 0

def polarize_pc1(pc1_values, threshold=MEAN_THRESHOLD):
    """Polarize PC1 values based on mean"""
    if pc1_values is None or len(pc1_values) == 0:
        return pc1_values
    mean_val = np.mean(pc1_values)
    if abs(mean_val) > threshold:
        return -pc1_values if mean_val > 0 else pc1_values
    return pc1_values

def windowed_PCA(vcf_path, window_size=50000, window_step=10000, min_variants=5):
    """Perform windowed PCA with improved processing"""
    # Load VCF
    print("Loading VCF file...")
    callset = allel.read_vcf(vcf_path, fields=['variants/CHROM', 'variants/POS', 
                                              'calldata/GT'])
    
    if callset is None:
        raise ValueError("Could not load VCF file")
        
    # Extract data
    chrom = np.array(callset['variants/CHROM'])
    pos = np.array(callset['variants/POS'])
    genotypes = allel.GenotypeArray(callset['calldata/GT'])
    
    # Results containers
    results = defaultdict(list)
    
    # Process each chromosome
    for current_chrom in np.unique(chrom):
        print(f"\nProcessing {current_chrom}")
        
        # Get chromosome data
        chrom_mask = chrom == current_chrom
        chrom_pos = pos[chrom_mask]
        chrom_genotypes = genotypes[chrom_mask]
        
        if len(chrom_pos) == 0:
            continue
            
        # Process windows
        for start in range(min(chrom_pos), max(chrom_pos), window_step):
            end = start + window_size
            window_mask = (chrom_pos >= start) & (chrom_pos < end)
            
            if np.sum(window_mask) >= min_variants:
                window_geno = chrom_genotypes[window_mask]
                
                # Perform PCA
                pc1, var_explained, n_variants = process_window_pca(window_geno)
                
                if pc1 is not None:
                    # Polarize PC1
                    pc1 = polarize_pc1(pc1)
                    
                    # Store results
                    results['pc1'].append(pc1)
                    results['positions'].append(start + window_size//2)
                    results['chromosomes'].append(current_chrom)
                    results['var_explained'].append(var_explained)
                    results['n_variants'].append(n_variants)
    
    return results

def plot_windowed_pca(results, output_dir):
    """Create enhanced PCA plot"""
    fig = go.Figure()
    
    # Plot PC1 for each chromosome
    unique_chroms = sorted(set(results['chromosomes']))
    colors = ['red', 'blue', 'green', 'orange', 'purple']  # Add more colors if needed
    
    for i, chrom in enumerate(unique_chroms):
        mask = [c == chrom for c in results['chromosomes']]
        pos = [p/1000000 for p, m in zip(results['positions'], mask) if m]  # Convert to Mb
        pc1_values = [pc1 for pc1, m in zip(results['pc1'], mask) if m]
        
        # Add trace for each chromosome
        fig.add_trace(go.Scatter(
            x=pos,
            y=pc1_values,
            mode='lines',
            name=chrom.split('#')[-1],
            line=dict(color=colors[i % len(colors)]),
        ))
    
    # Update layout
    fig.update_layout(
        title='Windowed PCA across chromosomes',
        xaxis_title='Position (Mb)',
        yaxis_title='PC1',
        showlegend=True,
        template='plotly_white',
        height=600,
        width=1200
    )
    
    # Save plot
    fig.write_html(os.path.join(output_dir, 'windowed_pca_plot.html'))
    fig.write_image(os.path.join(output_dir, 'windowed_pca_plot.png'))

if __name__ == "__main__":
    vcf_path = "/home/kjohnwill/yeast_PCA/reformatted_output.vcf"
    output_dir = "/home/kjohnwill/yeast_PCA"
    
    print("\nStarting PCA analysis...")
    try:
        # Run analysis
        results = windowed_PCA(
            vcf_path,
            window_size=50000,
            window_step=10000,
            min_variants=5
        )
        
        # Create plots
        if results and len(results['pc1']) > 0:
            plot_windowed_pca(results, output_dir)
            print("\nAnalysis complete! Check output directory for plots.")
        else:
            print("\nNo PCA results generated. Try adjusting parameters.")
            
    except Exception as e:
        print(f"Error during analysis: {e}")
