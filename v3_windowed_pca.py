# Dependencies
import os
import numpy as np
import pandas as pd
import allel
import plotly.graph_objects as go

# Set output directory
output_dir = "/home/kjohnwill/yeast_PCA"
os.makedirs(output_dir, exist_ok=True)

# Set file paths - update to your filtered variants file
vcf = "/home/kjohnwill/yeast_PCA/filtered_variants.vcf.gz"

def process_window_pca(genotypes, scaler='patterson'):
    """Perform PCA on a window of genotypes"""
    try:
        # Convert to alternate allele counts
        ac = genotypes.to_n_alt()
        
        # Skip windows with no variation
        if np.all(ac == ac[0]):
            return None, 0
            
        # Perform PCA
        coords, model = allel.pca(ac, n_components=2, scaler=scaler)
        
        return (coords[:, 0], # PC1
                ac.shape[0])  # number of variants
    except Exception as e:
        print(f"Error in PCA: {str(e)}")
        return None, 0

def plot_by_position(results, output_dir):
    """Create PCA plot with genomic positions on x-axis"""
    if not results or len(results['pc1']) == 0:
        print("No results to plot")
        return
        
    fig = go.Figure()
    
    # Updated sample names to match your VCF
    sample_names = ['DBVPG6044', 'DBVPG6765', 'S288C', 'SK1', 'UWOPS034614', 'Y12', 'YPS128']
    colors = ['red', 'blue', 'green', 'purple', 'orange', 'brown', 'pink']
    
    # Plot each sample
    for i in range(len(sample_names)):
        pc1_values = []
        positions = []
        
        # Collect data for this sample
        for j, pc1_window in enumerate(results['pc1']):
            if pc1_window is not None:
                pc1_values.append(pc1_window[i])
                positions.append(results['positions'][j] / 1000000)  # Convert to Mb
        
        # Add trace for this sample
        fig.add_trace(go.Scatter(
            x=positions,
            y=pc1_values,
            mode='lines',
            name=sample_names[i],
            line=dict(color=colors[i], width=1)
        ))
    
    # Update layout
    fig.update_layout(
        title='Windowed PCA by Position',
        xaxis_title='Position (Mb)',
        yaxis_title='PC1',
        showlegend=True,
        template='plotly_white',
        height=800,
        width=1500,
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="right",
            x=0.99
        )
    )
    
    # Save plot
    fig.write_html(os.path.join(output_dir, 'windowed_pca_by_position.html'))
    fig.write_image(os.path.join(output_dir, 'windowed_pca_by_position.png'))

def windowed_PCA(vcf, window_size=50000, window_step=10000, min_variants=5):
    """Perform windowed PCA analysis"""
    # Load VCF file with specific fields
    print("Loading VCF file...")
    callset = allel.read_vcf(vcf, fields=['variants/CHROM', 'variants/POS', 
                                         'calldata/GT', 'samples'])
    
    if callset is None:
        print("Error: Could not load VCF file")
        return None
        
    # Extract data
    chrom = np.array(callset['variants/CHROM'])
    pos = np.array(callset['variants/POS'])
    genotypes = allel.GenotypeArray(callset['calldata/GT'])
    samples = callset['samples']
    
    # Print diagnostic information
    print(f"\nData Summary:")
    print(f"Total variants: {len(pos)}")
    print(f"Genotype shape: {genotypes.shape}")
    print(f"Number of samples: {genotypes.shape[1]}")
    print(f"Sample names: {samples}")
    print(f"Unique chromosomes: {np.unique(chrom)}")
    
    # Initialize results
    pc1_values = []
    positions = []
    chromosomes = []
    n_variants = []
    
    # Process each chromosome
    for current_chrom in np.unique(chrom):
        print(f"\nProcessing {current_chrom}")
        
        # Get variants for this chromosome
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
                pc1, n_var = process_window_pca(window_geno)
                
                if pc1 is not None:
                    pc1_values.append(pc1)
                    positions.append(start + window_size//2)
                    chromosomes.append(current_chrom)
                    n_variants.append(n_var)
    
    return {
        'pc1': pc1_values,
        'positions': positions,
        'chromosomes': chromosomes,
        'n_variants': n_variants
    }

# Main workflow
if __name__ == "__main__":
    print("\nStarting PCA analysis...")
    try:
        # Run PCA analysis
        results = windowed_PCA(
            vcf, 
            window_size=50000,    # 50kb windows
            window_step=10000,    # 10kb steps
            min_variants=5        # Minimum 5 variants per window
        )
        
        # Create position-based plot
        if results and len(results['pc1']) > 0:
            plot_by_position(results, output_dir)
            print("\nAnalysis complete! Check output directory for plots:")
            print("1. windowed_pca_by_position.png")
        else:
            print("\nNo PCA results were generated. Try adjusting parameters.")
            
    except Exception as e:
        print(f"Error during analysis: {str(e)}")
