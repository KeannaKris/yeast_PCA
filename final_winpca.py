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
vcf = "/home/kjohnwill/yeast_PCA/merged_variants.vcf"

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
       title='Windowed PCA by Position (Nuclear Chromosomes)',
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
   fig.write_html(os.path.join(output_dir, 'windowed_pca_position_02.html'))
   fig.write_image(os.path.join(output_dir, 'PCA_windowed_position_02.png'))

def plot_combined_positions(results, output_dir):
   """Create PCA plot showing both chromosome and base pair positions"""
   if not results or len(results['pc1']) == 0:
       print("No results to plot")
       return
       
   fig = go.Figure()
   
   # Sample info
   sample_names = ['DBVPG6044', 'DBVPG6765', 'S288C', 'SK1', 'UWOPS034614', 'Y12', 'YPS128']
   colors = ['red', 'blue', 'green', 'purple', 'orange', 'brown', 'pink']
   
   # Create chromosome mapping - excludes chrMT
   chrom_lengths = {
       'SGDref#1#chrI': 230218,
       'SGDref#1#chrII': 813184,
       'SGDref#1#chrIII': 316620,
       'SGDref#1#chrIV': 1531933,
       'SGDref#1#chrV': 576874,
       'SGDref#1#chrVI': 270161,
       'SGDref#1#chrVII': 1090940,
       'SGDref#1#chrVIII': 562643,
       'SGDref#1#chrIX': 439888,
       'SGDref#1#chrX': 745751,
       'SGDref#1#chrXI': 666816,
       'SGDref#1#chrXII': 1078177,
       'SGDref#1#chrXIII': 924431,
       'SGDref#1#chrXIV': 784333,
       'SGDref#1#chrXV': 1091291,
       'SGDref#1#chrXVI': 948066
   }
   
   # Calculate cumulative positions
   cumulative_length = 0
   chrom_starts = {}
   for chrom in chrom_lengths:
       chrom_starts[chrom] = cumulative_length
       cumulative_length += chrom_lengths[chrom]
   
   # Plot each sample
   for i in range(len(sample_names)):
       pc1_values = []
       abs_positions = []
       chr_positions = []
       
       # Collect data for this sample
       for j, pc1_window in enumerate(results['pc1']):
           if pc1_window is not None:
               pc1_values.append(pc1_window[i])
               chrom = results['chromosomes'][j]
               rel_pos = results['positions'][j]
               abs_pos = chrom_starts[chrom] + rel_pos
               abs_positions.append(abs_pos / 1000000)  # Convert to Mb
               chr_positions.append(rel_pos / 1000000)  # Convert to Mb
       
       # Add traces
       fig.add_trace(go.Scatter(
           x=abs_positions,
           y=pc1_values,
           mode='lines',
           name=f'{sample_names[i]} (Genome)',
           line=dict(color=colors[i], width=1, dash='solid'),
           showlegend=True
       ))
   
   # Add chromosome boundaries and labels
   for chrom in chrom_lengths:
       start_pos = chrom_starts[chrom] / 1000000  # Convert to Mb
       fig.add_vline(x=start_pos, line_dash="dash", line_color="gray", opacity=0.5)
       fig.add_annotation(
           x=start_pos + (chrom_lengths[chrom]/2000000),  # Middle of chromosome
           y=1.05,
           text=chrom.split('#')[-1],
           textangle=45,
           xref="x",
           yref="paper",
           showarrow=False,
           font=dict(size=10)
       )
   
   # Update layout
   fig.update_layout(
       title='Windowed PCA by Genomic Position (Nuclear Chromosomes)',
       xaxis_title='Genomic Position (Mb)',
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
       ),
       margin=dict(t=100)
   )
   
   # Save plot
   fig.write_html(os.path.join(output_dir, 'windowed_pca_genomic_02.html'))
   fig.write_image(os.path.join(output_dir, 'PCA_windowed_genomic_02.png'))

def windowed_PCA(vcf, window_size=50000, window_step=10000, min_variants=5):
   """Perform windowed PCA analysis for nuclear chromosomes"""
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
   
   # Process each chromosome (excluding chrMT)
   for current_chrom in np.unique(chrom):
       if 'chrMT' not in current_chrom:
           print(f"\nProcessing {current_chrom}")
           
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

def mt_PCA(vcf, window_size=1000, window_step=500, min_variants=3):
   """Perform PCA analysis specifically for mitochondrial DNA"""
   print("Loading VCF file for MT analysis...")
   callset = allel.read_vcf(vcf, fields=['variants/CHROM', 'variants/POS', 
                                        'calldata/GT', 'samples'])
   
   if callset is None:
       print("Error: Could not load VCF file")
       return None
       
   # Extract data
   chrom = np.array(callset['variants/CHROM'])
   pos = np.array(callset['variants/POS'])
   genotypes = allel.GenotypeArray(callset['calldata/GT'])
   
   # Filter for only chrMT
   mt_mask = np.char.find(chrom.astype(str), 'chrMT') >= 0
   mt_pos = pos[mt_mask]
   mt_genotypes = genotypes[mt_mask]
   
   # Process windows for MT
   pc1_values = []
   positions = []
   
   if len(mt_pos) > 0:
       for start in range(min(mt_pos), max(mt_pos), window_step):
           end = start + window_size
           window_mask = (mt_pos >= start) & (mt_pos < end)
           
           if np.sum(window_mask) >= min_variants:
               window_geno = mt_genotypes[window_mask]
               pc1, _ = process_window_pca(window_geno)
               
               if pc1 is not None:
                   pc1_values.append(pc1)
                   positions.append(start + window_size//2)
   
   return {
       'pc1': pc1_values,
       'positions': positions
   }

def plot_mt_pca(results, output_dir):
   """Create PCA plot for mitochondrial DNA"""
   if not results or len(results['pc1']) == 0:
       print("No MT results to plot")
       return
       
   fig = go.Figure()
   
   sample_names = ['DBVPG6044', 'DBVPG6765', 'S288C', 'SK1', 'UWOPS034614', 'Y12', 'YPS128']
   colors = ['red', 'blue', 'green', 'purple', 'orange', 'brown', 'pink']
   
   for i in range(len(sample_names)):
       pc1_values = []
       positions = []
       
       for j, pc1_window in enumerate(results['pc1']):
           if pc1_window is not None:
               pc1_values.append(pc1_window[i])
               positions.append(results['positions'][j] / 1000)  # Convert to kb
       
       fig.add_trace(go.Scatter(
           x=positions,
           y=pc1_values,
           mode='lines',
           name=sample_names[i],
           line=dict(color=colors[i], width=1)
       ))
   
   fig.update_layout(
       title='Mitochondrial DNA Windowed PCA',
       xaxis_title='Position (kb)',
       yaxis_title='PC1',
       showlegend=True,
       template='plotly_white',
       height=600,
       width=1200,
       legend=dict(
           yanchor="top",
           y=0.99,
           xanchor="right",
           x=0.99
       )
   )
   
   fig.write_html(os.path.join(output_dir, 'mt_pca_02.html'))
   fig.write_image(os.path.join(output_dir, 'mt_pca_02.png'))

# Main workflow
if __name__ == "__main__":
   print("\nStarting PCA analysis...")
   try:
       # Run nuclear chromosomes PCA
       results = windowed_PCA(
           vcf, 
           window_size=50000,
           window_step=10000,
           min_variants=5
       )
       
       # Create nuclear plots
       if results and len(results['pc1']) > 0:
           plot_by_position(results, output_dir)
           plot_combined_positions(results, output_dir)
           print("\nNuclear analysis complete!")
           
       # Run MT analysis
       print("\nStarting MT analysis...")
       mt_results = mt_PCA(
           vcf,
           window_size=1000,  # Smaller windows for MT
           window_step=500,
           min_variants=3
       )
       
       # Create MT plot
       if mt_results and len(mt_results['pc1']) > 0:
           plot_mt_pca(mt_results, output_dir)
           print("MT analysis complete!")
           
       print("\nCheck output directory for plots:")
       print("1. PCA_windowed_position_02.png")
       print("2. PCA_windowed_genomic_02.png")
       print("3. mt_pca_02.png")
           
   except Exception as e:
       print(f"Error during analysis: {str(e)}")
