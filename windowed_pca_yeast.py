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
vcf = "/home/kjohnwill/yeast_PCA/variants03.vcf"

# Windowed PCA function and load VCF file
def windowed_PCA(vcf, window_size=1000000, window_step=10000, min_variants=10):
    # Load VCF file
    print("Loading VCF file...")
    callset = allel.read_vcf(vcf)

    # Extract chromosome and position data
    chrom = np.array(callset['variants/CHROM'])
    pos = np.array(callset['variants/POS'])

    # Extract genotype data
    print("Extracting genotype data...")
    genotypes = allel.GenotypeArray(callset['calldata/GT'])
    genotype_matrix = genotypes.to_n_alt().T  # Convert genotype data into allele counts

    # Initialize PC1 results list and window midpoints
    pc1_values = []
    window_midpoints = []
    chromosomes = []

    # Perform windowed PCA
    print("Performing windowed PCA...")
    for start in range(0, np.max(pos), window_step):
        end = start + window_size

        # Get indices of variants within window
        mask = (pos >= start) & (pos < end)
        if np.sum(mask) >= min_variants:
            windowed_genotypes = genotype_matrix[:, mask]

            # Replace NaN or Inf values with zeros
            windowed_genotypes = np.nan_to_num(windowed_genotypes, nan=0.0, posinf=0.0, neginf=0.0)

            # Check for invalid values
            if np.any(np.isnan(windowed_genotypes)) or np.any(np.isinf(windowed_genotypes)):
                print(f"Skipping window {start}-{end} due to invalid data.")
                continue

            # Perform PCA
            try:
                coords, model = allel.pca(windowed_genotypes, n_components=1)
                pc1_values.extend(coords[:, 0])
                window_mid = (start + end) // 2
                window_midpoints.extend([window_mid] * len(coords[:, 0]))
                chromosomes.extend([chrom[mask][0]] * len(coords[:, 0]))
            except ValueError as e:
                print(f"Error performing PCA on window {start}-{end}: {e}")
                continue
        else:
            print(f"Skipping window {start}-{end}: Insufficient variants.")
            pc1_values.extend([np.nan] * np.sum(mask))
            window_mid = (start + end) // 2
            window_midpoints.extend([window_mid] * np.sum(mask))
            chromosomes.extend(['NA'] * np.sum(mask))

    return np.array(window_midpoints), np.array(pc1_values), np.array(chromosomes)

# Visualize results across genome
def plot_windowed_pca(window_midpoints, pc1_values, chromosomes):
    print("Plotting PC1 values across genome...")

    # Plot PC1 values for each chromosome
    unique_chromosomes = np.unique(chromosomes)
    for chrom in unique_chromosomes:
        if chrom == 'NA':
            continue
        mask = chromosomes == chrom
        plt.plot(window_midpoints[mask], pc1_values[mask], label=f'Chr {chrom}', marker='o', markersize=3)
    plt.xlabel('Genomic position (Mb)')
    plt.ylabel('PC1')
    plt.title('Windowed PCA: PC1 across genomic position')
    plt.legend(loc='best')
    plt.grid(True)
    plt.show()

# Main workflow
if __name__ == "__main__":
    window_midpoints, pc1_values, chromosomes = windowed_PCA(vcf, window_size=1000000, window_step=10000, min_variants=10)
    plot_windowed_pca(window_midpoints, pc1_values, chromosomes)
