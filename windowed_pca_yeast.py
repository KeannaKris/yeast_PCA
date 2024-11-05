# Dependencies
import os
import numpy as np
import pandas as pd
import allel
import matplotlib.pyplot as plt

# Set output directory
output_dir = "/mnt/c/Projects/PCA/output/"
os.makedirs(output_dir, exist_ok=True)

# Set file paths
vcf = "/home/kjohnwill/yeast_PCA/output_variants.vcf"

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
    for start in range(0, pos[-1], window_step):
        end = start + window_size

        # Get indices of variants within window
        mask = (pos >= start) & (pos < end)
        if np.sum(mask) >= min_variants:
            windowed_genotypes = genotype_matrix[:, mask]

            # Perform PCA in window
            coords, model = allel.pca(windowed_genotypes, n_components=1)

            # Store PC1 and midpoint for window
            pc1_values.extend(coords[:, 0])
            window_mid = (start + end) // 2
            window_midpoints.extend([window_mid] * len(coords[:, 0]))
            chromosomes.extend([chrom[mask][0]] * len(coords[:, 0]))
        else:
            pc1_values.extend([np.nan] * np.sum(mask))
            window_mid = (start + end) // 2
            window_midpoints.extend([window_mid] * np.sum(mask))
            chromosomes.extend([chrom[mask][0]] * np.sum(mask))

    return np.array(window_midpoints), np.array(pc1_values), np.array(chromosomes)

# Visualize results across genome
def plot_windowed_pca(window_midpoints, pc1_values, chromosomes):
    print("Plotting PC1 values across genome...")

    # Plot PC1 values for each chromosome
    unique_chromosomes = np.unique(chromosomes)
    for chrom in unique_chromosomes:
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
