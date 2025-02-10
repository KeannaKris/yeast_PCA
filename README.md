# Windowed PCA Analysis for S. cerevisiae8 yeast
Data (scerevisiae8.fa): https://github.com/waveygang/wfmash/tree/main/data

The reference and samples (7) can be found in the same file.
  
  Reference: SDGref
  
  Samples:
    >DBVPG6044
    >DBVPG6765
    >S288C
    >SK1
    >UWOPS034614
    >Y12
    >YPS128
    
---
## Data preprocessing
#### Seperating reference and samples
Reference
```{shell}
samtools faidx data/scerevisiae8.fa $(grep ^SGDref data/scerevisiae8.fa.fai | cut -f 1) > data/SGDref.fa
```
Samples
```{shell}
samtools faidx data/scerevisiae8.fa $(grep ^SGDref -v data/scerevisiae8.fa.fai | cut -f 1) > data/genome.fa
#check that samples are in file
grep '^>' data/genome.fa | awk -F'#' '{print $1}' | sort -u
```
---
#### Using samtools and minimap2 to perform genome alignment
Perform on each sample
Extract the sample
```{shell}
samtools faidx data/genome.fa $(grep -o "^>S288C#1#[^ ]*" data/genome.fa | tr -d '>') > data/S288C.fa
```
Align the sequence to the reference using minimap2 (fast sequence aligner)
```{shell}
minimap2 -ax asm5 data/SGDref.fa S288C.fa > data/S288C.sam
```
Convert SAM to BAM
```{shell}
samtools view -bS data/S288C.sam > data/S288C.bam
```
Sorts BAM file by chromosomal coordinates
```{shell}
samtools sort data/S288C.bam -o data/S288C.sorted.bam
```
Create index file from sorted BAM
```{shell}
samtools index data/S288C.sorted.bam
```
---
#### Using bcftools to perform variant calling
created a loop to process each chromosome seperately
```{shell}
for chr in SGDref#1#chrI SGDref#1#chrII SGDref#1#chrIII SGDref#1#chrIV SGDref#1#chrV SGDref#1#chrVI SGDref#1#chrVII SGDref#1#chrVIII SGDref#1#chrIX SGDref#1#chrX SGDref#1#chrXI SGDref#1#chrXII SGDref#1#chrXIII SGDref#1#chrXIV SGDref#1#chrXV SGDref#1#chrXVI SGDref#1#chrMT; do
    echo "Processing $chr"
    bcftools mpileup -r "$chr" -Ou -f data/SGDref.fa \
    data/DBVPG6044.sorted.bam \
    data/DBVPG6765.sorted.bam \
    data/S288C.sorted.bam \
    data/SK1.sorted.bam \
    data/UWOPS034614.sorted.bam \
    data/Y12.sorted.bam \
    data/YPS128.sorted.bam | bcftools call --ploidy 1 -mv -Ou > "${chr}.bcf"
done
```
Then concatenate/merge
```{shell}
bcftools concat SGDref#1#chr*_mpileup.bcf -Ob -o merged_mpileup.bcf
rm SGDref#1#chrVI_mpileup.bcf
bcftools concat SGDref#1#chr*_mpileup.bcf -Ob -o merged_mpileup.bcf

```
Convert BCF file into VCF file
```{shell}
bcftools view merged_mpileup.bcf -Ov -o merged_mpileup.vcf
```
---
## Windowed PCA Analysis
Had to seperate mitochondrial DNA analysis due to different characteristics compared to the nuclear chromosomes.
```{shell}
git clone https://github.com/KeannaKris/yeast_PCA/blob/main/final_winpca.py
python3 ginal_winpca.py
```
The outputs:

mt_pca_02.html(kb)

mt_pca_02.png

windowed_pca_genomic_02.html(shows chr and mb) 

PCA_windowed_genomic_02.png

windowed_pca_position_02.html(this one looks horrible, like my original)

PCA_windowed_position_02.png
