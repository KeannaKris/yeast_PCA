import subprocess
import os

# Define output directory
output_dir = "/mnt/c/Projects/PCA/output"
os.makedirs(output_dir, exist_ok=True)

# Set file paths
reference = "/home/kjohnwill/yeast_PCA/reference.fa"
fasta = "/mnt/c/Projects/PCA/data/scerevisiae8.fa.gz/scerevisiae8.fa.gz"
aligned_bam = os.path.join(output_dir, "aligned_reads.bam")
sorted_bam = os.path.join(output_dir, "sorted_reads.bam")
vcf = os.path.join(output_dir, "variants.vcf.gz")
output_sam = os.path.join(output_dir, "aligned_reads.sam")

# Index the reference genome
def index_reference(reference):
  print("Indexing reference genome...")
  index_command = f"bwa index {reference}"
  subprocess.run(index_command, shell=True, check=True)
  
# Align sequence with BWA
def align_sequence(reference, fasta, output_sam):               # defining function
  print ("Aligning sequences...")                               # print status message
  bwa_command = f"bwa mem {reference} {fasta} > {output_sam}"   # constructing command
  subprocess.run(bwa_command, shell=True, check=True)           # run command

# Convert SAM to BAM and sort
def convert_sort_bam(sam_file, aligned_bam, sorted_bam):
  print("Converting SAM to BAM and sorting...")

  bam_command = f"samtools view -S -b {sam_file} > {aligned_bam}"
  subprocess.run(bam_command, shell=True, check=True)

  # Sort BAM file
  sort_command = f"samtools sort {aligned_bam} -o {sorted_bam}"
  subprocess.run(sort_command, shell=True, check=True)
  
# Index BAM file
def index_bam(bam_file):
  print("Indexing BAM file...")
  index_command = f"samtools index {bam_file}"
  subprocess.run(index_command, shell=True, check=True)

# Call variants with BCFtools
def call_variants(reference, sorted_bam, vcf):
  print("Calling variants...")
  mpileup_command = f"bcftools mpileup -f {reference} {sorted_bam} | bcftools call -mv -Oz -o {vcf}"
  subprocess.run(mpileup_command, shell=True, check=True)

# Main workflow
if __name__== "__main__":
  reference = "/mnt/c/Projects/PCA/data/reference.fa"
  fasta = "/mnt/c/Projects/PCA/data/scerevisiae8.fa.gz/scerevisiae8.fa"

  # Align sequences
  align_sequence(reference, fasta, output_sam)

  # Convert SAM to BAM and sort
  convert_sort_bam(output_sam, aligned_bam, sorted_bam)

  
  # Index BAM file
  index_bam(sorted_bam)

  # Call variants creatinf VCF file
  call_variants(reference, sorted_bam, vcf)
  
