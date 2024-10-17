# Microorganism SNPs

Microorganism SNPs is a Python package for generating Single Nucleotide Polymorphisms (SNPs) in microorganism genomes, including bacteria, DNA viruses and RNA viruses. It provides tools for VCF calling and modification, specifically tailored for microorganism genomic data.


## Basic Workflow

1. Data Input:
	- The process begins with a single FASTA file containing genomic sequences from multiple strains of a microorganism.

2. Reference Genome Selection and Sequence Preparation:
	- One strain is designated as the reference genome.
	- The genomic sequences of all strains are segmented into smaller fragments.

3. Alignment:
	- The fragmented sequences are aligned to the chosen reference genome.

4. Variant Calling:
	- Based on the alignment results, genetic variations (such as SNPs, insertions, and deletions) are identified and cataloged.
	- VCF calling using bcftools mpileup
	- Custom VCF modification and filtering

5. Data output:
	- bam
	- raw vcf file
	- filtered vcf file

## Installation

### Prerequisites

- Anaconda or Miniconda

### Steps

1. Clone the repository:
   ```
   git clone https://github.com/yourusername/microorganism_snps.git
   cd microorganism_snps
   ```

2. Create and activate the conda environment:
   ```
   conda env create -f environment.yml
   conda activate microorganism_snps
   ```

3. Install the package:
   ```
   pip install -e .
   ```

## Usage

Here's a basic example of how to use the microorganism_snps package:

```
usage: microorganism_snps [-h] [--inputFasta INPUTFASTA] [--refGenome REFGENOME] [--wkdir WKDIR] [--outputbase OUTPUTBASE]
                          [--minFragSize MINFRAGSIZE] [--minDepth MINDEPTH] [--minFreq MINFREQ] [--segmentSize SEGMENTSIZE]

Microorganism SNPs analysis

optional arguments:
  -h, --help            show this help message and exit
  --inputFasta INPUTFASTA
                        [required] input fasta sequence for various strains sequences.
  --refGenome REFGENOME
                        [required] reference fasta sequence ID.
  --wkdir WKDIR         [optional] working directory (default='./')
  --outputbase OUTPUTBASE
                        [optional] output name base (default='output').
  --minFragSize MINFRAGSIZE
                        [optional] minimum fasta sequence size (default=50)
  --minDepth MINDEPTH   [optional] minimum depth required for variants searching (default=10)
  --minFreq MINFREQ     [optional] min alternative allele frequency (default=0.05, range from 0 to 1)
  --segmentSize SEGMENTSIZE
                        [optional] segment size for long sequences (default=10000)
```


## Example usage

```
microorganism_snps \
    --inputFasta /path/to/microorganism_snps/tests/cmv_sequences.fasta \
    --refGenome NC_006273.2 \
    --outputbase test_output
```

## Visualization of example output
![SNP Distribution](images/snp_distribution.png)