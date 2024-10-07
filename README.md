
# microorganism SNPs


## Description

generate SNP info from microorganisms with multiple strains sequences

## Installation

```
git clone https://github.com/foreverycc/microorganism_snps.git
```

## Usage

```
usage: microorganism_snps.py [-h] [--wkdir WKDIR] [--outputbase OUTPUTBASE] [--refGenome REFGENOME]
                             [--inputFasta INPUTFASTA] [--date DATE] [--minFragSize MINFRAGSIZE]
                             [--minDepth MINDEPTH] [--minFreq MINFREQ] [--segmentSize SEGMENTSIZE]

get SNP info for bacteria or viruses with multiple strains.

optional arguments:
  -h, --help            show this help message and exit
  --wkdir WKDIR         working directory
  --outputbase OUTPUTBASE
                        output name base.
  --refGenome REFGENOME
                        reference fasta sequence for non-human genome build (default=/dev/null).
  --inputFasta INPUTFASTA
                        input fasta sequence for various strains sequences (default=/dev/null).
  --date DATE           date, MUST in the format of '20230408'
  --minFragSize MINFRAGSIZE
                        [optional] minimum fasta sequence size (default = 50)
  --minDepth MINDEPTH   [optional] minimum depth required for target searching (default = 50)
  --minFreq MINFREQ     [optional] minimum frequency required for target searching (default = 0)
  --segmentSize SEGMENTSIZE
                        [optional] segment size for long sequences (default = 10000) 
```

## Credits


This package was created with Cookiecutter and the `audreyr/cookiecutter-pypackage` project template.

Cookiecutter: https://github.com/audreyr/cookiecutter

`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
