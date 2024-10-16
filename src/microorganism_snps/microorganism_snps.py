import os
import argparse
from parameters import Param
from data_processing import filter_fasta, split_fasta, split_ref_fasta
from snp_analysis import index_ref_genome, alignment, vcf_call, vcf_mod
from utils import conf_logger
from Bio import SeqIO  # Add this import at the top of your file

def main(args=None):
    '''
    Main function to orchestrate the analysis pipeline.
    
    Parameters:
    - args: Command line arguments parsed from the user
    
    Returns:
    - None
    '''
    if args is None:
        # If no args provided, parse them from command line
        parser = argparse.ArgumentParser(description="Microorganism SNPs analysis")
        parser.add_argument("--inputFasta", help = "input fasta sequence for various strains sequences (default=/dev/null).", default = "/dev/null")
        parser.add_argument("--refGenome", help = "reference fasta sequence ID.", default = None)
        parser.add_argument("--wkdir", help = "working directory", default = "./")
        parser.add_argument("--outputbase", help = "output name base.", default = "output")
        parser.add_argument("--date", help = "date, MUST in the format of '20241008'", default = "20241008")
        parser.add_argument("--minFragSize", help = "[optional] minimum fasta sequence size (default = 50)", default = 50)
        parser.add_argument("--minDepth", help = "[optional] minimum depth required for variants searching (default = 10)", default = 10)
        parser.add_argument("--minFreq", help = "[optional] min alternative allele frequency (default = 0.05, range from 0 to 1)", default = 0.05)
        parser.add_argument("--segmentSize", help = "[optional] segment size for long sequences (default = 10000)", default = 10000)
        args = parser.parse_args()

    # Initialize parameters
    basic_param = Param(wkdir=args.wkdir, 
                        outputBase=args.outputbase, 
                        refGenome=args.refGenome,
                        inputFasta=args.inputFasta,
                        date=args.date,
                        minFragSize=args.minFragSize,
                        minDepth=args.minDepth,
                        minFreq=args.minFreq,
                        segmentSize=args.segmentSize)

    basic_param.add_absPath()
    
    # Configure logger
    logger = conf_logger(basic_param.logFile)
    logger.info("Basic settings.")
    logger.info(basic_param)

    # Step 1: Index reference genome
    logger.info("Index reference genome.")
    os.chdir(basic_param.wkdir)
    seq_dict = filter_fasta(basic_param.inputFasta, basic_param.minFragSize)
    split_ref_fasta(seq_dict, basic_param.refGenome)
    index_ref_genome(f"{basic_param.refGenome}.fa")

    # Step 2: Split fasta
    logger.info("Filter and split FASTA sequences.")
    split_fasta(seq_dict, basic_param.segmentSize, basic_param.inputFasta_update)

    # Step 3: Alignment
    logger.info("Perform alignment.")
    alignment(f"{basic_param.refGenome}.fa", basic_param.inputFasta_update, basic_param.outputBase)

    # Step 4: VCF calling and modification
    logger.info("Call variants and modify VCF.")
    vcf_call(f"{basic_param.refGenome}.fa", basic_param.outputBase, basic_param.minDepth)
    vcf_mod(basic_param.outputBase, basic_param.minFreq)

if __name__ == "__main__":
    main()
