import os
from .parameters import Param
from .data_processing import filter_fasta, split_fasta, split_ref_fasta
from .snp_analysis import index_ref_genome, alignment, vcf_call, vcf_mod
from .utils import conf_logger

def main(args):
    '''
    Main function to orchestrate the analysis pipeline.
    
    Parameters:
    - args: Command line arguments parsed from the user
    
    Returns:
    - None
    '''
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
    
    ## Check if refGenome is provided
    if basic_param.refGenome is None:
        logger.error("Reference genome is not provided. Please use the --refGenome argument.")
        return
    
    ## Check if the reference genome exists in seq_dict
    if basic_param.refGenome not in seq_dict:
        logger.error(f"Reference genome {basic_param.refGenome} not found in the input sequences.")
        return

    split_ref_fasta(seq_dict, basic_param.refGenome)
    index_ref_genome(f"{basic_param.refGenome}.fa")

    # Step 2: Split fasta
    logger.info("Filter and split FASTA sequences.")
    split_fasta(seq_dict, basic_param.segmentSize, basic_param.inputFasta_update)

    # Step 3: Perform lignment
    logger.info("Perform alignment.")
    alignment(f"{basic_param.refGenome}.fa", basic_param.inputFasta_update, basic_param.outputBase)

    # Step 4: VCF calling and processing
    logger.info("Call variants and modify VCF.")
    vcf_call(f"{basic_param.refGenome}.fa", basic_param.outputBase, basic_param.minDepth)
    vcf_mod(basic_param.outputBase, basic_param.minFreq)