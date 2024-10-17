import argparse
import sys
from .main import main as run_analysis

def parse_arguments():
    parser = argparse.ArgumentParser(description="Microorganism SNPs analysis")
    parser.add_argument("--inputFasta", help="input fasta sequence for various strains sequences (default =/dev/null).", default="/dev/null")
    parser.add_argument("--refGenome", help="reference fasta sequence ID.", default=None)
    parser.add_argument("--wkdir", help="working directory, default = ./", default="./")
    parser.add_argument("--outputbase", help="output name base, default=output.", default="output")
    parser.add_argument("--date", help="date, MUST in the format of '20241008'", default=None)
    parser.add_argument("--minFragSize", help="[optional] minimum fasta sequence size (default = 50)", type=int, default=50)
    parser.add_argument("--minDepth", help="[optional] minimum depth required for variants searching (default = 10)", type=int, default=10)
    parser.add_argument("--minFreq", help="[optional] min alternative allele frequency (default = 0.05, range from 0 to 1)", type=float, default=0.05)
    parser.add_argument("--segmentSize", help="[optional] segment size for long sequences (default = 10000)", type=int, default=10000)
    return parser.parse_args()

def run_app():
    """Entry point for the application."""
    args = parse_arguments()
    try:
        run_analysis(args)
    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    run_app()