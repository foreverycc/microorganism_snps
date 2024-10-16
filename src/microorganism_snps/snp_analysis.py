import subprocess
import os

def index_ref_genome(ref_genome):
    '''
    Index the reference genome using BWA.
    
    Parameters:
    - ref_genome: Path to the reference genome FASTA file
    
    Returns:
    - None
    '''
    cmd = f"bwa index {ref_genome}"
    subprocess.run(cmd, stderr=subprocess.STDOUT, shell=True, executable='/bin/bash')

def alignment(ref_genome, input_fasta_update, output_base):
    '''
    Align sequences from the input FASTA against the reference genome and create a sorted BAM file.
    
    Parameters:
    - ref_genome: Path to the reference genome
    - input_fasta_update: Path to the updated input FASTA file
    - output_base: Base name for the output BAM file
    
    Returns:
    - None
    '''
    cmd = f"bwa mem -t 8 {ref_genome} {input_fasta_update} | samtools sort -@8 -o {output_base}.bam -; samtools index {output_base}.bam"
    subprocess.run(cmd, stderr=subprocess.STDOUT, shell=True, executable='/bin/bash')

def vcf_call(ref_genome, output_base, min_depth):
    """
    Perform VCF calling using bcftools mpileup and process the output.

    Args:
    ref_genome (str): Path to the reference genome file.
    output_base (str): Base name for output files.
    min_depth (int): Minimum depth for considering a variant.
    """
    # Generate raw mpileup
    cmd1 = f"bcftools mpileup -B -f {ref_genome} -d 10000 {output_base}.bam > {output_base}.raw.mpileup"
    subprocess.run(cmd1, stderr=subprocess.STDOUT, shell=True, executable='/bin/bash')

    # Process the raw mpileup
    with open(f"{output_base}.raw.mpileup") as fh_raw, open(f"{output_base}.mpileup", "w") as fh_proc:
        counter = 1
        for line in fh_raw:
            if line.startswith("#"):
                fh_proc.write(line)
                continue
            
            line_lst = line.strip().split(";")[0].split()
            current_position = int(line_lst[1])
            
            # Fill non-coverage regions as SNPs
            while counter < current_position:
                fh_proc.write(f"{line_lst[0]}\t{counter}\t.\tN\t<*>\t0\t.\tDP=1;I16=5,0,0,0,1275,325125,0,0,0,0,0,0,92,1864,0,0;QS=1,0;MQ0F=1\tPL\t0,15,12\n")
                counter += 1

            counter = min(counter, current_position)
            counter += 1

            if "INDEL" in line:
                fh_proc.write(line)
                continue

            depth = int(line_lst[-1].replace("DP=", ""))
            if depth >= min_depth and 'QS=1,0' in line:
                continue

            fh_proc.write(line)

    # Generate VCF from processed mpileup
    cmd2 = f"cat {output_base}.mpileup | bcftools view -Ov - > {output_base}.raw.vcf"
    subprocess.run(cmd2, stderr=subprocess.STDOUT, shell=True, executable='/bin/bash')

def vcf_mod(output_base, min_freq):
    """
    Modify the raw VCF file and generate a processed VCF file.

    Args:
    output_base (str): Base name for input and output files.
    min_freq (float): Minimum frequency for considering a variant.
    """
    def get_freq_lst(string):
        if "QS=" not in string:
            return ["0.5", "0.5"]
        return [x for x in string.split(";") if "QS=" in x][0].replace("QS=", "").split(",")

    with open(f"{output_base}.raw.vcf") as fh_raw, \
         open(f"{output_base}.proc.vcf", "w") as fh_proc, \
         open(f"{output_base}.selVariants.tab", "w") as fh_final:
        
        counter = 1
        for line in fh_raw:
            if line.startswith("#"):
                fh_proc.write(line)
                continue

            snptype = "indel" if "INDEL" in line else "exact"
            line_lst = line.strip().split()
            freq_lst = get_freq_lst(line_lst[7])
            MAF = 1 - max(float(x) for x in freq_lst)

            if 0 < MAF < min_freq:
                continue

            fh_proc.write(line)
            
            end_position = int(line_lst[1]) + len(line_lst[3]) - 1
            line_out = f"rs{counter}\t{line_lst[0]}\t{line_lst[1]}\t{end_position}\t{snptype}\t{line_lst[3]},{line_lst[4]}\t{line_lst[-1]}\t{','.join(freq_lst)}\n"
            fh_final.write(line_out)
            counter += 1

# Usage example:
# vcf_call("reference.fa", "output", 10)
# vcf_mod("output", 0.01)
