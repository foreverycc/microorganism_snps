"""Main module."""
#!/usr/bin/env python

import os
import argparse
import shutil
import glob
import time
import subprocess
import logging
import pandas as pd
from termcolor import colored
from Bio import SeqIO

class param():
    '''
    Get basic settings
    '''
    def __init__(self, wkdir, outputBase, refGenome, inputFasta, date, minFragSize, minDepth, minFreq, segmentSize):

        self.wkdir = wkdir
        self.outputBase = outputBase
        self.refGenome = refGenome
        self.inputFasta = inputFasta
        self.date = date
        self.minFragSize = int(minFragSize)
        self.minDepth = int(minDepth)
        self.minFreq = float(minFreq)
        self.segmentSize = int(segmentSize)
        
    def add_absPath(self):
        if os.path.exists(self.wkdir):
            self.wkdir = os.path.abspath(self.wkdir)
            self.logFile = os.path.join(self.wkdir, "{}.{}.log.txt".format(self.outputBase, self.date))
        else:
            raise ValueError("The wkdir path is incorrect, please enter the right path.")

        if os.path.exists(self.refGenome) and os.path.exists(self.inputFasta):
            self.refGenome = os.path.abspath(self.refGenome)
            self.inputFasta = os.path.abspath(self.inputFasta)
            self.inputFasta_update = self.inputFasta.replace(".fasta", ".fa").replace(".fa", ".updated.fasta")
        else:
            raise ValueError("refFa path does not exist!")

#         self.indexDir = os.path.join(self.wkdir, "reference_genome")
#         os.system("mkdir {}".format(self.indexDir))
        
    def __str__(self):
        output_str = ""
        for key, value in self.__dict__.items():
            output_str += "{}: {}\n".format(key, value)
        return output_str


class microAln(param):
    '''
    perform alignment and vcf calling for microorganism genomes
    '''
    def __init__(self, param_obj):
        if isinstance(param_obj, param):
            self.__dict__ = {key:value for key, value in param_obj.__dict__.items()}
        else:
            raise TypeError("param_obj must be 'param' class")

    def indexRefGenome(self):
        cmd = "bwa index {}".format(self.refGenome)
#         print(cmd)
        subprocess.run(cmd, stderr=subprocess.STDOUT, shell=True, executable='/bin/bash')
    
    def filterFasta(self):
        self.seqDict = SeqIO.to_dict(SeqIO.parse(self.inputFasta, format = "fasta"))
        tmp = {seq:self.seqDict[seq] for seq in self.seqDict if len(self.seqDict[seq]) > self.minFragSize}
        self.seqDict = tmp
        
    def splitFasta(self):
        seqDict_new = {}
        for seq in self.seqDict:
            length = len(self.seqDict[seq])
            if  length < self.segmentSize:
                seqDict_new[seq] = self.seqDict[seq]
            else:
                n = 1
                while (length/n > self.segmentSize):
                    n += 1
                for i in range(n):
                    if i != (n - 1):
                        seqDict_new[seq + "_" + str(i)] = self.seqDict[seq][int(length * i/n) : int(length * (i+1)/n)]
                    elif i == (n-1):
                        seqDict_new[seq + "_" + str(i)] = self.seqDict[seq][int(length * i/n) : ]
        SeqIO.write(seqDict_new.values(), self.inputFasta_update, format = "fasta")
                    
        
    def alignment(self):
        cmd = "bwa mem -t 8 {} {} | samtools sort -@8 -o {}.bam -; samtools index {}.bam".format(self.refGenome, self.inputFasta_update, self.outputBase, self.outputBase)
        subprocess.run(cmd, stderr=subprocess.STDOUT, shell=True, executable='/bin/bash')

    def vcfCall(self):
        '''
        mpileup format
            NC_001664.4	26	.	C	<*>	0	.	DP=5;I16=5,0,0,0,1275,325125,0,0,0,0,0,0,92,1864,0,0;QS=1,0;MQ0F=1	PL	0,15,12
            NC_001664.4	27	.	C	<*>	0	.	DP=5;I16=5,0,0,0,1275,325125,0,0,0,0,0,0,95,1951,0,0;QS=1,0;MQ0F=1	PL	0,15,12
        '''
#         cmd = "bcftools mpileup -f {} {}.bam | bcftools call -mv -Ob -o {}.bcf; bcftools view {}.bcf | vcfutils.pl varFilter - > {}.vcf".format(self.refGenome, self.outputBase, self.outputBase, self.outputBase, self.outputBase)
        # cmd 1
        cmd1 = "bcftools mpileup -B -f {} -d 10000 {}.bam > {}.raw.mpileup".format(self.refGenome, self.outputBase, self.outputBase) 
        
        # generate vcf
        subprocess.run(cmd1, stderr=subprocess.STDOUT, shell=True, executable='/bin/bash')

        with open("{}.raw.mpileup".format(self.outputBase)) as fh_raw:
            with open("{}.mpileup".format(self.outputBase), "w") as fh_proc:
                counter = 1
                for line in fh_raw:
                    if ("#" in line):
                        fh_proc.write(line)
                        continue
                    line_lst = line.strip().split(";")[0].split()
                    # fill non-coverage regions as SNPs
#                     print ("counter: " + str(counter))
#                     print ("int(line_lst[1]): " + (line_lst[1]))
                    while counter < int(line_lst[1]):
                        line_out = "\t".join([line_lst[0], str(counter), ".", "N", "<*>", "0", ".", "DP=1;I16=5,0,0,0,1275,325125,0,0,0,0,0,0,92,1864,0,0;QS=1,0;MQ0F=1	PL	0,15,12"]) + "\n"
#                         print ("  counter: " + str(counter))
                        counter += 1
#                         print("  line_out: " + line_out)
                        fh_proc.write(line_out)

                    counter = min(counter, int(line_lst[1]))
                    counter += 1
                    
                    if ("INDEL" in line):
                        fh_proc.write(line)
                        continue
                        
                    depth = int(line_lst[-1].replace("DP=", ""))
                    if (depth >= self.minDepth) and ('QS=1,0' in line):
                        continue
                    
                    fh_proc.write(line)
        
        # cmd 2
        cmd2 = "cat {}.mpileup | bcftools view -Ov - > {}.raw.vcf".format(self.outputBase, self.outputBase)
        subprocess.run(cmd2, stderr=subprocess.STDOUT, shell=True, executable='/bin/bash')

    def vcfMod(self):
        '''
        Example modified VCF format:
            SNP1	chr1	10583	10583	exact	A,G,	0.144169,0.855831,
            SNP2	chr1	10611	10611	exact	C,G,	0.981175,0.018825,
            SNP3	chr1	13302	13304	indel	C,T,	0.885675,0.114325,
        VCF raw format:
            NC_006273.2	126	.	A	G,C,<*>	0	.	DP=96;I16=85,0,11,0,21675,5.52712e+06,2805,715275,0,0,0,0,1755,41401,256,6286;QS=0.885417,0.09375,0.0208333,0;VDB=0.938987;SGB=-0.676189;RPB=0.391867;MQB=1;BQB=1;MQ0F=1	PL	0,163,7,226,22,7,255,34,13,7
            NC_006273.2	132	.	A	C,<*>	0	.	DP=100;I16=98,0,2,0,24990,6.37245e+06,510,130050,0,0,0,0,2148,50608,13,145;QS=0.98,0.02,0;VDB=0.28;SGB=-0.453602;RPB=0.132653;MQB=1;BQB=1;MQ0F=1	PL	0,255,11,255,17,11
            NC_006273.2	136	.	A	C,<*>	0	.	DP=105;I16=104,0,1,0,26520,6.7626e+06,255,65025,0,0,0,0,2244,52742,25,625;QS=0.990476,0.00952381,0;SGB=-0.379885;RPB=1;MQB=1;BQB=1;MQ0F=1	PL	0,255,12,255,15,12        
        '''
        def getFreqLst(string):
            '''
            DP=100;I16=98,0,2,0,24990,6.37245e+06,510,130050,0,0,0,0,2148,50608,13,145;QS=0.98,0.02,0;VDB=0.28;SGB=-0.453602;RPB=0.132653;MQB=1;BQB=1;MQ0F=1
            '''
            if "QS=" not in string:
                return (["0.5","0.5"])
            str_lst = [x for x in string.split(";") if "QS=" in x][0].replace("QS=", "").split(",")
            return ([x for x in str_lst])
            
        with open("{}.raw.vcf".format(self.outputBase)) as fh_raw:
            with open("{}.proc.vcf".format(self.outputBase), "w") as fh_proc:
                with open("{}.vcf".format(self.outputBase), "w") as fh_final:
                    counter = 1
                    for line in fh_raw:
                        if "#" in line:
                            fh_proc.write(line)
                            continue
                        if "INDEL" in line:
                            snptype = "indel"
                        else:
                            snptype = "exact"

                        line_lst = line.strip().split()
                        freq_lst = getFreqLst(line_lst[7])
                        MAF = 1 - max([float(x) for x in freq_lst])
                        if (MAF > 0) and (MAF < self.minFreq):
                            continue
                        fh_proc.write(line)
                        line_out = "\t".join(["rs"+str(counter), line_lst[0], line_lst[1], str(int(line_lst[1]) + len(line_lst[3]) - 1), snptype, line_lst[3] + "," + line_lst[4], line_lst[-1],  ",".join(freq_lst)]) + "\n"
                        counter += 1
                        fh_final.write(line_out)
        
        
def conf_logger(basic_param = None):
    if isinstance(basic_param, param):
        logger = logging.getLogger()
        logger.setLevel(logging.INFO)
        lfh = logging.FileHandler(basic_param.logFile)
        lfh.setFormatter(logging.Formatter('%(asctime)s:%(levelname)s:%(name)s: %(message)s'))
        logger.addHandler(lfh)
        lsh = logging.StreamHandler()
        lsh.setFormatter(logging.Formatter('%(message)s'))
        logger.addHandler(lsh)

        return logger
    else:
        raise TypeError("basic_param must be 'param' class")


LOGGER = logging.getLogger()

def main(args = None):
    # assign basic parameters
    basic_param = param(wkdir=args.wkdir, 
                        outputBase=args.outputbase, 
                        refGenome=args.refGenome,
                        inputFasta=args.inputFasta,
                        date = args.date,
                        minFragSize = args.minFragSize,
                        minDepth = args.minDepth,
                        minFreq = args.minFreq,
                        segmentSize = args.segmentSize)

    basic_param.add_absPath()

    # set logger
    conf_logger(basic_param)
    LOGGER.info(colored("basic settings.", "red"))
    LOGGER.info(basic_param)

    # step 1: index reference genome
    LOGGER.info (colored("index reference genome.", "red"))
    os.chdir(basic_param.wkdir)
    microAln_obj = microAln(basic_param)
#     print(microAln_obj)
    microAln_obj.indexRefGenome()
    LOGGER.info(colored("index reference genome finished.\n", "red"))

    # step 2: filter and split fasta
    LOGGER.info(colored("filter and split fasta.", "red"))
    os.chdir(basic_param.wkdir)
    microAln_obj.filterFasta()
    microAln_obj.splitFasta()
    LOGGER.info(colored("filter and split fasta finished.\n", "red"))

    # step 3: alignment
    LOGGER.info (colored("perform alignment", "red"))
    microAln_obj.alignment()
    LOGGER.info(colored("Alignment finished.\n", "red"))

    # step 4: vcf calling
    LOGGER.info(colored("call variants.", "red"))
    microAln_obj.vcfCall()
    microAln_obj.vcfMod()
    LOGGER.info(colored("call variants finished.\n", "red"))
    
    
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = "get SNP file for bacteria or viruses with multiple strains.")
    parser.add_argument("--wkdir", help = "working directory")
    parser.add_argument("--outputbase", help = "output name base.")
    parser.add_argument("--refGenome", help = "reference fasta sequence for non-human genome build (default=/dev/null).", default = "/dev/null")
    parser.add_argument("--inputFasta", help = "input fasta sequence for various strains sequences (default=/dev/null).", default = "/dev/null")
    parser.add_argument("--date", help = "date, MUST in the format of '20230408'")
    parser.add_argument("--minFragSize", help = "[optional] minimum fasta sequence size (default = 50)", default = 50)
    parser.add_argument("--minDepth", help = "[optional] minimum depth required for target searching (default = 50)", default = 50)
    parser.add_argument("--minFreq", help = "[optional] minimum frequency required for target searching (default = 0)", default = 0)
    parser.add_argument("--segmentSize", help = "[optional] segment size for long sequences (default = 10000)", default = 10000)
    
    args = parser.parse_args()

    if args.wkdir == None or args.outputbase == None or args.refGenome == None or args.inputFasta == None:
        print ("Please fill in all the parameters.")
        print (os.system("microorganism_snps.py -h"))

    else:
        main(args)
