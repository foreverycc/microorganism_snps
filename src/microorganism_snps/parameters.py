import os
from datetime import date

class Param:
    '''
    Class to handle input parameters and set up basic settings for the analysis.
    '''
    def __init__(self, wkdir, outputBase, refGenome, inputFasta, minFragSize, minDepth, minFreq, segmentSize):
        # Initialize all parameters
        self.wkdir = wkdir
        self.outputBase = outputBase
        self.refGenome = refGenome
        self.inputFasta = inputFasta
        self.date = str(date.today())
        self.minFragSize = int(minFragSize)
        self.minDepth = int(minDepth)
        self.minFreq = float(minFreq)
        self.segmentSize = int(segmentSize)
        
    def add_absPath(self):
        '''
        Convert relative paths to absolute paths for working directory, reference genome, and input FASTA files.
        Raise an error if the paths are invalid.
        '''
        if os.path.exists(self.wkdir):
            self.wkdir = os.path.abspath(self.wkdir)
            self.logFile = os.path.join(self.wkdir, "{}.{}.log.txt".format(self.outputBase, self.date))
        else:
            raise ValueError("The wkdir path is incorrect, please enter the right path.")
        
        if os.path.exists(self.inputFasta):
            self.inputFasta = os.path.abspath(self.inputFasta)
            # Prepare updated fasta file paths
            self.inputFasta_update = self.inputFasta.replace(".fasta", ".fa").replace(".fa", ".updated.fasta")
        else:
            raise ValueError("fasta path does not exist!")

    def __str__(self):
        '''
        Override the string representation of the class to print all the attributes.
        '''
        output_str = ""
        for key, value in self.__dict__.items():
            output_str += "{}: {}\n".format(key, value)
        return output_str
