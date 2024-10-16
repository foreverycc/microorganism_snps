from Bio import SeqIO

def filter_fasta(input_fasta, min_frag_size):
    '''
    Filter the input FASTA sequences based on minimum fragment size.
    
    Parameters:
    - input_fasta: Path to the input FASTA file
    - min_frag_size: Minimum fragment size threshold
    
    Returns:
    - A dictionary containing sequences that meet the size threshold
    '''
    seq_dict = SeqIO.to_dict(SeqIO.parse(input_fasta, format="fasta"))
    filtered_seq_dict = {seq: seq_dict[seq] for seq in seq_dict if len(seq_dict[seq]) > min_frag_size}
    return filtered_seq_dict

def split_ref_fasta(seq_dict, ref_genome):
    '''
    separate the reference genome from the input fasta file
    '''
    if ref_genome in seq_dict:
        # Write the reference genome sequence to a new file
        SeqIO.write([seq_dict[ref_genome]], f"{ref_genome}.fa", "fasta")
    else:
        raise ValueError(f"Reference genome {ref_genome} not found in the sequence dictionary.")

def split_fasta(seq_dict, segment_size, input_fasta_update):
    '''
    Split long sequences into smaller fragments based on the segment size.
    
    Parameters:
    - seq_dict: Dictionary of sequences to split
    - segment_size: Maximum size for each segment
    - input_fasta_update: Path to save the updated FASTA file
    
    Returns:
    - None
    '''
    seq_dict_new = {}
    for seq in seq_dict:
        length = len(seq_dict[seq])
        if length < segment_size:
            seq_dict_new[seq] = seq_dict[seq]
        else:
            n = 1
            while (length / n > segment_size):
                n += 1
            for i in range(n):
                if i != (n - 1):
                    seq_dict_new[seq + "_" + str(i)] = seq_dict[seq][int(length * i / n): int(length * (i + 1) / n)]
                else:
                    seq_dict_new[seq + "_" + str(i)] = seq_dict[seq][int(length * i / n):]
    
    # Write the updated sequences to a new FASTA file
    SeqIO.write(seq_dict_new.values(), input_fasta_update, format="fasta")
