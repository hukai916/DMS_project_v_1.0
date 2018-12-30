"""
Kai's personal collection of little tools.
"""

from itertools import product
import subprocess
import pysam
import pandas as pd
from toolbox_sam import md2list # convert md label into reference based coordinates
from Bio import SeqIO
import json
from pathlib import Path

def wt_codon(wtfile, pos):
    """
    Retrieve the 3-letter wt codon give the wtfile and target pos.
    """
    wt_file = SeqIO.parse(wtfile.as_posix(), 'fasta')
    wt_seq  = next(wt_file).seq
    return(wt_seq[pos-1:pos+2])

def fastq_splitter_wrapper(raw_fastq_file, seq_id_file, out_fastq_file):
    """
    Split fastq files given seq_id_file and raw_fastq_file
    and expected out_fastq_file.
    """

    command = ' '.join(["seqtk subseq", raw_fastq_file.as_posix(), seq_id_file.as_posix(), ">", out_fastq_file.as_posix()])
    subprocess.run(command, shell=True)

def quality_control_library(wtfile, amplicon_range_list, mut_pos_list):
    """ For each mutation site, make 4*4*2 variants using WT sequence.
        Then split into corresponding amplicons.
        wtfile: WT fasta file
        amplicon_range_list: a list of (amplicon_start, amplicon_end)
        mut_pos_list: a list of mutation coordinates in WT seq f
        The quality score will be arbitrarily set to "J" for all nt.
        The default mutation is NNG (32 mutations per aa site).
    """
    wt_fa = open(wtfile).readlines()[1]
    library_list = []
    for i in range(len(amplicon_range_list)):
        amplicon_range = amplicon_range_list[i]
        mut_pos        = mut_pos_list[i]
        for site in mut_pos:
            seg1   = wt_fa[:site-1]
            seg2   = wt_fa[site+2:]
            lib_insert   = [''.join(x) for x in product('ATCG', 'ATCG', 'GC')]
            for insert in lib_insert:
                lib_seq      = ''.join([seg1,insert,seg2])
                lib_seq_fa   = lib_seq[amplicon_range[0]-1:amplicon_range[1]]
                lib_seq_name = '_'.join(['@Amplicon', str(i+1), str(site), insert])
                library_fq   = '\n'.join([lib_seq_name, lib_seq_fa, "+", "J"*len(lib_seq_fa)])
                library_list.append(library_fq)
    return(library_list)

def parse_bamfile(bamfile, wtfile):
    """
    Parse bamfile into a pd.df object that contains various information for easy access.
    """
    wt_name = wtfile.name
    bamfile = pysam.AlignmentFile(bamfile, 'rb')
    itr = list(bamfile.fetch(wt_name))[:] # Note the ref seq name must match its file name.

    target = [[item.cigartuples,
               item.query_name,
               item.query_sequence,
               item.query_qualities,
               item.get_reference_positions(), # Mapped positions of query on reference.
               item.get_tag("MD")]
             for item in itr]

    df = pd.DataFrame(target, columns=['Cigar',
                                       'Query_name',
                                       'Query_seq',
                                       'Query_quality',
                                       'Ref_pos',
                                       'Query_md'])

    df['Cigar_len'] = df.apply(lambda row:len(row['Cigar']), axis=1)
    df['Cigar_first'] = df.apply(lambda row:row['Cigar'][0][0], axis=1)
    df['Query_name']  = [item.query_name for item in itr]
    df['Mut_pos'] = df.apply(lambda row: md2list(row['Query_md'], row['Ref_pos'][0]), axis=1)

    return(df)
