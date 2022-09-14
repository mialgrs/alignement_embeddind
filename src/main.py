#! /usr/bin/python3
"""Script to make alignment between 2 sequences."""

import configparser
from bin import readfile, score
from bin import alignglobal, alignlocal, alignglocal


if __name__=="__main__":
    config = configparser.ConfigParser()
    config.read('embedding.cfg')
    #initialize embedding file for the 2 proteins
    emb1 = config['paths']['to_data'] + config['files']['prot_comp_emb']
    emb2 = config['paths']['to_data'] + config['files']['prot_int_emb']
    mat_emb1 = readfile.read_emb(emb1)
    mat_emb2 = readfile.read_emb(emb2)

    #compute and save dot prod between the 2 sequences
    file_dot = config['paths']['to_res'] + config['files']['dot']
    mat_dot = score.dot_product(mat_emb1, mat_emb2, file_dot)

    #load fasta into list of str
    prot1 = config['paths']['to_data'] + config['files']['prot_comp_fasta']
    prot2 = config['paths']['to_data'] + config['files']['prot_int_fasta']
    prot_fasta1 = readfile.read_fasta(prot1)
    prot_fasta2 = readfile.read_fasta(prot2)

    gap = int(config['alignment']['gap'])
    file_align = config['paths']['to_res'] + config['files']['align']

    #global alignment
    if config.getboolean('alignment','global'):
        file_result = config['paths']['to_res'] + "global_align.txt"
        mat_align = alignglobal.alignment(mat_dot, gap, file_align)
        glob_1, glob_2 = alignglobal.needleman_wunsch(mat_align, \
            prot_fasta1, prot_fasta2, file_result)

    #local alignment
    if config.getboolean('alignment','local'):
        file_result = config['paths']['to_res'] + "local_align.txt"
        mat_align = alignlocal.alignment(mat_dot, gap, file_align)
        loc_1, loc_2 = alignlocal.smith_waterman(mat_align, \
                prot_fasta1, prot_fasta2, file_result)

    #glocal alignment
    if config.getboolean('alignment','glocal'):
        file_result = config['paths']['to_res'] + "glocal_align.txt"
        mat_align = alignglobal.alignment(mat_dot, gap, file_align)
        gloc_1, gloc_2 = alignglocal.semi_global(mat_align, \
                prot_fasta1, prot_fasta2, file_result)
