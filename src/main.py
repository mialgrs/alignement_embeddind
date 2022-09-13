#! /usr/bin/python3

import configparser
import numpy as np
import readfile 
import alignglobal
import alignlocal
import alignglocal
import score


if __name__=="__main__":
    config = configparser.ConfigParser()   
    config.read('embedding.cfg')
    #init emb file for 2 prot
    emb1 = config['paths']['to_data'] + config['files']['prot1_emb']
    emb2 = config['paths']['to_data'] + config['files']['prot2_emb']
    mat_emb1 = readfile.read_emb(emb1)
    mat_emb2 = readfile.read_emb(emb2)

    #compute and save dot prod between the 2
    dot_file = config['paths']['to_res'] + config['files']['dot']
    mat_dot = score.dot_product(mat_emb1, mat_emb2, dot_file)

    prot1 = config['paths']['to_data'] + config['files']['prot1_fasta']
    prot2 = config['paths']['to_data'] + config['files']['prot2_fasta']
    prot_fasta1 = readfile.read_fasta(prot1)
    prot_fasta2 = readfile.read_fasta(prot2) 

    #global alignment
    if config.getboolean('alignment','global') == True :
        mat_align = alignglobal.alignment(mat_dot) 
        glob_bif, glob_adk = \
            alignglobal.needleman_wunsch(mat_align, prot_fasta1, prot_fasta2)
    
    #local alignment
    if config.getboolean('alignment','local') == True:
        mat_align = alignlocal.alignment(mat_dot)
        loc_bif, loc_adk = \
            alignlocal.smith_waterman(mat_align, prot_fasta1, prot_fasta2)

    #glocal alignment
    if config.getboolean('alignment','glocal') == True:
        gloc_bif, gloc_adk = \
            alignglocal.semi_global(mat_align, prot_fasta1, prot_fasta2)