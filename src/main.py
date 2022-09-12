#! /usr/bin/python3

import numpy as np
import readfile 
import alignglobal
import alignlocal
import score


if __name__=="__main__":
    emb1 = "../data/6PF2K_1bif.t5emb"

    emb2 = "../data/adk_2ak3a.t5emb"
    mat_emb1 = readfile.read_emb(emb1)
    mat_emb2 = readfile.read_emb(emb2)

    dot_file = "../results/dot_test.txt"
    mat_dot = score.dot_product(mat_emb1, mat_emb2, dot_file)

    prot1 = "../data/6PF2K_1BIF.fasta"
    prot2 = "../data/ADK_2AK3A.fasta"

    mat_align = alignglobal.alignment(mat_dot)
    #print(mat_align)
    prot_fasta1 = readfile.read_fasta(prot1)
    prot_fasta2 = readfile.read_fasta(prot2) 
    glob_bif, glob_adk = \
        alignglobal.needleman_wunsch(mat_align, prot_fasta1, prot_fasta2)
    print(glob_bif)
    print(f'{glob_adk}\n')
    
    loc_bif, loc_adk = \
        alignlocal.smith_waterman(mat_align, prot_fasta1, prot_fasta2)
    print(loc_bif)
    print(loc_adk)
    #with open ("../results/bif_adk.txt", "w") as file:
    #    file.write(f'{glob_bif}\n{glob_adk}\n\n')
    #    file.write(f'{loc_bif}\n{loc_adk}')
    # faire un file avec toutes les data à modif ?
        