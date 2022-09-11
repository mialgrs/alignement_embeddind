#! /usr/bin/python3

import numpy as np
import pandas as pd
import readfile 
import alignglobal
import score


if __name__=="__main__":

    emb1 = "../data/5_3_exonuclease_1bgxt.t5emb"

    emb2 = "../data/6PF2K_1bif.t5emb"
    mat_emb1 = readfile.read_emb(emb1)
    mat_emb2 = readfile.read_emb(emb2)

    dot_file = "../results/dot_test.txt"
    mat_dot = score.dot_product(mat_emb1, mat_emb2, dot_file)

    prot1 = "../data/5_3_EXONUCLEASE_1BGXT.fasta"
    prot2 = "../data/6PF2K_1BIF.fasta"

    mat_align = alignglobal.alignment(mat_dot)
    #print(mat_align)
    prot_fasta1 = readfile.read_fasta(prot1)
    prot_fasta2 = readfile.read_fasta(prot2) 
    align_exon, align_bif = \
        alignglobal.needleman_wunsch(mat_align, prot_fasta1, prot_fasta2)
    print(align_exon)
    print(align_bif)
    with open ("exon_bif.txt", "w") as file:
        file.write(f'{align_exon}\n')
        file.write(align_bif)


    # faire un file avec toutes les data à modif ?
        