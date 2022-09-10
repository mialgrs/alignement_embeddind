#! /usr/bin/python3

import numpy as np
import pandas as pd
import readfile 
import alignglobal
import score


if __name__=="__main__":

    emb1 = "../data/testeur.txt"

    emb2 = "../data/testeur2.txt"
    mat_emb1 = readfile.read_emb(emb1)
    mat_emb2 = readfile.read_emb(emb2)
    #print(mat_emb1[0,0])
    dot_file = "../results/dot_test.txt"
    mat_dot = score.dot_product(mat_emb1, mat_emb2, dot_file)
    print(len(mat_dot))
    print(len(mat_dot[0]))
    # dot prod entre les 2 seq fait 
    #np.amax == max()
    print(np.amax(mat_dot))
    #mtn faut faire la matrice align
    prot1 = "../data/test1.txt"
    prot2 = "../data/test2.txt"

    mat_align = alignglobal.alignment(mat_dot)
    print(mat_align)
    prot_fasta1 = readfile.read_fasta(prot1)
    prot_fasta2 = readfile.read_fasta(prot2)
    print(prot_fasta1)
    align_exon, align_bif = \
        alignglobal.needleman_wunsch(mat_align, prot_fasta1, prot_fasta2)
    print(align_exon)
    print(align_bif)
        