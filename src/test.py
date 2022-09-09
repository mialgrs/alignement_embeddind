#! /usr/bin/python3

import numpy as np
import pandas as pd
import read_fasta 

def read_emb(file):
    """function qui recup emb """
    with open (file, 'r') as emb:
        aaa = np.zeros((1,1024))
        print(aaa)
        print(aaa.shape)
        #liste = []
        for line in emb:
            vecteur = line.strip().split()
            vecteur = [float(i) for i in vecteur]
            aaa = np.concatenate((aaa, np.array([vecteur])))
            #liste.append(vecteur)
            hey = aaa[1:]
        #emb_arr = np.array(liste)
        return hey
        #return emb_arr

def dot_product(emb_arr1, emb_arr2, file_out):
    res = np.dot(emb_arr1, emb_arr2.T)
    # mettre option os si file existe ou non
    np.savetxt(file_out, res, delimiter="\t")
    return res

emb1 = "5_3_exonuclease_1bgxt.t5emb"

emb2 = "Cohesin_1aoha.t5emb"
emb_arr1 = read_emb(emb1)
emb_arr2 = read_emb(emb2)

tkt = dot_product(emb_arr1, emb_arr2, "wesh_cbon2.txt")
print(len(tkt[0]))
# dot prod entre les 2 seq fait 
print(np.amax(tkt))
#mtn faut faire la matrice align
prot1 = "5_3_EXONUCLEASE_1BGXT.fasta"
prot2 = "6PF2K_1BIF.fasta"

#prot_fasta1 = read_fasta(prot1)
#prot_fasta2 = read_fasta(prot2)

def alignment (dotprod):
    matrice = np.zeros((len(dotprod)+1, len(dotprod[0])+1))
    print(matrice.shape)
    for i in range(1, len(matrice)):
        for j in range(1, len(matrice[0])):
            diag = matrice[i-1,j-1] + dotprod[i-1,j-1]
            droite = matrice[i,j-1] + 0
            haut = matrice[i-1,j] + 0
            matrice[i,j] = max(diag,droite,haut)
    np.savetxt("mat_align2.txt", matrice, delimiter="\t")
    return matrice

alignment(tkt)

def needle_recurs(mat, fasta1, fasta2, i, j, seq1, seq2):
    # condition init
    if i == 0 and j == 0:
        return seq1[::-1], seq2[::-1]
    elif mat[i,j] == mat[i-1,j]:
        seq1 += "-"
        seq2 += fasta2[i-1]
        return needle_recurs(mat, fasta1, fasta2, i-1, j, seq1, seq2)
    elif mat[i,j] == mat[i,j-1]:
        seq1 += fasta1[j-1]
        seq2 += "-"
        return needle_recurs(mat, fasta1, fasta2, i, j-1, seq1, seq2)
    else:
        seq1 += fasta1[i-1]
        seq2 += fasta2[j-1] 
        return needle_recurs(mat, fasta1, fasta2, i-1, j-1, seq1, seq2)
    """while i > 0 and j > 0:
        if mat_align[i,j] == mat_align[i-1,j]:
            seq_align1 += "-"
            seq_align2 += fasta_seq2[i-1] #pas sur que ce soit i
            i -= 1
        elif mat_align[i,j] == mat_align[i,j-1]:
            seq_align1 += fasta_seq1[j-1] #pas sur que ce soit j
            seq_align2 += "-"
            j -= 1
        else:
            seq_align1 += fasta_seq1[i-1]
            seq_align2 += fasta_seq2[j-1]
            i -= 1
            j -= 1
    seq_align1 = seq_align1[::-1]
    seq_align2 = seq_align2[::-1]"""

# need to get sequence by recursivity
def needleman_wunsch(mat_align, fasta_seq1, fasta_seq2):
    # global
    seq_align1 = ""
    seq_align2 = ""
    i = len(mat_align)
    j = len(mat_align[0])

    needle_recurs(mat_align, fasta_seq1, fasta_seq2, i, j, seq_align1, seq_align2)

    return seq_align1, seq_align2


        