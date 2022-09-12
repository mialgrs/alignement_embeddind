#! /usr/bin/python3
"""Functions for local alignment."""

import numpy as np

def alignment(dotprod):
    """Get an alignment matrice between 2 sequences.
    
    It uses dot prod matrice to determine values of 
    the align matrice with Smith&Waterman algo.
    Parameters
    ----------
    dotprod : numpy array
        Dot product matrice between 2 sequences.
    
    Return
    ------
    numpy array
        Alignment matrice with Smith&Waterman algorithm."""
    matrice = np.zeros((len(dotprod)+1, len(dotprod[0])+1))
    for i in range(1, len(matrice)):
        for j in range(1, len(matrice[0])):
            diag = matrice[i-1,j-1] + dotprod[i-1,j-1]
            left = matrice[i,j-1] + 0
            up = matrice[i-1,j] + 0
            matrice[i,j] = max(diag,left,up, 0)
    #np.savetxt("mat_align2.txt", matrice, delimiter="\t")
    return matrice

def smith_recurs(mat, fasta1, fasta2, i, j, seq1, seq2):
    """align recurs
    
    Parameters
    ----------
    mat : numpy array
    
    fasta1 : list
    
    fasta2 : list
    
    i : int
    
    j : int
    
    seq1 : list
    
    seq2 : list
    
    
    Returns
    -------
    list, list 

    """
    if mat[i,j] == 0:
        return seq1, seq2
    elif max(mat[i-1,j],mat[i,j-1],mat[i-1,j-1]) == mat[i-1,j-1]:
        seq1.append(fasta1[i-1])
        seq2.append(fasta2[j-1])
        return smith_recurs(mat, fasta1, fasta2, i-1, j-1, seq1, seq2)
    elif max(mat[i-1,j],mat[i,j-1],mat[i-1,j-1]) == mat[i,j-1]:
        seq1.append("-") #pas sur que ce soit j
        seq2.append(fasta2[j-1])
        return smith_recurs(mat, fasta1, fasta2, i, j-1, seq1, seq2)
    else:
        seq1.append(fasta1[i-1])
        seq2.append("-") #pas sur que ce soit i
        return smith_recurs(mat, fasta1, fasta2, i-1, j, seq1, seq2)

def smith_waterman(mat_align, fasta_seq1, fasta_seq2):
    """fontion init align global
    
    Parameters
    ----------
    mat_align : numpy array
    
    fasta_seq1 : list

    fasta_seq2 : list
    
    Returns
    -------
    str, str

    """
    seq_align1 = {}
    seq_align2 = {}
    seq1 = []
    seq2 = []
    print(np.argwhere(mat_align == np.amax(mat_align)))
    ind = np.argwhere(mat_align == np.amax(mat_align))
    for i in range(len(ind)):
        seq_align1[i], seq_align2[i] = smith_recurs(mat_align, fasta_seq1, fasta_seq2, \
            ind[i,0], ind[i,1], seq1, seq2)
    for key in seq_align1.keys():
        seq_align1[key] = "".join(seq_align1[key][::-1])
        seq_align2[key] = "".join(seq_align2[key][::-1])
    return seq_align1, seq_align2

