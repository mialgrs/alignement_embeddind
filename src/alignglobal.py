#! /usr/bin/python3
"""functions for global alignment."""
import numpy as np

def alignment (dotprod):
    """
    dot prod
    
    Parameters
    ----------
    dotprod : np.array
    
    Returns
    -------
    np.array

    """
    matrice = np.zeros((len(dotprod)+1, len(dotprod[0])+1))
    #print(matrice.shape)
    for i in range(1, len(matrice)):
        for j in range(1, len(matrice[0])):
            diag = matrice[i-1,j-1] + dotprod[i-1,j-1]
            droite = matrice[i,j-1] + 0
            haut = matrice[i-1,j] + 0
            matrice[i,j] = max(diag,droite,haut)
    np.savetxt("mat_align2.txt", matrice, delimiter="\t")
    return matrice

def needle_recurs(mat, fasta1, fasta2, i, j, seq1, seq2):
    """align recurs
    
    Parameters
    ----------
    mat : 
    
    fasta1 :
    
    
    fasta2 :
    
    i :
    
    j :
    
    seq1 : 
    
    seq2 :
    
    
    Returns
    -------
    list, list 
    """
    #print(f"{i}, {j}")
    if i == 0 and j == 0:
        return seq1, seq2
    elif max(mat[i-1,j],mat[i,j-1],mat[i-1,j-1]) == mat[i-1,j-1]:
        seq1.append(fasta1[i-1])
        seq2.append(fasta2[j-1])
        return needle_recurs(mat, fasta1, fasta2, i-1, j-1, seq1, seq2)
    elif max(mat[i-1,j],mat[i,j-1],mat[i-1,j-1]) == mat[i,j-1]:
        seq1.append("-") #pas sur que ce soit j
        seq2.append(fasta2[j-1])
        return needle_recurs(mat, fasta1, fasta2, i, j-1, seq1, seq2)
    else:
        seq1.append(fasta1[i-1])
        seq2.append("-") #pas sur que ce soit i
        return needle_recurs(mat, fasta1, fasta2, i-1, j, seq1, seq2)
    """while i > 0 and j > 0:
        if mat[i,j] == mat[i-1,j]:
            seq1.append(fasta1[i-1])
            seq2.append("-") #pas sur que ce soit i
            i -= 1
        elif mat[i,j] == mat[i,j-1]:
            seq1.append("-") #pas sur que ce soit j
            seq2.append(fasta2[j-1])
            j -= 1
        else:
            seq1.append(fasta1[i-1])
            seq2.append(fasta2[j-1])
            i -= 1
            j -= 1
    while i > 0 and j > 0:
        if max(mat[i-1,j],mat[i,j-1],mat[i-1,j-1]) == mat[i-1,j-1]:
            seq1.append(fasta1[i-1])
            seq2.append(fasta2[j-1])
            i -= 1
            j -= 1
            
        elif max(mat[i-1,j],mat[i,j-1],mat[i-1,j-1]) == mat[i,j-1]:
            seq1.append("-") #pas sur que ce soit j
            seq2.append(fasta2[j-1])
            j -= 1
        else:
            seq1.append(fasta1[i-1])
            seq2.append("-") #pas sur que ce soit i
            i -= 1"""

# need to get sequence by recursivity
def needleman_wunsch(mat_align, fasta_seq1, fasta_seq2):
    """fontion init align global
    
    Parameters
    ----------
    mat_align :
    
    fasta_seq1 :

    fasta_seq2 :
    
    Returns
    -------
    str, str

    """
    # global
    seq_align1 = []
    seq_align2 = []
    i = len(mat_align)-1
    j = len(mat_align[0])-1

    needle_recurs(mat_align, fasta_seq1, fasta_seq2, \
        i, j, seq_align1, seq_align2)
    seq_align1 = "".join(seq_align1[::-1])
    seq_align2 = "".join(seq_align2[::-1])

    return seq_align1, seq_align2