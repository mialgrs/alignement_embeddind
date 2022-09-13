#! /usr/bin/python3
"""Functions for semi-global alignment.

The alignment matrix is computed the same way 
as a global alignment matrix so the function of the global align
is used for the semi-global align."""
import numpy as np 

def semi_recurs(mat, fasta1, fasta2, i, j, seq1, seq2):
    """Traceback the semi-global alignment by recursivity. 
    
    Parameters
    ----------
    mat : numpy array
        Alignment matrice between the 2 sequences.
    fasta1 : list
        Sequence of the first protein.
    fasta2 : list
        Sequence of the second protein.
    i : int
        The index where the traceback begin.
    j : int
        The column where the traceback begin.
    seq1 : list
        The aligned sequence of the first protein.
    seq2 : list
        The aligned sequence of the second protein.
    
    Returns
    -------
    list, list 
        The 2 aligned sequences.
    """
    if j == 0:
        return seq1, seq2
    elif max(mat[i-1,j],mat[i,j-1],mat[i-1,j-1]) == mat[i-1,j-1]:
        seq1.append(fasta1[i-1])
        seq2.append(fasta2[j-1])
        return semi_recurs(mat, fasta1, fasta2, i-1, j-1, seq1, seq2)
    elif max(mat[i-1,j],mat[i,j-1],mat[i-1,j-1]) == mat[i,j-1]:
        seq1.append("-")
        seq2.append(fasta2[j-1])
        return semi_recurs(mat, fasta1, fasta2, i, j-1, seq1, seq2)
    else:
        seq1.append(fasta1[i-1])
        seq2.append("-") 
        return semi_recurs(mat, fasta1, fasta2, i-1, j, seq1, seq2)

def semi_global(mat_align, fasta_seq1, fasta_seq2):
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
    seq_align1 = []
    seq_align2 = []
    j = len(mat_align[0])-1
    i = np.argmax(mat_align[:, j])
    semi_recurs(mat_align, fasta_seq1, fasta_seq2, \
        i, j, seq_align1, seq_align2)
    seq_align1 = "".join(seq_align1[::-1])
    seq_align2 = "".join(seq_align2[::-1])

    return seq_align1, seq_align2