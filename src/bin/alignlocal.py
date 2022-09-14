#! /usr/bin/python3
"""Functions for local alignment."""

import numpy as np

def alignment(dotprod, gap, file_out):
    """Get an alignment matrice between 2 sequences.

    It uses the dot product matrix and the Smith&Waterman algorithm  matrix
    to determine the score at each position in the alignment matrix.

    Parameters
    ----------
    dotprod : numpy array
        Dot product matrice between two sequences.
    gap : int
        The added value when you put a gap.
    file_out : str
        Name of the file to export the alignment matrix.

    Return
    ------
    numpy array
        Alignment matrix with Smith&Waterman algorithm.
    """
    matrice = np.zeros((len(dotprod)+1, len(dotprod[0])+1))
    for i in range(1, len(matrice)):
        for j in range(1, len(matrice[0])):
            diag = matrice[i-1,j-1] + dotprod[i-1,j-1]
            left = matrice[i,j-1] + gap
            up = matrice[i-1,j] + gap
            matrice[i,j] = max(diag,left,up, 0)
    np.savetxt(file_out, matrice, delimiter='\t')
    return matrice

def smith_recurs(mat, fasta1, fasta2, i, j, seq1, seq2):
    """Traceback the local alignment by recursivity.

    Parameters
    ----------
    mat : numpy array
        Alignment matrice between the two sequences.
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
        The two aligned sequences in reverse.
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

def smith_waterman(mat_align, fasta_seq1, fasta_seq2, file_out):
    """Recover the aligned sequences as a string.

    It initializes the variables to call the function smith_recurs(),
    puts the sequences in the right order and exports the result to a file.

    Parameters
    ----------
    mat_align : numpy array
        Alignment matrice between the two sequences.
    fasta_seq1 : list
        Sequence of the first protein.
    fasta_seq2 : list
        Sequence of the second protein.
    file_out : str
        Name of the file to export the alignment matrix.

    Returns
    -------
    dict, dict
        The aligned sequences, for a same key
        we have the sequences associated between them.
    """
    seq_align1 = {}
    seq_align2 = {}
    ind = np.argwhere(mat_align == np.amax(mat_align))
    for i in range(len(ind)):
        seq1 = []
        seq2 = []
        seq_align1[i], seq_align2[i] = \
            smith_recurs(mat_align, fasta_seq1, fasta_seq2, \
            ind[i,0], ind[i,1], seq1, seq2)
    for key in seq_align1.keys():
        seq_align1[key] = "".join(seq_align1[key][::-1])
        seq_align2[key] = "".join(seq_align2[key][::-1])
        with open (file_out, 'a') as file:
            file.write(f'{seq_align2[key]}\n{seq_align1[key]}\n\n')

    return seq_align1, seq_align2
