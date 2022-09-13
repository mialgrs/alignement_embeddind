#! /usr/bin/python3
"""Functions for global alignment."""
import numpy as np

def alignment (dotprod, gap, file_out):
    """Get an alignment matrice between 2 sequences.

    It uses dot prod matrice to determine values of
    the align matrice with Needleman&Wunsch algo.

    Parameters
    ----------
    dotprod : numpy array
        Dot product matrice between 2 sequences.
    gap : int

    file_out : str


    Returns
    -------
    numpy array
        Alignment matrice with Needleman&Wunsch algorithm.
    """
    matrice = np.zeros((len(dotprod)+1, len(dotprod[0])+1))
    for i in range(1, len(matrice)):
        for j in range(1, len(matrice[0])):
            diag = matrice[i-1,j-1] + dotprod[i-1,j-1]
            left = matrice[i,j-1] + gap
            up = matrice[i-1,j] + gap
            matrice[i,j] = max(diag,left,up)
    np.savetxt(file_out, matrice, delimiter='\t')
    return matrice

def needle_recurs(mat, fasta1, fasta2, i, j, seq1, seq2):
    """Traceback the global alignment by recursivity.

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
    if i == 0 and j == 0:
        return seq1, seq2

    elif max(mat[i-1,j],mat[i,j-1],mat[i-1,j-1]) == mat[i-1,j-1]:
        seq1.append(fasta1[i-1])
        seq2.append(fasta2[j-1])
        return needle_recurs(mat, fasta1, fasta2, i-1, j-1, seq1, seq2)

    elif max(mat[i-1,j],mat[i,j-1],mat[i-1,j-1]) == mat[i,j-1]:
        seq1.append("-")
        seq2.append(fasta2[j-1])
        return needle_recurs(mat, fasta1, fasta2, i, j-1, seq1, seq2)

    else:
        seq1.append(fasta1[i-1])
        seq2.append("-")
        return needle_recurs(mat, fasta1, fasta2, i-1, j, seq1, seq2)

def needleman_wunsch(mat_align, fasta_seq1, fasta_seq2, file_out):
    """Initialize and

    Parameters
    ----------
    mat_align : numpy array
        Alignment matrice between the 2 sequences.
    fasta_seq1 : list
        Sequence of the first protein.
    fasta_seq2 : list
        Sequence of the first protein.
    file_out : str


    Returns
    -------
    str, str
        The 2 aligned sequences.
    """
    seq_align1 = []
    seq_align2 = []
    i = len(mat_align)-1
    j = len(mat_align[0])-1

    needle_recurs(mat_align, fasta_seq1, fasta_seq2, \
        i, j, seq_align1, seq_align2)
    seq_align1 = "".join(seq_align1[::-1])
    seq_align2 = "".join(seq_align2[::-1])

    with open(file_out, 'w') as file:
        file.write(f'{seq_align1}\n{seq_align2}')

    return seq_align1, seq_align2
    