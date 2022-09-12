#! /usr/bin/python3
"""Functions for semi-global alignment.

The alignment matrix is computed the same way 
as a global alignment matrix so the function of the global align
is used for the semi-global align."""
def needle_recurs(mat, fasta1, fasta2, i, j, seq1, seq2):
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

# need to get sequence by recursivity
def needleman_wunsch(mat_align, fasta_seq1, fasta_seq2):
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