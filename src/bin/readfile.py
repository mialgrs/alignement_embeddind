#! /usr/bin/python3
"""Functions to read data from files."""

import numpy as np

def read_emb(file):
    """Put embedding data from a file to an array.

    Parameters
    ----------
    file : str
        Name of the embedding file.

    Returns
    -------
    numpy array
        Matrix of embedding values by lenght of the protein sequence."""
    with open (file, 'r') as emb:
        #aaa = np.zeros((1,1024))
        liste = []
        for line in emb:
            vecteur = line.strip().split()
            vecteur = [float(i) for i in vecteur]
            #aaa = np.concatenate((aaa, np.array([vecteur])))
            liste.append(vecteur)
            #hey = aaa[1:]
        emb_arr = np.array(liste)
        #return hey
        return emb_arr

def read_fasta(file):
    """Get a fasta file as entry to return the corresponding seq in a list.

    Parameters
    ----------
    file : str
        Name of fasta file.

    Returns
    -------
    list
        Sequence with each position separated.
    """
    line_seq = ""
    with open(file, "r") as fasta:
        for line in fasta:
            if not line.startswith(">"):
                line_seq += line.strip()
        seq = list(line_seq)
    return seq
