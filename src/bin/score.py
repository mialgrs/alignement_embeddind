#! /usr/bin/python3
"""Function to compute the dot product between to 2 sequences."""
import numpy as np

def dot_product(emb_arr1, emb_arr2, file_out):
    """Dot product between the vector of each position of the two sequences.

    Parameters
    ----------
    emb_arr1 : array
        Embedding matrix of the first protein.
    emb_arr2 : array
        Embedding matrix of the second protein.
    file_out : str
        Name of the file to export the dot product.

    Returns
    -------
    numpy array
        Two dimensional array of the lenghs of the two sequences.
    """
    res = np.dot(emb_arr1, emb_arr2.T)
    np.savetxt(file_out, res, delimiter="\t")

    return res
