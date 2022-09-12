#! /usr/bin/python3
"""Function to compute the dot product between to 2 sequences."""
import numpy as np
 
def dot_product(emb_arr1, emb_arr2, file_out):
    """on fait un dot product entre les vect des pos de chaque prot 2 a 2.
    
    Parameters
    ----------
    emb_arr1 : array
        mat emb de la prot1
    emb_arr2 : array
        mat emb de la prot2
    file_out : str
        Name of the file to export the dot product.
    
    Returns
    -------
    array
        array 2 dim de taille len(prot1, prot2) avec dot prod entre chaque pos.
        """
    res = np.dot(emb_arr1, emb_arr2.T)
    # mettre option os si file existe ou non
    #np.savetxt(file_out, res, delimiter="\t")
    return res

