#! /usr/bin/python3

import numpy as np
import pandas as pd

def read_emb(file):
    """function qui recup emb """
    with open (file, 'r') as emb:
        aaa = np.zeros((1,1024))
        print(aaa)
        print(aaa.shape)
        for line in emb:
            vecteur = line.strip().split()
            vecteur = [float(i) for i in vecteur]
            aaa = np.concatenate((aaa, np.array([vecteur])))
            hey = aaa[1:]
        return hey

file1 = "5_3_exonuclease_1bgxt.t5emb"

file2 = "6PF2K_1bif.t5emb"
aaa_arr = read_emb(file1)
#print(aaa_arr)
cohesin_arr = read_emb(file2)
#print(cohesin_arr)


a = np.array([[1, 2, 3, 4, 5], [2, 2, 3, 4, 5]])

b = np.array([[2, 3, 4, 5, 6], [3, 3, 4, 5, 6]])

c = np.dot(aaa_arr, cohesin_arr.T)
print(c)

