import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve
from Constant import *


b = np.zeros((8,n_harm))
Mps = 0

for i in range(1,n_harm+1):
    omega_n = i*np.pi/tau_k

    Mps = -2*Mp/(np.pi*i) * (np.cos(np.pi*i/(tau_k)*e/2) - np.cos(np.pi*i/(tau_k)*(e/2 - tau_k)))

    exp_pos_m = np.exp(omega_n*hm)
    exp_neg_m = np.exp(-omega_n*hm)
    exp_pos_a = np.exp(omega_n*h)
    exp_neg_a = np.exp(-omega_n*h)
    A = np.zeros((8,8))
    b_bis = np.zeros(8)

    # First row
    A[0][0] = exp_pos_m
    A[0][1] = exp_neg_m
    A[0][4] = -exp_pos_m
    A[0][5] = -exp_neg_m

    # Second row
    A[1][2] = exp_pos_m
    A[1][3] = exp_neg_m
    A[1][6] = -exp_pos_m
    A[1][7] = -exp_neg_m

    # Third row
    A[2][2] = mu2*exp_pos_m
    A[2][3] = -mu2*exp_neg_m
    A[2][6] = -mu1*exp_pos_m
    A[2][7] = mu1*exp_neg_m

    # Fourth row
    A[3][0] = mu2*exp_pos_m
    A[3][1] = -mu2*exp_neg_m
    A[3][4] = -mu1*exp_pos_m
    A[3][5] = mu1*exp_neg_m

    # Fifth row
    A[4][0] = 1
    A[4][1] = -1


    # Sixth row
    A[5][2] = 1
    A[5][3] = -1

    # seventh row

    A[6][4] = exp_pos_a
    A[6][5] = -exp_neg_a

    # eigth row
    A[7][6] = exp_pos_a
    A[7][7] = -exp_neg_a

    # Vector definition
    b_bis[0] = -mu0*Mps
    x = np.linalg.solve(A, b_bis)

    b[:,i-1] = x
"""
def write_matrix_to_file(matrix, file_name):
    with open(file_name, 'w') as file:
        for row in matrix:
            row_str = ' '.join(map(str, row))
            file.write(row_str + '\n')


file_name = 'matrice2.txt'
write_matrix_to_file(b, file_name)
print("La matrice a été écrite dans le fichier", file_name)
"""