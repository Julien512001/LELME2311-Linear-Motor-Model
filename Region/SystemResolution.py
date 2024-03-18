import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve
from Constant import *


b = np.zeros((8,n_harm))
Mps = 0

for i in range(1,n_harm+1):
    omega_n = i*np.pi/L

    Mps = Mps + 2*L/(np.pi*i)*(np.cos(i*np.pi) - np.cos(i*np.pi*(d-L)/L))

    exp_pos_m = np.exp(omega_n*hm)
    exp_neg_m = np.exp(-omega_n*hm)
    exp_pos_a = np.exp(omega_n*ha) 
    exp_neg_a = np.exp(-omega_n*ha)
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

    x = solve(A, b_bis)

    b[:,i-1] = x

