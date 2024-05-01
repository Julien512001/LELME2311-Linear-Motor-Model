import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve


class LinearMotor:

    n_harm = 10

    # Permeability definition
    # Air : mu
    mu0 = 4*np.pi*10**(-7)
    # Region 1 : mu
    mu1 = 1.05*mu0
    # Region2 : mu
    mu2 = 1*mu0


    def __init__(self, tau_k, Br, hm, ha):
        self.tau_k = tau_k
        self.Mp = Br/self.mu0
        self.hm = hm
        self.ha = ha
        self.h = ha + hm
        self.e = 1/3 * tau_k

    def systemResolution(self):
        b = np.zeros((8,self.n_harm))
        Mps = 0

        for i in range(1,self.n_harm+1):
            omega_n = i*np.pi/self.tau_k

            Mps = -2*self.Mp/(np.pi*i) * (np.cos(np.pi*i/(self.tau_k)*self.e/2) - np.cos(np.pi*i/(self.tau_k)*(self.e/2 - self.tau_k)))

            exp_pos_m = np.exp(omega_n*self.hm)
            exp_neg_m = np.exp(-omega_n*self.hm)
            exp_pos_a = np.exp(omega_n*self.h)
            exp_neg_a = np.exp(-omega_n*self.h)
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
            A[2][2] = self.mu2*exp_pos_m
            A[2][3] = -self.mu2*exp_neg_m
            A[2][6] = -self.mu1*exp_pos_m
            A[2][7] = self.mu1*exp_neg_m

            # Fourth row
            A[3][0] = self.mu2*exp_pos_m
            A[3][1] = -self.mu2*exp_neg_m
            A[3][4] = -self.mu1*exp_pos_m
            A[3][5] = self.mu1*exp_neg_m

            # Fifth row
            A[4][0] = 1.0
            A[4][1] = -1.0


            # Sixth row
            A[5][2] = 1.0
            A[5][3] = -1.0

            # seventh row

            A[6][4] = exp_pos_a
            A[6][5] = -exp_neg_a

            # eigth row
            A[7][6] = exp_pos_a
            A[7][7] = -exp_neg_a

            # Vector definition
            b_bis[0] = -self.mu0*Mps
            x = np.linalg.solve(A, b_bis)

            b[:,i-1] = x

        self.b = b
    
    def Region1(self):
        # Compute the coefficient of B for the region 1 with the magnets
        """
        Bq0_1 = np.zeros_like(P1)
        Bqc_1 = np.zeros_like(P1)
        Bqs_1 = np.zeros_like(P1)
        Bpc_1 = np.zeros_like(P1)
        Bps_1 = np.zeros_like(P1)

        Bq1 = np.zeros_like(P1)
        Bp1 = np.zeros_like(P1)


        for i in range(1, n_harm+1):

            Mps = -2*Mp/(np.pi*i) * (np.cos(np.pi*i/(tau_k)*e/2) - np.cos(np.pi*i/(tau_k)*(e/2 - tau_k)))

            omega_n = np.pi*i/tau_k
            an_1 = b[0,i-1]
            bn_1 = b[1,i-1]
            cn_1 = b[2,i-1]
            dn_1 = b[3,i-1]

            Bqc_1 = an_1*np.exp(omega_n*P1) - bn_1*np.exp(-omega_n*P1)
            Bqs_1 = cn_1*np.exp(omega_n*P1) - dn_1*np.exp(-omega_n*P1)

            Bpc_1 = -cn_1*np.exp(omega_n*P1) - dn_1*np.exp(-omega_n*P1)
            Bps_1 = an_1*np.exp(omega_n*P1) + bn_1*np.exp(-omega_n*P1) + mu0*Mps

            Bq1 = Bq1 + Bqs_1*np.sin(omega_n*Q1) + Bqc_1*np.cos(omega_n*Q1)
            Bp1 = Bp1 + Bpc_1*np.cos(omega_n*Q1) + Bps_1*np.sin(omega_n*Q1)
        """

    
    def myMotor(self):
        self.systemResolution()
    
p1 = LinearMotor(100e-3, 1.21, 10e-3, 0.25e-3)
p1.myMotor()


