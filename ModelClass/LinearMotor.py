import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve


class LinearMotor:

    n_harm = 100

    # Permeability definition
    # Air : mu
    mu0 = 4*np.pi*10**(-7)
    # Region 1 : mu
    mu1 = 1.05*mu0
    # Region2 : mu
    mu2 = 1*mu0

    # Section du câble en mm^2
    S = 10
    Nt = 10

    # Densité de courant en A/mm^2
    J = 5

    # Courant de la machine
    I = J * S /Nt



    def __init__(self, tau_k, Br, hm, ha, Lz):
        self.tau_k = tau_k
        self.Mp = Br/self.mu0
        self.hm = hm
        self.ha = ha
        self.Lz = Lz
        self.h = ha + hm
        self.e = 1/3 * tau_k
        self.y = self.h + 0.1e-3
        self.x = -self.tau_k

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

        Bq1 = 0
        Bp1 = 0

        for i in range(1, self.n_harm+1):

            Mps = -2*self.Mp/(np.pi*i) * (np.cos(np.pi*i/(self.tau_k)*self.e/2) - np.cos(np.pi*i/(self.tau_k)*(self.e/2 - self.tau_k)))

            omega_n = np.pi*i/self.tau_k
            an_1 = self.b[0,i-1]
            bn_1 = self.b[1,i-1]
            cn_1 = self.b[2,i-1]
            dn_1 = self.b[3,i-1]

            Bqc_1 = an_1*np.exp(omega_n*self.y) - bn_1*np.exp(-omega_n*self.y)
            Bqs_1 = cn_1*np.exp(omega_n*self.y) - dn_1*np.exp(-omega_n*self.y)

            Bpc_1 = -cn_1*np.exp(omega_n*self.y) - dn_1*np.exp(-omega_n*self.y)
            Bps_1 = an_1*np.exp(omega_n*self.y) + bn_1*np.exp(-omega_n*self.y) + self.mu0*Mps

            Bq1 = Bq1 + Bqs_1*np.sin(omega_n*self.x) + Bqc_1*np.cos(omega_n*self.x)
            Bp1 = Bp1 + Bpc_1*np.cos(omega_n*self.x) + Bps_1*np.sin(omega_n*self.x)
        
        self.Bq1 = Bq1
        self.Bp1 = Bp1
    
    def Region2(self):

        Bq2 = 0
        Bp2 = 0


        for i in range(1, self.n_harm+1):

            Mps = -2*self.Mp/(np.pi*i) * (np.cos(np.pi*i/(self.tau_k)*self.e/2) - np.cos(np.pi*i/(self.tau_k)*(self.e/2 - self.tau_k)))
            
            omega_n = np.pi*i/self.tau_k
            an_2 = self.b[4,i-1]
            bn_2 = self.b[5,i-1]
            cn_2 = self.b[6,i-1]
            dn_2 = self.b[7,i-1]

            Bqc_2 = an_2*np.exp(omega_n*self.y) - bn_2*np.exp(-omega_n*self.y)
            Bqs_2 = cn_2*np.exp(omega_n*self.y) - dn_2*np.exp(-omega_n*self.y)

            Bpc_2 = -cn_2*np.exp(omega_n*self.y) - dn_2*np.exp(-omega_n*self.y)
            Bps_2 = an_2*np.exp(omega_n*self.y) + bn_2*np.exp(-omega_n*self.y)

            Bq2 = Bq2 + Bqs_2*np.sin(omega_n*self.x) + Bqc_2*np.cos(omega_n*self.x)
            Bp2 = Bp2 + Bpc_2*np.cos(omega_n*self.x) + Bps_2*np.sin(omega_n*self.x)
        
        self.Bq2 = Bq2
        self.Bp2 = Bp2

    def get_F(self):
        
        Aps = 0
        Apc = 0
        A_harm = 0
        A_fond = 0

        phi_fond = 0
        phi_harm = 0

        F_active = 0
        F_ripple = 0
        F_THD = 0

        for i in range(1, self.n_harm+1):
            omega_n = np.pi*i/self.tau_k

            Mpc = 0
            Mps = -2*self.Mp/(np.pi*i) * (np.cos(np.pi*i/(self.tau_k)*self.e/2) - np.cos(np.pi*i/(self.tau_k)*(self.e/2 - self.tau_k)))

            A = self.b[4, i-1]/omega_n
            B = self.b[5, i-1]/omega_n
            C = self.b[6, i-1]/omega_n
            D = self.b[7, i-1]/omega_n

            Aps = C*np.exp(self.y*omega_n) + D*np.exp(-self.y*omega_n) - self.mu0*Mpc/omega_n
            Apc = A*np.exp(self.y*omega_n) + B*np.exp(-self.y*omega_n) + self.mu0*Mps/omega_n

            if (i == 1):
                A_fond = Aps*np.sin(omega_n*self.x) + Apc*np.cos(omega_n*self.x)

                phi_fond = self.Nt*2*self.Lz*A_fond

                F_active = 3/2 * omega_n * phi_fond * self.I

            else:
                A_int = Aps*np.sin(omega_n*self.x) + Apc*np.cos(omega_n*self.x)
                A_harm += A_int

                phi_int = self.Nt*2*self.Lz*A_int
                phi_harm += phi_int

                F_int = 3/2 * omega_n * phi_int * self.I
                F_ripple += F_int
                F_THD += F_int**2


        F = F_active + F_ripple
        self.F = F
        self.F_active = F_active
        self.F_ripple = F_ripple
        self.THD = np.sqrt(F_THD)/np.abs(F_active)

    def myMotor(self):
        self.systemResolution()
        self.Region1()
        self.Region2()
        self.get_F()
        print("F = {}".format(self.F))
        print("F_active = {}".format(self.F_active))
        print("F_ripple = {}".format(self.F_ripple))
        print("THD = {}".format(self.THD))

# tau_k, Br, hm, ha, Lz
p1 = LinearMotor(100e-3, 1.21, 10e-3, 0.25e-3, 45e-3)
p1.myMotor()


