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

    Nt = 1

    # Densité de courant en A/m^2
    J = 5*1e6
    eta = 0.5



    def __init__(self, tau_k, e, Br, hm, ha, Lz, lp, lq):
        self.tau_k = tau_k
        self.e = e
        self.Mp = Br/self.mu0
        self.hm = hm
        self.ha = ha
        self.Lz = Lz
        self.h = ha + hm
        self.y = self.h
        self.x = -self.tau_k
        self.Sbob = lp*lq
    
    def constant(self):
        # Courant de la machine
        self.I = self.J * self.Sbob  * self.eta/self.Nt

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
        return b

    def get_F_active(self):
        self.constant()
        b = self.systemResolution()
    
        omega_n = np.pi/self.tau_k
        Mps = -2*self.Mp/(np.pi) * (np.cos(np.pi/(self.tau_k)*self.e/2) - np.cos(np.pi/(self.tau_k)*(self.e/2 - self.tau_k)))

        A = b[4, 0]/omega_n
        B = b[5, 0]/omega_n
        C = b[6, 0]/omega_n
        D = b[7, 0]/omega_n

        Aps = C*np.exp(self.y*omega_n) + D*np.exp(-self.y*omega_n)
        Apc = A*np.exp(self.y*omega_n) + B*np.exp(-self.y*omega_n) + self.mu0*Mps/omega_n

        q = np.linspace(-self.tau_k, self.tau_k, 1000)

        A_fond = np.max(Aps*np.sin(omega_n*q) + Apc*np.cos(omega_n*q))

        psi_fond = self.Nt*2*self.Lz*A_fond

        F_active = 3/2 * omega_n * psi_fond * self.I

        self.F_active = F_active

        return F_active

    def get_THD(self):
        self.constant()
        b = self.systemResolution()
        q = np.linspace(-self.tau_k, self.tau_k, 1000)

        omega_n = np.pi/self.tau_k

        Mps = -2*self.Mp/(np.pi) * (np.cos(np.pi/(self.tau_k)*self.e/2) - np.cos(np.pi/(self.tau_k)*(self.e/2 - self.tau_k)))

        A = b[4, 0]/omega_n
        B = b[5, 0]/omega_n
        C = b[6, 0]/omega_n
        D = b[7, 0]/omega_n

        Aps = C*np.exp(self.y*omega_n) + D*np.exp(-self.y*omega_n)
        Apc = A*np.exp(self.y*omega_n) + B*np.exp(-self.y*omega_n) + self.mu0*Mps/omega_n

        A_fond = np.max(Aps*np.sin(omega_n*q) + Apc*np.cos(omega_n*q))
        psi_fond = self.Nt*2*self.Lz*A_fond
        
        psi_THD = 0

        for i in range(3, self.n_harm+1):
            omega_n = np.pi*i/self.tau_k

            Mpc = 0
            Mps = -2*self.Mp/(np.pi*i) * (np.cos(np.pi*i/(self.tau_k)*self.e/2) - np.cos(np.pi*i/(self.tau_k)*(self.e/2 - self.tau_k)))

            A = self.b[4, i-1]/omega_n
            B = self.b[5, i-1]/omega_n
            C = self.b[6, i-1]/omega_n
            D = self.b[7, i-1]/omega_n


            Aps = C*np.exp(self.y*omega_n) + D*np.exp(-self.y*omega_n) - self.mu0*Mpc/omega_n
            Apc = A*np.exp(self.y*omega_n) + B*np.exp(-self.y*omega_n) + self.mu0*Mps/omega_n

            A_int = np.max(Aps*np.sin(omega_n*q) + Apc*np.cos(omega_n*q))
            psi_int = self.Nt*2*self.Lz*A_int
            psi_THD += psi_int*psi_int

        THD = np.sqrt(psi_THD)/psi_fond
        self.THD = THD

        return THD

    def myMotor(self):
        self.get_F_active()
        self.get_THD()
        print("F_active = {}".format(self.F_active))
        print("THD = {}".format(self.THD))

# tau_k, e, Br, hm, ha, Lz, lp, lq

p1 = LinearMotor(20e-3, 5e-3, 1.21, 7e-3, 3e-3, 3e-3, 1e-3, 3e-3)
p1.myMotor()


