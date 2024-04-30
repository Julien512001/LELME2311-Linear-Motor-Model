from SystemResolution import *


P = np.linspace(hm, h, n_mesh) # y
Q = np.linspace(-L+d/2, L-d/2, n_mesh) # x
p,q = np.meshgrid(P,Q)

Jl0 = 0
B0 = 0
C0 = 0

Al0 = np.zeros_like(p)
Aps = np.zeros_like(p)
Apc = np.zeros_like(p)
A_harm = np.zeros_like(p)
A_fond = np.zeros_like(p)
Al = np.zeros_like(p)

Al0 = - mu1 * Jl0 * p**2 + B0 * p + C0


for i in range(1, n_harm+1):
    omega_n = np.pi*i/L

    Mpc = 0
    Mps = 2*Mp/(np.pi*i)*(np.cos(i*np.pi) - np.cos(i*np.pi*(d-L)/L))
    Jls = 5
    Jlc = 0

    A = b[4, i-1]/omega_n
    B = b[5, i-1]/omega_n
    C = b[6, i-1]/omega_n
    D = b[7, i-1]/omega_n

    Aps = C*np.exp(p*omega_n) + D*np.exp(-p*omega_n) - mu0*Mpc/omega_n + mu1*Jls/(omega_n**2)
    Apc = A*np.exp(p*omega_n) + B*np.exp(-p*omega_n) + mu0*Mps/omega_n + mu1*Jlc/(omega_n**2)

    if (i == 1):
        A_fond = Aps*np.sin(omega_n*q) + Apc*np.cos(omega_n*q)
    else:
        A_harm += Aps*np.sin(omega_n*q) + Apc*np.cos(omega_n*q)
    



Al = Al0 + A_fond + A_harm