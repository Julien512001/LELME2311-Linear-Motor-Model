from SystemResolution import *


# We compute the potential vector field for the air-gap (The second region)

Jl0 = 0
B0 = 0
C0 = 0

Aps = np.zeros_like(P2)
Apc = np.zeros_like(P2)
A_harm = np.zeros_like(P2)
A_fond = np.zeros_like(P2)
Al = np.zeros_like(P2)

phi_fond = np.zeros_like(P2)
phi_harm = np.zeros_like(P2)

F_active = np.zeros_like(P2)
F_ripple = np.zeros_like(P2)


for i in range(1, n_harm+1):
    omega_n = np.pi*i/tau_k

    Mpc = 0
    Mps = -2*Mp/(np.pi*i) * (np.cos(np.pi*i/(tau_k)*e/2) - np.cos(np.pi*i/(tau_k)*(e/2 - tau_k)))
    # Jl = 0, car ici on veut calculer l'effet du champ sur le bobinage et pas l'effet du bobinage sur le champ
    Jls = 0
    Jlc = 0

    A = b[4, i-1]/omega_n
    B = b[5, i-1]/omega_n
    C = b[6, i-1]/omega_n
    D = b[7, i-1]/omega_n

    Aps = C*np.exp(P2*omega_n) + D*np.exp(-P2*omega_n) - mu0*Mpc/omega_n + mu1*Jls/(omega_n**2)
    Apc = A*np.exp(P2*omega_n) + B*np.exp(-P2*omega_n) + mu0*Mps/omega_n + mu1*Jlc/(omega_n**2)



    if (i == 1):
        A_fond = Aps*np.sin(omega_n*q) + Apc*np.cos(omega_n*q)
        phi_fond = Nt*2*Lz*A_fond[250][250]
        F_active = 3/2 * omega_n * phi_fond * I
        print("Fondamentale")
        print(F_active)
    else:
        A_int = Aps*np.sin(omega_n*q) + Apc*np.cos(omega_n*q)
        A_harm += A_int
        phi_int = Nt*2*Lz*A_int[250][250]
        phi_harm += phi_int
        F_int = 3/2 * omega_n * phi_int * I
        F_ripple += F_int
        
        print("----")
        print(i)
        print(F_int)
        


Al = A_fond + A_harm
phi = phi_fond + phi_harm
F = F_active + F_ripple
