import numpy as np
import matplotlib.pyplot as plt


# Problem dimensions in [mm]
tau_k = 100e-3           # half period of the problem : 42 mm
hm = 10e-3          # PM height (region 1) : 45mm
ha = 0.25e-3          # air-gap height (region 2) : 15mm
h = hm + ha         # Total height
e = 1/3 * tau_k             # PM air-gap lenght --> to be defined in term of tau_k
q_femm = tau_k/2


# Linspace and mesh definition
n_mesh = 500

q = np.linspace(-tau_k, tau_k, n_mesh)

# For the region 1
p1 = np.linspace(0, hm, n_mesh)
Q1, P1 = np.meshgrid(q, p1)

# For the region 2
p2 = np.linspace(hm, h, n_mesh)
Q2, P2 = np.meshgrid(q, p2)

# Permeability definition
# Air : mu
mu0 = 4*np.pi*10**(-7)
# Region 1 : mu
mu1 = 1.05*mu0
# Region2 : mu
mu2 = 1*mu0

# PM definition
Br = 1.21
Mp = Br/mu0

# Number of harmonics
n_harm = 100


# Define the function
def magnetization_profile(q, e, Mp):
    profile = np.zeros_like(q)
    profile[(q >= -tau_k) & (q < -tau_k + e/2)] = 0
    profile[(q >= -tau_k + e/2) & (q < -e/2)] = Mp
    profile[(q >= -e/2) & (q <= e/2)] = 0
    profile[(q >= e/2) & (q < tau_k - e/2)] = -Mp
    profile[(q >= tau_k - e/2) & (q <= tau_k)] = 0
    return profile

f = np.zeros_like(q)


for i in range(1, n_harm+1):
    Mps = -2*Mp/(np.pi*i) * (np.cos(np.pi*i/(tau_k)*e/2) - np.cos(np.pi*i/(tau_k)*(e/2 - tau_k)))
    f += Mps*np.sin(np.pi*i/tau_k * q)


# Plotting
plt.figure(figsize=(8, 6))
plt.plot(q, magnetization_profile(q, e, Mp), color='b', label='Magnetization Profile')
plt.plot(q, f, color='r', label='Fourier')
plt.xlabel('q [mm]')
plt.ylabel('$M_p$')
plt.title('Magnetization Profile')
plt.grid(True)
plt.legend()
plt.show()