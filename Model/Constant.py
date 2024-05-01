import numpy as np

# Problem dimensions in [mm]
tau_k = 42e-3           # half period of the problem : 42 mm
hm = 45e-3          # PM height (region 1) : 45mm
ha = 15e-3          # air-gap height (region 2) : 15mm
h = hm + ha         # Total height
e = 1/3 * tau_k             # PM air-gap lenght --> to be defined in term of tau_k
q_femm = tau_k/2


# Problem dimensions in [mm]
tau_k = 60e-3           # half period of the problem : 100 mm
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


# Section du câble en mm^2
S = 10
Nt = 10

# Densité de courant en A/mm^2
J = 5

# Courant de la machine
I = J * S



# Winding dimension
Lz = 45e-3                          # Lz est la profondeur de la machine
y_median = hm + ha + 0.1e-3         # y_median est la hauteur du milieu de la bobine

