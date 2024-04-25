import numpy as np


Br = 1.21

# Air : mu
mu0 = 4*np.pi*10**(-7)

Mp = Br/mu0



# Définition de L comme étant la longueur d'un aimant plus un airgap
L = 0.042
# d = 0.038
# d/2 = 0.019
d = 2.0*L/3.0

n_harm = 50
n_mesh = 500


# Section du câble en mm^2
S = 10

# Densité de courant en A/mm^2
Jl0 = 5

# Courant de la machine
I = Jl0 * S


L2 = 1

# Region 1 height (magnet)
#hm = 0.45
hm = 45e-3
# Region 2 height (airgap)
# ha = 0.15
ha = 15e-3

# h = 0.60
h = hm + ha

# Region 1 : mu
mu1 = 1.05*mu0
# Region2 : mu
mu2 = 1*mu0
