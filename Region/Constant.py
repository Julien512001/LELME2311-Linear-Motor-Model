import numpy as np


Br = 1.21

# Air : mu
mu0 = 4*np.pi*10**(-7)

Mp = Br/mu0
#L = (4e-3 + 0.5e-3)/2.0
print(Mp)
#L = 1e-3 + 1e-3 + 2.5e-3 + 1e-3 + 1e-3
#d = 5e-3
# d/2 = 20 mm


# Définition de L comme étant la longueur d'un aimant plus un airgap
L = 0.042
# d = 0.038
# d/2 = 0.019
d = 2.0*L/3.0
#d = 0.02*2
#L = d + 4e-3/2.0
print(L)
print(d)


n_harm = 50
n_mesh = 500


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
