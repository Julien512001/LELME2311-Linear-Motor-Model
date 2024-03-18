import numpy as np


Br = 1.21

Mp = Br/(4*np.pi*10**(-7))
print(Mp)

#L = (4e-3 + 0.5e-3)/2.0
L = 1e-3 + 1e-3 + 2.5e-3 + 1e-3 + 1e-3
d = 5e-3

n_harm = 50
n_mesh = 500


# Region 1 height (magnet)
hm = 1e-3
# Region 2 height (airgap)
ha = 1e-3

h = hm + ha

# Region 1 : mu
mu1 = 1.05
# Region2 : mu
mu2 = 1
# Air : mu
mu0 = 1