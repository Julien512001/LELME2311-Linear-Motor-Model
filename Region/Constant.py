import numpy as np




Mp = 12600
L = 1.0
d = 0.8

n_harm = 100
n_mesh = 500


# Region 1 height (magnet)
hm = 0.1
# Region 2 height (airgap)
ha = 0.2

h = hm + ha

# Region 1 : mu
mu1 = 1.05
# Region2 : mu
mu2 = 1
# Air : mu
mu0 = 1