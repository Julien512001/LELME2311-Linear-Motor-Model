from Potentiel import *
from Constant import *


P = np.linspace(hm, h, n_mesh) # y
Q = np.linspace(-L+d/2, L-d/2, n_mesh) # x
p,q = np.meshgrid(P,Q)


phi = np.zeros_like(p)

phi0 = 2*L2*Al0
phi_harm = 2*L2*A_harm



