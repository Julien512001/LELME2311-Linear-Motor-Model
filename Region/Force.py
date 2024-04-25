from Potentiel import *
from flux import *
from Constant import *


P = np.linspace(hm, h, n_mesh) # y
Q = np.linspace(-L+d/2, L-d/2, n_mesh) # x
p,q = np.meshgrid(P,Q)


omega_1 = np.pi/L
F_active = 3/2 * omega_1 * phi0 * I