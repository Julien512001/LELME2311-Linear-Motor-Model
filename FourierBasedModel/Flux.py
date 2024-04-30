from Potentiel import *
from Constant import *



phi = np.zeros_like(P2)

phi_fond = Nt*2*Lz*A_fond
phi_harm = Nt*2*Lz*A_harm

phi = phi_fond + phi_harm

print(phi.shape)
