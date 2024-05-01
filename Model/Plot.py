import numpy as np
import matplotlib.pyplot as plt
from Region1 import *
from Region2 import *
from Potentiel import *
#from Flux import *


fig = plt.figure()
plt.title("Bq1 interface")
plt.xlim(-q_femm, q_femm)
plt.plot(q, Bq1[-1, :])

plt.figure()
plt.title("Bp1 interface")
plt.xlim(-q_femm, q_femm)
plt.plot(q, Bp1[-1, :])

plt.figure()
plt.title("Bq2 interface")
plt.xlim(-q_femm, q_femm)
plt.plot(q, Bq2[-1,:])

plt.figure()
plt.title("Bp2 interface")
plt.xlim(-q_femm, q_femm)
plt.plot(q, Bp2[-1,:])

"""
# Magnetic potential plot
plt.figure()
plt.title("A")
plt.contourf(Q2,P2, Al, levels = 60)
plt.xlim(-q_femm, q_femm)
ax = plt.gca()
#ax.set_aspect(1)
plt.colorbar()

plt.figure()
plt.title("A")
plt.contour(Q2,P2, Al, levels = 60)
plt.xlim(-q_femm, q_femm)
ax = plt.gca()
#ax.set_aspect(1)
plt.colorbar()
"""
"""
# flux plot
plt.figure()
plt.title("$\phi$")
plt.contour(Q2,P2, phi, levels = 30)
#plt.xlim(-q_femm, q_femm)
ax = plt.gca()
#ax.set_aspect(1)
plt.colorbar()
"""
plt.show()