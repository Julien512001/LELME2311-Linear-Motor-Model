import numpy as np
import matplotlib.pyplot as plt
from Region1 import *
from Region2 import *
from Potentiel import *
#from Flux import *



"""
plt.figure(figsize=(8, 6))
plt.plot(q, Bq1[-1, :], color='b', label='$B_{q1}$')
plt.xlim(-q_femm, q_femm)
plt.xlabel('$q [mm]$')
plt.ylabel('$B_{q1} [T]$')
plt.title('Magnetic field')
plt.grid(True)
plt.legend()

plt.figure(figsize=(8, 6))
plt.plot(q, Bp1[-1, :], color='b', label='$B_{p1}$')
plt.xlim(-q_femm, q_femm)
plt.xlabel('$q [mm]$')
plt.ylabel('$B_{p1} [T]$')
plt.title('Magnetic field')
plt.grid(True)
plt.legend()
"""
"""
plt.figure(figsize=(8, 6))
plt.plot(q, Bp1[-1, :], color='b', label='$B_{p1}$')
plt.plot(q, Bq1[-1, :], color='r', label='$B_{q1}$')
plt.xlim(-q_femm, q_femm)
plt.xlabel('$q [mm]$')
plt.ylabel('$B_{p1}$ and $B_{q1}$')
plt.title('Magnetic field')
plt.grid(True)
plt.legend()
"""
"""
plt.figure(figsize=(8, 6))
plt.plot(q, Bq2[0, :], color='b', label='$B_{q2}$')
plt.xlim(-q_femm, q_femm)
plt.xlabel('$q [mm]$')
plt.ylabel('$B_{q2} [T]$')
plt.title('Magnetic field')
plt.grid(True)
plt.legend()

plt.figure(figsize=(8, 6))
plt.plot(q, Bp2[0, :], color='b', label='$B_{p2}$')
plt.xlim(-q_femm, q_femm)
plt.xlabel('$q [mm]$')
plt.ylabel('$B_{p2} [T]$')
plt.title('Magnetic field')
plt.grid(True)
plt.legend()
"""

plt.figure(figsize=(8, 6))
plt.plot(q, Bp1[0, :], color='b', label='$B_{p1}$')
plt.plot(q, Bq1[0, :], color='r', label='$B_{q1}$')
plt.xlim(-q_femm, q_femm)
plt.xlabel('$q [mm]$')
plt.ylabel('$B_{p1}$ and $B_{q1}$')
plt.title('Magnetic field')
plt.grid(True)
plt.legend()

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