import numpy as np
import matplotlib.pyplot as plt
from Region1 import *
from Region2 import *

P1 = np.linspace(0, hm, n_mesh) # y
Q1 = np.linspace(-L+d/2, L-d/2, n_mesh) # x
p1,q1 = np.meshgrid(P1,Q1)

P2 = np.linspace(hm, h, n_mesh) # y
Q2 = np.linspace(-L+d/2, L-d/2, n_mesh) # x
p2,q2 = np.meshgrid(P2,Q2)
"""
fig = plt.figure()
plt.contour(q1,p1,B1)
plt.colorbar()
ax = plt.gca()
ax.set_aspect(1)

fig = plt.figure()
plt.contour(q2,p2,B2)
plt.colorbar()
ax = plt.gca()
ax.set_aspect(1)
"""

"""
fig = plt.figure()
plt.title("Bq1 interface")
plt.plot(Q, Bq1[:,0])

plt.figure()
plt.title("Bp1 interface")
plt.plot(Q, Bp1[:,0])
"""
"""
plt.figure()
plt.title("fourier")
plt.plot(Q,f)
plt.axvline(-L+d/2, c='r')
plt.axvline(L-d/2, c='r')
"""




plt.figure()
plt.title("Bq2 interface")
plt.plot(Q, Bq2[:,450])

plt.figure()
plt.title("Bp2 interface")
plt.plot(Q, Bp2[:,450])



"""
B1 = np.sqrt(Bq1**2 + Bp1**2)

plt.figure()
plt.title("B1")
plt.plot(Q, B1[:,0])
"""
plt.show()