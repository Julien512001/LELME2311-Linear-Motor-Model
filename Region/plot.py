import numpy as np
import matplotlib.pyplot as plt
from Region1 import *
from Region2 import *


'''
fig = plt.figure()
plt.contourf(q,p,B1)
plt.colorbar()
'''


fig = plt.figure()
plt.title("Bq1 interface")
plt.plot(Q, Bq1[:,0])

plt.figure()
plt.title("Bp1 interface")
plt.plot(Q, Bp1[:,0])

plt.figure()
plt.plot(Q,f)
plt.axvline(-L+d/2, c='r')
plt.axvline(L-d/2, c='r')

'''
print(B2.shape)
fig = plt.figure()
plt.contourf(q,p,B2)
plt.colorbar()
'''

plt.figure()
plt.title("Bq2 interface")
plt.plot(Q, Bq2[:,0])

plt.figure()
plt.title("Bp2 interface")
plt.plot(Q, Bp2[:,0])




B1 = np.sqrt(Bq1**2 + Bp1**2)

plt.figure()
plt.title("B1")
plt.plot(Q, B1[:,0])

plt.show()