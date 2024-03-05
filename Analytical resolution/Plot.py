import numpy as np
import matplotlib.pyplot as plt
from B_Field import *

from mpl_toolkits import mplot3d

print(Bp_1.shape)
print(x.shape)
print(y.shape)
'''
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(x, y, Bp_1)


plt.figure()
plt.contourf(x, y, Bq_1)
'''

# q in x direction (tangential)
# p in y direction (perpendicular)

plt.figure()
plt.plot(x, Bq_1)

plt.figure()
plt.plot(Bp_1, y)




plt.show()