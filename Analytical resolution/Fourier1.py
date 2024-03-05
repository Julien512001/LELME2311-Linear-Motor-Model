import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap

plt.rcParams['figure.figsize'] = [8,8]
plt.rcParams.update({'font.size':18})

# Define domain
dx = .001
L = np.pi
x = L * np.arange(-1+dx, dx+1, dx)
n = len(x)
nquart = int(np.floor(n/4))

# Define the hat function
f = np.zeros_like(x)
f[0:nquart] = 1.0
f[nquart:3*nquart] = 0.0
f[3*nquart:4*nquart] = -1.0

fig, ax = plt.subplots()
ax.plot(x,f,'-',color='k',LineWidth = 2)


# We compute Fourier series
name = "Accent"
cmap = get_cmap('tab10')
colors = cmap.colors
ax.set_prop_cycle(color=colors)

A0 = np.sum(f * np.ones_like(x))*dx
fFS = A0/2

n = 50
A = np.zeros(n)
B = np.zeros(n)
for k in range(n):
    A[k] = np.sum(f * np.cos(np.pi*(k+1)*x/L))*dx
    B[k] = np.sum(f * np.sin(np.pi*(k+1)*x/L))*dx
    fFS = fFS + A[k]*np.cos(np.pi*(k+1)*x/L) + B[k]*np.sin(np.pi*(k+1)*x/L)
    if k==n-1:
        ax.plot(x, fFS, '-')
plt.show()






