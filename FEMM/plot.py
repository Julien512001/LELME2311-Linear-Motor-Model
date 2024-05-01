import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

df1 = pd.read_csv('Data/Bp1_0.txt', delimiter='\t')
df2 = pd.read_csv('Data/Bq1_0.txt', delimiter='\t')

df3 = pd.read_csv('Data/Bp1_1.txt', delimiter='\t')
df4 = pd.read_csv('Data/Bq1_1.txt', delimiter='\t')

df5 = pd.read_csv('Data/Bp2_1.txt', delimiter='\t')
df6 = pd.read_csv('Data/Bq2_1.txt', delimiter='\t')

#df1 = pd.read_csv("Position/odometry.txt", delimiter=',')


# Extraction des valeurs des deux colonnes dans des tableaux NumPy

tau_k = 60e-3

Bp1_0 = df1.iloc[:, 0].values.astype(float)
x = np.linspace(-tau_k, tau_k, len(Bp1_0))

Bq1_0 = df2.iloc[:, 1].values.astype(float)
print(Bq1_0)

plt.figure(figsize=(8, 6))
plt.plot(x, Bp1_0, color='b', label='$B_{p1}$')
plt.plot(x, Bq1_0, color='r', label='$B_{q1}$')
plt.xlabel('$q [mm]$')
plt.ylabel('$B_{p1}$ and $B_{q1}$')
plt.title('Magnetic field')
plt.grid(True)
plt.legend()


Bp1_1 = df3.iloc[:, 1].values.astype(float)
x = np.linspace(-tau_k, tau_k, len(Bp1_0))

Bq1_1 = df4.iloc[:, 1].values.astype(float)
print(Bq1_0)

plt.figure(figsize=(8, 6))
plt.plot(x, Bp1_1, color='b', label='$B_{p1}$')
plt.plot(x, Bq1_1, color='r', label='$B_{q1}$')
plt.xlabel('$q [mm]$')
plt.ylabel('$B_{p1}$ and $B_{q1}$')
plt.title('Magnetic field')
plt.grid(True)
plt.legend()

Bp2_1 = df5.iloc[:, 1].values.astype(float)
x = np.linspace(-tau_k, tau_k, len(Bp1_0))

Bq2_1 = df6.iloc[:, 1].values.astype(float)
print(Bq1_0)

plt.figure(figsize=(8, 6))
plt.plot(x, Bp2_1, color='b', label='$B_{p2}$')
plt.plot(x, Bq2_1, color='r', label='$B_{q2}$')
plt.xlabel('$q [mm]$')
plt.ylabel('$B_{p2}$ and $B_{q2}$')
plt.title('Magnetic field')
plt.grid(True)
plt.legend()


plt.show()