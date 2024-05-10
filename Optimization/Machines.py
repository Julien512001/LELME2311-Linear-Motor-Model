import pandas as pd
import matplotlib.pyplot as plt
from LinearMotor import *

# Charger les données depuis le fichier CSV
data = pd.read_csv('Optimization/optimisation3.csv')

# Extraire les colonnes Force et THD dans des listes
machine = data["Machine"].tolist()
tau_k = data['tau_k'].tolist()
e = data['e'].tolist()
hm = data['hm'].tolist()
ha = data['ha'].tolist()
Lz = data['Lz'].tolist()
lq = data['lq'].tolist()
F_active = data['F_active'].tolist()
THD = data['THD'].tolist()
mass = data["mass"].tolist()
time = data['time'].tolist()

# Tracer la force en fonction de THD
plt.figure(figsize=(10, 6))
plt.scatter(THD, F_active, color='blue', alpha=0.5)
plt.title('Force en fonction de THD')
plt.xlabel('THD')
plt.ylabel('Force')
plt.grid(True)

plt.figure(figsize=(10, 6))
plt.scatter(THD, time, color='blue', alpha=0.5)
plt.title('Force en fonction de THD')
plt.xlabel('THD')
plt.ylabel('Time')
plt.grid(True)

i = 1

p = LinearMotor(tau_k[i], e[i], hm[i], ha[i], Lz[i], lq[i])

print("tau_k = {}".format(tau_k[i]))
print("e = {}".format(e[i]))
lmagnet = tau_k[i] - e[i]
print("lmagnet = {}".format(lmagnet))



F_active = p.get_F_active()
print("F_active = {}".format(F_active))

THD = p.get_THD()
print("THD = {}".format(THD))

mass = p.get_totalMass()
print("mass = {}".format(mass*1000))

vmax = p.get_vmax()
print("vmax = {}".format(vmax))

t = p.get_time()
print("t = {}".format(t))

Sbob = p.get_Sbob()
print("Sbob = {}".format(Sbob))

Uspire = p.get_Uspire()
print("Uspire = {}".format(Uspire))

Nspire = p.get_Nspire()
print("Nspire = {}".format(Nspire))

Sfil = p.get_Sfil()
print("Sfil = {}".format(Sfil))

I = p.get_current()
print("I = {}".format(I))

lspire = p.get_lSpire()
print("lspire = {}".format(lspire))

lbob = p.get_lBobinage()
print("lbob = {}".format(lbob))

opti_I = p.get_optiCurrent()
print("opti_I = {}".format(opti_I))

opti_mass = p.get_optiMass()
print("opti_mass = {}".format(opti_mass))

opti_vmax = p.get_optiVmax()
print("opti_vmax = {}".format(opti_vmax))

opti_time = p.get_optiTime()
print("opti_time = {}".format(opti_time))

lq, lp = p.get_lqlp()
print("lq = {}".format(lq))
print("lp = {}".format(lp))


# Define the function
def dimension(q, e, tau_k):
    profile = np.zeros_like(q)
    profile[(q >= -tau_k) & (q < -tau_k + e/2)] = 0
    profile[(q >= -tau_k + e/2) & (q < -e/2)] = 1
    profile[(q >= -e/2) & (q <= e/2)] = 0
    profile[(q >= e/2) & (q < tau_k - e/2)] = 1
    profile[(q >= tau_k - e/2) & (q <= tau_k)] = 0
    return profile

n_mesh = 500

q = np.linspace(-tau_k[i] - 0.1, tau_k[i] + 0.1, n_mesh)


# Plotting
plt.figure(figsize=(8, 6))
plt.plot(q, dimension(q, e[i], tau_k[i]), color='b', label='Magnetization Profile')
plt.xlabel('q [mm]')
plt.ylabel('$M_p$')
plt.title('Magnetization Profile')
plt.grid(True)
plt.legend()
plt.show()

plt.show()