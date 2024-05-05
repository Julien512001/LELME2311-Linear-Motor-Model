import pandas as pd
import matplotlib.pyplot as plt
from LinearMotor import *

# Charger les données depuis le fichier CSV
data = pd.read_csv('Optimization/OptimalMachines.csv')

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


lq, lp = p.get_lqlp()
print("lq = {}".format(lq))
print("lp = {}".format(lp))

plt.show()