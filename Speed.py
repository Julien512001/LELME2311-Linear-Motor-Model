import numpy as np
import matplotlib.pyplot as plt

# Distance à parcourir en mètres
distance = 0.05

# Temps total
t_total = 10
# Intervalle de temps (dt)
dt = 0.01

# Création du vecteur temps
t = np.arange(0, t_total, dt)

# Calcul de l'accélération nécessaire pour atteindre la distance donnée
a_max = 5

# Création du vecteur d'accélération
a = np.zeros_like(t)
a[:int(len(t)/2)] = a_max
a[int(len(t)/2):] = -a_max

# Calcul de la vitesse en intégrant l'accélération
v = np.cumsum(a) * dt

# Correction de la vitesse pour que la distance parcourue soit exactement celle souhaitée
v += (distance - np.sum(v) * dt) / t_total

# Calcul de la position en intégrant la vitesse
x = np.cumsum(v) * dt

# Tracé des profils d'accélération, de vitesse et de position
plt.figure(figsize=(10, 8))
plt.plot(t, a, label="accélération")
plt.title('Profil d\'accélération')
plt.grid(True)
plt.legend()
plt.xlabel('Temps [s]')
plt.ylabel('$a_x [m/s^2]$')

plt.figure(figsize=(10, 8))
plt.plot(t, v, label="vitesse")
plt.title('Profil de vitesse')
plt.grid(True)
plt.legend()
plt.xlabel('Temps [s]')
plt.ylabel('$v_x [m/s]$')

plt.figure(figsize=(10, 8))
plt.plot(t, x, label="position")
plt.title('Profil de position')
plt.grid(True)
plt.legend()
plt.xlabel('Temps [s]')
plt.ylabel('$x [m]$')

plt.show()
