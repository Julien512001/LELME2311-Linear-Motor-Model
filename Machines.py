import pandas as pd
import matplotlib.pyplot as plt

# Charger les données depuis le fichier CSV
data = pd.read_csv('OptimalMachines.csv')

# Extraire les colonnes Force et THD dans des listes
machine = data["Machine"].tolist()
force_list = data['F_active'].tolist()
thd_list = data['THD'].tolist()
mass = data["mass"].tolist()

# Tracer la force en fonction de THD
plt.figure(figsize=(10, 6))
plt.scatter(thd_list, force_list, color='blue', alpha=0.5)
plt.title('Force en fonction de THD')
plt.xlabel('THD')
plt.ylabel('Force')
plt.grid(True)


# Ajouter le numéro de la machine à côté de chaque point
for i, machine in enumerate(machine):
    plt.text(thd_list[i], force_list[i], str(machine), fontsize=8, va='bottom', ha='right')


plt.show()