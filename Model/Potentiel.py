from SystemResolution import *


# We compute the potential vector field for the air-gap (The second region)
def get_F(x, y):

    Aps = 0
    Apc = 0
    A_harm = 0
    A_fond = 0
    Al = 0

    phi_fond = 0
    phi_harm = 0

    F_active = 0
    F_ripple = 0
    F_THD = 0


    for i in range(1, n_harm+1):
        omega_n = np.pi*i/tau_k

        Mpc = 0
        Mps = -2*Mp/(np.pi*i) * (np.cos(np.pi*i/(tau_k)*e/2) - np.cos(np.pi*i/(tau_k)*(e/2 - tau_k)))

        A = b[4, i-1]/omega_n
        B = b[5, i-1]/omega_n
        C = b[6, i-1]/omega_n
        D = b[7, i-1]/omega_n

        Aps = C*np.exp(y*omega_n) + D*np.exp(-y*omega_n) - mu0*Mpc/omega_n
        Apc = A*np.exp(y*omega_n) + B*np.exp(-y*omega_n) + mu0*Mps/omega_n

        if (i == 1):
            A_fond = Aps*np.sin(omega_n*x) + Apc*np.cos(omega_n*x)

            phi_fond = Nt*2*Lz*A_fond

            F_active = 3/2 * omega_n * phi_fond * I

        else:
            A_int = Aps*np.sin(omega_n*x) + Apc*np.cos(omega_n*x)
            A_harm += A_int

            phi_int = Nt*2*Lz*A_int
            phi_harm += phi_int

            F_int = 3/2 * omega_n * phi_int * I
            F_ripple += F_int
            F_THD += F_int**2


    F = F_active + F_ripple
    THD = np.sqrt(F_THD)/np.abs(F_active)
    return F, F_active, F_ripple, THD

F, F_active, F_ripple, THD = get_F(25e-3, y_median)
print(F_active)
print(THD)

"""
step = 0.001
X = np.arange(-tau_k, tau_k, step)
F = np.zeros_like(X)
F_active = np.zeros_like(X)
F_ripple = np.zeros_like(X)

i = 0
for x in X:
    F[i], F_active[i], F_ripple[i] = get_F(x, y_median)
    i+=1

plt.figure()
#plt.plot(X, F, label = "F")
plt.plot(X, F_active, label="F_active")
#plt.plot(X, F_ripple, label="F_ripple")
plt.legend()

plt.show()
"""