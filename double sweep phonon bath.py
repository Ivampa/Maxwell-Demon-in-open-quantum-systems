"""
    CODE MUST BE CLEANED: NOT A FINAL VERSION!!!!

Simulation of an open quantum system formed by a serial double quantum dot
capacitively coupled to a lateral double dot (A–G) and to a phonon bath.

The code performs a double sweep over the chemical potentials μ_L and μ_R
for different values of E₂, building the total Liouvillian (electronic + phononic),
computing the steady-state populations, and evaluating:

- Electronic currents (I_L, I_R)
- Electronic heat currents (J_L, J_R)
- Phonon heat current (J_ph)
- Entropy flows (Ṡ_L, Ṡ_R, Ṡ_ph)
- Work and efficiency (engine / refrigerator / heat pump regimes)

The aim is to study non-equilibrium thermodynamics, entropy production
and performance regimes (engine, refrigerator and heat pump) in a
quantum-dot-based Maxwell-demon-like setup.
"""

import numpy as np
from scipy.linalg import eigh, null_space
import matplotlib.pyplot as plt

# --- 1. Constantes ---
# --- 1. Constantes ---
Kb=1
kB=1



E1 = 1.2

EA = 0.2993
EG = 0.2994
EC = 0.35
Gamma_12 = 0.2
Gamma_AG = 0.05
T_L = 1.0
T_R = 1.0

phi=1

Te_ph=0.1


E2=E1


mu_L=0
mu_R=0






# --- 2. Hamiltoniano en base local ---
H = np.array([
    [EA + E1,        Gamma_AG,       Gamma_12,    0,            0,          0],
    [Gamma_AG,       EG + E1,          0,         Gamma_12,     0,          0],
    [Gamma_12,       0,             EA+E2+EC,     Gamma_AG,     0,          0],
    [0,              Gamma_12,      Gamma_AG,     EG + E2,      0,          0],
    [0,              0,               0,          0,            EA,    Gamma_AG],
    [0,              0,               0,          0,      Gamma_AG,       EG]
], dtype=complex)

# --- 3. Diagonalización ---
energies, eigenvectors = eigh(H)
eigenvectors = eigenvectors / np.linalg.norm(eigenvectors, axis=0)

# --- 4. Reordenamiento ---
indice_0mas, indice_0menos = None, None
for i in range(6):
    if np.all(np.abs(eigenvectors[:4, i]) < 1e-10):
        if indice_0mas is None:
            indice_0mas = i
        else:
            indice_0menos = i

if energies[indice_0mas] < energies[indice_0menos]:
    indice_0mas, indice_0menos = indice_0menos, indice_0mas

indices_ABCD = [i for i in range(6) if i not in [indice_0mas, indice_0menos]]
orden_ABCD = sorted(indices_ABCD, key=lambda i: energies[i])
orden_final = orden_ABCD + [indice_0mas, indice_0menos]

eigenvectors = eigenvectors[:, orden_final]
energies = energies[orden_final]




# --- 9. Matrices de transición electrónicas ---

def fermi(E, mu, T):
    return 1 / (1 + np.exp((E - mu) / (kB*T)))

def fermibarrido(E, mu, T):
    return 1 / (1 + np.exp((E - mu-mu_R) / (kB*T)))


def n_ph(delta):
    
        exp_arg =(float(np.real(delta)) / (kB * Te_ph))
        return 1 / (np.exp(exp_arg) - 1)

P_L, P_R = np.zeros((6, 6)), np.zeros((6, 6))
fL, fR = np.zeros((6, 6)), np.zeros((6, 6))

# --- P_L ---

P_L[0, 4] = np.abs(np.conj(eigenvectors[0, 0]) * eigenvectors[4, 4] + np.conj(eigenvectors[1, 0]) * eigenvectors[5, 4])**2
P_L[0, 5] = np.abs(np.conj(eigenvectors[0, 0]) * eigenvectors[4, 5] + np.conj(eigenvectors[1, 0]) * eigenvectors[5, 5])**2

P_L[1, 4] = np.abs(np.conj(eigenvectors[0, 1]) * eigenvectors[4, 4] + np.conj(eigenvectors[1, 1]) * eigenvectors[5, 4])**2
P_L[1, 5] = np.abs(np.conj(eigenvectors[0, 1]) * eigenvectors[4, 5] + np.conj(eigenvectors[1, 1]) * eigenvectors[5, 5])**2

P_L[2, 4] = np.abs(np.conj(eigenvectors[0, 2]) * eigenvectors[4, 4] + np.conj(eigenvectors[1, 2]) * eigenvectors[5, 4])**2
P_L[2, 5] = np.abs(np.conj(eigenvectors[0, 2]) * eigenvectors[4, 5] + np.conj(eigenvectors[1, 2]) * eigenvectors[5, 5])**2

P_L[3, 4] = np.abs(np.conj(eigenvectors[0, 3]) * eigenvectors[4, 4] + np.conj(eigenvectors[1, 3]) * eigenvectors[5, 4])**2
P_L[3, 5] = np.abs(np.conj(eigenvectors[0, 3]) * eigenvectors[4, 5] + np.conj(eigenvectors[1, 3]) * eigenvectors[5, 5])**2

P_L[4, 0] = np.abs(np.conj(eigenvectors[4, 4]) * eigenvectors[0, 0] + np.conj(eigenvectors[5, 4]) * eigenvectors[1, 0])**2
P_L[4, 1] = np.abs(np.conj(eigenvectors[4, 4]) * eigenvectors[0, 1] + np.conj(eigenvectors[5, 4]) * eigenvectors[1, 1])**2
P_L[4, 2] = np.abs(np.conj(eigenvectors[4, 4]) * eigenvectors[0, 2] + np.conj(eigenvectors[5, 4]) * eigenvectors[1, 2])**2
P_L[4, 3] = np.abs(np.conj(eigenvectors[4, 4]) * eigenvectors[0, 3] + np.conj(eigenvectors[5, 4]) * eigenvectors[1, 3])**2

P_L[5, 0] = np.abs(np.conj(eigenvectors[4, 5]) * eigenvectors[0, 0] + np.conj(eigenvectors[5, 5]) * eigenvectors[1, 0])**2
P_L[5, 1] = np.abs(np.conj(eigenvectors[4, 5]) * eigenvectors[0, 1] + np.conj(eigenvectors[5, 5]) * eigenvectors[1, 1])**2
P_L[5, 2] = np.abs(np.conj(eigenvectors[4, 5]) * eigenvectors[0, 2] + np.conj(eigenvectors[5, 5]) * eigenvectors[1, 2])**2
P_L[5, 3] = np.abs(np.conj(eigenvectors[4, 5]) * eigenvectors[0, 3] + np.conj(eigenvectors[5, 5]) * eigenvectors[1, 3])**2

# --- P_R ---

P_R[0, 4] = np.abs(np.conj(eigenvectors[3, 0]) * eigenvectors[5, 4])**2
P_R[0, 5] = np.abs(np.conj(eigenvectors[3, 0]) * eigenvectors[5, 5])**2

P_R[1, 4] = np.abs(np.conj(eigenvectors[3, 1]) * eigenvectors[5, 4])**2
P_R[1, 5] = np.abs(np.conj(eigenvectors[3, 1]) * eigenvectors[5, 5])**2

P_R[2, 4] = np.abs(np.conj(eigenvectors[3, 2]) * eigenvectors[5, 4])**2
P_R[2, 5] = np.abs(np.conj(eigenvectors[3, 2]) * eigenvectors[5, 5])**2

P_R[3, 4] = np.abs(np.conj(eigenvectors[3, 3]) * eigenvectors[5, 4])**2
P_R[3, 5] = np.abs(np.conj(eigenvectors[3, 3]) * eigenvectors[5, 5])**2

P_R[4, 0] = np.abs(np.conj(eigenvectors[5, 4]) * eigenvectors[3, 0])**2
P_R[4, 1] = np.abs(np.conj(eigenvectors[5, 4]) * eigenvectors[3, 1])**2
P_R[4, 2] = np.abs(np.conj(eigenvectors[5, 4]) * eigenvectors[3, 2])**2
P_R[4, 3] = np.abs(np.conj(eigenvectors[5, 4]) * eigenvectors[3, 3])**2

P_R[5, 0] = np.abs(np.conj(eigenvectors[5, 5]) * eigenvectors[3, 0])**2
P_R[5, 1] = np.abs(np.conj(eigenvectors[5, 5]) * eigenvectors[3, 1])**2
P_R[5, 2] = np.abs(np.conj(eigenvectors[5, 5]) * eigenvectors[3, 2])**2
P_R[5, 3] = np.abs(np.conj(eigenvectors[5, 5]) * eigenvectors[3, 3])**2

for i in range(6):
    fL[i, 4] = fermi(energies[i] - energies[4], mu_L, T_L)#aquí es donde tienes que cambiar las transiciones
    fL[i, 5] = fermi(energies[i] - energies[5], mu_L, T_L)#si
    fR[i, 4] = fermi(energies[i] - energies[4], mu_R, T_R)#aqui no
    fR[i, 5] = fermi(energies[i] - energies[5], mu_R, T_R)#no

    fL[4, i] = 1 - fL[i, 4]#si
    fL[5, i] = 1 - fL[i, 5]#si
    fR[4, i] = 1 - fR[i, 4]#no
    fR[5, i] = 1 - fR[i, 5]#no

Gamma = P_L * fL + P_R * fR
GammaL= P_L * fL
GammaR= P_R * fR




M = np.array([
    [-(Gamma[4,0] + Gamma[5,0]),     0,                     0,                     0,       Gamma[0,4], Gamma[0,5]],
    [0, -(Gamma[4,1] + Gamma[5,1]),  0,                     0,       Gamma[1,4], Gamma[1,5]],
    [0, 0, -(Gamma[4,2] + Gamma[5,2]), 0,       Gamma[2,4], Gamma[2,5]],
    [0, 0, 0, -(Gamma[4,3] + Gamma[5,3]),       Gamma[3,4], Gamma[3,5]],
    [Gamma[4,0], Gamma[4,1], Gamma[4,2], Gamma[4,3], -(Gamma[0,4] + Gamma[1,4] + Gamma[2,4] + Gamma[3,4]), 0],
    [Gamma[5,0], Gamma[5,1], Gamma[5,2], Gamma[5,3], 0, -(Gamma[0,5] + Gamma[1,5] + Gamma[2,5] + Gamma[3,5])]
], dtype=float)

# M[j, i] = tasa de transición de i → j (por túnel electrónico)
# M[i, i] -= suma de tasas de salida de i



# --- 5. Ángulos de mezcla ---
theta1 = np.arctan(2 * Gamma_AG**2 / (EA - EG))
theta2 = np.arctan(2 * Gamma_AG**2 / (EA - EG + EC))
theta0 = np.arctan(2 * Gamma_AG**2 / (EA - EG))

# --- 6. Definición de estados mezclados ---
ket_1p = [np.cos(theta1 / 2),  np.sin(theta1 / 2)]
ket_1m = [-np.sin(theta1 / 2), np.cos(theta1 / 2)]

ket_2p = [np.cos(theta2 / 2),  np.sin(theta2 / 2)]
ket_2m = [-np.sin(theta2 / 2), np.cos(theta2 / 2)]

ket_0p = [np.cos(theta0 / 2),  np.sin(theta0 / 2)]
ket_0m = [-np.sin(theta0 / 2), np.cos(theta0 / 2)]

# --- 7. Matriz de cambio de base: local → (1±, 2±, 0±) ---
# Orden de la base local: [1A, 1G, 2A, 2G, 0A, 0G]
# Orden nueva base:       [1+, 1-, 2+, 2-, 0+, 0-]

transform_matrix = np.zeros((6, 6), dtype=complex)

# 1A y 1G → 1±
transform_matrix[0, 0] = ket_1p[0]  # 1A → 1+
transform_matrix[0, 1] = ket_1m[0]  # 1A → 1-
transform_matrix[1, 0] = ket_1p[1]  # 1G → 1+
transform_matrix[1, 1] = ket_1m[1]  # 1G → 1-

# 2A y 2G → 2±
transform_matrix[2, 2] = ket_2p[0]
transform_matrix[2, 3] = ket_2m[0]
transform_matrix[3, 2] = ket_2p[1]
transform_matrix[3, 3] = ket_2m[1]

# 0A y 0G → 0±
transform_matrix[4, 4] = ket_0p[0]
transform_matrix[4, 5] = ket_0m[0]
transform_matrix[5, 4] = ket_0p[1]
transform_matrix[5, 5] = ket_0m[1]

# transform_matrix = transform_matrix.T

# --- 7. Representación de autovectores ABCD en nueva base ---
autovectores_transformados = transform_matrix @ eigenvectors

# --- 8. Impresión por consola ---
base_labels = ['|1+>', '|1->', '|2+>', '|2->', '|0,+>', '|0,->']
print("Autovectores ABCD en la nueva base (1±, 2±, 0A, 0G):\n")
for idx in range(6):
    print(f"Estado {chr(65+idx)}:")
    for j in range(6):
        val = autovectores_transformados[j, idx]
        print(f"  {base_labels[j]}: {val.real:.4f}")
    print()

epsilon = 0  # umbral para filtrar transiciones pequeñas
T = np.zeros((6, 6), dtype=float)

# --- 1. Operador fonónico O_fonon en base (1±,2±,0±) ---
O_fonon = np.zeros((6, 6), dtype=complex)
O_fonon[0, 1] = O_fonon[1, 0] = 1  # 1+ <-> 1-
O_fonon[2, 3] = O_fonon[3, 2] = 1  # 2+ <-> 2-
O_fonon[4, 5] = O_fonon[5, 4] = 1 # 0+ <-> 0-

# --- 2. Proyección del operador fonónico a la base de autoestados ABCDFG ---
# (usando que eigenvectors está expresado en base local, y transform_matrix es base nueva en local)
O_ph = eigenvectors.conj().T @ transform_matrix.conj().T @ O_fonon @ transform_matrix @ eigenvectors
# O_ph = autovectores_transformados.conj().T @ O_fonon @ autovectores_transformados
O_ph=O_ph.T
T_ph = np.abs(O_ph)**2

# --- 3. Construcción de la matriz T de transiciones fonónicas ---
for i in range(6):
    for j in range(6):
        if i != j and T_ph[i, j] > epsilon:
            deltaE = abs(energies[j] - energies[i])
            if energies[j] > energies[i]:
                T[j, i] = phi*T_ph[j, i] * n_ph(deltaE)          # absorción
                
            else:
                T[j, i] = phi*T_ph[j, i] * (1 + n_ph(deltaE))    # emisión
                

for i in range(6):
    T[i,i] = -np.sum(T[:,i]) 

                
                
                
                
                


# --- 4. Impresión de resultados ---
print("\nMatriz de transiciones fonónicas (base ABCDFG):")
for i in range(6):
    for j in range(6):
        print(f"T[{i},{j}] = {T[i,j]:.6f}")
    print()







# --- Preparamos el barrido --- 
mu_L_values = np.linspace(2, -2, 200)  # mu_L baja
mu_R_values = np.linspace(-2, 2, 200)  # mu_R sube

I_L_values = []
I_R_values = []
JL_values = []
JR_values = []
Jph_values = []
mu_diff_values = []  # Para guardar mu = mu_L - mu_R en cada paso
WL_values = []
WR_values = []


# --- Loop sobre mu_L y mu_R pareados ---
for mu_L, mu_R in zip(mu_L_values, mu_R_values):
    mu = mu_L - mu_R
    mu_diff_values.append(mu)

    # recalcular matrices de transiciones fL y fR
    fL = np.zeros((6, 6))
    fR = np.zeros((6, 6))
    for i in range(6):
        fL[i, 4] = fermi(energies[i] - energies[4], mu_L, T_L)
        fL[i, 5] = fermi(energies[i] - energies[5], mu_L, T_L)
        fR[i, 4] = fermi(energies[i] - energies[4], mu_R, T_R)
        fR[i, 5] = fermi(energies[i] - energies[5], mu_R, T_R)

        fL[4, i] = 1 - fL[i, 4]
        fL[5, i] = 1 - fL[i, 5]
        fR[4, i] = 1 - fR[i, 4]
        fR[5, i] = 1 - fR[i, 5]

    Gamma = P_L * fL + P_R * fR
    GammaL = P_L * fL
    GammaR = P_R * fR

    M = np.array([
        [-(Gamma[4,0] + Gamma[5,0]),     0,                     0,                     0,       Gamma[0,4], Gamma[0,5]],
        [0, -(Gamma[4,1] + Gamma[5,1]),  0,                     0,       Gamma[1,4], Gamma[1,5]],
        [0, 0, -(Gamma[4,2] + Gamma[5,2]), 0,       Gamma[2,4], Gamma[2,5]],
        [0, 0, 0, -(Gamma[4,3] + Gamma[5,3]),       Gamma[3,4], Gamma[3,5]],
        [Gamma[4,0], Gamma[4,1], Gamma[4,2], Gamma[4,3], -(Gamma[0,4] + Gamma[1,4] + Gamma[2,4] + Gamma[3,4]), 0],
        [Gamma[5,0], Gamma[5,1], Gamma[5,2], Gamma[5,3], 0, -(Gamma[0,5] + Gamma[1,5] + Gamma[2,5] + Gamma[3,5])]
    ], dtype=float)

    # --- Matriz total de transiciones ---
    P = M + T
    # --- Expansión de P a 36x36 ---
    P_extendida = np.zeros((36, 36), dtype=float)
    P_extendida[:6, :6] = P
    
    # --- Estado estacionario ---
    nullvecs = null_space(P)
    rho_ss = nullvecs[:, 0].real
    rho_ss /= np.sum(rho_ss)
    print("rho_ss =", rho_ss)
    evals = np.linalg.eigvals(P_extendida)
    print("Autovalores de L_total:", np.sort(np.real(evals)))


    # # --- Estado estacionario ---
    # nullvecs = null_space(P)
    # rho_ss = nullvecs[:, 0].real
    # rho_ss /= np.sum(rho_ss)
    # evals = np.linalg.eigvals(P)
    # print("Autovalores de L_total:", np.sort(np.real(evals)))

    # --- Cálculo de corrientes electrónicas ---
    I_L_ss, I_R_ss = 0, 0
    for X in range(4):
        I_L_ss += -GammaL[X,4] * rho_ss[4]
        I_L_ss += -GammaL[X,5] * rho_ss[5]
        I_L_ss +=  GammaL[4,X] * rho_ss[X]
        I_L_ss +=  GammaL[5,X] * rho_ss[X]

    for X in range(4):
        I_R_ss += -GammaR[X,4] * rho_ss[4]
        I_R_ss += -GammaR[X,5] * rho_ss[5]
        I_R_ss +=  GammaR[4,X] * rho_ss[X]
        I_R_ss +=  GammaR[5,X] * rho_ss[X]

    I_L_ss *= -1
    I_R_ss *= -1
    I_L_values.append(I_L_ss)
    I_R_values.append(I_R_ss)
    ## --- Cálculo de potencias W_L y W_R ---
    W_L_ss, W_R_ss = 0, 0

    for X in range(4):  # índices A,B,C,D
     for W in range(4,6):  # índices F,G
        W_L_ss += mu_L * (GammaL[W,X] * rho_ss[X] - GammaL[X,W] * rho_ss[W])
        W_R_ss += mu_R * (GammaR[W,X] * rho_ss[X] - GammaR[X,W] * rho_ss[W])

     # Guardar las potencias
    WL_values.append(W_L_ss)
    WR_values.append(W_R_ss)
    
   

    


    # --- Corrientes de calor electrónicas ---
    J_L_ss, J_R_ss = 0, 0
    for X in range(4):
        for W in range(4,6):
            J_L_ss += (GammaL[W,X] * rho_ss[X] - GammaL[X,W] * rho_ss[W]) * (energies[X] - energies[W] - mu_L)
            J_R_ss += (GammaR[W,X] * rho_ss[X] - GammaR[X,W] * rho_ss[W]) * (energies[X] - energies[W] - mu_R)

    JL_values.append(J_L_ss)
    JR_values.append(J_R_ss)

    # --- Corriente de calor fonónica ---
    Jph = 0
    for j in range(6):
        for i in range(6):
            if i != j:
                Jph += (T[i,j] * rho_ss[j] - T[j,i] * rho_ss[i]) * (energies[j] - energies[i])
    Jph_values.append(Jph/2)



# --- Gráfica de corrientes electrónicas I_L y I_R ---
plt.figure(figsize=(7,4))

# Asegurar que los arrays sean NumPy para operaciones vectoriales
mu_array = np.array(mu_L_values) - np.array(mu_R_values)
IL_array = np.array(I_L_values)
IR_array = np.array(I_R_values)

plt.plot(mu_diff_values, IL_array, label='I_L (q units)', color='purple')
plt.plot(mu_diff_values, IR_array, label='I_R (q units)', color='red')

# Máscara para sombreado donde hay potencia útil
mask = (IL_array > 0) & (mu_array < 0)
plt.fill_between(mu_array, 0, IL_array, where=mask, color='green', alpha=0.3)

plt.xlabel(r'$\mu = (\mu_L - \mu_R)/k_BT$')
plt.ylabel(r'Corriente $I/q$ [1/tiempo]')

plt.title('Corrientes electrónicas vs diferencia de potencial (adimensional)')
plt.grid(True, linestyle='--', linewidth=0.5)
plt.legend()

# Recuadro discreto con parámetros del sistema dentro del gráfico
param_text = (
    r"$E_1=1$, $E_2=1$, $E_C=0.35$"
    r"$\Gamma_{12}=0.2$, $\Gamma_{AG}=0.05$"
    r"$T_L=1.0$, $T_R=1.0$"
)
plt.gca().text(0.02, 0.95, param_text, transform=plt.gca().transAxes,
               fontsize=9, va='top', ha='left',
               bbox=dict(boxstyle='round,pad=0.4', facecolor='whitesmoke', edgecolor='gray', alpha=0.8))


plt.tight_layout()
plt.axhline(0, color='black', linewidth=1.5)
plt.show()



# --- Calcular la suma de corrientes de calor ---
J_total_values = -np.array(JL_values) - np.array(JR_values) - np.array(Jph_values)
W_total_values = np.array(WL_values) + np.array(WR_values)

# --- Comprobar si existe producción de trabajo ---
if np.any(J_total_values > 0):
    print("Existe una región donde se produce trabajo (J_total > 0).")
else:
    print("No hay regiones donde se produzca trabajo (J_total <= 0 en todo el dominio).")

# --- Gráfica de corrientes de calor ---
plt.figure(figsize=(7,4))
plt.plot(mu_diff_values, JL_values, label='J_L', color='blue')
plt.plot(mu_diff_values, JR_values, label='J_R', color='orange')
plt.plot(mu_diff_values, Jph_values, label='J_ph', color='green')
plt.plot(mu_diff_values, J_total_values, label='J_total', color='black', linestyle='--')  # curva suma

# --- Sombrear donde J_total > 0 ---
mu_array = np.array(mu_diff_values)
Jtot_array = np.array(J_total_values)
mask = (Jtot_array > 0)
plt.fill_between(mu_array, 0, Jtot_array, where=mask, color='orange', alpha=0.3)

plt.xlabel('μ = μ_L - μ_R')
plt.ylabel('Corriente de calor')
plt.title('Corrientes de calor vs diferencia de potencial en baño +-')
plt.grid(True, linestyle='--', linewidth=0.5)
plt.legend()
plt.tight_layout()
plt.axhline(0, color='black', linewidth=1.5)  #

plt.show()


# --- Gráfica de potencias W_L, W_R y W_total ---
plt.figure(figsize=(7,4))
plt.plot(mu_diff_values, WL_values, label='W_L', color='blue')
plt.plot(mu_diff_values, WR_values, label='W_R', color='orange')
plt.plot(mu_diff_values, W_total_values, label='W_total', color='black', linestyle='--')

plt.xlabel('μ = μ_L - μ_R')
plt.ylabel('Potencia')
plt.title('Potencias vs diferencia de potencial')
plt.grid(True, linestyle='--', linewidth=0.5)
plt.legend()
plt.tight_layout()
plt.axhline(0, color='black', linewidth=1.5)

plt.show()

# --- Buscar dónde la potencia total es positiva ---
mu_potencia_positiva = [mu for mu, W in zip(mu_diff_values, W_total_values) if W > 0]

if mu_potencia_positiva:
    print("\nPotencia total positiva en los siguientes valores de mu:")
    print(mu_potencia_positiva)
else:
    print("\nNo hay regiones donde la potencia total sea positiva.")



# Convertir a arrays
JL = np.array(JL_values)
JR = np.array(JR_values)
W_total = np.array(W_total_values)
mu_diff = np.array(mu_diff_values)

# # Suma que debería igualar J_ph
# J_check = -JL - JR - W_total

# # Gráfica
# plt.figure(figsize=(8, 5))
# plt.plot(mu_diff, JL, label='J_L', color='blue')
# plt.plot(mu_diff, JR, label='J_R', color='orange')
# plt.plot(mu_diff, W_total, label='W_total', color='black', linestyle='--')
# plt.plot(mu_diff, J_check, label='J_L + J_R + W_total', color='green')



# plt.xlabel('μ = μ_L - μ_R')
# plt.ylabel('Energía (corriente de calor / potencia)')
# plt.title('Chequeo de conservación de energía ')
# plt.grid(True, linestyle='--', linewidth=0.5)
# plt.axhline(0, color='black', linewidth=1.2)
# plt.legend()
# plt.tight_layout()
# plt.show()


# --- Cálculo de flujos de entropía ---
Sdot_L = np.array(JL_values) / T_L
Sdot_R = np.array(JR_values) / T_R
Sdot_ph = np.array(Jph_values) / Te_ph
Sdot_total = Sdot_L + Sdot_R + Sdot_ph
Sdot_sumLR = Sdot_L + Sdot_R

# --- Gráfica de flujos de entropía ---
plt.figure(figsize=(8, 5))
plt.plot(mu_diff_values, Sdot_L, label=r'$\dot{S}_L$', color='hotpink')
plt.plot(mu_diff_values, Sdot_R, label=r'$\dot{S}_R$', color='mediumorchid')
plt.plot(mu_diff_values, Sdot_sumLR, label=r'$\dot{S}_L + \dot{S}_R$', color='deepskyblue')
plt.plot(mu_diff_values, Sdot_ph, label=r'$\dot{S}_{ph}$', color='green')
plt.plot(mu_diff_values, Sdot_total, label=r'$\dot{S}_{total}$', color='black', linestyle='--')

plt.xlabel(r'$\mu = \mu_L - \mu_R$')
plt.ylabel(r'Flujo de entropía $\dot{S}$ / Kb')
plt.title('Flujos de entropía ')
plt.grid(True, linestyle='--', linewidth=0.5)
plt.axhline(0, color='black', linewidth=1.2)
plt.legend()
plt.tight_layout()
plt.show()





# --- Conversión de unidades a kB*T (T referencial: T_L o T_R) ---
T_ref = T_L  # puedes ajustar si prefieres usar T_R u otro

# Normalizar todas las cantidades relevantes
mu_diff_values = np.array(mu_diff_values) / (kB * T_ref)
mu_L_values = np.array(mu_L_values) / (kB * T_ref)
mu_R_values = np.array(mu_R_values) / (kB * T_ref)

I_L_values = np.array(I_L_values)  # Corrientes electrónicas ya están en unidades arbitrarias (q)
I_R_values = np.array(I_R_values)

JL_values = np.array(JL_values) / (kB * T_ref)
JR_values = np.array(JR_values) / (kB * T_ref)
Jph_values = np.array(Jph_values) / (kB * T_ref)

WL_values = np.array(WL_values) / (kB * T_ref)
WR_values = np.array(WR_values) / (kB * T_ref)

# --- Gráfica de corrientes electrónicas I_L y I_R ---
plt.figure(figsize=(7,4))
plt.plot(mu_diff_values, I_L_values, label='I_L (q units)', color='purple')
plt.plot(mu_diff_values, I_R_values, label='I_R (q units)', color='red')
mu_array = mu_L_values - mu_R_values
IL_array = I_L_values
mask = (IL_array > 0) & (mu_array < 0)
plt.fill_between(mu_array, 0, IL_array, where=mask, color='green', alpha=0.3)
plt.xlabel(r'$\mu = (\mu_L - \mu_R)/k_BT$')
plt.ylabel('Corriente I/q')
plt.title('Corrientes electrónicas vs diferencia de potencial (adimensional)')
plt.grid(True, linestyle='--', linewidth=0.5)
plt.legend()
plt.tight_layout()
plt.axhline(0, color='black', linewidth=1.5)
plt.show()

# --- Gráfica de corrientes de calor ---
plt.figure(figsize=(7,4))
plt.plot(mu_diff_values, JL_values, label='J_L', color='blue')
plt.plot(mu_diff_values, JR_values, label='J_R', color='orange')
plt.plot(mu_diff_values, Jph_values, label='J_ph', color='green')
plt.xlabel(r'$\mu = (\mu_L - \mu_R)/k_BT$')
plt.ylabel(r'Corriente de calor $J / (k_B T)$')
plt.title('Corrientes de calor (normalizadas) vs diferencia de potencial')
plt.grid(True, linestyle='--', linewidth=0.5)
plt.legend()
plt.tight_layout()
plt.axhline(0, color='black', linewidth=1.5)
plt.show()

# --- Cálculo de J_total y trabajo
J_total_values = -JL_values - JR_values - Jph_values
W_total_values = WL_values + WR_values

if np.any(J_total_values > 0):
    print("Existe una región donde se produce trabajo (J_total > 0).")
else:
    print("No hay regiones donde se produzca trabajo (J_total <= 0 en todo el dominio).")

plt.figure(figsize=(7,4))

# Aplicar máscara de rango
mask_range = (mu_diff_values >= -0.5) & (mu_diff_values <= 0.25)

# Filtrar los datos
mu_filtered = mu_diff_values[mask_range]
JL_filtered = JL_values[mask_range]
JR_filtered = JR_values[mask_range]
Jph_filtered = Jph_values[mask_range]
Jtotal_filtered = J_total_values[mask_range]

# Máscara de valores positivos para sombreado
mask_positive = Jtotal_filtered > 0

# Gráficas
plt.plot(mu_filtered, JL_filtered, label='J_L', color='blue')
plt.plot(mu_filtered, JR_filtered, label='J_R', color='orange')
plt.plot(mu_filtered, Jph_filtered, label='J_ph', color='green')
plt.plot(mu_filtered, Jtotal_filtered, label='J_total', color='black', linestyle='--')
plt.fill_between(mu_filtered, 0, Jtotal_filtered, where=mask_positive, color='orange', alpha=0.3)

# Estética
plt.xlabel(r'$\mu = (\mu_L - \mu_R)/k_BT$')
plt.ylabel(r'Corriente de calor $J / (k_B T)$ [1/tiempo]')
plt.title('Corriente de calor total vs diferencia de potencial')
plt.grid(True, linestyle='--', linewidth=0.5)
plt.legend()
plt.tight_layout()
plt.axhline(0, color='black', linewidth=1.5)
plt.show()

plt.figure(figsize=(7,4))

# Aplicar máscara de rango en mu
mask_range = (mu_diff_values >= -0.5) & (mu_diff_values <= 0.25)

# Filtrar los datos
mu_filtered = mu_diff_values[mask_range]
JL_filtered = JL_values[mask_range]
JR_filtered = JR_values[mask_range]
Jph_filtered = Jph_values[mask_range]
Jtotal_filtered = J_total_values[mask_range]

# Máscara para sombrear donde J_total > 0
mask_positive = Jtotal_filtered > 0

# Graficar
plt.plot(mu_filtered, JL_filtered, label='J_L', color='blue')
plt.plot(mu_filtered, JR_filtered, label='J_R', color='orange')
plt.plot(mu_filtered, Jph_filtered, label='J_ph', color='green')
plt.plot(mu_filtered, Jtotal_filtered, label='J_total', color='black', linestyle='--')
plt.fill_between(mu_filtered, 0, Jtotal_filtered, where=mask_positive, color='orange', alpha=0.3)

# Estética y límites
plt.xlabel(r'$\mu = (\mu_L - \mu_R)/k_BT$')
plt.ylabel(r'Corriente de calor $J / (k_B T)$ [1/tiempo]')
plt.title('Corriente de calor total vs diferencia de potencial')
plt.ylim(-0.003, 0.003)
plt.grid(True, linestyle='--', linewidth=0.5)
plt.axhline(0, color='black', linewidth=1.5)
plt.legend()
plt.tight_layout()
plt.show()

# --- Gráfica de potencias ---
plt.figure(figsize=(7,4))
plt.plot(mu_diff_values, WL_values, label='W_L', color='blue')
plt.plot(mu_diff_values, WR_values, label='W_R', color='orange')
plt.plot(mu_diff_values, W_total_values, label='W_total', color='black', linestyle='--')
plt.xlabel(r'$\mu = (\mu_L - \mu_R)/k_BT$')
plt.ylabel(r'Trabajo $W / (k_B T)$ [1/tiempo]')

plt.title('Flujos de trabajo')
plt.grid(True, linestyle='--', linewidth=0.5)
plt.legend()
plt.tight_layout()
plt.axhline(0, color='black', linewidth=1.5)
plt.show()

mu_potencia_positiva = [mu for mu, W in zip(mu_diff_values, W_total_values) if W > 0]
if mu_potencia_positiva:
    print("\nPotencia total positiva en los siguientes valores de mu (adimensionales):")
    print(mu_potencia_positiva)
else:
    print("\nNo hay regiones donde la potencia total sea positiva.")

# --- Flujos de entropía ---
Sdot_L = JL_values / (T_L / T_ref)
Sdot_R = JR_values / (T_R / T_ref)
Sdot_ph = Jph_values / (Te_ph / T_ref)
Sdot_total = Sdot_L + Sdot_R + Sdot_ph
Sdot_sumLR = Sdot_L + Sdot_R

plt.figure(figsize=(8, 5))
plt.plot(mu_diff_values, Sdot_L, label=r'$\dot{S}_L$', color='hotpink')
plt.plot(mu_diff_values, Sdot_R, label=r'$\dot{S}_R$', color='mediumorchid')
plt.plot(mu_diff_values, Sdot_sumLR, label=r'$\dot{S}_L + \dot{S}_R$', color='deepskyblue')
plt.plot(mu_diff_values, Sdot_ph, label=r'$\dot{S}_{ph}$', color='green')
plt.plot(mu_diff_values, Sdot_total, label=r'$\dot{S}_{total}$', color='black', linestyle='--')
plt.xlabel(r'$\mu = (\mu_L - \mu_R)/k_BT$')
plt.ylabel(r'Flujo de entropía $\dot{S}$ [$k_B$/tiempo]')

plt.title('Flujos de entropía (normalizados)')
plt.grid(True, linestyle='--', linewidth=0.5)
plt.axhline(0, color='black', linewidth=1.2)
plt.legend()
plt.tight_layout()
plt.show()


import matplotlib.pyplot as plt

# Diccionario de parámetros
parametros = {
    'E₁': E1,
    'E₂': E2,
    'E_A': EA,
    'E_G': EG,
    'E_C': EC,
    'Γ₁₂': Gamma_12,
    'Γ_AG': Gamma_AG,
    'T_L': T_L,
    'T_R': T_R,
    'T_ph': Te_ph,
    
    'φ': phi
}

# Crear figura
fig, ax = plt.subplots(figsize=(7, 3))
ax.axis('off')

# Mostrar parámetros como texto
text = '\n'.join([f'{key} = {value}' for key, value in parametros.items()])
plt.text(0.05, 0.95, text, fontsize=13, verticalalignment='top')

plt.title("Parámetros del sistema", fontsize=15, weight='bold')
plt.tight_layout()
plt.show()
