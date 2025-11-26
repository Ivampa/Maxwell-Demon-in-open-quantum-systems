"""
 CODE MUST BE CLEANED: NOT A FINAL VERSION!!!!
 Future implementation in qutip is a WIP

Steady-state analysis of an interacting double quantum dot system
with electronic reservoirs and a phononic bath.

The Hamiltonian is diagonalized and expressed in a mixed eigenbasis
(1±, 2±, 0±), allowing the construction of both electronic and phononic
transition matrices. A Liouvillian formalism is used to compute the
non-equilibrium steady state as a function of the chemical bias μ_L − μ_R.

The code evaluates:
- Population steady state from the null space of the Liouvillian
- Electronic and phononic heat currents
- Power and entropy production
- Regions of work extraction and thermodynamic consistency

The model is inspired by quantum thermodynamic engines and
Maxwell-demon-like behavior in nanoscale quantum systems.
"""




import numpy as np
from scipy.linalg import eigh, null_space
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# --- 1. Constantes ---
Kb=1
kB=1



E1 = 1

EA = 0.2993
EG = 0.2994
EC = 0.35
Gamma_12 = 0.2
Gamma_AG = 0.05
T_L = 1.0
T_R = 1.0
Te_ph = 1
phi=1



# E1 = 1.2

# EA = 0.2994
# EG = 0.2995
# EC = 0.35
# Gamma_12 = 0.35
# Gamma_AG = 0.35
   

# kB = 1.0
# T_L = 4.0
# T_R = 4.0
# Te_ph = 0.1




# # --- Barrido de E2 ---
# E2_values = [
#     E1,
#     E1 - EC / 4,
#     E1 - EC / 2,
#     E1 - 0.75 * EC,
#     E1 - EC,
#     E1 - 1.5 * EC,
#     E1 - 2.5 * EC,
#     E1 - 4 * EC
# ]


# --- Barrido de E2 ---
E2_values = [
    E1,
    E1 - EC / 4,
    E1 - EC / 2,
    E1 - 0.75 * EC,
    E1 - EC,
    E1 - 1.5 * EC,
    E1 - 2.5 * EC,
    E1 - 4 * EC
]

# Para almacenar resultados por cada E2
all_Jph = []
all_JL = []
all_JR = []
all_Sph = []
all_SL = []
all_SR = []
all_SLR = []
all_Winput = []





def fermi(E, mu, T):
    return 1 / (1 + np.exp((E - mu) / (kB*T)))



# --- Bucle externo sobre E2 ---
for E2 in E2_values:
    
    # Redefinir Hamiltoniano con nuevo E2
    H = np.array([
        [EA + E1,        Gamma_AG,       Gamma_12,    0,            0,          0],
        [Gamma_AG,       EG + E1,          0,         Gamma_12,     0,          0],
        [Gamma_12,       0,             EA+E2+EC,     Gamma_AG,     0,          0],
        [0,              Gamma_12,      Gamma_AG,     EG + E2,      0,          0],
        [0,              0,               0,          0,            EA,    Gamma_AG],
        [0,              0,               0,          0,      Gamma_AG,       EG]
    ], dtype=complex)

    # --- [Aquí va tu bloque actual de diagonalización y todo el código hasta el final del barrido mu] ---


    # --- 3. Diagonalización ---
    energies, eigenvectors = eigh(H)
    eigenvectors = eigenvectors / np.linalg.norm(eigenvectors, axis=0)
    
    # --- 4. Reordenamiento de autovectores ---
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
    
    energies = energies[orden_final]
    eigenvectors = eigenvectors[:, orden_final]
    
    # --- 5. Ángulos de mezcla ---
    theta1 = np.arctan(2 * Gamma_AG**2 / (EA - EG))
    theta2 = np.arctan(2 * Gamma_AG**2 / (EA - EG + EC))
    theta0 = np.arctan(2 * Gamma_AG**2 / (EA - EG))
    
    ket_1p = [np.cos(theta1/2),  np.sin(theta1/2)]
    ket_1m = [-np.sin(theta1/2), np.cos(theta1/2)]
    ket_2p = [np.cos(theta2/2),  np.sin(theta2/2)]
    ket_2m = [-np.sin(theta2/2), np.cos(theta2/2)]
    ket_0p = [np.cos(theta0/2),  np.sin(theta0/2)]
    ket_0m = [-np.sin(theta0/2), np.cos(theta0/2)]
    
    # --- 6. Matriz de transformación en base mixta → base (1±,2±,0±) ---
    Tmat = np.zeros((6,6), dtype=complex)
    Tmat[0,0], Tmat[0,1] = ket_1p[0], ket_1m[0]
    Tmat[1,0], Tmat[1,1] = ket_1p[1], ket_1m[1]
    Tmat[2,2], Tmat[2,3] = ket_2p[0], ket_2m[0]
    Tmat[3,2], Tmat[3,3] = ket_2p[1], ket_2m[1]
    Tmat[4,4], Tmat[4,5] = ket_0p[0], ket_0m[0]
    Tmat[5,4], Tmat[5,5] = ket_0p[1], ket_0m[1]
    
    # --- 7. Operador fonónico en base transformada ---
    O_fonon = np.zeros((6, 6), dtype=complex)
    O_fonon[0, 1] = O_fonon[1, 0] = 1
    O_fonon[2, 3] = O_fonon[3, 2] = 1
    O_fonon[4, 5] = O_fonon[5, 4] = 1
    
    O_ph = eigenvectors.conj().T @ Tmat.conj().T @ O_fonon @ Tmat @ eigenvectors
    
    # --- 8. Matriz de tasas fonónicas ---
    def n_ph(E):
        return 1 / (np.exp(E / (kB * Te_ph)) - 1)
    
    W = np.zeros((6, 6), dtype=float)
    for i in range(6):
        for j in range(6):
            if i != j and np.abs(O_ph[i, j]) > 0:
                dE = energies[j] - energies[i]
                if dE > 0:
                    W[i, j] = phi * np.abs(O_ph[i, j])**2 * n_ph(dE)
                else:
                    W[i, j] = phi * np.abs(O_ph[i, j])**2 * (1 + n_ph(-dE)) #aqui sí - porque no es abs
    
    # --- 9. Operadores de salto L_ij ---
    L_ops = []
    for i in range(6):
        for j in range(6):
            if i != j and W[i, j] > 0:
                L = np.zeros((6, 6), dtype=complex)
                L[i, j] = np.sqrt(W[j, i])
                L_ops.append(L)
    
    # --- 10. Superoperador disipativo D ---
    D = np.zeros((36, 36), dtype=complex)
    for m in range(6):
        for n in range(6):
            idx_mn = 6*m + n
            for p in range(6):
                for q in range(6):
                    idx_pq = 6*p + q
                    contrib = 0
                    for L in L_ops:
                        contrib += L[m, p] * np.conj(L[n, q])
                        if n == q:
                            contrib -= 0.5 * sum(np.conj(L[k, m]) * L[k, p] for k in range(6))
                        if m == p:
                            contrib -= 0.5 * sum(np.conj(L[k, n]) * L[k, q] for k in range(6))
                    D[idx_mn, idx_pq] = contrib
    
    # --- 11. Reordenamiento (diagonal primero) ---
    reorder = [6*i + i for i in range(6)] + [6*i + j for i in range(6) for j in range(6) if i != j]
    D_diagcoh = D[reorder, :][:, reorder]
    np.savetxt("D_diagcoh.txt", D_diagcoh.real, fmt="%.6e", delimiter="\t")
    
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
    
    # --- 12. Barrido doble de mu_L y mu_R ---
    
    mu_L_values = np.linspace(2, -2, 200)
    mu_R_values = np.linspace(-2,  2 , 200)
    
    
    
    rho_diag_collection = []
    I_L_values = []
    I_R_values = []
    JL_values = []
    JR_values = []
    Jph_values = []
    mu_diff_values = []  # Para guardar mu = mu_L - mu_R en cada paso
    WL_values = []
    WR_values = []
    e1_values=[]
    e2_values=[]
    

    
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
        P = M 
        #+ T
        # --- Expansión de P a 36x36 ---
        P_extendida = np.zeros((36, 36), dtype=float)
        P_extendida[:6, :6] = P
        M_extendida = np.zeros((36, 36), dtype=float)
        M_extendida[:6, :6] = M
        
        # --- Estado estacionario ---
        # nullvecs = null_space(P_extendida)
        # rho_ss = nullvecs[:, 0].real
        # rho_ss /= np.sum(rho_ss)
        # evals = np.linalg.eigvals(P_extendida)
        # print("Autovalores de L_total:", np.sort(np.real(evals)))
    
    
        
        
    
    
        
    
        # --- Superoperador total ---
        L_total = P_extendida + D_diagcoh
        L_reducida= np.zeros((8, 8), dtype=float)
        L_reducida[:8, :8] = L_total[:8, :8]
            # --- Imprimir todas las filas de L_total ---
        # print("Filas de la matriz L_total:")
        # for i, fila in enumerate(L_total):
        #     print(f"Fila {i}: {fila}")
        # print("Fila 0:", L_total[0])
        # print("Fila 1:", L_total[1])
        # print("Fila 2:", L_total[2])
        # print("Fila 3:", L_total[3])
        # print("Fila 4:", L_total[4])
        # print("Fila 5:", L_total[5])
        # print("Fila 6:", L_total[6])
        # print("Fila 7:", L_total[7])
        # print("Fila 8:", L_total[8])
        # print("Fila 9:", L_total[9])
    
    
    
        # # --- Cálculo del estado estacionario ---
        # rho_ss = null_space(L_total)
        # rho_ss = rho_ss[:, 0].real
        # rho_ss /= np.sum(rho_ss[:6])  # normalización de las poblaciones
        # coherencias_total = np.sum(np.abs(rho_ss[6:]))
        
        # --- Cálculo del estado estacionario ---
        L_duro=M_extendida + D_diagcoh
        nullvecs = null_space(L_duro)
        rho_ss = nullvecs[:, 0].real
        print("rho_ss =", rho_ss)
        rho_ss /= np.sum(rho_ss)  # normalización de las poblaciones
        
    
        coherencias_total = np.sum(np.abs(rho_ss[6:]))
        
        # print("Norma de coherencias:", coherencias_total)
        # evals = np.linalg.eigvals(L_total)
        # print("Autovalores de L_total:", np.sort(np.real(evals)))
        
        
    
    
        # --- Guardar resultado (solo parte diagonal si quieres)
        rho_diag_collection.append(rho_ss[:6])
        
        
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
        
        # --- Cálculo de potencias W_L y W_R ---
        W_L_ss, W_R_ss = 0, 0
        for X in range(4):
            for Z in range(4,6):
                W_L_ss += mu_L * (GammaL[Z,X] * rho_ss[X] - GammaL[X,Z] * rho_ss[Z])
                W_R_ss += mu_R * (GammaR[Z,X] * rho_ss[X] - GammaR[X,Z] * rho_ss[Z])
        
        WL_values.append(W_L_ss)
        WR_values.append(W_R_ss)
        
        # --- Corrientes de calor electrónicas ---
        J_L_ss, J_R_ss = 0, 0
        for X in range(4):
            for Z in range(4,6):
                J_L_ss += (GammaL[Z,X] * rho_ss[X] - GammaL[X,Z] * rho_ss[Z]) * (energies[X] - energies[Z] - mu_L)
                J_R_ss += (GammaR[Z,X] * rho_ss[X] - GammaR[X,Z] * rho_ss[Z]) * (energies[X] - energies[Z] - mu_R)
        
        JL_values.append(J_L_ss)
        JR_values.append(J_R_ss)
        
        # --- Corriente de calor fonónica ---
        Jph = 0
        for j in range(6):
            for i in range(6):
                if i != j:
                    Jph += (D_diagcoh[i,j] * rho_ss[j] - D_diagcoh[j,i] * rho_ss[i]) * (energies[j] - energies[i])
        Jph_values.append(Jph/2)
    
        
        el, er = 0, 0
        for j in range(6):
            for i in range(6):
                if i != j:
                    el += (GammaL[i,j] * rho_ss[j] - GammaL[j,i] * rho_ss[i]) * (energies[j] - energies[i])
                    er += (GammaR[i,j] * rho_ss[j] - GammaR[j,i] * rho_ss[i]) * (energies[j] - energies[i])
        e1_values.append(el)
        e2_values.append(er)
        
    T_ref=T_L
    
    
    # --- Guardar resultados normalizados para cada E2 ---
    all_Jph.append(np.array(Jph_values) / (kB * T_ref))
    
    all_JL.append(np.array(JL_values) / (kB * T_ref))
    all_JR.append(np.array(JR_values) / (kB * T_ref))
        
    all_Sph.append(np.array(Jph_values) / Te_ph)
    all_SL.append(np.array(JL_values) / T_L)
    all_SR.append(np.array(JR_values) / T_R)
    all_SLR.append(np.array(JL_values) / T_L + np.array(JR_values) / T_R)
    all_Winput.append((np.array(WL_values) + np.array(WR_values)) / (kB * T_ref))






# --- Etiquetas cualitativas del barrido E₂ ---
labels = [
    "E₂ = E₁",
    "E₂ = E₁ − E_C/4",
    "E₂ = E₁ − E_C/2",
    "E₂ = E₁ − 0.75E_C",
    "E₂ = E₁ − E_C",
    "E₂ = E₁ − 1.5E_C",
    "E₂ = E₁ − 2.5E_C",
    "E₂ = E₁ − 4E_C"
]

colormaps = {
    'Jph': cm.cool,
    'JL': cm.Blues,
    'JR': cm.Reds,
    'Sph': cm.cool,      # verdes más claros
    'SL': cm.Purples,       # rositas claros (válido)
    'SR': cm.spring,      # morados oscuros
    'SLR': cm.winter      # cian-dorado neutro
}

# Variable común en el eje x
x = mu_diff_values  # ya está normalizado

# --- Función para trazar una familia de curvas ---
# --- Función para trazar una familia de curvas ---
def plot_family(data_list, ylabel, title, cmap_name, mu_range=None):
    plt.figure(figsize=(7, 4))
    cmap = colormaps[cmap_name]
    colores = cmap(np.linspace(0.2, 0.9, len(data_list)))

    # Aplicar filtro de mu solo si se especifica un rango
    if mu_range is not None:
        x_filtered = [mu for mu in x if mu_range[0] <= mu <= mu_range[1]]
        indices = [i for i, mu in enumerate(x) if mu_range[0] <= mu <= mu_range[1]]
        data_list = [data[indices] for data in data_list]
    elif cmap_name in ['JR', 'SR']:
        x_filtered = [mu for mu in x if -4 <= mu <= 1]
        indices = [i for i, mu in enumerate(x) if -4 <= mu <= 1]
        data_list = [data[indices] for data in data_list]
    else:
        x_filtered = x

    for i, data in enumerate(data_list):
        alpha = 0.35 if i in [5, 7] else 1.0
        plt.plot(x_filtered, data, color=colores[i], label=labels[i], alpha=alpha)

    if cmap_name in ['Sph', 'SL', 'SR', 'SLR']:
        unidades = r' [$k_B/\mathrm{tiempo}$]'
    else:
        unidades = r' [$\mathrm{tiempo}^{-1}$]'

    plt.xlabel(r'$\mu = (\mu_L - \mu_R)/k_BT$')
    plt.ylabel(ylabel + unidades)
    plt.title(title)
    plt.legend(fontsize=9)
    plt.grid(True, linestyle='--', linewidth=0.5)
    plt.tight_layout()
    plt.axhline(0, color='black', linewidth=1.2)
    plt.show()


# Gráficas con llamada normal
plot_family(all_Jph, r'$J_{ph} / (k_B T)$', r'$J_{ph}$', 'Jph')
plot_family(all_JL, r'$J_L / (k_B T)$', r'$J_L$', 'JL')
plot_family(all_JR, r'$J_R / (k_B T)$', r'$J_R$', 'JR')
plot_family(all_Sph, r'$\dot{S}_{ph}$', r'$\dot{S}_{ph}$', 'Sph')
plot_family(all_SL, r'$\dot{S}_L$', r'$\dot{S}_L$', 'SL')
plot_family(all_SR, r'$\dot{S}_R$', r'$\dot{S}_R$', 'SR')

# Gráfica restringida solo para SL+SR
plot_family(all_SLR, r'$\dot{S}_L + \dot{S}_R$', r'$\dot{S}_L + \dot{S}_R$', 'SLR', mu_range=(-0.5, 0.8))




# --- Análisis de cocientes para ciertos E₂ ---

# Índices correspondientes a E₂ = E₁ (0) y E₂ = E₁ - E_C/2 (2)
indices_objetivo = [0, 2]

# --- Cocientes de energía ---
ratios_J = []
for i in indices_objetivo:
    numerador = np.abs(np.array(all_JL[i]) + np.array(all_JR[i]))
    denominador = np.array(all_Jph[i])
    ratios_J.append(numerador / denominador)

# --- Cocientes de entropía total electrónica / fonónica ---
ratios_S = []
for i in indices_objetivo:
    numerador = np.abs(np.array(all_SL[i]) + np.array(all_SR[i]))
    denominador = np.array(all_Sph[i])
    ratios_S.append(numerador / denominador)

# --- Cocientes individuales SR/Sph para todos los E₂ ---
ratios_SR_Sph = []
for i in range(len(all_SR)):
    ratios_SR_Sph.append(np.abs(np.array(all_SR[i])) / np.array(all_Sph[i]))

# --- Etiquetas y colores ---
labels_casos = [labels[i] for i in indices_objetivo]
colores_r = cm.RdPu(np.linspace(0.3, 0.8, len(ratios_SR_Sph)))

# --- Filtro de mu entre -2 y 1 ---
mu_filtrada = [mu for mu in mu_diff_values if -2 <= mu <= 1]
indices_mu = [i for i, mu in enumerate(mu_diff_values) if -2 <= mu <= 1]

# --- Gráfica 1: |J_L + J_R| / J_ph ---
colores_j = [cm.plasma(0.4), cm.plasma(0.7)]
plt.figure(figsize=(7, 4))
for i, r in enumerate(ratios_J):
    datos_filtrados = r[indices_mu]
    plt.plot(mu_filtrada, datos_filtrados, label=labels_casos[i], color=colores_j[i])
plt.xlabel(r'$\mu = (\mu_L - \mu_R)/k_BT$')
plt.ylabel(r'$| \!J_L + J_R| / J_{ph}$')
plt.legend()
plt.grid(True, linestyle='--', linewidth=0.5)
plt.title(r'$| \!J_L + J_R| / J_{ph}$')
plt.tight_layout()
plt.show()

# --- Gráfica 2: |S_L + S_R| / S_ph ---
colores_s = [cm.viridis(0.4), cm.viridis(0.7)]
plt.figure(figsize=(7, 4))
for i, r in enumerate(ratios_S):
    datos_filtrados = r[indices_mu]
    plt.plot(mu_filtrada, datos_filtrados, label=labels_casos[i], color=colores_s[i])
plt.xlabel(r'$\mu = (\mu_L - \mu_R)/k_BT$')
plt.ylabel(r'$|\dot{S}_L + \dot{S}_R| / \dot{S}_{ph}$')
plt.legend()
plt.title(r'$|\dot{S}_L + \dot{S}_R| / \dot{S}_{ph}$')
plt.grid(True, linestyle='--', linewidth=0.5)
plt.tight_layout()
plt.show()

# --- Gráfica 3: |S_R| / S_ph para todos los E₂ ---
colores_sr = cm.RdPu(np.linspace(0.3, 0.9, len(ratios_SR_Sph)))
plt.figure(figsize=(7, 4))
for i, r in enumerate(ratios_SR_Sph):
    datos_filtrados = r[indices_mu]
    plt.plot(mu_filtrada, datos_filtrados, label=labels[i], color=colores_sr[i])
plt.xlabel(r'$\mu = (\mu_L - \mu_R)/k_BT$')
plt.ylabel(r'$|\dot{S}_R| / \dot{S}_{ph}$')
plt.title(r'$|\dot{S}_R| / \dot{S}_{ph}$')
plt.legend(fontsize=9)
plt.grid(True, linestyle='--', linewidth=0.5)
plt.tight_layout()
plt.show()

print("len(all_Winput) =", len(all_Winput))
print("len(all_Jph) =", len(all_Jph))


# --- Índices E₂ = E₁ y E₁ - E_C/2 ---
indices_objetivo = [0, 2]
colores = [cm.Reds(0.6), cm.Reds(0.85)]

# --- Eficiencia bomba térmica ---
plt.figure(figsize=(7, 4))
for i, idx in enumerate(indices_objetivo):
    eta_b = np.abs(all_Jph[idx]) / np.abs(all_Winput[idx])
    plt.plot(mu_filtrada, eta_b[indices_mu], label=labels[idx], color=colores[i])
plt.xlabel(r'$\mu = (\mu_L - \mu_R)/k_BT$')
plt.ylabel(r'$\eta_{\mathrm{bomba}}$')
plt.title(r'$|J_{ph}| / |W_L + W_R|$')
plt.legend()
plt.grid(True, linestyle='--', linewidth=0.5)
plt.axhline(0, color='black')
plt.tight_layout()
plt.show()

# --- Eficiencia frigorífica ---
colores = [cm.Blues(0.6), cm.Blues(0.85)]
plt.figure(figsize=(7, 4))
for i, idx in enumerate(indices_objetivo):
    eta_f = np.abs(all_JL[idx]) / np.abs(all_Winput[idx])
    plt.plot(mu_filtrada, eta_f[indices_mu], label=labels[idx], color=colores[i])
plt.xlabel(r'$\mu = (\mu_L - \mu_R)/k_BT$')
plt.ylabel(r'$\eta_{\mathrm{frigo}}$')
plt.title(r'$|J_L| / |W_L + W_R|$')
plt.legend()
plt.grid(True, linestyle='--', linewidth=0.5)
plt.axhline(0, color='black')
plt.tight_layout()
plt.show()



import matplotlib.patches as mpatches

plt.figure(figsize=(7, 4))

# --- Eficiencia bomba térmica ---
colores_bomba = [cm.Blues(0.3), cm.Blues(0.85)]
for i, idx in enumerate(indices_objetivo):
    eta_b = np.abs(all_Jph[idx]) / np.abs(all_Winput[idx])
    plt.plot(mu_filtrada, eta_b[indices_mu],
             label=f'{labels[idx]} (frigo)', color=colores_bomba[i])

# --- Eficiencia frigorífica ---
colores_frigo = [cm.Reds(0.3), cm.Reds(0.85)]
for i, idx in enumerate(indices_objetivo):
    eta_f = np.abs(all_JR[idx]) / np.abs(all_Winput[idx])
    plt.plot(mu_filtrada, eta_f[indices_mu],
             label=f'{labels[idx]} (bomba)', color=colores_frigo[i])

# --- Leyenda adicional con definiciones ---
bomba_patch = mpatches.Patch(color='white', label=r'Frigo: $\eta = J_{ph} / |W_L + W_R|$')
frigo_patch = mpatches.Patch(color='white', label=r'Bomba: $\eta = J_R / |W_L + W_R|$')

plt.xlabel(r'$\mu = (\mu_L - \mu_R)/k_BT$')
plt.ylabel(r'Eficiencia normalizada')
plt.title(r'Eficiencia isoterma (T_L=T_R) como "bomba" y como "frigorífico"')
plt.legend(handles=[bomba_patch, frigo_patch] + plt.gca().get_legend_handles_labels()[0],
           fontsize=9, loc='upper left', frameon=True)
plt.grid(True, linestyle='--', linewidth=0.5)
plt.axhline(0, color='black')
plt.tight_layout()
plt.show()

# Cocientes para E₂ = E₁ − E_C/2
idx = 2
mu_filtrada = [mu for mu in mu_diff_values if -2 <= mu <= 1]
indices_mu = [i for i, mu in enumerate(mu_diff_values) if -2 <= mu <= 1]

# --- Gráfica: |J_L + J_R| / J_ph ---
plt.figure(figsize=(7, 4))
plt.plot(mu_filtrada, (((all_JL[idx] + all_JR[idx]) / all_Jph[idx]))[indices_mu])
plt.xlabel(r'$\mu = (\mu_L - \mu_R)/k_BT$')
plt.ylabel(r'$|\!J_L + J_R| / J_{ph}$')
plt.title(r'Para $E_2 = E_1 - E_C/2$')
plt.grid(True, linestyle='--', linewidth=0.5)
plt.tight_layout()
plt.show()

# --- Gráfica: |S_L + S_R| / S_ph ---
plt.figure(figsize=(7, 4))
plt.plot(mu_filtrada, (np.abs(all_SL[idx] + all_SR[idx]) / np.abs(all_Sph[idx]))[indices_mu])
plt.xlabel(r'$\mu = (\mu_L - \mu_R)/k_BT$')
plt.ylabel(r'$|\dot{S}_L + \dot{S}_R| / \dot{S}_{ph}$')
plt.title(r'Para $E_2 = E_1 - E_C/2$')
plt.grid(True, linestyle='--', linewidth=0.5)
plt.tight_layout()
plt.show()

# --- Gráfica: |S_R| / S_ph ---
plt.figure(figsize=(7, 4))
plt.plot(mu_filtrada, (np.abs(all_SR[idx]) / np.abs(all_Sph[idx]))[indices_mu])
plt.xlabel(r'$\mu = (\mu_L - \mu_R)/k_BT$')
plt.ylabel(r'$|\dot{S}_R| / \dot{S}_{ph}$')
plt.title(r'Para $E_2 = E_1 - E_C/2$')
plt.grid(True, linestyle='--', linewidth=0.5)
plt.tight_layout()
plt.show()


print("La entropía mínima de SL+SR ocurre en:", labels[np.argmin([np.min(slr) for slr in all_SLR])])

sorted_indices = np.argsort([np.min(slr) for slr in all_SLR])
print("La segunda entropía mínima de SL+SR ocurre en:", labels[sorted_indices[1]])

minimos = sorted([np.min(slr) for slr in all_SLR])
print("Diferencia entre los dos mínimos más bajos de SL+SR:", minimos[1] - minimos[0])

print("\n⎯⎯⎯⎯⎯⎯⎯⎯ Cocientes Jph / SLR_min para cada resonancia ⎯⎯⎯⎯⎯⎯⎯⎯")
for i in range(len(all_Jph)):
    slr = all_SLR[i]
    jph = all_Jph[i]
    idx_min = np.argmin(slr)
    cociente = np.abs(jph[idx_min] / slr[idx_min]) if slr[idx_min] != 0 else np.nan
    print(f"{labels[i]} → Jph / (SL+SR)_min = {cociente:.4f}")

