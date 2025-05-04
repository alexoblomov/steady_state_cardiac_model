# @author : Alanna Kennard
# exploration of liminal values for pthorax defining the different cases
# g tol wrt arterial compliance given
# maximal reserve volume compensation and heart rate
# compensation

import numpy as np
import matplotlib.pyplot as plt
from parameters import *

G = np.linspace(g_earth,1.7*980, 3)

P_thorax = np.linspace(- 4 * 1333,3 * 1333,8)

#P_thorax = -4*1333;
P_RA = P_thorax + dP_RA

cases = np.zeros((len(P_thorax),len(G)), dtype=np.uint8)

Vd_total_vec = np.zeros(len(G))
Q_vec = np.zeros(len(G))
F_vec = np.zeros(len(G))
Ppa_vec = np.zeros(len(G))

sol_Vd_Pthorax_G = np.zeros((len(P_thorax),len(G)))
sol_Q_Pthorax_G = np.zeros((len(P_thorax),len(G)))
sol_F_Pthorax_G = np.zeros((len(P_thorax),len(G)))
sol_Ppa_Pthorax_G = np.zeros((len(P_thorax),len(G)))

Hu = Hu_patient
Hl = Hl_patient

for j in range(len(P_thorax)):
    dP_RA = P_RA[j] - P_thorax[j]
    for i in range(len(G)):
        if P_thorax[j] <= - dP_RA:
            Vd_total = Vtotal - Cp * (C_RVD / C_LVD) * (dP_RA) \
                       - (Tp * Gs + Csa) * Psa_u_star \
                       - (Tp * Gs_l + Csa_l) \
                       * rho * G[i] * Hu - Cs_l * rho * G[i] * (- Hl)
            Q = Gs * Psa_u_star + rho * G[i] * Hu / Rs_l
            F = Q / (C_RVD * (P_RA[j] - P_thorax[j]))
            Ppv = P_thorax[j] + (C_RVD / C_LVD) * dP_RA 
            Ppa = Ppv + Q * Rp

            cases[j, i] = 1
        elif P_thorax[j] > - dP_RA and P_thorax[j] < rho * G[i] * Hu - dP_RA:
            Vd_total = Vtotal - Cp * (C_RVD / C_LVD) * dP_RA \
                      - (Tp * Gs + Csa)* Psa_u_star \
                      - (Tp * Gs_l + Csa_l) * (rho * G[i] * Hu - Cs_l *
                                               rho * G[i] * (- Hl)) \
                      - (Csv_l - Tp * Gs_l) * (P_thorax[j] + dP_RA)
            Psv_l = - rho * G[i] * Hl + P_thorax[j] + dP_RA
            Psv_u = 0
            Psa_l = Psa_u_star + rho * G[i] * (Hu - Hl)

            Qs_u = Psa_u_star * 1 / Rs_u 
            Qs_l = (Psa_u_star + rho * G[i] * Hu \
                    - P_thorax[j]-dP_RA) * 1/Rs_l
            Q = Qs_u + Qs_l
            F = Q / (C_RVD * (dP_RA))
            F_also = (Gs * Psa_u_star +
                      (1/Rs_l)* (rho*G[i]*Hu-P_thorax[i] - dP_RA)) / (
                      C_RVD*dP_RA)


            Ppv = P_thorax[j] + (C_RVD / C_LVD) * dP_RA 
            Ppa = Ppv + Q * Rp

            cases[j, i] = 2
        elif P_thorax[j] >= rho * G[i] * Hu - dP_RA:

            Vd_total = Vtotal - Cp * (C_RVD / C_LVD) * dP_RA - \
                    (Tp*Gs + Csa_l+Csa_u)*Psa_u_star \
                    +(Tp*Gs + Csa_u - Csv_l) * rho * G[i] * Hu \
                    + (Csa_l+Csv_l) * rho * G[i] * (-Hl) \
                    + (Csv_l+Csv_u - Tp*Gs)* (P_thorax[j] - dP_RA)

            Psv_l = P_thorax[j] + dP_RA + rho * G[i] * (- Hl)
            Psv_u = P_thorax[j] + dP_RA - rho * G[i] * Hu
            Psa_l = Psa_u_star + rho * G[i] * (Hu - Hl)
            Qs_u = (Psa_u_star - Psv_u) / Rs_u


            Qs_l = (Psa_l - Psv_l) / Rs_l
            Q = Qs_u + Qs_l
            F = Q / (C_RVD * (dP_RA))
            Ppv = P_thorax[j] + (C_RVD / C_LVD) * dP_RA
            Ppa = Ppv + Q * Rp
            cases[j, i] = 3
        if Vd_total > 0:
            Vd_total_vec[i] = Vd_total
            Q_vec[i] = Q
            F_vec[i] = F
            Ppa_vec[i] = Ppa
        else:
            Vd_total_vec[i] = float('nan')
            Q_vec[i] = float('nan')
            F_vec[i] = float('nan')
            Ppa_vec[i] = float('nan')
    sol_Vd_Pthorax_G[j,:] = Vd_total_vec
    sol_Q_Pthorax_G[j,:] = Q_vec
    sol_F_Pthorax_G[j,:] = F_vec
    sol_Ppa_Pthorax_G[j,:] = Ppa_vec

#conversions:
G = G / 100 / (g_earth / 100)
sol_Q_Pthorax_G = sol_Q_Pthorax_G * 60 / 1000
sol_F_Pthorax_G = sol_F_Pthorax_G * 60
sol_Ppa_Pthorax_G = sol_Ppa_Pthorax_G / 1333
sol_Vd_Pthorax_G = sol_Vd_Pthorax_G / 1000

#family of plots for Reserve Volume vs. G for different Pthorax values:

plt.figure(figsize=(15, 12))
plt.subplots_adjust(hspace=0.5)
plt.suptitle("Reserve volume vs g", fontsize=18, y=0.95)
pthorax_titles = ["P thorax " + str(pth / 1333) + " mmHg" for pth in P_thorax]

for n, plt_title in enumerate(pthorax_titles):
    ax = plt.subplot(4, 2, n + 1)

    idx_case_1 = cases == 1
    idx_case_2 = cases == 2
    idx_case_3 = cases == 3
    ax.plot(G,sol_Vd_Pthorax_G[n,:])

    ax.set_title(plt_title)
    ax.set_xlabel("g multiple")
    ax.set_ylabel("VT0")
plt.savefig("cleaned_fixed_VT0_vs_g_V0_w_height_factor.png")