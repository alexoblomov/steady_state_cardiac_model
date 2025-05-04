# @author : Alanna Kennard
# sensitivity of g tol wrt arterial compliance given
# maximal reserve volume compensation and heart rate
# compensation

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_theme()
from matplotlib.ticker import FormatStrFormatter

from parameters import *

# pressures
P_thorax = -4*1333
P_RA = P_thorax + dP_RA

# for a given hypothetical maximum HR
max_F = 180 # maximum allowed heart rate

# discretization
dx = 100
Csa_u_range = np.linspace(0.5*Csa_u, 1.5*Csa_u, dx)
Csa_l_range = np.linspace(0.5*Csa_l, 1.5*Csa_l, dx)

gtol_vs_CsaU_vs_CsaL= np.zeros((len(Csa_u_range), len(Csa_l_range)))

F = max_F
Hu = Hu_patient
Hl = Hl_patient

# Setting Vd_total to zero to find g tol at max Heart rate for case I

for i in range(dx):
    for j in range(dx):
            Csa_u = Csa_u_range[i]
            Csa_l = Csa_l_range[j]
            Csa = Csa_l + Csa_u
            Cs_l = Csa_l + Csv_l
            g_tol = (Vtotal - Cp*(C_RVD/C_LVD)*(dP_RA) -
                                (Tp*Gs+Csa)*Psa_u_star)/((Tp*Gs_l+ Csa_l)
                                                         *rho*Hu-Cs_l*rho*Hl)
            gtol_vs_CsaU_vs_CsaL[i,j] = g_tol

# conversions:
gtol_vs_CsaU_vs_CsaL = gtol_vs_CsaU_vs_CsaL / 100 / (g_earth / 100)
Csa_u_range = Csa_u_range * 1333 / 1000 # ml/mmHg
Csa_l_range = Csa_l_range * 1333 / 1000 # ml/mmHg
Csa_u_range = np.around(Csa_u_range, 6) # ml/mmHg
Csa_l_range = np.around(Csa_l_range, 6) # ml/mmHg

# heatmap plot

num_ticks = 10
# the index of the position of yticks
yticks_ = np.linspace(0, dx - 1, num_ticks, dtype=int)
# the content of labels of these yticks
yticklabels_ = [str(Csa_u_range[idx]) for idx in yticks_]

xticks_ = np.linspace(0, dx - 1, num_ticks, dtype=int)
# the content of labels of these yticks
xticklabels_ = [str(Csa_l_range[idx]) for idx in xticks_]

fig, ax = plt.subplots()
im = ax.imshow(gtol_vs_CsaU_vs_CsaL)
fig.colorbar(im, ax=ax)

# Show all ticks and label them with the respective list entries
ax.set_xticks(xticks_, labels=xticklabels_)
ax.set_yticks(yticks_, labels=yticklabels_)

ax.grid(False)
# Rotate the tick labels and set their alignment.
plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")
fig.tight_layout(pad=2  )

plt.ylabel('upper compliance (ml/mmHg)', fontsize = 12)
plt.xlabel('lower compliance (ml/mmHg)', fontsize = 12)
plt.savefig('compliance_gtol_height_factor.png')
