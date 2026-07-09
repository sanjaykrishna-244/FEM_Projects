import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.linalg import eigh
from shipdeckframe.doubleframe import K_red, M_red, details, xy, free

N = xy.shape[0]
# __ STATIC LOAD ANALYSIS _____________
alldf = np.arange(0, N)
eles = len(alldf) - 1
F = np.zeros(3*N)
F[184] = -10000
u = np.zeros_like(F)
u[free] = np.linalg.solve(K_red, F[np.ix_(free)])

'''# __ EXPORTING DETAILS ________________
dfu = pd.DataFrame({
    'x' : xy[:, 0],
    'y' : xy[:, 1],
    'u' : u[np.arange(0, 3*N, 3)],
    'v' : u[np.arange(1, 3*N, 3)],
    'theta' : u[np.arange(2, 3*N, 3)],
})
with pd.ExcelWriter('shipdeckframe/doubleframe.xlsx', mode = 'a', engine = 'openpyxl', if_sheet_exists='replace') as writer:
    dfu.to_excel(writer, sheet_name='Deflection', index= False)
'''
# __ PLOT _____________________________
plt.figure(figsize = (10, 5))
for j in details:
    xj = details[j][3]
    yj = details[j][4]
    plt.plot(xj, yj, color = (0, 0, 0), linestyle = ':')
    plt.plot(xj + 1 * u[np.ix_(details[j][0]*3)], yj + 1 * u[np.ix_(details[j][0]*3 + 1)], color = (1/8, 1 - 1/8, 1/10))
plt.title("Deflection under 100000 N load of Ship-deck like double frame")
plt.grid()
plt.axis('equal')
plt.xlabel("X-in m")
plt.ylabel("Y-in m")
plt.savefig('shipdeckframe/staticload.png', dpi = 300)
plt.show()
