import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.linalg import eigh
from shipdeckframe.doubleframe import K_red, M_red, details, xy, free

N = xy.shape[0]
# __ MODAL ANALYSIS ___________________
modes = np.zeros((3*N, 3*N))
freqs, modes[np.ix_(free, free)] = eigh(K_red, M_red)
freqs = np.sqrt(np.abs(freqs)) / (2 * np.pi)

# __ EXPORTING DETAILS ________________
dff = pd.DataFrame({
    'S.No' : np.arange(1, 31),
    'Frequency (Hz)' : freqs[:30],
})
with pd.ExcelWriter('shipdeckframe/doubleframe.xlsx', mode = 'a', engine = 'openpyxl', if_sheet_exists='replace') as writer:
    dff.to_excel(writer, sheet_name='Frequencies', index= False)

# __ PLOT _____________________________
for j in details:
    xj = details[j][3]
    yj = details[j][4]
    plt.plot(xj, yj, color = (0, 0, 0), linestyle = ':')
    for i in range(3, 6):
        if j == 1:
            plt.plot(xj + 5 * modes[:, i][np.ix_(details[j][0]*3)], yj + 5 * modes[:, i][np.ix_(details[j][0]*3 + 1)], color = (i/8, 1 - i/8, i/10), label = f"{i-2}th mode shape")
        else:
            plt.plot(xj + 5 * modes[:, i][np.ix_(details[j][0]*3)], yj + 5 * modes[:, i][np.ix_(details[j][0]*3 + 1)], color = (i/8, 1 - i/8, i/10))
plt.title("Mode shapes of first few modes of Ship-deck like double frame")
plt.legend()
plt.grid()
plt.xlabel("X-in m")
plt.ylabel("Y-in m")
plt.show()