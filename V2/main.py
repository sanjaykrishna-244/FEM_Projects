import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import V2.geometry.crosssection as sec
import V2.geometry.definition as define
import V2.material as mat

import V2.mesh as mesh

steel = mat.Material()
rect = sec.Rectangle()
beam1 = define.Beam_element((0, 0), (1, 0), steel, rect)
beam2 = define.Beam_element((1, 0), (1, 1), steel, rect)
beam3 = define.Beam_element((1, 1), (0, 0), steel, rect)

print(beam3.sin())

frame = np.array([(beam1, 20), (beam2, 20), (beam3, 30)])

nodes, nodemap, XY = mesh.OneDMeshing(frame)

plt.figure(figsize=(8,4))

plt.scatter(XY[:,0], XY[:,1], color='k')

for i, (x, y) in enumerate(XY):
    plt.text(x, y, str(i),
             fontsize=8,
             color='red')

plt.axis("equal")
plt.grid(True)
plt.show()