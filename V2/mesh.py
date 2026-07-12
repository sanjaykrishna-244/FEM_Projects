import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def OneDMeshing(frame):
    X = []
    Y = []
    N = 0
    for beam in frame:
        n = beam[1] + 1
        N += n
        details = beam[0]
        x = np.linspace(details.start[0], details.end[0], n)
        x = x.tolist()
        X += x
        y = np.linspace(details.start[1], details.end[1], n)
        y = y.tolist()
        Y += y

    XY = np.array([X, Y]).T
    XY = np.round(XY, 5)
    _, idx = np.unique(XY, axis = 0, return_index = True)
    XY = XY[np.sort(idx)]
    nodes = {i : tuple(xy) for i, xy in enumerate(XY)}
    nodemap = {tuple(xy) : i for i, xy in enumerate(XY)}
    return nodes, nodemap, XY
