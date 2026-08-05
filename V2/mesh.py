import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import V2.entities as ent


def OneDMeshing(frame):
    """
    frame is defined as a list of 2 membered tuples
        tuple[0] -> the defined beam oject
        tuple[1] -> number of elements in the beam object

    helper objects:
        XY_beams -> stores the beam, nparray(coords of the nodes in the element)
                 -> helps retrieving the element connectivity
        
        XY -> stores the global node coordinates 
           -> helps to duplicate things and sort order.
        
        nodesmap -> stores ths mapping from coordinates to nodes

    outputs:
    1) nodes -> numpy array of node objects
    2) elements -> numpy array of element objects
    """    
    X = []
    Y = []
    N = 0
    XY_beams = []
    for beam in frame:
        n = beam[1] + 1
        N += n
        beamobject = beam[0]
        beamtype = beamobject.type
        if beamtype == "Straight":
            x = np.round(np.linspace(beamobject.start[0], beamobject.end[0], n), 5)
            y = np.round(np.linspace(beamobject.start[1], beamobject.end[1], n), 5)
        elif beamtype == "Curved":
            theta1 = beamobject.start_angle()
            theta2 = beamobject.end_angle()
            
        xy = np.array([x, y]).T
        XY_beams.append((beam, xy))
        x = x.tolist(); X += x
        y = y.tolist(); Y += y

    XY = np.array([X, Y]).T
    _, idx = np.unique(XY, axis = 0, return_index = True)
    XY = XY[np.sort(idx)]
    nodesmap = {tuple(xy) : ent.node(i, tuple(xy)) for i, xy in enumerate(XY)}
    nodes = list(nodesmap.values())

    elements = []
    i = 0
    for beam in XY_beams:
        xy = beam[1]
        for j in range(len(xy) - 1):
            elem_i = ent.element1D(i, nodesmap[tuple(xy[j])], nodesmap[tuple(xy[j+1])], beam[0])
            elements.append(elem_i)
            i += 1
    

    return nodes, elements, nodesmap
