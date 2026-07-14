import numpy as np

def assemble(nodes, elements):
    alldofs = len(nodes) * 3
    K = np.zeros((alldofs, alldofs))
    M = np.zeros((alldofs, alldofs))

    for element in elements:
        ke = element.Ke()
        me = element.Me()

        dofs = element.dof()

        K[np.ix_(dofs, dofs)] += ke
        M[np.ix_(dofs, dofs)] += me

    return K, M