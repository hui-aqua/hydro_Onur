from . import aster
from . import one_dimensional
from . import two_dimensional
from . import wake_effect

import numpy as np



def distribute_force(elementList: list, number_of_node:int)->np.array:
    """
    Transfer the forces on elements to their corresponding nodes.\n
    :return: [np.array].shape=(N,3) Unit [N]. The hydrodynamic forces on nodes.
    """
    force_on_nodes=np.zeros((number_of_node,3))
    for each in elementList:
        for n in each.element:
            force_on_nodes[n]+=each.hydro_total_forces/len(each.element)
    return force_on_nodes
