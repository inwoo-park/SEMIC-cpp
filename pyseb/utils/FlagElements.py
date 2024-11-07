import numpy as np
import os

def FlagElements(md, region):
    """
    FLAGELEMENTS - flag the elements in an region

       The region can be given with an exp file, a list of elements or vertices

       Usage:
          flag = FlagElements(md, region)

       Example:
          flag = FlagElements(md, 'all')
          flag = FlagElements(md, '')
          flag = FlagElements(md, 'Domain.exp')
          flag = FlagElements(md, '~Domain.exp')
    """
    if isinstance(region, np.ndarray) or isinstance(region, bool):
        if np.size(region, 0) == md.mesh.numberofelements:
            flag = region
        elif np.size(region, 0) == md.mesh.numberofvertices:
            flag = (np.sum(region[md.mesh.elements - 1] > 0, axis=1) == np.size(md.mesh.elements, 1))
        else:
            raise TypeError("Flaglist for region must be of same size as number of elements in model.")

    else:
        raise TypeError("Invalid region option")

    return flag
