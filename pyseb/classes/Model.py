#!/usr/bin/env python3
import numpy as np
from ..utils.FlagElements import FlagElements

class Model(object):
    def __init__(self) -> None: # {{{
        # initialize structure
        self.mesh = Mesh2d()
        pass
        # }}}

    def extract(self, area): # {{{
        """EXTRACT - extract a model according to an Argus contour or flag list

        This routine extracts a submodel from a bigger model with respect to a given contour
        md must be followed by the corresponding exp file or flags list
        It can either be a domain file (argus type, .exp extension), or an array of element flags.
        If user wants every element outside the domain to be
        extract2d, add '~' to the name of the domain file (ex: '~HO.exp')
        an empty string '' will be considered as an empty domain
        a string 'all' will be considered as the entire domain

        Usage:
            md2 = extract(md, area)

        Examples:
            md2 = extract(md, 'Domain.exp')

        See also: EXTRUDE, COLLAPSE
        """

        #copy model
        md1 = copy.deepcopy(self)

        #get elements that are inside area
        flag_elem = FlagElements(md1, area)
        if not np.any(flag_elem):
            raise RuntimeError("extracted model is empty")

        #kick out all elements with 3 dirichlets
        spc_elem = np.nonzero(np.logical_not(flag_elem))[0]
        spc_node = np.unique(md1.mesh.elements[spc_elem, :]) - 1
        flag = np.ones(md1.mesh.numberofvertices)
        flag[spc_node] = 0
        pos = np.nonzero(np.logical_not(np.sum(flag[md1.mesh.elements - 1], axis=1)))[0]
        flag_elem[pos] = 0

        #extracted elements and nodes lists
        pos_elem = np.nonzero(flag_elem)[0]
        pos_node = np.unique(md1.mesh.elements[pos_elem, :]) - 1

        #keep track of some fields
        numberofvertices1 = md1.mesh.numberofvertices
        numberofelements1 = md1.mesh.numberofelements
        numberofvertices2 = np.size(pos_node)
        numberofelements2 = np.size(pos_elem)
        flag_node = np.zeros(numberofvertices1)
        flag_node[pos_node] = 1

        #Create Pelem and Pnode (transform old nodes in new nodes and same thing for the elements)
        Pelem = np.zeros(numberofelements1, int)
        Pelem[pos_elem] = np.arange(1, numberofelements2 + 1)
        Pnode = np.zeros(numberofvertices1, int)
        Pnode[pos_node] = np.arange(1, numberofvertices2 + 1)

        #renumber the elements (some node won't exist anymore)
        elements_1 = copy.deepcopy(md1.mesh.elements)
        elements_2 = elements_1[pos_elem, :]
        elements_2[:, 0] = Pnode[elements_2[:, 0] - 1]
        elements_2[:, 1] = Pnode[elements_2[:, 1] - 1]
        elements_2[:, 2] = Pnode[elements_2[:, 2] - 1]
        if md1.mesh.__class__.__name__ == 'mesh3dprisms':
            elements_2[:, 3] = Pnode[elements_2[:, 3] - 1]
            elements_2[:, 4] = Pnode[elements_2[:, 4] - 1]
            elements_2[:, 5] = Pnode[elements_2[:, 5] - 1]

        #OK, now create the new model!

        #take every field from model
        md2 = copy.deepcopy(md1)

        #modify some specific fields
        #Mesh
        md2.mesh.numberofelements = numberofelements2
        md2.mesh.numberofvertices = numberofvertices2
        md2.mesh.elements = elements_2

        flag_elem_2d = flag_elem[np.arange(0, md1.mesh.numberofelements2d)]
        pos_elem_2d = np.nonzero(flag_elem_2d)[0]
        flag_node_2d = flag_node[np.arange(0, md1.mesh.numberofvertices2d)]
        pos_node_2d = np.nonzero(flag_node_2d)[0]

        md2.mesh.numberofelements2d = np.size(pos_elem_2d)
        md2.mesh.numberofvertices2d = np.size(pos_node_2d)
        md2.mesh.elements2d = md1.mesh.elements2d[pos_elem_2d, :]
        md2.mesh.elements2d[:, 0] = Pnode[md2.mesh.elements2d[:, 0] - 1]
        md2.mesh.elements2d[:, 1] = Pnode[md2.mesh.elements2d[:, 1] - 1]
        md2.mesh.elements2d[:, 2] = Pnode[md2.mesh.elements2d[:, 2] - 1]

        md2.mesh.x2d = md1.mesh.x[pos_node_2d]
        md2.mesh.y2d = md1.mesh.y[pos_node_2d]

        return md2
        # }}}

class Mesh2d(object):
    def __init__(self) -> None: # {{{
        pass
        # }}}
