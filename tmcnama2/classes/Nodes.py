import sys
sys.path.append("..")
import numpy as np
from itertools import count
from lib.stamping_functions import stamp_y_sparse

class Nodes:
    
    _node_counter = count(0)
    _index_counter = count(0)
    node_key_ = {}
    all_node_key_ = {}
    ground_node_ = "Unassigned"

    def __init__(self, name, v_nom, phase):
        self.name = name
        self.v_nom = v_nom
        self.phase = phase
        if phase == 'N':
            self.ground = True
        else:
            self.ground = False
        self.id = self._node_counter.__next__()
        self.index = None
        Nodes.node_key_[self.name] = self.id

    def assign_indexes(self, use_sparse):
        self.index = self._index_counter.__next__()

