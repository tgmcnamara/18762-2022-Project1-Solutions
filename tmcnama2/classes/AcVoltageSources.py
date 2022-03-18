import numpy as np
from itertools import count
from classes.Nodes import Nodes
from lib.stamping_functions import stamp_y_sparse, stamp_j_sparse

class AcVoltageSources:
    _acvs_counter = count(0)
    acvs_key_ = {}

    def __init__(self, name, phase, from_node, to_node, amplitude, angle, frequency):
        self.name = name
        self.phase = phase
        self.amplitude = amplitude[0]*np.sqrt(2/3)
        self.frequency = frequency[0]
        self.angle_rad = angle[0]*np.pi/180.0
        self.from_node = from_node
        self.to_node = to_node
        self.from_node_index = None
        self.to_node_index = None
        self.id = self._acvs_counter.__next__()
        AcVoltageSources.acvs_key_[(self.from_node, self.to_node, self.name)] = self.id

    def assign_indexes(self, node_list):
        from_node_obj = node_list[Nodes.node_key_[self.from_node]]
        self.from_node_index = from_node_obj.index
        to_node_obj = node_list[Nodes.node_key_[self.to_node]]
        self.to_node_index = to_node_obj.index
        self.current_index = Nodes._index_counter.__next__()

    def stamp_sparse(self, t, y_row, y_col, y_val, y_counter, j_row, j_val, j_counter):
        V = self.amplitude*np.sin(2*np.pi*self.frequency*t + self.angle_rad)
        # current out of from node
        y_counter = stamp_y_sparse(self.from_node_index, self.current_index, 1, y_row, y_col, y_val, y_counter)
        # current into to node
        y_counter = stamp_y_sparse(self.to_node_index, self.current_index, -1, y_row, y_col, y_val, y_counter)
        # current_index row: V_from - V_to = V
        y_counter = stamp_y_sparse(self.current_index, self.from_node_index, 1, y_row, y_col, y_val, y_counter)
        y_counter = stamp_y_sparse(self.current_index, self.to_node_index, -1, y_row, y_col, y_val, y_counter)
        j_counter = stamp_j_sparse(self.current_index, V, j_row, j_val, j_counter)
        return (y_counter, j_counter)

    def stamp_dense(self, t, Y, J):
        V = self.amplitude*np.sin(2*np.pi*self.frequency*t + self.angle_rad)
        # current out of to node
        Y[self.to_node_index, self.current_index] += 1
        # current into from node
        Y[self.from_node_index, self.current_index] += -1
        # current_index row: V_from - V_to = V
        Y[self.current_index, self.from_node_index] += 1
        Y[self.current_index, self.to_node_index] += -1
        J[self.current_index, 0] += V
        
