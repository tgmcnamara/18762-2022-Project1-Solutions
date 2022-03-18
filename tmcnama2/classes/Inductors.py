import sys
sys.path.append("..")
import numpy as np
from itertools import count
from classes.Nodes import Nodes
from lib.stamping_functions import stamp_y_sparse, stamp_j_sparse


class Inductors:
    _inductor_counter = count(0)
    inductor_key_ = {}

    def __init__(self, name, phase, from_node, to_node, l):
        self.name = name
        self.phase = phase
        self.l = l[0]
        self.from_node = from_node
        self.to_node = to_node
        self.from_node_index = None
        self.to_node_index = None
        self.id = self._inductor_counter.__next__()
        Inductors.inductor_key_[self.name] = self.id
        self.companion_node_index = None
        self.companion_current_index = None
        self.most_recent_voltage = 0.0
        self.most_recent_current = 0.0
        self.most_recent_t = 0.0

    def assign_indexes(self, node_list):
        from_node_obj = node_list[Nodes.node_key_[self.from_node]]
        self.from_node_index = from_node_obj.index
        to_node_obj = node_list[Nodes.node_key_[self.to_node]]
        self.to_node_index = to_node_obj.index
        self.companion_node_index = Nodes._index_counter.__next__()
        self.companion_current_index = Nodes._index_counter.__next__()

    def stamp_sparse(self, V, dt, y_row, y_col, y_val, y_counter, j_row, j_val, j_counter):
        companion_g = dt/(2*self.l)
        i_source = V[self.companion_current_index,0] + companion_g*(V[self.from_node_index,0] - V[self.to_node_index,0])
        # from_node row: companion_current out
        y_counter = stamp_y_sparse(self.from_node_index, self.companion_current_index, 1, y_row, y_col, y_val, y_counter)
        
        # companion_current_node row: v_from_node - v_companion_node = 0
        y_counter = stamp_y_sparse(self.companion_current_index, self.from_node_index, 1, y_row, y_col, y_val, y_counter)
        y_counter = stamp_y_sparse(self.companion_current_index, self.companion_node_index, -1, y_row, y_col, y_val, y_counter)

        # companion_node_row: -i_companion + (v_companion_node - v_to_node)*dt/(2*L) = -i_source
        y_counter = stamp_y_sparse(self.companion_node_index, self.companion_current_index, -1, y_row, y_col, y_val, y_counter)
        y_counter = stamp_y_sparse(self.companion_node_index, self.companion_node_index, companion_g, y_row, y_col, y_val, y_counter)
        y_counter = stamp_y_sparse(self.companion_node_index, self.to_node_index, -companion_g, y_row, y_col, y_val, y_counter)
        j_counter = stamp_j_sparse(self.companion_node_index, -i_source, j_row, j_val, j_counter)

        # to_node_row: add dt/2L*v_to_node to Y, +i_source to J
        y_counter = stamp_y_sparse(self.to_node_index, self.to_node_index, companion_g, y_row, y_col, y_val, y_counter)        
        y_counter = stamp_y_sparse(self.to_node_index, self.companion_node_index, -companion_g, y_row, y_col, y_val, y_counter)        
        j_counter = stamp_j_sparse(self.to_node_index, i_source, j_row, j_val, j_counter)

        return (y_counter, j_counter)

    def stamp_dense(self, V, dt, Y, J):
        companion_g = dt/(2*self.l)
        i_source = V[self.companion_current_index,0] + companion_g*(V[self.from_node_index,0] - V[self.to_node_index,0])
        # from_node row: companion_current out
        Y[self.from_node_index, self.companion_current_index] += 1
        # companion_current_node row: v_from_node - v_companion_node = 0
        Y[self.companion_current_index, self.from_node_index] += 1
        Y[self.companion_current_index, self.companion_node_index] += -1
        # companion_node_row: -i_companion + (v_companion_node - v_to_node)*dt/(2*L) = i_source
        Y[self.companion_node_index, self.companion_current_index] += -1
        Y[self.companion_node_index, self.companion_node_index] += companion_g
        Y[self.companion_node_index, self.to_node_index] += -companion_g
        J[self.companion_node_index,0] += -i_source
        
        # to_node_row: add dt/2L*v_to_node to Y, -i_source to J
        Y[self.to_node_index, self.companion_node_index] += -companion_g
        Y[self.to_node_index, self.to_node_index] += companion_g
        J[self.to_node_index,0] += i_source

    def stamp_shorted_sparse(self, y_row, y_col, y_val, y_counter):
        # companion_current row: v_from_node - v_to_node = 0
        y_counter = stamp_y_sparse(self.companion_current_index, self.from_node_index, 1, y_row, y_col, y_val, y_counter)
        y_counter = stamp_y_sparse(self.companion_current_index, self.to_node_index, -1, y_row, y_col, y_val, y_counter)
        # v_from_row += i_companion
        y_counter = stamp_y_sparse(self.from_node_index, self.companion_current_index, 1, y_row, y_col, y_val, y_counter)
        # v_to_row += -i_companion
        y_counter = stamp_y_sparse(self.to_node_index, self.companion_current_index, -1, y_row, y_col, y_val, y_counter)
        # v_companion_row: v_from_node - v_companion_node = 0
        y_counter = stamp_y_sparse(self.companion_node_index, self.from_node_index, 1, y_row, y_col, y_val, y_counter)
        y_counter = stamp_y_sparse(self.companion_node_index, self.companion_node_index, -1, y_row, y_col, y_val, y_counter)
        return y_counter

    def stamp_shorted_dense(self, Y):
        # companion_current row: v_from_node - v_to_node = 0
        Y[self.companion_current_index, self.from_node_index] = 1
        Y[self.companion_current_index, self.to_node_index] = -1
        # v_from_row += i_companion
        Y[self.from_node_index, self.companion_current_index] = 1
        # v_to_row += -i_companion
        Y[self.to_node_index, self.companion_current_index] = -1
        # v_companion_row: v_from_node - v_companion_node = 0
        Y[self.companion_node_index, self.from_node_index] = 1
        Y[self.companion_node_index, self.companion_node_index] = -1