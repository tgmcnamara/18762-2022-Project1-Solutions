import numpy as np
from itertools import count
from classes.Nodes import Nodes
from lib.stamping_functions import stamp_y_sparse, stamp_j_sparse


class Capacitors:
    _capacitor_counter = count(0)
    capacitor_key_ = {}

    def __init__(self, name, phase, from_node, to_node, c):
        self.name = name
        self.phase = phase
        self.c = c[0]
        self.from_node = from_node
        self.to_node = to_node
        self.from_node_index = None
        self.to_node_index = None
        self.id = self._capacitor_counter.__next__()
        Capacitors.capacitor_key_[self.name] = self.id
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
        companion_r = dt/(2*self.c)
        companion_v = (V[self.from_node_index,0] - V[self.to_node_index,0]) + companion_r*V[self.companion_current_index,0]

        # +1/r in from_node row
        y_counter = stamp_y_sparse(self.from_node_index, self.from_node_index, 1/companion_r, y_row, y_col, y_val, y_counter)
        y_counter = stamp_y_sparse(self.from_node_index, self.companion_node_index, -1/companion_r, y_row, y_col, y_val, y_counter)
        # (V[companion_node]-V[from_node])/R_companion + i_companion= 0 in companion_node row
        y_counter = stamp_y_sparse(self.companion_node_index, self.from_node_index, -1/companion_r, y_row, y_col, y_val, y_counter)
        y_counter = stamp_y_sparse(self.companion_node_index, self.companion_node_index, 1/companion_r, y_row, y_col, y_val, y_counter)
        y_counter = stamp_y_sparse(self.companion_node_index, self.companion_current_index, 1, y_row, y_col, y_val, y_counter)
        # V[companion_node] - V[to_node] = companion_v
        y_counter = stamp_y_sparse(self.companion_current_index, self.companion_node_index, 1, y_row, y_col, y_val, y_counter)
        y_counter = stamp_y_sparse(self.companion_current_index, self.to_node_index, -1, y_row, y_col, y_val, y_counter)
        j_counter = stamp_j_sparse(self.companion_current_index, companion_v, j_row, j_val, j_counter)
        # i_companion into to_node row
        y_counter = stamp_y_sparse(self.to_node_index, self.companion_current_index, -1, y_row, y_col, y_val, y_counter)

        return (y_counter, j_counter)

    def stamp_dense(self, V, dt, Y, J):
        companion_r = dt/(2*self.c)
        companion_v = (V[self.from_node_index,0] - V[self.to_node_index,0]) + companion_r*V[self.companion_current_index,0]

        # from_node row: add (v_from - v_companion)*1/companion_r
        Y[self.from_node_index, self.from_node_index] += 1/companion_r
        Y[self.from_node_index, self.companion_node_index] += -1/companion_r
        # companion_node row: (V[companion_node]-V[from_node])*1/companion_r + i_companion= 0
        Y[self.companion_node_index, self.from_node_index] += -1/companion_r
        Y[self.companion_node_index, self.companion_node_index] += 1/companion_r
        Y[self.companion_node_index, self.companion_current_index] += 1
        # V[companion_node] - V[to_node] = companion_v
        Y[self.companion_current_index, self.companion_node_index] += 1
        Y[self.companion_current_index, self.to_node_index] += -1
        J[self.companion_current_index] += companion_v
        # i_companion into to_node row
        Y[self.to_node_index, self.companion_current_index] += -1

    def stamp_open_sparse(self, y_row, y_col, y_val, y_counter):
        # companion_node row: v_from_node - v_companion_node = 0
        y_counter = stamp_y_sparse(self.companion_node_index, self.from_node_index, 1, y_row, y_col, y_val, y_counter)
        y_counter = stamp_y_sparse(self.companion_node_index, self.companion_node_index, -1, y_row, y_col, y_val, y_counter)
        # companion_current row: companion_current = 0
        y_counter = stamp_y_sparse(self.companion_current_index, self.companion_current_index, 1, y_row, y_col, y_val, y_counter)
        return y_counter

    def stamp_open_dense(self, Y):
        # companion_node row: v_from_node - v_companion_node = 0
        Y[self.companion_node_index, self.from_node_index] += 1
        Y[self.companion_node_index, self.companion_node_index] += -1        
        # companion_current row: companion_current = 0
        Y[self.companion_current_index, self.companion_current_index] += 1        

    def calculate_voltage_across(self, V):
        return (V[self.from_node_index,0] - V[self.to_node_index,0])

    def calculate_current_across(self, V):
        return V[self.companion_current_index,0]