import numpy as np
from itertools import count
from classes.Nodes import Nodes
from lib.stamping_functions import stamp_y_sparse, stamp_j_sparse


class Switches:
    _switch_counter = count(0)
    switch_key_ = {}

    def __init__(self, name, phase, from_node, to_node, t_open, t_close):
        self.name = name
        self.phase = phase
        self.from_node = from_node
        self.to_node = to_node
        self.from_node_index = None
        self.to_node_index = None
        self.t_open = t_open[0]
        self.t_close = t_close[0]
        self.id = self._switch_counter.__next__()
        Switches.switch_key_[(self.from_node, self.to_node, self.name)] = self.id
        

    def assign_indexes(self, node_list):
        from_node_obj = node_list[Nodes.node_key_[self.from_node]]
        self.from_node_index = from_node_obj.index
        to_node_obj = node_list[Nodes.node_key_[self.to_node]]
        self.to_node_index = to_node_obj.index
        self.switch_current_index = Nodes._index_counter.__next__()

    def stamp_constant_sparse(self, y_row, y_col, y_val, y_counter):
        # switch current out of from node
        y_counter = stamp_y_sparse(self.from_node_index, self.switch_current_index, 1, y_row, y_col, y_val, y_counter)
        # switch current into to node
        y_counter = stamp_y_sparse(self.to_node_index, self.switch_current_index, -1, y_row, y_col, y_val, y_counter)
        return y_counter

    def stamp_constant_dense(self, Y):
        # switch current out of from node
        Y[self.from_node_index, self.switch_current_index] += 1
        # switch current into to node
        Y[self.to_node_index, self.switch_current_index] += -1

    def stamp_time_dependent_sparse(self, t, y_row, y_col, y_val, y_counter):
        if t >= self.t_open and t < self.t_close:
            # switch open: switch current = 0
            y_counter = stamp_y_sparse(self.switch_current_index, self.switch_current_index, 1, y_row, y_col, y_val, y_counter)
        else:
            # switch closed: v_from - v_to = 0
            y_counter = stamp_y_sparse(self.switch_current_index, self.from_node_index, 1, y_row, y_col, y_val, y_counter)
            y_counter = stamp_y_sparse(self.switch_current_index, self.to_node_index, -1, y_row, y_col, y_val, y_counter)
        return y_counter

    def stamp_time_dependent_dense(self, t, Y):
        if t >= self.t_open and t < self.t_close:
            # switch open: switch current = 0
            Y[self.switch_current_index, self.switch_current_index] += 1
        else:
            # switch closed: v_from - v_to = 0
            Y[self.switch_current_index, self.from_node_index] += 1
            Y[self.switch_current_index, self.to_node_index] += -1

