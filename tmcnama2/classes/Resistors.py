
import numpy as np
from itertools import count
from classes.Nodes import Nodes
from lib.stamping_functions import stamp_y_sparse, stamp_j_sparse


class Resistors:
    _resistor_counter = count(0)
    resistor_key_ = {}

    def __init__(self, name, phase, from_node, to_node, r):
        self.name = name
        self.phase = phase
        self.r = r[0]
        self.from_node = from_node
        self.to_node = to_node
        self.from_node_index = None
        self.to_node_index = None
        self.id = self._resistor_counter.__next__()
        Resistors.resistor_key_[self.name] = self.id

    def assign_indexes(self, node_list):
        from_node_obj = node_list[Nodes.node_key_[self.from_node]]
        self.from_node_index = from_node_obj.index
        to_node_obj = node_list[Nodes.node_key_[self.to_node]]
        self.to_node_index = to_node_obj.index

    def stamp_sparse(self, ylin_row, ylin_col, ylin_val, ylin_counter):
        ylin_counter = stamp_y_sparse(self.from_node_index, self.from_node_index, 1/self.r, ylin_row, ylin_col, ylin_val, ylin_counter)
        ylin_counter = stamp_y_sparse(self.from_node_index, self.to_node_index, -1/self.r, ylin_row, ylin_col, ylin_val, ylin_counter)
        ylin_counter = stamp_y_sparse(self.to_node_index, self.from_node_index, -1/self.r, ylin_row, ylin_col, ylin_val, ylin_counter)
        ylin_counter = stamp_y_sparse(self.to_node_index, self.to_node_index, 1/self.r, ylin_row, ylin_col, ylin_val, ylin_counter)
        return(ylin_counter)
        
    def stamp_dense(self, Y):
        Y[self.from_node_index, self.from_node_index] += 1/self.r
        Y[self.from_node_index, self.to_node_index] += -1/self.r
        Y[self.to_node_index, self.from_node_index] += -1/self.r
        Y[self.to_node_index, self.to_node_index] += 1/self.r

    def calculate_voltage_across(self, V):
        return V[self.from_node_index] - V[self.to_node_index]

    def calculate_current_across(self, V):
        v = self.calculate_voltage_across(V)
        return (v/self.r)