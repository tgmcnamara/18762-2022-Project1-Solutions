import numpy as np
from lib.create_sparse_matrix import create_sparse_matrix

def stamp_t0_dense(devices, size_Y):
    Y = np.zeros((size_Y, size_Y), dtype=np.float)
    J = np.zeros((size_Y, 1), dtype=np.float)

    # short inductors
    for ele in devices['inductors']:
        ele.stamp_shorted_dense(Y)
    # open capacitors
    for ele in devices['capacitors']:
        ele.stamp_open_dense(Y)
    # switches
    for ele in devices['switches']:
        ele.stamp_time_dependent_dense(0, Y)
    # constant sources
    for ele in devices['ac_voltage_sources']:
        ele.stamp_dense(0, Y, J)
    # induction motors
    for ele in devices['induction_motors']:
        ele.stamp_t0_dense(Y)

    return(Y, J)

def stamp_t0_sparse(devices, size_Y):
    num_Y_entries = size_Y * 100
    y_val = np.zeros(num_Y_entries)
    y_row = np.zeros(num_Y_entries)
    y_col = np.zeros(num_Y_entries)
    y_entry_counter = 0

    j_val = np.zeros(num_Y_entries)
    j_row = np.zeros(num_Y_entries)
    j_entry_counter = 0

    # short inductors
    for ele in devices['inductors']:
        (y_entry_counter) = ele.stamp_shorted_sparse(y_row, y_col, y_val, y_entry_counter)
    # open capacitors
    for ele in devices['capacitors']:
        (y_entry_counter) = ele.stamp_open_sparse(y_row, y_col, y_val, y_entry_counter)
    # switches at t=0
    for ele in devices['switches']:
        (y_entry_counter) = ele.stamp_time_dependent_sparse(0, y_row, y_col, y_val, y_entry_counter)
    # ac voltage sources at t=0
    for ele in devices['ac_voltage_sources']:
        (y_entry_counter, j_entry_counter) = ele.stamp_sparse(0, y_row, y_col, y_val, y_entry_counter, j_row, j_val, j_entry_counter)
    # induction motors do nothing at t=0, but stamp 1's on the diagonal to keep matrix non-singular
    for ele in devices['induction_motors']:
        (y_entry_counter) = ele.stamp_t0_sparse(y_row, y_col, y_val, y_entry_counter)

    y_val = y_val[:y_entry_counter]
    y_row = y_row[:y_entry_counter]
    y_col = y_col[:y_entry_counter]
    
    j_val = j_val[:j_entry_counter]
    j_row = j_row[:j_entry_counter]
    
    ground_ind = devices['ground_node'].index

    # delete all stamp entries in ground row
    non_ground_y_rows = np.nonzero(y_row != ground_ind)
    y_val = y_val[non_ground_y_rows]
    y_row = y_row[non_ground_y_rows]
    y_col = y_col[non_ground_y_rows]

    non_ground_j_rows = np.nonzero(j_row != ground_ind)
    j_val = j_val[non_ground_j_rows]
    j_row = j_row[non_ground_j_rows]
    # delete all stamp entries in ground row
    non_ground_y_cols = np.nonzero(y_col != ground_ind)
    y_val = y_val[non_ground_y_cols]
    y_row = y_row[non_ground_y_cols]
    y_col = y_col[non_ground_y_cols]

    (Y, J) = create_sparse_matrix(y_row, y_col, y_val, j_row, j_val, size_Y)

    return(Y, J)