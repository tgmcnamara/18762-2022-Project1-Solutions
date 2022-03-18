import numpy as np
from lib.create_sparse_matrix import create_sparse_matrix

def stamp_time_dependent_dense(V, t, dt, devices, size_Y):
    Y = np.zeros((size_Y, size_Y), dtype=np.float)
    J = np.zeros((size_Y, 1), dtype=np.float)

    # capacitors
    for ele in devices['capacitors']:
        ele.stamp_dense(V, dt, Y, J)
    # inductors
    for ele in devices['inductors']:
        ele.stamp_dense(V, dt, Y, J)
    # switches
    for ele in devices['switches']:
        ele.stamp_time_dependent_dense(t, Y)
    # sources
    for ele in devices['ac_voltage_sources']:
        ele.stamp_dense(t, Y, J)

    return(Y, J)

def stamp_time_dependent_sparse(V, t, dt, devices, size_Y):
    num_Y_entries = size_Y*100

    y_val = np.zeros((num_Y_entries))
    y_row = np.zeros((num_Y_entries))
    y_col = np.zeros((num_Y_entries))
    y_counter = 0

    j_val = np.zeros((num_Y_entries))
    j_row = np.zeros((num_Y_entries))
    j_counter = 0

    # capacitors
    for ele in devices['capacitors']:
        (y_counter, j_counter) = ele.stamp_sparse(V, dt, y_row, y_col, y_val, y_counter, j_row, j_val, j_counter)
    # inductors
    for ele in devices['inductors']:
        (y_counter, j_counter) = ele.stamp_sparse(V, dt, y_row, y_col, y_val, y_counter, j_row, j_val, j_counter)
    # switches
    for ele in devices['switches']:
        (y_counter) = ele.stamp_time_dependent_sparse(t, y_row, y_col, y_val, y_counter)
    # stamp sources
    for ele in devices['ac_voltage_sources']:
        (y_counter, j_counter) = ele.stamp_sparse(t, y_row, y_col, y_val, y_counter, j_row, j_val, j_counter)

    y_val = y_val[:y_counter]
    y_row = y_row[:y_counter]
    y_col = y_col[:y_counter]
    
    j_val = j_val[:j_counter]
    j_row = j_row[:j_counter]

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