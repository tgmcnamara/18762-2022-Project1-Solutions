import numpy as np
from lib.create_sparse_matrix import create_sparse_matrix

def stamp_constant_dense(devices, size_Y):
    Y = np.zeros((size_Y, size_Y), dtype=np.float)
    J = np.zeros((size_Y, 1), dtype=np.float)

    # resistors
    for ele in devices['resistors']:
        ele.stamp_dense(Y)
    # constant sources

    for ele in devices['switches']:
        ele.stamp_constant_dense(Y)

    ground_ind = devices['ground_node'].index
    # remove ground node stamps
    Y[ground_ind,:] = np.zeros((1,size_Y))
    Y[:,ground_ind] = np.zeros((size_Y,))
    Y[ground_ind, ground_ind] = 1
    J[ground_ind, 0] = 0

    return(Y, J)

def stamp_constant_sparse(devices, size_Y):
    num_Y_entries = size_Y*100
    y_val = np.zeros(num_Y_entries)
    y_row = np.zeros(num_Y_entries)
    y_col = np.zeros(num_Y_entries)
    y_entry_counter = 0

    j_val = np.zeros(num_Y_entries)
    j_row = np.zeros(num_Y_entries)
    j_entry_counter = 0

    # resistors
    for ele in devices['resistors']:
        y_entry_counter = ele.stamp_sparse(y_row, y_col, y_val, y_entry_counter)
   
    # constant sources

    # switches
    for ele in devices['switches']:
        y_entry_counter = ele.stamp_constant_sparse(y_row, y_col, y_val, y_entry_counter)

    y_val = y_val[:y_entry_counter]
    y_row = y_row[:y_entry_counter]
    y_col = y_col[:y_entry_counter]
    
    j_val = j_val[:j_entry_counter]
    j_row = j_row[:j_entry_counter]

    # delete all stamp entries in ground row
    ground_ind = devices['ground_node'].index
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

    # always stamp that ground node = 0
    y_val = np.append(y_val, 1)
    y_row = np.append(y_row, ground_ind)
    y_col = np.append(y_col, ground_ind)

    (Y, J) = create_sparse_matrix(y_row, y_col, y_val, j_row, j_val, size_Y)

    return(Y, J)