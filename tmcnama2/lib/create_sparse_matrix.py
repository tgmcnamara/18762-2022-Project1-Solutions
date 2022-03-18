from scipy.sparse import csc_matrix
import numpy as np

def create_sparse_matrix(y_row, 
                         y_col, 
                         y_val,
                         j_row,
                         j_val, 
                         size_Y):
    # generate a CSC matrix given lists of row indices, column indices, and values
    # y_row, y_col, y_val must all be same size
    # j_row, j_val must be same size
    # size_Y declares the dimension of the matrix/vector overall
    Y = csc_matrix((y_val, (y_row, y_col)), shape=(size_Y, size_Y), dtype=np.float64)
    J = csc_matrix((j_val, (j_row, np.zeros((j_row.shape[0],)))), shape=(size_Y, 1), dtype=np.float64)

    return((Y, J))
