def stamp_y_sparse(row_index, col_index, value, y_row, y_col, y_val, y_counter):
    y_val[y_counter] = value
    y_row[y_counter] = row_index
    y_col[y_counter] = col_index
    y_counter += 1
    return (y_counter)

def stamp_j_sparse(row_index, value, j_row, j_val, j_counter):
    j_val[j_counter] = value
    j_row[j_counter] = row_index
    j_counter += 1
    return (j_counter)
