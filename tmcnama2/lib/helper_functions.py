# Some helper functions

def almost_equal(d1, d2, epsilon=1e-7):
    # used to avoid numerical float issues that can happen
    # with direct comparison in python sometimes
    return abs(d2 - d1) < epsilon

def check_error_type(err_max_ind, devices):
    # debugging function
    # given the index of the largest difference between V and V_sol,
    # this will tell you what device/node the index corresponds to
    error_info = [("Node %s" % ele.name) for ele in devices['nodes'] if ele.index == err_max_ind]
    if len(error_info) != 0:
        return error_info[0]
    
    error_info = [("Capacitor %s companion node" % ele.name) for ele in devices['capacitors'] if ele.companion_node_index == err_max_ind]
    if len(error_info) != 0:
        return error_info[0]

    error_info = [("Capacitor %s current" % ele.name) for ele in devices['capacitors'] if ele.companion_current_index == err_max_ind]
    if len(error_info) != 0:
        return error_info[0]

    error_info = [("Inductor %s companion node" % ele.name) for ele in devices['inductors'] if ele.companion_node_index == err_max_ind]
    if len(error_info) != 0:
        return error_info[0]

    error_info = [("Inductor %s current" % ele.name) for ele in devices['inductors'] if ele.companion_current_index == err_max_ind]
    if len(error_info) != 0:
        return error_info[0]

    error_info = [("Voltage source %s current" % ele.name) for ele in devices['ac_voltage_sources'] if ele.current_index == err_max_ind]
    if len(error_info) != 0:
        return error_info[0]

    error_info = [("Switch %s current" % ele.name) for ele in devices['switches'] if ele.switch_current_index == err_max_ind]
    if len(error_info) != 0:
        return error_info[0]

    error_info = [("IM %s ids" % ele.name) for ele in devices['induction_motors'] if ele.ids_index == err_max_ind]
    if len(error_info) != 0:
        return error_info[0]

    error_info = [("IM %s iqs" % ele.name) for ele in devices['induction_motors'] if ele.iqs_index == err_max_ind]
    if len(error_info) != 0:
        return error_info[0]

    error_info = [("IM %s idr" % ele.name) for ele in devices['induction_motors'] if ele.idr_index == err_max_ind]
    if len(error_info) != 0:
        return error_info[0]
    
    error_info = [("IM %s iqr" % ele.name) for ele in devices['induction_motors'] if ele.iqr_index == err_max_ind]
    if len(error_info) != 0:
        return error_info[0]

    error_info = [("IM %s omega_r" % ele.name) for ele in devices['induction_motors'] if ele.omega_r_index == err_max_ind]
    if len(error_info) != 0:
        return error_info[0]

    return ("Unable to find index")
