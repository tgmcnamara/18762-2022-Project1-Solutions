import numpy as np
import scipy
from scikits.umfpack import spsolve
from scipy.linalg import solve
from copy import deepcopy
from time import time
from lib.process_case import process_case
from lib.initialize import initialize
from lib.stamp_constant_terms import stamp_constant_sparse, stamp_constant_dense
from lib.stamp_time_dep_terms import stamp_time_dependent_dense, stamp_time_dependent_sparse
from lib.stamp_nonlinear_terms import stamp_nonlinear_dense, stamp_nonlinear_sparse
from lib.stamp_t0_terms import stamp_t0_sparse, stamp_t0_dense
from lib.helper_functions import check_error_type, almost_equal
from classes.Nodes import Nodes
from lib.plot_results import plot_results_linear_network, plot_results_tpim_network

def main(network_file, dt, tolerance, use_sparse):
    # parse file and create devices
    devices = process_case(network_file, use_sparse)

    # initialize settings for simulation
    t_start = 0
    iteration_count = 0
    inner_loop_count = 0
    time_at_start = time()
    contains_nonlinear_devices = (len(devices['induction_motors']) > 0)
    # run for .2 s for RL and RLC networks, 1s for TPIM networks
    if contains_nonlinear_devices:
        t_final = 1.0-dt
    else:
        t_final = 0.2-dt

    # Turn on for debugging 
    print_error_info = False
    
    # initialize Y, V, and J
    size_Y = Nodes._index_counter.__next__()
    V = initialize(devices, size_Y)
    # generate constant matrices
    if use_sparse:
        (Y_constant, J_constant) = stamp_constant_sparse(devices, size_Y)
    else:
        (Y_constant, J_constant) = stamp_constant_dense(devices, size_Y)
    V_hist = np.empty(0, dtype=np.float)

    t = t_start
    # solve for t=0
    if use_sparse:
        (Y_t0, J_t0) = stamp_t0_sparse(devices, size_Y)
    else:
        (Y_t0, J_t0) = stamp_t0_dense(devices, size_Y)

    Y = Y_constant + Y_t0
    J = J_constant + J_t0

    if use_sparse:
        V_sol = spsolve(Y, J).reshape(-1,1)
    else:
        V_sol = solve(Y, J).reshape(-1,1)

    t_hist = [0]
    V_hist = V_sol
    V = V_sol
    V_prev = deepcopy(V)
    t_backlash_log = []
    sparsity_log = []

    while True:
        # increase t
        t += dt

        # update time-dependent devices
        if use_sparse:
            (Y_time_dep, J_time_dep) = stamp_time_dependent_sparse(V_prev, t, dt, devices, size_Y)
        else:
            (Y_time_dep, J_time_dep) = stamp_time_dependent_dense(V_prev, t, dt, devices, size_Y)

        # solve circuit
        if contains_nonlinear_devices:
            inner_loop_count = 0
            #V = V_prev
            max_error = tolerance*1.1
            while True:
                inner_loop_count += 1
                if use_sparse:
                    (Y_nonlinear, J_nonlinear) = stamp_nonlinear_sparse(V, V_prev, dt, devices, size_Y)
                    sparsity_log.append(Y.nnz/(size_Y**2))
                else:
                    (Y_nonlinear, J_nonlinear) = stamp_nonlinear_dense(V, V_prev, dt, devices, size_Y)
                    sparsity_log.append(np.count_nonzero(Y)/size_Y**2)
                Y = Y_constant + Y_time_dep + Y_nonlinear
                J = J_constant + J_time_dep + J_nonlinear
                t_backslash = time()
                if use_sparse:
                    V_sol = spsolve(Y, J).reshape(-1,1)
                else:
                    V_sol = solve(Y, J).reshape(-1,1)
                t_backlash_log.append((time() - t_backslash))
                max_error = np.amax(np.abs(V_sol-V))
                if print_error_info:
                    max_error_index = np.argmax(np.abs(V_sol - V))
                    err_info = check_error_type(max_error_index, devices)
                    print("Max error: %f, %s" % (max_error, err_info))
                V = deepcopy(V_sol)
                if max_error < tolerance:
                    devices['induction_motors'][0].check_exact(V_sol, V_prev, dt)
                    print("NR converged at t=%f in %d iterations" % (t, inner_loop_count))
                    break
            iteration_count += inner_loop_count
        else:
            Y = Y_constant + Y_time_dep
            J = J_constant + J_time_dep
            
            if use_sparse:
                sparsity_log.append(Y.nnz/(size_Y**2))
                t_backslash = time()
                V_sol = spsolve(Y, J).reshape(-1,1)
            else:
                sparsity_log.append(np.count_nonzero(Y)/size_Y**2)
                t_backslash = time()
                V_sol = solve(Y, J).reshape(-1,1)
            iteration_count += 1
            t_backlash_log.append((time() - t_backslash))
                            
        # log info
        t_hist.append(t)
        V_hist = np.concatenate((V_hist, V_sol), axis=1)
        V_prev = V_sol
        if almost_equal(t, t_final):
            break

    time_at_finish = time()
    print("Simulation took %f s, %d iterations" % ((time_at_finish-time_at_start), iteration_count))
    print("Use sparse: %s, #nodes: %d, avg backslash operation: %f s, avg sparsity: %.2f percent" % 
            (str(use_sparse),size_Y, np.mean(t_backlash_log), 100*np.mean(sparsity_log)))
    network_name = network_file.split('/')[-1]
    if network_name in ['IM_circuit.mat', 'IM_switch_circuit.mat']:
        plot_results_tpim_network(network_name, devices, V_hist, t_hist, dt)
    elif network_name in ['RL_circuit.mat', 'RLC_circuit.mat']:
        plot_results_linear_network(network_name, devices, V_hist, t_hist)
    else:
        print("Plotting not available for new network.")
    
