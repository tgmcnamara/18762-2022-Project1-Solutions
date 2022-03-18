import numpy as np
from matplotlib import pyplot as plt
from classes.Nodes import Nodes
from classes.Inductors import Inductors
from classes.Capacitors import Capacitors
from scipy.io import loadmat
from scipy.interpolate import interp1d

def plot_results_linear_network(network_name, devices, V_hist, t_hist):
    if network_name == 'RL_circuit.mat':
        simulink_results_file = 'testcases/project1_rl_circuit_simulink_output.mat'
        voltages_to_plot = ['n3_a', 'n3_b', 'n3_c']
        currents_to_plot = ['L2a', 'L2b', 'L2c']
    elif network_name == 'RLC_circuit.mat':
        simulink_results_file = 'testcases/project1_bonus_rlc_circuit_simulink_output.mat'
        voltages_to_plot = ['n4_a', 'n4_b', 'n4_c']
        # The current exiting the C1's is equivalent to the sum of the currents going through
        # the load branches on each phase
        currents_to_plot = ['C1a', 'C1b', 'C1c']
    
    # Load simulink data for direct comparison of results
    simulink_data = loadmat(simulink_results_file)
    simulink_t = simulink_data['t_out'][0,:]
    interp_v = np.zeros((3,V_hist.shape[1]))
    interp_i = np.zeros((3,V_hist.shape[1]))
    v_error = np.zeros((3,V_hist.shape[1]))
    i_error = np.zeros((3,V_hist.shape[1]))

    # interpolate to our time domain so the timeseries are comparable
    for ind in range(3):
        simulink_v = simulink_data['voltages_out'][ind,:]
        simulink_i = simulink_data['currents_out'][ind,:]
        interp_v_func = interp1d(simulink_t, simulink_v)
        interp_v[ind,:] = interp_v_func(t_hist)
        interp_v_func = interp1d(simulink_t, simulink_i)
        interp_i[ind,:] = interp_v_func(t_hist)


    # plot voltages and calculate errors
    plot1 = plt.figure(1)
    for ind in range(3):
        node_name = voltages_to_plot[ind]
        vp_index = devices['nodes'][Nodes.node_key_[node_name]].index
        v_plot = V_hist[vp_index,:]
        i_error[ind,:] = np.abs(interp_v[ind,:] - v_plot)/np.max(np.abs(interp_v[ind,:]))
        plt.plot(t_hist, v_plot)

    #plt.plot(simulink_t, simulink_data['voltages_out'][0,:], simulink_t, simulink_data['voltages_out'][1,:], simulink_t, simulink_data['voltages_out'][2,:])
    plt.grid()
    plt.legend(voltages_to_plot)
    plt.xlabel('time (s)')
    plt.ylabel('V')
    plt.title("Three-phase voltage response for network %s" % network_name)
    
    # plot currents
    plot2 = plt.figure(2)
    for ind in range(3):
        current_device_name = currents_to_plot[ind]
        if current_device_name[0] == 'L':
            i_index = devices['inductors'][Inductors.inductor_key_[current_device_name]].companion_current_index
        else:
            i_index = devices['inductors'][Capacitors.capacitor_key_[current_device_name]].companion_current_index
        i_plot = V_hist[i_index,:]
        i_error[ind,:] = np.abs(interp_i[ind,:] - i_plot)/np.max(np.abs(interp_i[ind,:]))
        plt.plot(t_hist, i_plot)

    #plt.plot(simulink_t, simulink_data['currents_out'][0,:], simulink_t, simulink_data['currents_out'][1,:], simulink_t, simulink_data['currents_out'][2,:])
    plt.grid()
    plt.legend(['Phase A load current', 'Phase B load current', 'Phase C load current'])
    plt.xlabel('time (s)')
    plt.ylabel('A')
    plt.title("Three-phase current response for network %s" % network_name)
    
    # # plot simulink error
    plot3 = plt.figure(3)
    plt.plot(t_hist, v_error[0,:], t_hist, v_error[1,:], t_hist, v_error[2,:])
    plt.plot(t_hist, i_error[0,:], t_hist, i_error[1,:], t_hist, i_error[2,:])
    plt.xlabel('time (s)')
    plt.ylim([0, .1])
    plt.title("Error of generated signals compared to Simulink output for %s" % network_name)
    plt.ylabel('relative error compared to Simulink values')
    plt.legend(['Va', 'Vb', 'Vc', 'Ia', 'Ib', 'Ic'])
    plt.grid()
    plt.show()


def plot_results_tpim_network(network_name, devices, V_hist, t_hist, dt):

    # extract timeseries we need to plot
    tpim_hist = devices['induction_motors'][0].get_info_hist(V_hist, dt)

    # for the switch circuit, we were asked to only display the disturbance and response
    if network_name == 'IM_switch_circuit.mat':
        t0 = .65
        tf = 1.0
        simulink_results_file = 'testcases/im_switch_circuit_signals.mat'
    else:
        t0 = 0
        tf = 1.5
        simulink_results_file = 'testcases/im_simple_circuit_signals.mat'
    t = t_hist

    plot_simulink = True
    # load simulink data
    simulink_data = loadmat(simulink_results_file)
    simulink_t = simulink_data['t_out'][0,:]
    simulink_ids = simulink_data['out_signals'][1,:]
    simulink_iqs = simulink_data['out_signals'][0,:]
    simulink_idr = simulink_data['out_signals'][3,:]
    simulink_iqr = simulink_data['out_signals'][2,:]
    # matlab outputs mechanical speed (which is omega_r/pole pairs)
    simulink_omega = simulink_data['out_signals'][4,:]*(devices['induction_motors'][0].n_poles/2)
    simulink_te = simulink_data['out_signals'][5,:]
    simulink_vds = simulink_data['out_signals'][7,:]
    simulink_vqs = simulink_data['out_signals'][6,:]

    # plot stator currents
    plot3 = plt.figure(3)
    plt.plot(t, tpim_hist['ids_hist'])
    plt.plot(t, tpim_hist['iqs_hist'])
    if plot_simulink:
        plt.plot(simulink_t, simulink_ids)
        plt.plot(simulink_t, simulink_iqs)
        plt.legend([r'$i_{ds}$ local', r'$i_{qs}$ local', r'$i_{ds}$ simulink', r'$i_{qs}$ simulink'], loc='upper right')
    else:
        plt.legend([r'$i_{ds}$', r'$i_{qs}$'], loc='upper right')
    plt.xlabel('time (s)')
    plt.ylabel('Transformed stator currents (A)')
    plt.xlim([t0, tf])
    plt.grid()
    plt.title("TPIM stator currents for network %s" % network_name)
    # plot rotor currents
    plot4 = plt.figure(4)
    plt.plot(t, tpim_hist['idr_hist'])
    plt.plot(t, tpim_hist['iqr_hist'])
    if plot_simulink:
        plt.plot(simulink_t, simulink_idr)
        plt.plot(simulink_t, simulink_iqr)
        plt.legend([r'$i_{dr}$ local', r'$i_{qr}$ local', r'$i_{dr}$ simulink', r'$i_{qr}$ simulink'], loc='upper right')
    else:
        plt.legend([r'$i_{dr}$', r'$i_{qr}$'], loc='upper right')
    plt.xlim([t0, tf])
    plt.grid()
    plt.autoscale(enable=True, axis='y')
    plt.xlabel('time (s)')
    plt.ylabel('Transformed rotor currents (A)')
    plt.legend([r'$I_{dr}$', r'$I_{qr}$'], loc='upper right')
    plt.title("TPIM rotor currents for network %s" % network_name)
    plt.legend(['idr calc', 'iqr calc', 'idr simulink', 'iqr simulink'])
    # plot omega
    plot5 = plt.figure(5)
    plt.plot(t, (tpim_hist['omega_r_hist']))
    if plot_simulink:
        plt.plot(simulink_t, simulink_omega)
        plt.legend([r'$\omega_r$ local', r'$\omega_r$ simulink'])
    plt.xlabel('time (s)')
    plt.xlim([t0, tf])
    plt.grid()
    plt.autoscale(enable=True, axis='y')
    plt.ylabel(r'$\omega_r$ (rad/s))')
    plt.title("TPIM rotor frequency for network %s" % network_name)
    # plot torque
    plot6 = plt.figure(6)
    plt.plot(t, tpim_hist['te_hist'])
    if plot_simulink:
        plt.plot(simulink_t, simulink_te)
        plt.legend([r'$T_e$ local', r'$T_e$ simulink'])
    plt.xlabel('time (s)')
    plt.xlim([t0, tf])
    plt.grid()
    plt.autoscale(enable=True, axis='y')
    plt.ylabel(r'$T_e$ (N m)')
    plt.title("TPIM electrical torque for network %s" % network_name)
    plt.show()