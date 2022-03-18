import numpy as np
from itertools import count
from classes.Nodes import Nodes
from lib.stamping_functions import stamp_y_sparse, stamp_j_sparse

class InductionMotors:
    _motor_counter = count(0)
    motor_key_ = {}

    def __init__(self,
                 name,
                 phase_a_node,
                 phase_b_node,
                 phase_c_node,
                 neutral_node,
                 power_nom,
                 v_nom,
                 frequency,
                 lm,
                 rs,
                 rr,
                 lls,
                 llr,
                 j_im,
                 tm,
                 d_fric):
        self.name = name
        self.phase_a_node = phase_a_node
        self.phase_b_node = phase_b_node
        self.phase_c_node = phase_c_node
        self.neutral_node = neutral_node
        self.phase_a_node_index = None
        self.phase_b_node_index = None
        self.phase_c_node_index = None
        self.neutral_node_index = None
        self.id = self._motor_counter.__next__()
        InductionMotors.motor_key_[(self.phase_a_node, self.phase_b_node, self.phase_c_node, self.name)] = self.id
        self.power_nom = power_nom[0]
        self.v_nom = v_nom[0]
        self.motor_freq = frequency[0]
        self.lm = lm[0]
        self.rs = rs[0]
        self.rr = rr[0]
        self.lls = lls[0]
        self.llr = llr[0]
        self.lss = self.lls + self.lm
        self.lrr = self.llr + self.lm
        self.j_im = j_im[0]
        self.tm = tm[0]
        self.n_poles = 4
        # To match the Simulink model, the friction coefficient needs to be scaled
        # by 1/n_poles (F is applied to the mechanical omega, which is 1/p*omega_r)
        self.d_fric = d_fric[0]/self.n_poles
        # Slide deck said 3/4*Lm*n_poles but Simulink model uses
        # Te = 3/2*Lm*n_poles and that causes much faster
        # settling time for me, which matches Simulink better
        self.Te_alpha = 3/2*self.lm*self.n_poles
        self.vds_index = None
        self.vqs_index = None
        self.ids_index = None
        self.idr_index = None
        self.iqs_index = None
        self.iqr_index = None
        self.omega_r_index = None
        lambd = np.pi * 2/3
        # we're in the stationary frame
        self.ref_theta = 0
        self.p_transform = 2/3*np.array([[0.5, 0.5, 0.5],                                                                       # 0 component
                                         [np.cos(self.ref_theta), np.cos(self.ref_theta-lambd), np.cos(self.ref_theta+lambd)],  # d component
                                         [np.sin(self.ref_theta), np.sin(self.ref_theta-lambd), np.sin(self.ref_theta+lambd)]]) # q component
        # This is the transform simulink uses, where the q current is aligned with phase a
        # self.p_transform = 2/3*np.array([[0.5, 0.5, 0.5], 
        #                                  [np.sin(self.ref_theta), np.sin(self.ref_theta-lambd), np.sin(self.ref_theta+lambd)],
        #                                  [np.cos(self.ref_theta), np.cos(self.ref_theta-lambd), np.cos(self.ref_theta+lambd)]])
        self.p_inv = np.linalg.inv(self.p_transform)

    def assign_indexes(self, node_list):
        phase_a_node_obj = node_list[Nodes.node_key_[self.phase_a_node]]
        self.phase_a_node_index = phase_a_node_obj.index
        phase_b_node_obj = node_list[Nodes.node_key_[self.phase_b_node]]
        self.phase_b_node_index = phase_b_node_obj.index
        phase_c_node_obj = node_list[Nodes.node_key_[self.phase_c_node]]
        self.phase_c_node_index = phase_c_node_obj.index
        neutral_node_obj = node_list[Nodes.node_key_[self.neutral_node]]
        self.neutral_node_index = neutral_node_obj.index
        
        self.vds_index = Nodes._index_counter.__next__()
        self.vqs_index = Nodes._index_counter.__next__()
        self.ids_index = Nodes._index_counter.__next__()
        self.idr_index = Nodes._index_counter.__next__()
        self.iqs_index = Nodes._index_counter.__next__()
        self.iqr_index = Nodes._index_counter.__next__()
        self.omega_r_index = Nodes._index_counter.__next__()

    def check_exact(self, V, V_prev, dt):
        # (vds, vqs) = self.calc_V0dq(V)
        vds = V[self.vds_index,0]
        vqs = V[self.vqs_index,0]
        ids = V[self.ids_index,0]
        idr = V[self.idr_index,0]
        iqs = V[self.iqs_index,0]
        iqr = V[self.iqr_index,0]
        omega_r = V[self.omega_r_index,0]
        
        #(vds_prev, vqs_prev) = self.calc_V0dq(V_prev)
        vds_prev = V_prev[self.vds_index,0]
        vqs_prev = V_prev[self.vqs_index,0]
        ids_prev = V_prev[self.ids_index,0]
        idr_prev = V_prev[self.idr_index,0]
        iqs_prev = V_prev[self.iqs_index,0]
        iqr_prev = V_prev[self.iqr_index,0]
        omega_r_prev = V_prev[self.omega_r_index,0]

        vds_eq = -vds + (2*self.lss/dt + self.rs)*ids + 2*self.lm/dt*idr + (-vds_prev + (self.rs-2*self.lss/dt)*ids_prev - 2*self.lm/dt*idr_prev)
        vdr_eq = 2*self.lm/dt*ids + (2*self.lrr/dt + self.rr)*idr - self.lm*omega_r*iqs - self.lrr*omega_r*iqr + (-2*self.lm/dt*ids_prev + (self.rr-2*self.lrr/dt)*idr_prev - omega_r_prev*(self.lm*iqs_prev + self.lrr*iqr_prev))
        vqs_eq = -vqs + (2*self.lss/dt + self.rs)*iqs + 2*self.lm/dt*iqr + (-vqs_prev + (self.rs - 2*self.lss/dt)*iqs_prev - 2*self.lm/dt*iqr_prev)
        vqr_eq = (self.rr + 2*self.lrr/dt)*iqr + 2*self.lm/dt*iqs + omega_r*(self.lm*ids + self.lrr*idr) + ((self.rr - 2*self.lrr/dt)*iqr_prev - 2*self.lm/dt*iqs_prev + omega_r_prev*((self.lm*ids_prev + self.lrr*idr_prev)))

        g_eq = 2*self.j_im/dt
        Te = self.Te_alpha*(idr*iqs - iqr*ids)
        Te_prev = self.Te_alpha*(idr_prev*iqs_prev - iqr_prev*ids_prev)
        omega_eq = (self.d_fric + g_eq)*omega_r + Te + ((self.d_fric - g_eq)*omega_r_prev + Te_prev + 2*self.tm)
        if np.amax(np.abs([vds_eq, vdr_eq, vqs_eq, vqr_eq, omega_eq])) > 1e-8:
            print("WARNING: one of the TPIM nonlinear equations did not converge to less than 1e-8")

    def calc_tpim_stamps(self, V, V_prev, dt):
        ids = V[self.ids_index,0]
        idr = V[self.idr_index,0]
        iqs = V[self.iqs_index,0]
        iqr = V[self.iqr_index,0]
        omega_r = V[self.omega_r_index,0]
        
        vds_prev = V_prev[self.vds_index,0]
        vqs_prev = V_prev[self.vqs_index,0]
        ids_prev = V_prev[self.ids_index,0]
        idr_prev = V_prev[self.idr_index,0]
        iqs_prev = V_prev[self.iqs_index,0]
        iqr_prev = V_prev[self.iqr_index,0]
        omega_r_prev = V_prev[self.omega_r_index,0]

        stamp_dict = {}

        # vds equation
        stamp_dict['dfds_dids'] = self.rs + 2*self.lss/dt
        stamp_dict['dfds_didr'] = 2*self.lm/dt
        # includes -vds --> second row of P*Vabc
        stamp_dict['dfds_dvds'] = -1
        vds_hist = (-vds_prev + (self.rs - 2*self.lss/dt)*ids_prev - 2*self.lm/dt*idr_prev)
        #print("Vds_hist: %.10f" % vds_hist)
        stamp_dict['fds_jstamp'] = -vds_hist

        # vdr equation
        stamp_dict['dfdr_dids'] = 2*self.lm/dt
        stamp_dict['dfdr_didr'] = self.rr + 2*self.lrr/dt
        stamp_dict['dfdr_diqs'] = -self.lm*omega_r
        stamp_dict['dfdr_diqr'] = -self.lrr*omega_r
        stamp_dict['dfdr_dfomega'] = -self.lm*iqs - self.lrr*iqr
        vdr_hist = (-2*self.lm/dt*ids_prev + (self.rr-2*self.lrr/dt)*idr_prev - omega_r_prev*(self.lm*iqs_prev + self.lrr*iqr_prev))
        #print("Vdr_hist: %.10f" % vdr_hist)
        stamp_dict['fdr_jstamp'] = -vdr_hist - omega_r*(self.lm*iqs + self.lrr*iqr)
        
        # vqs equation
        stamp_dict['dfqs_diqs'] = self.rs + 2*self.lss/dt
        stamp_dict['dfqs_diqr'] = 2*self.lm/dt
        stamp_dict['dfqs_dvqs'] = -1
        vqs_hist = (-vqs_prev + (self.rs - 2*self.lss/dt)*iqs_prev - 2*self.lm/dt*iqr_prev)
        #print("Vqs_hist: %.10f" % vqs_hist)
        stamp_dict['fqs_jstamp'] = -vqs_hist
        
        # vqr equation
        stamp_dict['dfqr_dids'] = self.lm*omega_r
        stamp_dict['dfqr_didr'] = self.lrr*omega_r
        stamp_dict['dfqr_diqs'] = 2*self.lm/dt
        stamp_dict['dfqr_diqr'] = self.rr + 2*self.lrr/dt
        stamp_dict['dfqr_dfomega'] = self.lm*ids + self.lrr*idr
        vqr_hist = ((self.rr - 2*self.lrr/dt)*iqr_prev - 2*self.lm/dt*iqs_prev + omega_r_prev*((self.lm*ids_prev + self.lrr*idr_prev)))
        #print("Vqr_hist: %.10f" % vqr_hist)
        #print('-'*10)
        stamp_dict['fqr_jstamp'] = -vqr_hist + omega_r*(self.lm*ids + self.lrr*idr)

        # omega equation
        # using norton equivalent
        # g_eq = 2*self.j_im/dt
        # # Te = self.Te_alpha*(idr*iqs - iqr*ids)
        # Te_prev = self.Te_alpha*(idr_prev*iqs_prev - iqr_prev*ids_prev)
        # stamp_dict['dfomega_dids'] = -self.Te_alpha*iqr
        # stamp_dict['dfomega_didr'] = self.Te_alpha*iqs
        # stamp_dict['dfomega_diqs'] = self.Te_alpha*idr
        # stamp_dict['dfomega_diqr'] = -self.Te_alpha*ids
        # stamp_dict['dfomega_dfomega'] = self.d_fric + g_eq
        # omega_hist = (self.d_fric - g_eq)*omega_r_prev - Te_prev + 2*self.tm
        # stamp_dict['fomega_jstamp'] = -omega_hist - self.Te_alpha*(idr*iqs - iqr*ids)
        g_eq = 2*self.j_im/dt
        # Te = self.Te_alpha*(idr*iqs - iqr*ids)
        Te_prev = self.Te_alpha*(idr_prev*iqs_prev - iqr_prev*ids_prev)
        stamp_dict['dfomega_dids'] = -self.Te_alpha*iqr
        stamp_dict['dfomega_didr'] = self.Te_alpha*iqs
        stamp_dict['dfomega_diqs'] = self.Te_alpha*idr
        stamp_dict['dfomega_diqr'] = -self.Te_alpha*ids
        stamp_dict['dfomega_dfomega'] = self.d_fric + g_eq
        omega_hist = (self.d_fric - g_eq)*omega_r_prev + Te_prev + 2*self.tm
        stamp_dict['fomega_jstamp'] = -omega_hist + self.Te_alpha*(idr*iqs - iqr*ids)

        # voltage controlled voltage source equations
        stamp_dict['dvds_dva'] = self.p_transform[1,0]
        stamp_dict['dvds_dvb'] = self.p_transform[1,1]
        stamp_dict['dvds_dvc'] = self.p_transform[1,2]

        stamp_dict['dvqs_dva'] = self.p_transform[2,0]
        stamp_dict['dvqs_dvb'] = self.p_transform[2,1]
        stamp_dict['dvqs_dvc'] = self.p_transform[2,2]

        # current controlled current source equations
        # apply inverse park transform to I0dq_stator to get i_a, i_b, i_c
        # v_phase_a row: i_tpim_a = ids 
        stamp_dict['dia_dids'] = self.p_inv[0,1] # 1
        stamp_dict['dia_diqs'] = self.p_inv[0,2] # 0
        # v_phase_b row: i_tpim_b = -.5*ids - sqrt(3)/2*iqs
        stamp_dict['dib_dids'] = self.p_inv[1,1]#-.5
        stamp_dict['dib_diqs'] = self.p_inv[1,2]#-np.sqrt(3)/2
        # v_phase_c row: i_tpim_c = -.5*ids + sqrt(3)/2*iqs
        stamp_dict['dic_dids'] = self.p_inv[2,1] #-.5
        stamp_dict['dic_diqs'] = self.p_inv[2,2] #np.sqrt(3)/2
        return stamp_dict

    def stamp_sparse(self, V, V_prev, dt, y_row, y_col, y_val, y_entry_counter, j_row, j_val, j_entry_counter):
        stamp_dict = self.calc_tpim_stamps(V, V_prev, dt)
        
        # ids row - fds equation
        y_entry_counter = stamp_y_sparse(self.ids_index, self.ids_index, stamp_dict['dfds_dids'], y_row, y_col, y_val, y_entry_counter)
        y_entry_counter = stamp_y_sparse(self.ids_index, self.idr_index, stamp_dict['dfds_didr'], y_row, y_col, y_val, y_entry_counter)
        y_entry_counter = stamp_y_sparse(self.ids_index, self.vds_index, stamp_dict['dfds_dvds'], y_row, y_col, y_val, y_entry_counter)
        j_entry_counter = stamp_j_sparse(self.ids_index, stamp_dict['fds_jstamp'], j_row, j_val, j_entry_counter)

        # idr row - fdr equation
        y_entry_counter = stamp_y_sparse(self.idr_index, self.ids_index, stamp_dict['dfdr_dids'], y_row, y_col, y_val, y_entry_counter)
        y_entry_counter = stamp_y_sparse(self.idr_index, self.idr_index, stamp_dict['dfdr_didr'], y_row, y_col, y_val, y_entry_counter)
        y_entry_counter = stamp_y_sparse(self.idr_index, self.iqs_index, stamp_dict['dfdr_diqs'], y_row, y_col, y_val, y_entry_counter)
        y_entry_counter = stamp_y_sparse(self.idr_index, self.iqr_index, stamp_dict['dfdr_diqr'], y_row, y_col, y_val, y_entry_counter)
        y_entry_counter = stamp_y_sparse(self.idr_index, self.omega_r_index, stamp_dict['dfdr_dfomega'], y_row, y_col, y_val, y_entry_counter)
        j_entry_counter = stamp_j_sparse(self.idr_index, stamp_dict['fdr_jstamp'], j_row, j_val, j_entry_counter)

        # iqs row - fqs equation
        y_entry_counter = stamp_y_sparse(self.iqs_index, self.iqs_index, stamp_dict['dfqs_diqs'], y_row, y_col, y_val, y_entry_counter)
        y_entry_counter = stamp_y_sparse(self.iqs_index, self.iqr_index, stamp_dict['dfqs_diqr'], y_row, y_col, y_val, y_entry_counter)
        y_entry_counter = stamp_y_sparse(self.iqs_index, self.vqs_index, stamp_dict['dfqs_dvqs'], y_row, y_col, y_val, y_entry_counter)
        j_entry_counter = stamp_j_sparse(self.iqs_index, stamp_dict['fqs_jstamp'], j_row, j_val, j_entry_counter)
        
        # iqr row - fqr equation
        y_entry_counter = stamp_y_sparse(self.iqr_index, self.ids_index, stamp_dict['dfqr_dids'], y_row, y_col, y_val, y_entry_counter)
        y_entry_counter = stamp_y_sparse(self.iqr_index, self.idr_index, stamp_dict['dfqr_didr'], y_row, y_col, y_val, y_entry_counter)
        y_entry_counter = stamp_y_sparse(self.iqr_index, self.iqs_index, stamp_dict['dfqr_diqs'], y_row, y_col, y_val, y_entry_counter)
        y_entry_counter = stamp_y_sparse(self.iqr_index, self.iqr_index, stamp_dict['dfqr_diqr'], y_row, y_col, y_val, y_entry_counter)
        y_entry_counter = stamp_y_sparse(self.iqr_index, self.omega_r_index, stamp_dict['dfqr_dfomega'], y_row, y_col, y_val, y_entry_counter)
        j_entry_counter = stamp_j_sparse(self.iqr_index, stamp_dict['fqr_jstamp'], j_row, j_val, j_entry_counter)

        # omega_r row - f_omega equation
        y_entry_counter = stamp_y_sparse(self.omega_r_index, self.ids_index, stamp_dict['dfomega_dids'], y_row, y_col, y_val, y_entry_counter)
        y_entry_counter = stamp_y_sparse(self.omega_r_index, self.idr_index, stamp_dict['dfomega_didr'], y_row, y_col, y_val, y_entry_counter)
        y_entry_counter = stamp_y_sparse(self.omega_r_index, self.iqs_index, stamp_dict['dfomega_diqs'], y_row, y_col, y_val, y_entry_counter)
        y_entry_counter = stamp_y_sparse(self.omega_r_index, self.iqr_index, stamp_dict['dfomega_diqr'], y_row, y_col, y_val, y_entry_counter)
        y_entry_counter = stamp_y_sparse(self.omega_r_index, self.omega_r_index, stamp_dict['dfomega_dfomega'], y_row, y_col, y_val, y_entry_counter)
        j_entry_counter = stamp_j_sparse(self.omega_r_index, stamp_dict['fomega_jstamp'], j_row, j_val, j_entry_counter)

        # voltage controlled voltage sources
        #     dvds_dva*va + dvds_dvb*vb + dvds_dvc*vc = vds - v_neutral
        # --> dvds_dva*va + dvds_dvb*vb + dvds_dvc*vc - vds + v_neutral = 0
        y_entry_counter = stamp_y_sparse(self.vds_index, self.phase_a_node_index, stamp_dict['dvds_dva'], y_row, y_col, y_val, y_entry_counter)
        y_entry_counter = stamp_y_sparse(self.vds_index, self.phase_b_node_index, stamp_dict['dvds_dvb'], y_row, y_col, y_val, y_entry_counter)
        y_entry_counter = stamp_y_sparse(self.vds_index, self.phase_c_node_index, stamp_dict['dvds_dvc'], y_row, y_col, y_val, y_entry_counter)
        y_entry_counter = stamp_y_sparse(self.vds_index, self.vds_index, -1, y_row, y_col, y_val, y_entry_counter)
        y_entry_counter = stamp_y_sparse(self.vds_index, self.neutral_node_index, 1, y_row, y_col, y_val, y_entry_counter)
        #     dvqs_dva*va + dvqs_dvb*vb + dvqs_dvc*vc = vqs - v_neutral
        # --> dvqs_dva*va + dvqs_dvb*vb + dvqs_dvc*vc - vqs + v_neutral = 0
        y_entry_counter = stamp_y_sparse(self.vqs_index, self.phase_a_node_index, stamp_dict['dvqs_dva'], y_row, y_col, y_val, y_entry_counter)
        y_entry_counter = stamp_y_sparse(self.vqs_index, self.phase_b_node_index, stamp_dict['dvqs_dvb'], y_row, y_col, y_val, y_entry_counter)
        y_entry_counter = stamp_y_sparse(self.vqs_index, self.phase_c_node_index, stamp_dict['dvqs_dvc'], y_row, y_col, y_val, y_entry_counter)
        y_entry_counter = stamp_y_sparse(self.vqs_index, self.vqs_index, -1, y_row, y_col, y_val, y_entry_counter)
        y_entry_counter = stamp_y_sparse(self.vqs_index, self.neutral_node_index, 1, y_row, y_col, y_val, y_entry_counter)

        # current controlled current sources - all currents are LEAVING phase node (i.e. positive)
        #                                      and going INTO ground node
        y_entry_counter = stamp_y_sparse(self.phase_a_node_index, self.ids_index, stamp_dict['dia_dids'], y_row, y_col, y_val, y_entry_counter)
        y_entry_counter = stamp_y_sparse(self.neutral_node_index, self.ids_index, -stamp_dict['dia_dids'], y_row, y_col, y_val, y_entry_counter)
        y_entry_counter = stamp_y_sparse(self.phase_a_node_index, self.iqs_index, stamp_dict['dia_diqs'], y_row, y_col, y_val, y_entry_counter)
        y_entry_counter = stamp_y_sparse(self.neutral_node_index, self.iqs_index, -stamp_dict['dia_diqs'], y_row, y_col, y_val, y_entry_counter)
        y_entry_counter = stamp_y_sparse(self.phase_b_node_index, self.ids_index, stamp_dict['dib_dids'], y_row, y_col, y_val, y_entry_counter)
        y_entry_counter = stamp_y_sparse(self.phase_b_node_index, self.iqs_index, stamp_dict['dib_diqs'], y_row, y_col, y_val, y_entry_counter)
        y_entry_counter = stamp_y_sparse(self.neutral_node_index, self.ids_index, -stamp_dict['dib_dids'], y_row, y_col, y_val, y_entry_counter) 
        y_entry_counter = stamp_y_sparse(self.neutral_node_index, self.iqs_index, -stamp_dict['dib_diqs'], y_row, y_col, y_val, y_entry_counter) 
        y_entry_counter = stamp_y_sparse(self.phase_c_node_index, self.ids_index, stamp_dict['dic_dids'], y_row, y_col, y_val, y_entry_counter)
        y_entry_counter = stamp_y_sparse(self.phase_c_node_index, self.iqs_index, stamp_dict['dic_diqs'], y_row, y_col, y_val, y_entry_counter) 
        y_entry_counter = stamp_y_sparse(self.neutral_node_index, self.ids_index, -stamp_dict['dic_dids'], y_row, y_col, y_val, y_entry_counter)
        y_entry_counter = stamp_y_sparse(self.neutral_node_index, self.iqs_index, -stamp_dict['dic_diqs'], y_row, y_col, y_val, y_entry_counter)
        return (y_entry_counter, j_entry_counter)

    def stamp_dense(self, V, V_prev, dt, Y, J):
        stamp_dict = self.calc_tpim_stamps(V, V_prev, dt)
        
        # ids row - vds equation
        Y[self.ids_index, self.ids_index] += stamp_dict['dfds_dids']
        Y[self.ids_index, self.idr_index] += stamp_dict['dfds_didr']
        Y[self.ids_index, self.vds_index] += stamp_dict['dfds_dvds']
        J[self.ids_index, 0] += stamp_dict['fds_jstamp']

        # idr row - vdr equation
        Y[self.idr_index, self.ids_index] += stamp_dict['dfdr_dids']
        Y[self.idr_index, self.idr_index] += stamp_dict['dfdr_didr']
        Y[self.idr_index, self.iqs_index] += stamp_dict['dfdr_diqs']
        Y[self.idr_index, self.iqr_index] += stamp_dict['dfdr_diqr']
        Y[self.idr_index, self.omega_r_index] += stamp_dict['dfdr_dfomega']
        J[self.idr_index, 0] += stamp_dict['fdr_jstamp']

        # iqs row - vqs
        Y[self.iqs_index, self.iqs_index] += stamp_dict['dfqs_diqs']
        Y[self.iqs_index, self.iqr_index] += stamp_dict['dfqs_diqr']
        Y[self.iqs_index, self.vqs_index] += stamp_dict['dfqs_dvqs']
        J[self.iqs_index, 0] += stamp_dict['fqs_jstamp']
        
        # iqr row - vqr
        Y[self.iqr_index, self.ids_index] += stamp_dict['dfqr_dids']
        Y[self.iqr_index, self.idr_index] += stamp_dict['dfqr_didr']
        Y[self.iqr_index, self.iqs_index] += stamp_dict['dfqr_diqs']
        Y[self.iqr_index, self.iqr_index] += stamp_dict['dfqr_diqr']
        Y[self.iqr_index, self.omega_r_index] += stamp_dict['dfqr_dfomega'] 
        J[self.iqr_index, 0] += stamp_dict['fqr_jstamp']

        # omega_r row - omega_r
        Y[self.omega_r_index, self.ids_index] += stamp_dict['dfomega_dids']
        Y[self.omega_r_index, self.idr_index] += stamp_dict['dfomega_didr']
        Y[self.omega_r_index, self.iqs_index] += stamp_dict['dfomega_diqs']
        Y[self.omega_r_index, self.iqr_index] += stamp_dict['dfomega_diqr']
        Y[self.omega_r_index, self.omega_r_index] += stamp_dict['dfomega_dfomega']
        J[self.omega_r_index, 0] += stamp_dict['fomega_jstamp']

        # voltage controlled voltage sources
        #     dvds_dva*va + dvds_dvb*vb + dvds_dvc*vc = vds - v_neutral
        # --> dvds_dva*va + dvds_dvb*vb + dvds_dvc*vc - vds + v_neutral = 0
        Y[self.vds_index, self.phase_a_node_index] += stamp_dict['dvds_dva']
        Y[self.vds_index, self.phase_b_node_index] += stamp_dict['dvds_dvb']
        Y[self.vds_index, self.phase_c_node_index] += stamp_dict['dvds_dvc']
        Y[self.vds_index, self.vds_index] += -1
        Y[self.vds_index, self.neutral_node_index] += 1
        #     dvqs_dva*va + dvqs_dvb*vb + dvqs_dvc*vc = vqs - v_neutral
        # --> dvqs_dva*va + dvqs_dvb*vb + dvqs_dvc*vc - vqs + v_neutral = 0
        Y[self.vqs_index, self.phase_a_node_index] += stamp_dict['dvqs_dva']
        Y[self.vqs_index, self.phase_b_node_index] += stamp_dict['dvqs_dvb']
        Y[self.vqs_index, self.phase_c_node_index] += stamp_dict['dvqs_dvc']
        Y[self.vqs_index, self.vqs_index] += -1
        Y[self.vqs_index, self.neutral_node_index] += 1
        
        # current controlled current sources - all currents are LEAVING phase node (i.e. positive)
        #                                      and going INTO ground node (i.e. negative)
        Y[self.phase_a_node_index, self.ids_index] += stamp_dict['dia_dids']
        Y[self.phase_a_node_index, self.iqs_index] += stamp_dict['dia_diqs']
        Y[self.neutral_node_index, self.ids_index] += -stamp_dict['dia_dids']   
        Y[self.neutral_node_index, self.iqs_index] += -stamp_dict['dia_diqs']        
        Y[self.phase_b_node_index, self.ids_index] += stamp_dict['dib_dids'] 
        Y[self.phase_b_node_index, self.iqs_index] += stamp_dict['dib_diqs']
        Y[self.neutral_node_index, self.ids_index] += -stamp_dict['dib_dids'] 
        Y[self.neutral_node_index, self.iqs_index] += -stamp_dict['dib_diqs'] 
        Y[self.phase_c_node_index, self.ids_index] += stamp_dict['dic_dids']
        Y[self.phase_c_node_index, self.iqs_index] += stamp_dict['dic_diqs'] 
        Y[self.neutral_node_index, self.ids_index] += -stamp_dict['dic_dids']
        Y[self.neutral_node_index, self.iqs_index] += -stamp_dict['dic_diqs']

    def stamp_t0_sparse(self, y_row, y_col, y_val, y_entry_counter):
        # initial conditions of TPIM are all currents=0, omega_r=0
        # so just stamp 1's on the diagonals and leave 0 in the J
        # to avoid matrix singularity
        y_entry_counter = stamp_y_sparse(self.idr_index, self.idr_index, 1, y_row, y_col, y_val, y_entry_counter)
        y_entry_counter = stamp_y_sparse(self.ids_index, self.ids_index, 1, y_row, y_col, y_val, y_entry_counter)
        y_entry_counter = stamp_y_sparse(self.iqr_index, self.iqr_index, 1, y_row, y_col, y_val, y_entry_counter)
        y_entry_counter = stamp_y_sparse(self.iqs_index, self.iqs_index, 1, y_row, y_col, y_val, y_entry_counter)
        y_entry_counter = stamp_y_sparse(self.omega_r_index, self.omega_r_index, 1, y_row, y_col, y_val, y_entry_counter)
        return (y_entry_counter)

    def stamp_t0_dense(self, Y):
        # initial conditions of TPIM
        # voltages still coupled
        #     dvds_dva*va + dvds_dvb*vb + dvds_dvc*vc = vds - v_neutral
        # --> dvds_dva*va + dvds_dvb*vb + dvds_dvc*vc - vds + v_neutral = 0
        Y[self.vds_index, self.phase_a_node_index] += self.p_transform[1,0]
        Y[self.vds_index, self.phase_b_node_index] += self.p_transform[1,1]
        Y[self.vds_index, self.phase_c_node_index] += self.p_transform[1,2]
        Y[self.vds_index, self.vds_index] += -1
        Y[self.vds_index, self.neutral_node_index] += 1
        #     dvqs_dva*va + dvqs_dvb*vb + dvqs_dvc*vc = vqs - v_neutral
        # --> dvqs_dva*va + dvqs_dvb*vb + dvqs_dvc*vc - vqs + v_neutral = 0
        Y[self.vqs_index, self.phase_a_node_index] += self.p_transform[2,0]
        Y[self.vqs_index, self.phase_b_node_index] += self.p_transform[2,1]
        Y[self.vqs_index, self.phase_c_node_index] += self.p_transform[1,2]
        Y[self.vqs_index, self.vqs_index] += -1
        Y[self.vqs_index, self.neutral_node_index] += 1
        
        # all currents=0, omega_r=0
        # so just stamp 1's on the diagonals and leave 0 in the J
        # to avoid matrix singularity
        Y[self.idr_index, self.idr_index] += 1
        Y[self.ids_index, self.ids_index] += 1
        Y[self.iqr_index, self.iqr_index] += 1
        Y[self.iqs_index, self.iqs_index] += 1
        Y[self.omega_r_index, self.omega_r_index] += 1
        

    def calc_V0dq(self, V):
        V_abc = np.array([V[self.phase_a_node_index, 0], V[self.phase_b_node_index, 0], V[self.phase_c_node_index, 0]]).reshape(-1,1)
        V_0dq = np.matmul(self.p_transform, V_abc)
        return (V_0dq[1,0], V_0dq[2,0])

    def get_info_hist(self, V_hist, dt):
        info_hist = {}
        ids_hist = V_hist[self.ids_index, :]
        idr_hist = V_hist[self.idr_index, :]
        iqs_hist = V_hist[self.iqs_index, :]
        iqr_hist = V_hist[self.iqr_index, :]
        info_hist['te_hist'] = -self.Te_alpha*(idr_hist*iqs_hist - iqr_hist*ids_hist)
        info_hist['ids_hist'] = ids_hist
        info_hist['idr_hist'] = idr_hist
        info_hist['iqr_hist'] = iqr_hist
        info_hist['iqs_hist'] = iqs_hist
        info_hist['omega_r_hist'] = V_hist[self.omega_r_index, :]
        info_hist['vds_hist'] = V_hist[self.vds_index, :]
        info_hist['vqs_hist'] = V_hist[self.vqs_index, :]
        info_hist['d_omega_true'] = np.gradient(V_hist[self.omega_r_index, :])/dt
        info_hist['d_omega_exp'] = 1/self.j_im*(info_hist['te_hist'] - self.tm - self.d_fric*info_hist['omega_r_hist'])
        return info_hist