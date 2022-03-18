"""Read testcase files.

Author(s): Naeem Turner-Bandele
Created Date: 01-13-2021
Updated Date: 01-15-2021
Email: nturnerb@cmu.edu
Status: Development

Reads an input mat file and then passes the outputs to a dictionary for use.

Usage:
	r = Reader(input_file='input.mat')   # instantiate a reader object by providing a .mat case file
	r.read()  # read the file and assign all NumPy array object to a dictionary
	casedata = r.all_objects

"""

from scipy.io import loadmat
import numpy as np


class Reader:
    all_objects = {}  # stores all read objects from the .mat file

    def __init__(self, **kwargs):
        self.input_file = kwargs.get("input_file",
                                     "../testcases/RL_circuit.mat")
        self._nodes = None
        self._resistors = None
        self._inductors = None
        self._capacitors = None
        self._sources = None
        self._inductionMotors = None
        self._switches = None

    def read_nodes(self, node_data):
        """Read and clean up the input node data.

        Args:
          node_data (ndarray): An array of arrays with node circuit data.

        Returns:
          None

        """
        dt = np.dtype([('Name', np.unicode_, 16), ('Vnom', np.float64, (1,)),
                       ('Phases', np.unicode_, 4)])
        nodes = [
            (ele[0].item(), ele[1].item(), ele[2].item()) for ele in node_data
        ]
        self._nodes = np.array(nodes, dtype=dt)

    def read_sources(self, source_data):
        """Read and clean up the input voltage source data.

        Args:
          source_data (ndarray): An array of arrays with voltage source circuit data.

        Returns:
          None

        """
        dt = np.dtype([('Name', np.unicode_, 16), ('Phases', np.unicode_, 4),
                       ('FromNode', np.unicode_, 16), ('ToNode', np.unicode_,
                                                       16),
                       ('VoltageAmplitude', np.float64, (1,)),
                       ('Angle', np.float64, (1,)),
                       ('Frequency', np.float64, (1,))])
        sources = [(ele[0].item(), ele[1].item(), ele[2].item(), ele[3].item(),
                    ele[4].item(), ele[5].item(), ele[6].item())
                   for ele in source_data]
        self._sources = np.array(sources, dtype=dt)

    def read_resistors(self, resistor_data):
        """Read and clean up the input resistor data.

        Args:
          resistor_data (ndarray): An array of arrays with the resistor circuit data.

        Returns:
          None

        """
        dt = np.dtype([('Name', np.unicode_, 16), ('Phases', np.unicode_, 4),
                       ('FromNode', np.unicode_, 16), ('ToNode', np.unicode_,
                                                       16),
                       ('R', np.float64, (1,))])
        resistors = [(ele[0].item(), ele[1].item(), ele[2].item(),
                      ele[3].item(), ele[4].item()) for ele in resistor_data]
        self._resistors = np.array(resistors, dtype=dt)

    def read_inductors(self, inductor_data):
        """Read and clean up the input inductor data.

        Args:
          inductor_data (ndarray): An array of arrays with the inductor circuit data.

        Returns:
          None

        """
        dt = np.dtype([('Name', np.unicode_, 16), ('Phases', np.unicode_, 4),
                       ('FromNode', np.unicode_, 16), ('ToNode', np.unicode_,
                                                       16),
                       ('L', np.float64, (1,))])
        inductors = [(ele[0].item(), ele[1].item(), ele[2].item(),
                      ele[3].item(), ele[4].item()) for ele in inductor_data]
        self._inductors = np.array(inductors, dtype=dt)

    def read_capacitors(self, capacitor_data):
        """Read and clean up the input capacitor data.

        Args:
          capacitor_data (ndarray): An array of arrays with the capacitor circuit data.

        Returns:
          None

        """
        dt = np.dtype([('Name', np.unicode_, 16), ('Phases', np.unicode_, 4),
                       ('FromNode', np.unicode_, 16), ('ToNode', np.unicode_,
                                                       16),
                       ('C', np.float64, (1,))])
        capacitors = [(ele[0].item(), ele[1].item(), ele[2].item(),
                       ele[3].item(), ele[4].item()) for ele in capacitor_data]
        self._capacitors = np.array(capacitors, dtype=dt)

    def read_inductionMotors(self, inductionMotor_data):
        """Read and clean up the input induction motor data.

        Args:
          inductionMotor_data (ndarray): An array of arrays with the induction motor circuit data.

        Returns:
          None

        """
        dt = np.dtype([('Name', np.unicode_, 16),
                       ('connectedNodeA', np.unicode_, 16),
                       ('connectedNodeB', np.unicode_, 16),
                       ('connectedNodeC', np.unicode_, 16),
                       ('Pn', np.float64, (1,)), ('Vn', np.float64, (1,)),
                       ('Fn', np.float64, (1,)), ('Lm', np.float64, (1,)),
                       ('Rs', np.float64, (1,)), ('Rr', np.float64, (1,)),
                       ('Lls', np.float64, (1,)), ('Llr', np.float64, (1,)),
                       ('J', np.float64, (1,)), ('Tm', np.float64, (1,)),
                       ('D', np.float64, (1,))])
        inductionMotors = [
            (ele[0].item(), ele[1].item(), ele[2].item(), ele[3].item(),
             ele[4].item(), ele[5].item(), ele[6].item(), ele[7].item(),
             ele[8].item(), ele[9].item(), ele[10].item(), ele[11].item(),
             ele[12].item(), ele[13].item(), ele[14].item())
            for ele in inductionMotor_data
        ]
        self._inductionMotors = np.array(inductionMotors, dtype=dt)

    def read_switches(self, switch_data):
        """Read and clean up the input switch data.

        Args:
          switch_data (ndarray): An array of arrays with the switch circuit data.

        Returns:
          None

        """
        dt = np.dtype([('Name', np.unicode_, 16), ('Phases', np.unicode_, 4),
                       ('FromNode', np.unicode_, 16), ('ToNode', np.unicode_,
                                                       16),
                       ('TimeOpen', np.float64, (1,)),
                       ('TimeClose', np.float64, (1,))])
        switches = [(ele[0].item(), ele[1].item(), ele[2].item(), ele[3].item(),
                     ele[4].item(), ele[5].item()) for ele in switch_data]
        self._switches = np.array(switches, dtype=dt)

    def read(self):
        """
        Load an input .mat file with the circuit data, clean up the data, and then add the circuit array objects to
        a dictionary.
        Returns:
          None
        """
        input_data = loadmat(self.input_file)

        self.read_nodes(input_data['Nodes'])
        if 'Sources' in input_data:
            self.read_sources(input_data['Sources'])
        if 'Resistors' in input_data:
            self.read_resistors(input_data['Resistors'])
        if 'Inductors' in input_data:
            self.read_inductors(input_data['Inductors'])
        if 'Capacitors' in input_data:
            self.read_capacitors(input_data['Capacitors'])
        if 'InductionMotors' in input_data:
            self.read_inductionMotors(input_data['InductionMotors'])
        if 'Switches' in input_data:
            self.read_switches(input_data['Switches'])

        self.all_objects['nodes'] = self._nodes
        self.all_objects['sources'] = self._sources
        self.all_objects['resistors'] = self._resistors
        self.all_objects['inductors'] = self._inductors
        self.all_objects['capacitors'] = self._capacitors
        self.all_objects['inductionMotors'] = self._inductionMotors
        self.all_objects['switches'] = self._switches
