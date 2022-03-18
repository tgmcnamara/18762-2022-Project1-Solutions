import sys
sys.path.append("..")
from classes.Reader import Reader
from classes.Nodes import Nodes
from classes.Resistors import Resistors
from classes.Capacitors import Capacitors
from classes.Inductors import Inductors
from classes.Switches import Switches
from classes.AcVoltageSources import AcVoltageSources
from classes.InductionMotors import InductionMotors


def process_case(input_file, use_sparse):
    # instantiate a reader object by providing a .mat case file
    r = Reader(input_file=input_file)
    # read the file and assign all NumPy array object to a dictionary
    r.read()

    nodes = []
    resistors = []
    capacitors = []
    inductors = []
    induction_motors = []
    switches = []
    ac_voltage_sources = []
    devices = {}

    # create device objects
    for ele in r.all_objects['nodes']:
        nodes.append(Nodes(ele[0], ele[1], ele[2]))
    
    ground_node = [ele for ele in nodes if ele.ground]
    if len(ground_node) > 1:
        raise Exception("More than one ground node listed in input file.")
    elif len(ground_node) == 0:
        raise Exception("Input file contains no ground node to use as reference")
    else:
        ground_node = ground_node[0]
        Nodes.ground_node_ = ground_node.name
        devices['ground_node'] = ground_node

    if type(r.all_objects['resistors']) != type(None):
        for ele in r.all_objects['resistors']:
            resistors.append(Resistors(ele[0], ele[1], ele[2], ele[3], ele[4])) 
    
    if type(r.all_objects['capacitors']) != type(None):
        for ele in r.all_objects['capacitors']:
            capacitors.append(Capacitors(ele[0], ele[1], ele[2], ele[3], ele[4]))

    if type(r.all_objects['inductors']) != type(None):
        for ele in r.all_objects['inductors']:
            inductors.append(Inductors(ele[0], ele[1], ele[2], ele[3], ele[4]))
    
    if type(r.all_objects['switches']) != type(None):
        for ele in r.all_objects['switches']:
            switches.append(Switches(ele[0], ele[1], ele[2], ele[3], ele[4], ele[5]))

    if type(r.all_objects['inductionMotors']) != type(None):
        for ele in r.all_objects['inductionMotors']:
            induction_motors.append(InductionMotors(ele[0], ele[1], ele[2], ele[3], ground_node.name, ele[4], ele[5], ele[6],
                                                    ele[7], ele[8], ele[9], ele[10], ele[11], ele[12], ele[13], ele[14]))

    if type(r.all_objects['sources']) != type(None):
        for ele in r.all_objects['sources']:
            ac_voltage_sources.append(AcVoltageSources(ele[0], ele[1], ele[2], ele[3], ele[4], ele[5], ele[6]))
    # assign indexes
    for ele in nodes:
        ele.assign_indexes(use_sparse)
    for ele in resistors:
        ele.assign_indexes(nodes)
    for ele in capacitors:
        ele.assign_indexes(nodes)
    for ele in inductors:
        ele.assign_indexes(nodes)
    for ele in switches:
        ele.assign_indexes(nodes)
    for ele in induction_motors:
        ele.assign_indexes(nodes)
    for ele in ac_voltage_sources:
        ele.assign_indexes(nodes)

    devices['nodes'] = nodes
    devices['resistors'] = resistors
    devices['capacitors'] = capacitors
    devices['inductors'] = inductors
    devices['switches'] = switches
    devices['induction_motors'] = induction_motors
    devices['ac_voltage_sources'] = ac_voltage_sources

    return devices