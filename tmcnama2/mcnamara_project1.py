import argparse
import os
import sys
from run_time_domain_circuit import main

def str_to_bool(s):
    return s.lower() == "true"

parser = argparse.ArgumentParser()
parser.add_argument("network_file", type=str, 
                    help="Path to circuit network's .mat file, if applicable")
parser.add_argument("--dt", type=float, 
                    help="Time step for simulation",
                    default=1e-4)
parser.add_argument("--tol", type=float, 
                    help="Newton-Raphson convergence tolerance",
                    default=1e-8)
parser.add_argument("--use_sparse", type=str, 
                    help="True: use sparse matrices; False: use dense matrices",
                    choices=["True", "False"],
                    default='True')
# Dynamic time steps
# Solver

args = parser.parse_args()
args.use_sparse = str_to_bool(args.use_sparse)

main(args.network_file, args.dt, args.tol, args.use_sparse)