import numpy as np

def initialize(devices, size_Y):
    # instantiate the state-variable vector
    # this function is boring for now, but we could make it fancier
    # if we wanted to
    V = np.zeros((size_Y,1), dtype=np.float)
    return V