#fast_polarity boundscheck=False, wraparound=False, nonecheck=False

import numpy as np

def polarity_edge_finder_optimised(double[:] input_array):
    cdef long[:] output_array = np.zeros(len(input_array), dtype=np.int64)

    cdef int i
    i = 0

    while i < (len(input_array)-1):
        if input_array[i] * input_array[i+1] <= 0:
            output_array[i] = 1
        else:
            output_array[i] = 0
        i += 1

    return output_array
