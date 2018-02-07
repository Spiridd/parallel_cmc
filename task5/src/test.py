"""
Loads 3 matrices from files, calculates the result
and compares it with what was in the third file
"""
import numpy as np
import argparse
from math import sqrt
import time
import sys


def compare(M, W, eps=1e-3):
    errors = np.sum(np.divide(abs(M-W), M)>eps)
    shape = M.shape
    print("TEST: matrix size is {}x{}, there are {} errors".format(shape[0], shape[1], errors))


def get_np_matrix(filename):
    try:
        f = open(filename)
    except FileNotFoundError:
        print("TEST: Cannot open file " + filename)
        sys.exit(1)
    aux_values = np.fromfile(f, count=3, dtype=np.dtype('i'))
    data_t = 'f' if aux_values[0] == 0 else 'd'
    dim = aux_values[1:]
    print("dim = {}x{}".format(dim[0], dim[1]))
    ret = np.fromfile(f, dtype=np.dtype(data_t)).reshape(dim)
    f.close()
    return ret


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("file_A", type=str)
    parser.add_argument("file_B", type=str)
    parser.add_argument("file_C", type=str)
    args = parser.parse_args()

    A = get_np_matrix(args.file_A)
    B = get_np_matrix(args.file_B)
    C = get_np_matrix(args.file_C)

    m, n = A.shape
    h = B.shape[1]

    begin = time.time()
    C0 = A.dot(B)
    end = time.time()
    """
    print(A)
    print(B)
    print(C)
    print(C0)
    """
    print("TEST: numpy time: {:.5f} sec".format(end-begin))
    
    compare(C, C0)

if __name__ == '__main__':
    main()

