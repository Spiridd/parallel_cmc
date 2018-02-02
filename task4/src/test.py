"""
Loads 3 matrices from files, calculates the result
and compares it with what was in the third file
"""
import numpy as np
import argparse
from math import sqrt
import time


def compare(M, W, eps=1e-3):
    errors = np.sum(np.divide(abs(M-W), M)>eps)
    shape = M.shape
    print("matrix size is {}x{}, there are {} errors".format(shape[0], shape[1], errors))


def get_np_matrix(filename):
    f = open(filename)
    aux_values = np.fromfile(f, count=3, dtype=np.dtype('i'))
    data_t = 'f' if aux_values[0] == 0 else 'd'
    dim = aux_values[1:]
    ret = np.fromfile(f, dtype=np.dtype(data_t)).reshape(dim)
    f.close()
    return ret


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("file_A", type=str)
    parser.add_argument("file_b", type=str)
    parser.add_argument("file_c", type=str)
    args = parser.parse_args()

    A = get_np_matrix(args.file_A)
    b = get_np_matrix(args.file_b)
    c = get_np_matrix(args.file_c)

    m, n = A.shape
    h = b.shape[1]

    begin = time.time()
    c0 = A.dot(b)
    end = time.time()
    print("numpy time: {:.5f} sec".format(end-begin))
    
    compare(c, c0)

if __name__ == '__main__':
    main()

