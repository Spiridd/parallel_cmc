import numpy as np
import argparse
from math import sqrt


def compare(M, W, eps=1e-3):
    errors = np.sum(np.divide(abs(M-W), W)>eps)
    shape = M.shape
    print("matrix size is {}x{}, theres is {} errors".format(shape[0], shape[1], errors))


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
    parser.add_argument("file_B", type=str)
    parser.add_argument("file_C", type=str)
    args = parser.parse_args()

    A = get_np_matrix(args.file_A)
    B = get_np_matrix(args.file_B)
    C = get_np_matrix(args.file_C)

    n, m = A.shape
    h = B.shape[1]
    C0 = A.dot(B)

    compare(C, C0)

if __name__ == '__main__':
    main()

