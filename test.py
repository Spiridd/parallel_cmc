import numpy as np
import argparse
from math import sqrt


def compare(M, W, eps=1e-6):
    errors = np.sum(abs(M-W)>eps)
    print("there are {} positions and {} errors".format(len(M)**2, errors))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("file_A", type=str)
    parser.add_argument("file_B", type=str)
    parser.add_argument("file_C", type=str)
    args = parser.parse_args()

    A = np.fromfile(args.file_A, dtype=np.dtype('d'))
    B = np.fromfile(args.file_B, dtype=np.dtype('d'))
    C = np.fromfile(args.file_C, dtype=np.dtype('d'))
    size = int(sqrt(len(A)))
    shape = (size, size)
    A = A.reshape(shape)
    B = B.reshape(shape)
    C = C.reshape(shape)
    C0 = A.dot(B)

    compare(C, C0)

if __name__ == '__main__':
    main()

