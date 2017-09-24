import numpy as np
import argparse
from matplotlib import pyplot as plt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("times", type=str)
    args = parser.parse_args()

    times = np.fromfile(args.times, dtype=np.dtype('d'))
    x = tuple(range(len(times)))
    plt.plot(x, times)
    modes = ("ijk", "jik", "ikj", "kij", "kji", "jki")
    plt.xticks(x, modes)
    plt.ylabel("time, sec")
    plt.title("matrix multiplication")
    plt.grid()
    plt.show()

if __name__ == '__main__':
    main()

