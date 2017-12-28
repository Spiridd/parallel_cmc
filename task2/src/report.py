"""
This script reads data from a file, combines it in a new file and plots it:
* L1 cache misses on size of matrix for ijk, ikj indices and bsz=32, 36
* same for L2
* total cycles on ...
* TLB misses on ...
* time on ...
* flops on ...
"""

from matplotlib import pyplot as plt
import csv


def combine_data(filename1, filename2):
    new_filename = '../res/data.csv'
    fout = open(new_filename, 'w')

    part1 = open(filename1).readlines()
    part2 = open(filename2).readlines()

    for line1, line2 in zip(part1, part2):
        fout.write(line1.rstrip('\n'))
        # skip indices, size, bsz and time
        fields = line2.rstrip('\n').split(',')[3:-1]
        [fout.write(',' + field) for field in fields]
        fout.write('\n')
    fout.close()
    return new_filename


def add_line(line, fieldnames, dest):
    for value, name in zip(line, fieldnames):
        if name == 'bsz' or name == 'indices':
            continue
        v = int(value) if name != 'time' else float(value)
        try:
            dest[name].append(v)
        except KeyError:
            dest[name] = [v]


def convert_PAPI(string):
    buf = ['']
    for c in string:
        if c == 'D':
            buf.append('data')
        elif c == 'I':
            buf.append('instruction')
        elif c == 'C':
            buf.append('cache')
        elif c == 'M':
            buf.append('misses')
        else:
            print('error in interpreting: ' + string)
    return ' '.join(buf)


def get_title(field):
    if field.startswith('PAPI'):
        name = field[5:]
        if name == 'TOT_CYC':
            return 'total cpu cycles'
        elif name == 'FP_OPS':
            return 'total floating point operations'
        else:
            index = name.index('_')
            return name[:index] + convert_PAPI(name[index+1:])
    else:
        return field


def get_label(field):
    return 'time, sec' if field == 'time' else field


def plot_data(csv_filename):
    f = open(csv_filename, 'r')
    csv_reader = csv.reader(f)
    fieldnames = next(csv_reader)

    ijk_32 = {}
    ikj_32 = {}
    ijk_36 = {}
    for line in csv_reader:
        if line[2] == '32':
            if line[0] == 'ijk':
                add_line(line, fieldnames, ijk_32)
            else:
                add_line(line, fieldnames, ikj_32)
        else:
            if line[0] == 'ijk':
                add_line(line, fieldnames, ijk_36)

    for field in fieldnames:
        if field in ('indices', 'bsz', 'size'):
            continue
        plt.figure()
        plt.xlabel('size')
        plt.title(get_title(field))
        plt.ylabel(get_label(field))
        plt.plot(ijk_32['size'], ijk_32[field], 'b', label='ijk 32')
        plt.plot(ijk_36['size'], ijk_36[field], 'r', label='ijk 36')
        plt.plot(ikj_32['size'], ikj_32[field], 'k', label='ikj 32')
        plt.legend()
        plt.grid()
    plt.show()


def main():
    filename = combine_data('../res/results1.txt', '../res/results2.txt')
    plot_data(filename)


if __name__ == '__main__':
    main()

