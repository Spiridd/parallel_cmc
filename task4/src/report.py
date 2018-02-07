from matplotlib import pyplot as plt
import os


def get_filenames(dir_name):
    filepath = os.path.abspath(__file__)
    cur_dir = os.path.dirname(filepath)
    dir_name_abs = os.path.join(cur_dir, dir_name)
    return [f for f in os.listdir(dir_name_abs) if not f.startswith('desc')]


def read_data_from_files(filenames, files_dir):
    data = []
    for filename in filenames:
        with open(os.path.join(files_dir, filename), 'r') as f:
            record = {}
            lines = f.readlines()
            record['n_proc'] = int(lines[0].split()[1])
            record['size'] = lines[1].split()[1]
            record['time'] = float(lines[3].split()[-1])
            point_index = filename.index('.')
            record['mapping'] = filename[:point_index]
            record['mode'] = filename[point_index+1:filename.index('.', point_index+1)]
            data.append(record)
    return data


def plot_data(data):
    # add effectiveness, speedup
    modes = set([d['mode'] for d in data])
    for mode in modes:
        plt.figure()
        sizes = set([d['size'] for d in data])
        for size in sizes:
            some_data = [(d['time'], d['n_proc']) for d in data if d['mode']==mode and d['size']==size]
            some_data = sorted(some_data, key=lambda x: x[1])
            one_proc_time = some_data[0][0]
            keys = [x[1] for x in some_data]
            values = [x[0] for x in some_data]
            plt.semilogx(keys, values, ':o', label=size, basex=2)
        plt.title(mode + ' time')
        plt.xlabel('number of processes')
        plt.ylabel('time')
        plt.legend()
        plt.grid()
    plt.show()


def main():
    files_dir = '../out'
    filenames = get_filenames(files_dir)
    data = read_data_from_files(filenames, files_dir)
    plot_data(data)


if __name__ == '__main__':
    main()

