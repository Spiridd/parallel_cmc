from matplotlib import pyplot as plt
import os


def get_files_starting_with(word, dir_name):
    filepath = os.path.abspath(__file__)
    cur_dir = os.path.dirname(filepath)
    dir_name_abs = os.path.join(cur_dir, dir_name)
    return [f for f in os.listdir(dir_name_abs) if f.startswith(word)]


def read_data_from_files(filenames, files_dir):
    times = {}
    for filename in filenames:
        key = filename.lstrip('data_').rstrip('.txt') 
        with open(os.path.join(files_dir, filename), 'r') as f:
            temp = []
            for line in f:
                temp.append(float(line.split(',')[2]))
            times[key] = (min(temp), max(temp), sum(temp))
    return times


def plot_data(data):
    simple_data = []
    for key, value in data.items():
        simple_data.append((key, value[0], value[1], value[2]))
    simple_data = sorted(simple_data, key=lambda x: x[0])
    keys = [x[0] for x in simple_data]
    values_1 = [x[1] for x in simple_data]
    values_2 = [x[2] for x in simple_data]
    values_3 = [x[3] for x in simple_data]
    plt.plot(keys, values_1, 'r', label='min')
    plt.plot(keys, values_2, 'm', label='max')
    plt.plot(keys, values_3, 'y', label='overall')
    plt.legend()
    plt.grid()
    plt.show()


def main():
    files_dir = 'out'
    files = get_files_starting_with('data', files_dir)
    data = read_data_from_files(files, files_dir)
    plot_data(data)


if __name__ == '__main__':
    main()

