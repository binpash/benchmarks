#!/usr/bin/env python3

from sys import stderr
from pathlib import Path
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import sys, os
# Add the parent directory to sys.path
sys.path.insert(1, os.path.join(sys.path[0], '..'))
from all_scripts import get_all_scripts
from project_root import get_project_root

root = get_project_root()
data_path = root / 'infrastructure/target/dynamic_analysis.jsonl'
input_size_path = root / 'infrastructure/data/size_inputs.jsonl'

def get_input_sizes_df(df):
    sizes_df = pd.read_json(input_size_path, lines=True)
    def find_input_size(row):
        total = 0
        for file in row['all_input_files']:
            file_path = file
            file = Path(file).relative_to(row['category'])
            file_rest = str(Path(*file.parts[1:]))
            relevant = sizes_df[(sizes_df['path'] == file_rest) & (sizes_df['category'] == row['category'])]
            if relevant.empty:
                if 'sort' in file_path: # one of them sort files
                    continue
                if 'riker/input' in file_path: 
                    # these should all be in the file, 
                    # everything else is in an intermediate
                    continue 
                if 'aurpkg/outputs' in file_path:
                    # these are intemediate files
                    continue
                if 'bio/outputs' in file_path:
                    # these are intemediate files
                    continue
                if 'test_result' in file_path:
                    # these are final files
                    continue
                if 'tmp' in file_path:
                    # these are intermediate files
                    continue
                print('could not find input size for', file, row['script'], file=stderr)
                continue
            size, = relevant['size_bytes']
            total += size
        return total
    df['input_size'] = df.apply(find_input_size, axis=1)
    return df

def get_map_df():
    items = [
        (str(script.relative_to(root)), benchmark_name)
        for benchmark_name, scripts in get_all_scripts().items()
        for script in scripts
    ]
    return pd.DataFrame(items, columns=['script', 'benchmark'])

def plot_benchmark_times_split(df):
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(10, 6))
    sns.barplot(x='benchmark', y='user_time', data=df, color='blue', label='User time')
    sns.barplot(x='benchmark', y='system_time', data=df, color='red', label='System time')
    plt.xticks(rotation=60, ha='right')
    plt.subplots_adjust(bottom=0.25)
    plt.yscale('symlog', linthresh=0.1)
    plt.legend()
    plt.show()

def plot_benchmark_times(df,
                         ticks = ([0, 0.1, 1, 10, 100, 1000, 10000], 
                                  ['0', '0.1s', '1s', '10s', '100s', '1,000s', '10,000s']),
                         ylabel='CPU time',
                         linthresh=0.1):
    sns.set(style="whitegrid")
    plt.figure(figsize=(10, 6))
    sns.barplot(x='benchmark', y='time_in_commands', data=df, color='blue', label='Commands')
    sns.barplot(x='benchmark', y='time_in_shell', data=df, color='green', label='Shell')
    plt.xticks(rotation=60, ha='right')
    plt.xlabel('')
    plt.subplots_adjust(bottom=0.25)
    plt.yscale('symlog', linthresh=linthresh)
    plt.yticks(*ticks)
    plt.ylabel(ylabel)
    plt.legend()
    plt.show()

def plot_io(df,
            ticks=([0, 100000000, 1000000000, 10000000000, 100000000000, 1000000000000], 
                   ['0', '100MB', '1GB', '10GB', '100GB', '1TB']),
            ylabel='IO bytes',
            linthresh=100000000):
    sns.set(style="whitegrid")
    plt.figure(figsize=(10, 6))
    sns.barplot(x='benchmark', y='io_chars', data=df, color='green', label=None)
    plt.yscale('symlog', linthresh=linthresh)
    plt.yticks(*ticks)
    plt.ylabel(ylabel)
    plt.xticks(rotation=60, ha='right')
    plt.xlabel('')
    plt.subplots_adjust(bottom=0.25)
    plt.show()

def plot_memory(df,
                ticks=([0, 1000000, 10000000, 100000000, 1000000000], 
                       ['0', '1MB', '10MB', '100MB', '1GB']),
                ylabel='Memory high water mark (bytes)',
                linthresh=1000000):
    sns.set(style="whitegrid")
    plt.figure(figsize=(10, 6))
    sns.barplot(x='benchmark', y='max_unique_set_size', data=df, color='purple', label=None)
    plt.xticks(rotation=60, ha='right')
    plt.xlabel('')
    plt.subplots_adjust(bottom=0.25)
    plt.yscale('symlog', linthresh=linthresh)
    plt.yticks(*ticks)
    plt.ylabel(ylabel)
    plt.show()

def plot_time_vs_wall(df):
    # what fraction of the real (wall) runtime of the process is user or system time?
    df['time_occupied'] = (df['user_time'] + df['system_time']) / df['wall_time']
    sns.set(style="whitegrid")
    plt.figure(figsize=(10, 6))
    sns.barplot(x='benchmark', y='time_occupied', data=df, color='blue')
    plt.subplots_adjust(bottom=0.25)
    plt.xticks(rotation=60, ha='right')
    plt.xlabel('')
    plt.yticks(np.linspace(0, 6, 10))
    plt.ylabel('Proporion of CPU time to wall time')
    plt.show()

dynamic_analysis_script_translations = {
    "riker/scripts/vim/run.sh": "riker/scripts/vim/build.sh",
    "riker/scripts/xz/run.sh": "riker/scripts/xz/build.sh",
    "riker/scripts/redis/run.sh": "riker/scripts/redis/build.sh",
    "riker/scripts/xz-clang/run.sh": "riker/scripts/xz-clang/build.sh",
    "riker/scripts/lua/run.sh": "riker/scripts/lua/build.sh",
    "riker/scripts/memcached/run.sh": "riker/scripts/memcached/build.sh",
    "riker/scripts/sqlite/run.sh": "riker/scripts/sqlite/build.sh"
}

def read_data():
    df = pd.read_json(data_path, lines=True)
    for col in list(df.columns[1:7]) + list(df.columns[9:10]):
        df[col] = df[col].astype(float)

    df = get_input_sizes_df(df)

    # translate any script names that need it
    # unfortunately this is a bit of a hack to cover up some inconsistency between what the dynamic analysis sees as the "main script" run by a benchmark, 
    # and what really is the main script that the syntax analysis considers
    # TODO: move this mapping above into script_globs.json to make this discrepancy a little more explicit
    df['script'] = df['script'].apply(lambda x: dynamic_analysis_script_translations.get(x, x))

    # aggregate by benchmark
    map_df = get_map_df()
    df = df.merge(map_df, on='script')

    # report any benchmarks where the wall time is not approximately equal to the sum of user and system time
    for _, row in df.iterrows():
        if abs(row['wall_time'] - (row['user_time'] + row['system_time'])) > 0.1:
            print(f"Wall time for benchmark {row['benchmark']} maybe suspicious: {row['wall_time']} vs (u{row['user_time']} + s{row['system_time']})", file=stderr)

    # merge the read and write_chars
    df['io_chars'] = df['read_chars'] + df['write_chars']
    # calculate time in shell and time in commands
    df['time'] = df['user_time'] + df['system_time']
    df['time_in_shell'] = df['user_time_in_shell'] + df['system_time_in_shell']
    df['time_in_commands'] = df['time'] - df['time_in_shell']

    # sum all times
    bench_df = df.groupby('benchmark').agg({'user_time': 'sum', 
                                      'system_time': 'sum',
                                      'max_unique_set_size': 'max',
                                      'read_chars': 'sum', 
                                      'write_chars': 'sum', 
                                      'user_time_in_shell': 'sum', 
                                      'system_time_in_shell': 'sum', 
                                      'input_size': 'sum',
                                      'io_chars': 'sum',
                                      'time': 'sum',
                                      'time_in_shell': 'sum',
                                      'time_in_commands': 'sum',
                                      'wall_time': 'sum'}).reset_index()

    return df, bench_df

def main():
    _, df = read_data()

    
    
#     # this metric is bad in the cases that we don't have the data for the inputs
#     # the benchmark has a specific input on disk that it processed. how much of this input did it process per second?
#     df['bytes_per_second_input'] = df['input_size'] / (df['wall_time'])
#     sns.set(style="whitegrid")
#     plt.figure(figsize=(10, 6))
#     sns.barplot(x='benchmark', y='bytes_per_second_input', data=df, color='blue', label='Commands')
#     plt.xticks(rotation=60, ha='right')
#     plt.title('missing input')
#     plt.yscale('symlog', linthresh=1e1)
#     plt.yticks(np.logspace(1, 12, 12))
#     plt.ylabel('bytes (input file size) per second (wall time)')
#     plt.legend()
#     plt.show()

    # # fraction of the actual scheduled (user + system) time is in the shell. how much cpu work is in the launnched shell
    # df['time_in_shell_frac'] = (df['user_time_in_shell'] + df['system_time_in_shell']) / (df['user_time'] + df['system_time'])
    # sns.set(style="whitegrid")
    # plt.figure(figsize=(10, 6))
    # sns.barplot(x='benchmark', y='time_in_shell_frac', data=df, color='blue', label='Commands')
    # plt.xticks(rotation=60, ha='right')
    # plt.subplots_adjust(bottom=0.25)
    # plt.yticks(np.linspace(0, 1, 10))
    # plt.ylabel('fraction of time (utime + stime) in the shell')
    # plt.legend()
    # plt.show()

    
    plot_time_vs_wall(df)
    plot_benchmark_times(df)
    plot_io(df)
    plot_memory(df)

    df_rel_to_input = df.copy()
    df_rel_to_input['io_chars'] = df_rel_to_input['io_chars'] / df_rel_to_input['input_size']
    df_rel_to_input['max_unique_set_size'] = df_rel_to_input['max_unique_set_size'] / df_rel_to_input['input_size']
    df_rel_to_input['time_in_shell'] = df_rel_to_input['time_in_shell'] / df_rel_to_input['input_size']
    df_rel_to_input['time_in_commands'] = df_rel_to_input['time_in_commands'] / df_rel_to_input['input_size']

    plot_benchmark_times(df_rel_to_input, ylabel='CPU time per input byte',
                         ticks=([0, 0.00000001, 0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01], 
                                ['0', '10ns', '100ns',    '1us',    '10us', '100us',  '1ms', '10ms']),
                                linthresh=0.00000001)
    plot_io(df_rel_to_input, ylabel='IO per input byte',
            ticks=([0,    1,   10,     100,    1000,  10000,  100000, 1000000], 
                   ['0', '1B', '10B', '100B', '1KB', '10KB', '100KB', '1MB']),
                   linthresh=1)
    plot_memory(df_rel_to_input, ylabel='Memory per input byte',
                ticks=([0,   0.001,   0.01,     0.1,  1,    10,    100,   1000,  10000], 
                       ['0', '0.001B', '0.01B', '0.1B', '1B', '10B', '100B', '1KB', '10KB']),
                linthresh=0.001)
    
    # the benchmark did a certain amount of io. how many bytes per second was this?
    df_rel_to_wall = df.copy()
    df_rel_to_wall['io_chars'] = (df_rel_to_wall['read_chars'] + df_rel_to_wall['write_chars']) / (df_rel_to_wall['wall_time'])
    plot_io(df_rel_to_wall, ylabel='IO bytes per second wall time',
            ticks=([0, 1000000, 10000000, 100000000, 1000000000, 10000000000], 
                   ['0', '1MB', '10MB', '100MB', '1GB', '10GB']),
            linthresh=1000000)

if __name__ == '__main__':
    main()
