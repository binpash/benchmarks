#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from all_scripts import get_all_scripts
from project_root import get_project_root

# Data format example: (column titles not in the file, all times are seconds)
# benchmark_script, user_time, system_time, max_unique_set_size, read_chars, write_chars, user_time_in_shell, system_time_in_shell, all_input_files, wall_time
# covid-mts/scripts/1.sh,17.32,4.08,13512704,10364829524,7538064195,0.0,0.0,covid-mts/input/in.csv,22.41
# file-enc/scripts/encrypt_files.sh,1.55,0.54,2998272,925500199,923880800,0.0,0.01,file-enc/input/pcaps/http-download101c.pcapng;file-enc/input/pcaps/ftp-download101.pcapng;file-enc/input/pcaps/tr-twohosts.pcapng;file-enc/input/pcaps/split250_00004_20160704110759.pcapng;file-enc/input/pcaps/http-pcaprnet101.pcapng;file-enc/input/pcaps/sec-suspicious101.pcapng;file-enc/input/pcaps/challenge101-8.pcapng;file-enc/input/pcaps/challenge101-3.pcapng;file-enc/input/pcaps/challenge101-6.pcapng;file-enc/input/pcaps/http-openoffice101a.pcapng;file-enc/input/pcaps/tr-winsize.pcapng;file-enc/input/pcaps/ftp-crack101.pcapng;file-enc/input/pcaps/http-misctraffic101.pcapng;file-enc/input/pcaps/http-google101.pcapng;file-enc/input/pcaps/http-wiresharkdownload101.pcapng;file-enc/input/pcaps/general101d.pcapng;file-enc/input/pcaps/http-download101.pcapng;file-enc/input/pcaps/http-chappellu101.pcapng;file-enc/input/pcaps/http-college101.pcapng;file-enc/input/pcaps/http-download101d.pcapng;file-enc/input/pcaps/net-lost-route.pcapng;file-enc/input/pcaps/general101.pcapng;file-enc/input/pcaps/http-sfgate101.pcapng;file-enc/input/pcaps;file-enc/input/pcaps/http-browse101b.pcapng;file-enc/input/pcaps/split250_00001_20160704110759.pcapng;file-enc/input/pcaps/http-download-a.pcapng,2.09
# log-analysis/scripts/nginx.sh,0.57,0.17,15114240,417362164,232406881,0.0,0.0,log-analysis/input/nginx-logs;log-analysis/input/nginx-logs/log5;log-analysis/input/nginx-logs/log7,0.75

root = get_project_root()
data_path = root / 'infrastructure/target/dynamic_analysis.csv'
input_size_path = root / 'infrastructure/data/size-inputs.csv'

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
    plt.xticks(rotation=90)
    plt.yscale('symlog', linthresh=0.1)
    plt.legend()
    plt.show()

def plot_benchmark_times(df,
                         ticks = ([0, 0.001, 0.01, 0.1, 1, 10, 100, 1000], 
                                  ['0', '1ms', '10ms', '100ms', '1s', '10s', '100s', '1000s']),
                         ylabel='Time (s)',
                         linthresh=0.001):
    sns.set(style="whitegrid")
    plt.figure(figsize=(10, 6))
    sns.barplot(x='benchmark', y='time_in_commands', data=df, color='blue', label='Commands')
    sns.barplot(x='benchmark', y='time_in_shell', data=df, color='green', label='Shell')
    plt.xticks(rotation=90)
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
    sns.barplot(x='benchmark', y='io_chars', data=df, color='green', label='IO bytes')
    plt.yscale('symlog', linthresh=linthresh)
    plt.yticks(*ticks)
    plt.ylabel(ylabel)
    plt.xticks(rotation=90)
    plt.legend()
    plt.show()

def plot_memory(df,
                ticks=([0, 1000000, 10000000, 100000000, 1000000000, 10000000000, 100000000000], 
                       ['0', '1MB', '10MB', '100MB', '1GB', '10GB', '100GB']),
                ylabel='Memory (bytes)',
                linthresh=1000000):
    sns.set(style="whitegrid")
    plt.figure(figsize=(10, 6))
    sns.barplot(x='benchmark', y='max_unique_set_size', data=df, color='purple', label='Max unique set size')
    plt.xticks(rotation=90)
    plt.yscale('symlog', linthresh=linthresh)
    plt.yticks(*ticks)
    plt.ylabel(ylabel)
    plt.legend()
    plt.show()

def main(data_path):
    df = pd.read_csv(data_path, header=None)
    df.columns = ['script', 'user_time', 'system_time', 'max_unique_set_size', 'read_chars', 'write_chars', 'user_time_in_shell', 'system_time_in_shell', 'all_input_files', 'wall_time', 'start_time']
    for col in list(df.columns[1:7]) + list(df.columns[9:10]):
        df[col] = df[col].astype(float)
    df['all_input_files'] = df['all_input_files'].apply(lambda x: str(x).split(';'))

    # aggregate by benchmark
    map_df = get_map_df()
    df = df.merge(map_df, on='script')
    # sum all times
    df = df.groupby('benchmark').agg({'user_time': 'sum', 
                                      'system_time': 'sum',
                                      'max_unique_set_size': 'sum',
                                      'read_chars': 'sum', 
                                      'write_chars': 'sum', 
                                      'user_time_in_shell': 'sum', 
                                      'system_time_in_shell': 'sum', 
                                      'all_input_files': 'sum',
                                      'wall_time': 'sum'}).reset_index()

    # merge the read and write_chars
    df['io_chars'] = df['read_chars'] + df['write_chars']
    df = df.drop(columns=['read_chars', 'write_chars'])

    # calculate time in shell and time in commands
    df['time'] = df['user_time'] + df['system_time']
    df['time_in_shell'] = df['user_time_in_shell'] + df['system_time_in_shell']
    df['time_in_commands'] = df['time'] - df['time_in_shell']

    # report any benchmarks where the wall time is not approximately equal to the sum of user and system time
    for _, row in df.iterrows():
        if abs(row['wall_time'] - (row['user_time'] + row['system_time'])) > 0.1:
            print(f"Wall time for benchmark {row['benchmark']} maybe suspicious: {row['wall_time']} vs (u{row['user_time']} + s{row['system_time']})")

#     # relative numbers to input size
#     input_sizes = pd.read_csv(input_size_path, header=None)
#     input_sizes.columns = ['input_size', # bytes
#                            'input_file']
#     input_sizes['input_size'] = input_sizes['input_size'].apply(lambda x: int(x))
#     input_sizes['benchmark'] = input_sizes['input_file'].apply(lambda x: str(x).split('/')[0])
#     input_sizes = input_sizes.groupby('benchmark').agg({'input_size': 'sum'}).reset_index()
# 
#     df_rel_to_input = df.merge(input_sizes, on='benchmark')
#     df_rel_to_input['io_chars'] = df_rel_to_input['io_chars'] / df_rel_to_input['input_size']
#     df_rel_to_input['max_unique_set_size'] = df_rel_to_input['max_unique_set_size'] / df_rel_to_input['input_size']
#     df_rel_to_input['time_in_shell'] = df_rel_to_input['time_in_shell'] / df_rel_to_input['input_size']
#     df_rel_to_input['time_in_commands'] = df_rel_to_input['time_in_commands'] / df_rel_to_input['input_size']

    plot_benchmark_times(df)
    plot_io(df)
    plot_memory(df)

#     plot_benchmark_times(df_rel_to_input, ylabel='Time per input byte',
#                          ticks=([0, 0.00000001, 0.0000001, 0.000001, 0.00001, 0.0001], 
#                                 ['0', '10ns', '100ns',    '1us',    '10us', '100us']),
#                                 linthresh=0.00000001)
#     plot_io(df_rel_to_input, ylabel='IO per input byte',
#             ticks=([0,    1,   10,     100,    1000], 
#                    ['0', '1B', '10B', '100B', '1KB']),
#                    linthresh=1)
#     plot_memory(df_rel_to_input, ylabel='Memory per input byte',
#                 ticks=([0,   0.001,   0.01,     0.1,  1,    10,    100,   1000,  10000], 
#                        ['0', '0.001B', '0.01B', '0.1B', '1B', '10B', '100B', '1KB', '10KB']),
#                 linthresh=0.001)


if __name__ == '__main__':
    main(data_path)
