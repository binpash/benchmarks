#!/usr/bin/env python3

import pandas as pd
from do_pca import perform_pca_and_plot
import ast
from project_root import get_project_root
import sys
from pathlib import Path
import viz.syntax as stx
import viz.dynamic as dyn

root = get_project_root()
loc_data_path = root / 'infrastructure/target/lines_of_code.csv'
data_path = root / 'infrastructure/target/dynamic_analysis.jsonl'
OUTPUT = sys.argv[1] if len(sys.argv) > 1 else None

# make OUTPUT absolute path
if OUTPUT:
    OUTPUT = Path(OUTPUT).absolute()
else:
    OUTPUT = Path(__file__).parent / 'pca-row-plot.pdf'

def read_sys_results():
    syscall_data_path = root / 'infrastructure/data/no_of_syscalls.csv'
    df = pd.read_csv(syscall_data_path)
    df.rename(columns={'Benchmark-(small)': 'benchmark', 'System Calls': 'sys_calls'}, inplace=True)
    return df

def count_constructs(series):
    return len(set(series))

def count_unique_cmds(series):
    return len({node for node in series if 'command(' in node})

def read_loc_data():
    loc_data = pd.read_csv(loc_data_path, header=None)
    loc_data.columns = ['script', 'loc']
    map_df = stx.get_map_df()
    loc_data = loc_data.merge(map_df, on='script')
    loc_data_bench = loc_data.groupby('benchmark').agg({
        'loc': 'sum',
        'script': 'count'
    }).reset_index()
    loc_data_bench.rename(columns={'script': 'number_of_scripts'}, inplace=True)
    return loc_data, loc_data_bench

syntax_script, syntax_bench = stx.read_data(True)
dyn_script,    dyn_bench = dyn.read_data()
loc_data_script, loc_data_bench = read_loc_data()
syntax_script_all_cmds, syntax_bench_all_cmds = stx.read_data(False)
syntax_bench['constructs'] = syntax_bench['nodes'].apply(count_constructs)
sys_results = read_sys_results()

syntax_bench_all_cmds['unique_cmds'] = syntax_bench_all_cmds['nodes'].apply(count_unique_cmds)

embedding_df = pd.read_csv(root / 'infrastructure/data/embeddings.csv')
embedding_df['embedding'] = embedding_df['embedding'].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)

# Embedding is a list of numbers, turn them into columns
embedding_df = pd.concat([embedding_df['benchmark'], embedding_df['embedding'].apply(pd.Series)], axis=1)

big_bench = syntax_bench.merge(dyn_bench, on='benchmark')\
    .merge(loc_data_bench, on='benchmark')\
    .merge(syntax_bench_all_cmds[['benchmark', 'unique_cmds']], on='benchmark')\
    .merge(sys_results, on='benchmark')

perform_pca_and_plot(big_bench, embedding_df, str(OUTPUT))
