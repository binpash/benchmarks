import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Data format example:
# covid-mts/scripts/1.sh,command(cat):1;quoted_control:5;command(sed):1;command(cut):2;command(sort):2;command(uniq):1;command(awk):1;pipeline:1
# file-enc/scripts/encrypt_files.sh,command(mkdir):1;variable_use:4;command(openssl):1;quoted_control:2;function_command:1;command(export):1;assignment:1;command(cat):1;command(pure_func):1;file_redirection:1;pipeline:1;for_command:1
# ...

data_path = 'target/nodes_in_scripts.csv'
benchmark_mapping_path = 'target/scripts_to_benchmark.csv'


special_commands = ['eval', 'alias']

node_types = """
pipeline
background
subshell_command
while_command
for_command
case_command
if_command
and_command
or_command
negate_command
function_command
assignment
variable_use
file_redirection
dup_redirection
heredoc_redirection
home_tilde_control
dollar_paren_shell_control
dollar_paren_paren_arith_control
quoted_control
""".strip().split("\n") + special_commands
# Omitted these because they don't seem to be useful:
# ---
# redirection
# raw_command
# escaped_char

def normalize_node_name(node):
    # remove "_command" suffix if present
    return node.replace('_command', '')

def node_heatmap(df):
    # todo which of these are missing entirely?
    #unique_node_names = list(set(df['nodes'].apply(lambda x: [x for x in x.keys()]).sum()))

    heatmap_data = pd.DataFrame(index=list(map(normalize_node_name, node_types)), columns=df['benchmark'])
    for _, row in df.iterrows():
        for node, count in row['nodes'].items():
            heatmap_data.at[normalize_node_name(node), row['benchmark']] = count

    heatmap_data = heatmap_data.fillna(0)
    limit = 5
    heatmap_data = heatmap_data.applymap(lambda x: min(x, limit))
    annot_data = heatmap_data.applymap(lambda x: \
        #'>{}'.format(limit) \
        '*' \
            if x == limit else '')
    
    plt.figure(figsize=(50, 8))
    sns.heatmap(heatmap_data, cmap='Reds', annot=annot_data, fmt='', cbar_kws={'label': 'Node Count'})
    plt.xlabel('Benchmark')
    plt.ylabel('Node Names')
    plt.title('Node Heatmap')
    plt.show()

def extract_special_command(node):
    for sc in special_commands:
        if node == f'command({sc})':
            return sc
    return node

def merge_node_counts(series):
    merged_dict = {}
    for d in series:
        for k, v in d.items():
            merged_dict[k] = merged_dict.get(k, 0) + v
    return merged_dict

def main(data_path):
    df = pd.read_csv(data_path, header=None)
    df.columns = ['script', 'nodes']
    # Unpack node counts
    df['nodes'] = df['nodes'].apply(lambda x: dict([tuple(i.split(':')) for i in x.split(';')]) if isinstance(x, str) else {})
    # Transform nodes entries for 'command(eval)' and 'command(alias)' into 'eval' and 'alias'
    df['nodes'] = df['nodes'].apply(lambda x: {extract_special_command(k): v for k, v in x.items()})
    # Merge all the "command" nodes, we don't care about the individual commands here
    df['nodes'] = df['nodes'].apply(lambda x: {k: int(v) for k, v in x.items() if 'command(' not in k} | {'command': sum([int(v) for k, v in x.items() if 'command(' in k])})
    
    # Aggregate by benchmark
    map_df = pd.read_csv(benchmark_mapping_path, header=None)
    map_df.columns = ['script', 'benchmark']
    df = df.merge(map_df, on='script')
    df = df.groupby('benchmark').agg({'nodes': merge_node_counts}).reset_index()

    node_heatmap(df)

if __name__ == '__main__':
    main(data_path)