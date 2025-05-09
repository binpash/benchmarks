import os
import matplotlib.pyplot as plt
import numpy as np


# Set the plotting style if desired
# plt.style.use('ggplot')  # Example: ggplot style

# Creates a bar chart on the given axis.
def create_bar_chart(ax, labels, data, bar_width, position, color, label):
    return ax.bar([p + position for p in range(len(data))], data, bar_width, label=label, color=color)

# Configures the axes properties.
def setup_ax(ax, xlabel, ylabel, title, xticks, xticklabels):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels, rotation=45, ha='right')

# Saves the current plot to the specified directory.
def save_plot(output_dir, filename):
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{filename}.pdf"))

# Plots a comparison of execution times for Bash and hs.
def plot_benchmark_times_combined(benchmarks, bash_times, orch_times, output_dir, filename):
    fig, ax = plt.subplots(figsize=(10, 6))
    bar_width = 0.35

    create_bar_chart(ax, benchmarks, bash_times, bar_width, 0, 'b', 'Bash')
    create_bar_chart(ax, benchmarks, orch_times, bar_width, bar_width, 'r', 'hs')

    setup_ax(ax, 'Benchmarks', 'Execution Time (s)', 'Execution Time Comparison: Bash vs hs',
             [i + bar_width / 2 for i in range(len(benchmarks))], benchmarks)
    ax.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{filename}.pdf"))


# Plots individual comparison charts for each benchmark.
def plot_benchmark_times_individual(benchmarks, bash_times, orch_times, output_dir, filename):
    num_benchmarks = len(benchmarks)
    fig, axes = plt.subplots(num_benchmarks, 1, figsize=(10, 6 * num_benchmarks), squeeze=False)
    
    for ax, benchmark, bash_time, orch_time in zip(axes.flatten(), benchmarks, bash_times, orch_times):
        labels = ['Bash', 'hs']
        times = [bash_time, orch_time]
        create_bar_chart(ax, labels, times, 0.2, 0, 'b', 'Bash')
        create_bar_chart(ax, labels, times, 0.2, 0.2, 'r', 'hs')

        setup_ax(ax, '', 'Execution Time (s)', f'Execution Time Comparison for {benchmark}', range(len(labels)), labels)

    save_plot(output_dir, filename)

def sort_node_ids(node_ids):
    def parse_id(node_id):
        parts = node_id.split('+')
        concrete_id = int(parts[0])
        iter_ids = tuple(int(iter_id) for iter_id in parts[1].split('-')) if len(parts) > 1 else ()
        return (concrete_id,) + iter_ids

    sorted_ids = sorted(node_ids, key=parse_id, reverse=True)
    return sorted_ids

def plot_prog_blocks(prog_blocks, output_dir, filename):
    # Define colors for different statuses
    colors = {
        'READY': 'red',
        'EXE': 'orange',
        'SPEC_E': 'blue',
        'SPEC_F': 'lightblue',
        'COMMIT': 'green',
        'UNSAFE': 'purple',
        'INIT': 'grey'
    }

    first_time = prog_blocks[0][0]
    times = [(block[0] - first_time).total_seconds() for block in prog_blocks]

    unsorted_node_ids = {node[0] for block in prog_blocks for node in block[1]}
    node_ids = sort_node_ids(unsorted_node_ids)  # Sort the node IDs using the custom sorting function

    statuses = {node_id: [] for node_id in node_ids}

    for block_time, nodes in prog_blocks:
        elapsed_time = (block_time - first_time).total_seconds()
        for node_id, status in nodes:
            statuses[node_id].append((elapsed_time, status))

    fig_height = 0.5 * len(node_ids)
    fig_width = 12  # Fixed width
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    status_legend_handles = {}

    for node_id in node_ids:
        y_pos = node_ids.index(node_id)
        for i, (start_time, status) in enumerate(statuses[node_id]):
            end_time = times[-1] if i == len(statuses[node_id]) - 1 else statuses[node_id][i + 1][0]
            color = colors.get(status, 'grey')
            ax.broken_barh([(start_time, end_time - start_time)], (y_pos - 0.4, 0.8), facecolors=color)
            if status not in status_legend_handles:
                status_legend_handles[status] = plt.Rectangle((0, 0), 1, 1, fc=color)

    ax.set_xlabel("Time since first tick (seconds)")
    ax.set_ylabel("Node ID")
    ax.set_title("Node Status Over Time")
    ax.set_yticks(np.arange(len(node_ids)))
    ax.set_yticklabels(node_ids)
    ax.grid(True)
    ax.legend(status_legend_handles.values(), status_legend_handles.keys(), loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=3, frameon=True)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{filename}.pdf"))
