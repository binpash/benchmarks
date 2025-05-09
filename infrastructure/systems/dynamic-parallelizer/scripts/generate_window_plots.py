import os
import glob
import sys
import pandas as pd
import argparse
import matplotlib.pyplot as plt

# Set global font properties
plt.rcParams.update({'font.size': 22, 'font.family': 'serif'})

# Parse command line arguments
parser = argparse.ArgumentParser(description='Generate window plots')
parser.add_argument('input_directory', type=str, help='Path to the input directory')
parser.add_argument('output_file', type=str, help='Path to the output file')
parser.add_argument('--benchmarks', type=str, nargs='+', help='List of benchmarks to plot')
args = parser.parse_args()

csv_files = glob.glob(os.path.join(sys.argv[1], "*.csv"))
csv_files = sorted(csv_files, key=lambda x: int(x.split("_w")[1].split(".")[0]))

benchmark_results = []

for csv_file in csv_files:
    df_temp = pd.read_csv(csv_file)
    window = int(csv_file.split("_w")[1].split(".")[0])  # Extracting window number from filename
    for _, row in df_temp.iterrows():
        if args.benchmarks is None or row["benchmark_name"] in args.benchmarks:
            benchmark_results.append({
                "benchmark_name": row["benchmark_name"],
                "window": window,
            "hs_time": row["hs_time"],
            "base_time": row["sh_time"] if "w0" in csv_file else None  # Base time is sh_time from results_w0.csv
        })

df = pd.DataFrame(benchmark_results)

# Fill in base_time for all records based on benchmark_name matching
df['base_time'] = df.groupby('benchmark_name')['base_time'].transform(lambda x: x.ffill().bfill())
# Calculate hs_time / base_time ratio
df['hs_base_ratio'] = df['hs_time'] / df['base_time']
# Calculate relative speedup as base_time / hs_time
df['relative_speedup'] = df['base_time'] / df['hs_time']



plt.figure(figsize=(10, 6))
markers = ['o', 'v', '^', '<', '>', "s"]
colors = ['#FFB6C1', '#ADD8E6', '#90EE90', '#FFDAB9', '#87CEFA', '#98FB98']

# Relative speedup
for i, benchmark_name in enumerate(df['benchmark_name'].unique()):
    benchmark_df = df[df['benchmark_name'] == benchmark_name]
    plt.plot(benchmark_df['window'], benchmark_df['relative_speedup'], label=benchmark_name, linestyle="-", marker=markers[i % len(markers)], color=colors[i % len(colors)], linewidth=3, markersize=8)

# Update y-axis label to reflect the change to relative speedup
plt.ylabel('Relative Performance', fontsize=21, fontfamily='serif')

# Set axis labels with specific font properties
plt.xlabel('Window', fontsize=21, fontfamily='serif')

plt.legend(title='', fontsize='13.4')
plt.grid(True, which='both', linestyle='--', linewidth=0.3)

# Log scale on x-axis (messes up 0 value)
# plt.xscale('log', base=2)  # Corrected log2 scale setting
# plt.gca().xaxis.set_major_formatter(plt.FormatStrFormatter('%d'))
# Set x-axis and y-axis tick label properties

# Get current y-axis tick labels
current_yticks = plt.gca().get_yticks()

# Format y-axis tick labels
formatted_yticklabels = [f'{label:.0f}x' for label in current_yticks]

# Set formatted y-axis tick labels with specific font properties and tilt them
plt.yticks(ticks=current_yticks, labels=formatted_yticklabels, fontsize=15, fontfamily='serif')
# set_major_formatter(plt.FuncFormatter(lambda x, _: f'' if x % 1 == 0.5 else f'{x:.0f}x'))

# Set formatted y-axis tick labels with specific font properties
plt.yticks(ticks=current_yticks, labels=formatted_yticklabels, fontsize=20, fontfamily='serif')


plt.xticks(fontsize=17, fontfamily='serif')
plt.yticks(fontsize=17, fontfamily='serif')
plt.xticks(df['window'].unique())
plt.tight_layout()

plt.savefig(sys.argv[2], bbox_inches='tight', pad_inches=0.05)
