import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import numpy as np
from matplotlib.gridspec import GridSpec

sns.set_theme(style='ticks')

# Parse command line arguments
parser = argparse.ArgumentParser(description='Generate box and bar plots for benchmarks.')
parser.add_argument('input_file', type=str, help='Path to the input directory containing the CSV file')
parser.add_argument('output_file', type=str, help='Path to the output PDF file')
args = parser.parse_args()

# Read the data
data = pd.read_csv(args.input_file)

# Ensure numeric types for 'sh' and 'hs'
data['sh'] = pd.to_numeric(data['sh'], errors='coerce')
data['hs'] = pd.to_numeric(data['hs'], errors='coerce')

# Drop rows with NaNs in 'sh'
data.dropna(subset=['sh'], inplace=True)

# Save 'hs' execution time before computing Speedup vs sh
data['hs_exec_time'] = data['hs']

# Calculate the Speedup vs sh for 'hs' compared to 'sh'
data['hs'] = data['sh'] / data['hs']

# Melt the dataframe to long format for plotting, including 'sh' and 'hs_exec_time'
value_vars_to_melt = ['hs']
data_long = data.melt(
    id_vars=['benchmark', 'sub-benchmark', 'sh', 'hs_exec_time'],
    value_vars=value_vars_to_melt,
    var_name='Measurement',
    value_name='Speedup vs sh'
)

# Set global font properties
plt.rcParams.update({'font.size': 8, 'font.family': 'sans-serif'})

# Plotting
fig = plt.figure(figsize=(7, 4.7))
# Create a GridSpec with width ratios 1:2
gs = GridSpec(1, 3, width_ratios=[1, 1.5, 0.0])  # The last 0.1 is to adjust spacing
ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])

# Filter data for boxplot and barplot based on 'sub-benchmark'
grouped_data = data_long[data_long['sub-benchmark'] != "BAR"]
bar_data = data_long[data_long['sub-benchmark'] == "BAR"]

# Split the data based on 'sh' value
data_sh_ge_10 = grouped_data[grouped_data['sh'] >= 10]  # 'sh' >= 10
data_sh_lt_10 = grouped_data[grouped_data['sh'] < 10]   # 'sh' < 10

# Boxplot for grouped benchmarks
sns.boxplot(
    x='sub-benchmark',
    y='Speedup vs sh',
    hue='Measurement',
    data=grouped_data,
    ax=ax1,
    palette='pastel',
    showfliers=False
)
# Scatter plot for 'sh' >= 10 (blue dots)
sns.stripplot(
    x='sub-benchmark',
    y='Speedup vs sh',
    hue='Measurement',
    data=data_sh_ge_10,
    palette=['blue'],
    ax=ax1,
    dodge=True,
    size=7,
    alpha=1,
    jitter=True,
    legend=False
)

# Remove duplicate legend entries
ax1.legend_.remove()

# Overlay red dots for 'sh' < 10 with jitter
# Compute x positions
sub_benchmark_order = grouped_data['sub-benchmark'].unique()
sub_benchmark_map = {name: i for i, name in enumerate(sub_benchmark_order)}
measurement_order = grouped_data['Measurement'].unique()
n_measurements = len(measurement_order)
dodge_amount = 0.8 / n_measurements
measurement_offsets = np.linspace(
    -0.4 + dodge_amount / 2,
    0.4 - dodge_amount / 2,
    n_measurements
)
measurement_offset_map = {
    measurement: offset for measurement, offset in zip(measurement_order, measurement_offsets)
}

# Calculate x positions for red dots
x_positions = data_sh_lt_10.apply(
    lambda row: sub_benchmark_map[row['sub-benchmark']] +
    measurement_offset_map[row['Measurement']],
    axis=1
)

# Add jitter to x positions
jitter_amount = dodge_amount / 5  # Adjust the jitter amount as needed
jitter = np.random.uniform(-jitter_amount, jitter_amount, size=len(x_positions))
x_positions += jitter

y_positions = data_sh_lt_10['Speedup vs sh']

# Plot red dots with jitter
ax1.scatter(x_positions, y_positions, color='red', s=30, zorder=10)

# Group 'bar_data' to get one data point per bar
bar_data_grouped = bar_data.groupby(['benchmark', 'Measurement']).agg({
    'Speedup vs sh': 'mean',
    'hs_exec_time': 'mean'
}).reset_index()

# Barplot for "BAR" benchmarks
sns.barplot(
    x='benchmark',
    y='Speedup vs sh',
    hue='Measurement',
    data=bar_data_grouped,
    ax=ax2,
    palette='pastel',
    ci=None
)

# Annotate bars with hs execution time
measurement_order = bar_data_grouped['Measurement'].unique()

# for container, measurement in zip(ax2.containers, measurement_order):
#     # Get the data for this measurement
#     data_subset = bar_data_grouped[bar_data_grouped['Measurement'] == measurement]
#     for bar, (_, row) in zip(container, data_subset.iterrows()):
#         x = bar.get_x() + bar.get_width() / 2
#         y = bar.get_height()
#         hs_time = row['hs_exec_time']
#         ax2.text(x, y, f'{hs_time:.0f}s', ha='center', va='bottom', fontsize=8, rotation=90)

# Customize plot elements
y_ticks = np.array([0.125, 0.25, 0.5, 1, 2, 4, 8])
y_tick_labels = ['0.125×', '0.25×', '0.5×', '1×', '2×', '4×', '8×']

# Set y-scale and ticks for both axes
for ax in [ax1, ax2]:
    ax.set_yscale('log', base=2)
    ax.set_ylim([0.125, 8])  # Adjust the limits to include all ticks
    ax.axhline(y=1, color='black', linestyle='--', label='sh (Baseline)')
    ax.set_xlabel('')
    ax.tick_params(axis='x', rotation=45)
    # Add light dashed lines at specified ticks (excluding y=1)
    for y in y_ticks:
        if y != 1:
            ax.axhline(y=y, color='gray', linestyle='--', linewidth=0.5, zorder=0)

# Set y-ticks and labels for left plot (ax1)
ax1.set_yticks(y_ticks)
ax1.set_yticklabels(y_tick_labels)
ax1.set_ylabel('Relative performance vs sh', fontsize=16)

# Remove y-ticks and labels from right plot (ax2)
ax2.set_yticks(y_ticks)
ax2.set_yticklabels([])
ax2.tick_params(axis='y', which='both', left=False, labelleft=False)
ax2.set_ylabel('')

# Remove duplicate legend
ax2.get_legend().remove()

# Align all x-axis labels at the right-side of the text
for ax in [ax1, ax2]:
    for label in ax.get_xticklabels():
        label.set_horizontalalignment('right')

plt.tight_layout(rect=[0, 0.05, 1, 1])
plt.savefig(args.output_file, bbox_inches='tight', pad_inches=0.05)
