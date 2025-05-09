import pandas as pd
import sys

# Re-import necessary libraries and recreate DataFrame after reset
import matplotlib.pyplot as plt

# Check if input and output files are provided as command line arguments
if len(sys.argv) < 3:
    print("Please provide the input and output files as command line arguments.")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

# Read the input file into a DataFrame
df_single = pd.read_csv(input_file)

# Adjusting the plot to show relative execution time compared to SH, in greyscale
plt.figure(figsize=(12, 8))

# Calculate and plot the relative speedup compared to SH
x_values = df_single['benchmark_name'].unique()
x_ticks = range(len(x_values))

for i, benchmark in enumerate(x_values):
    benchmark_df = df_single[df_single['benchmark_name'] == benchmark]
    avg_hs_time = benchmark_df['hs_time'].mean()
    avg_sh_time = benchmark_df['sh_time'].mean()
    relative_speedup_hs = avg_hs_time / avg_sh_time
    plt.bar(i, relative_speedup_hs, color='grey', label='hS' if i == 0 else "")

plt.axhline(y=1, color='black', linestyle='--', linewidth=2, label='sh')

plt.ylabel('Relative execution time compared to sh', fontsize=14)
plt.xticks(x_ticks, x_values, rotation=45, ha="right", fontsize=12)
plt.legend()
plt.tight_layout()
plt.savefig(output_file)
