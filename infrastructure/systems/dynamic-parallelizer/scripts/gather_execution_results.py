import os
import csv
import sys

# Window size is taken from the command line arguments
if len(sys.argv) != 3:
    print("Usage: python parse_results.py <window>")
    sys.exit(1)

window = sys.argv[1]

output_dir = os.path.normpath(sys.argv[2])

# Combined results directory
combined_results_dir = os.path.join(output_dir)
os.makedirs(combined_results_dir, exist_ok=True)

# CSV file path
csv_file_path = os.path.join(combined_results_dir, f'results_w{window}.csv')

# Header for the CSV file
csv_header = ['benchmark_name', 'benchmark_dir', 'hs_time', 'sh_time']

# Walk the output directory and parse the times
with open(csv_file_path, 'w', newline='') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow(csv_header)

    for root, dirs, files in os.walk(output_dir):
        if 'hs_time' in files or 'sh_time' in files:
            benchmark_name = os.path.basename(root)
            benchmark_dir = os.path.relpath(root, output_dir)
            
            # Parse hs_time
            hs_time = None
            hs_time_file = os.path.join(root, 'hs_time')
            if os.path.exists(hs_time_file):
                with open(hs_time_file, 'r') as file:
                    hs_time = file.read().strip()

            # Parse sh_time
            sh_time = None
            sh_time_file = os.path.join(root, 'sh_time')
            if os.path.exists(sh_time_file):
                with open(sh_time_file, 'r') as file:
                    sh_time = file.read().strip()

            writer.writerow([benchmark_name, benchmark_dir, hs_time, sh_time])

print(f"Results have been written to {csv_file_path}")
