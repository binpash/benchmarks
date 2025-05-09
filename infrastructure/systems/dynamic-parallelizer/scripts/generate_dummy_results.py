import os
import random
import argparse

# Set up the argument parser
parser = argparse.ArgumentParser(description='Generate dummy benchmark outputs.')
parser.add_argument('output_dir', help='Path to the output directory where files will be created.')
parser.add_argument('--min', type=float, default=0.5, help='Minimum time in seconds.')
parser.add_argument('--max', type=float, default=2.0, help='Maximum time in seconds.')
parser.add_argument('--reverse', action='store_true', help='Reverse the timing making hs slower than sh.')

# Parse command line arguments
args = parser.parse_args()

# Function to write dummy times to files
def write_dummy_times(directory, min_time, max_time, reverse):
    hs_time = random.uniform(min_time, max_time)
    sh_time = random.uniform(min_time, max_time)

    # Ensure hs is slower than sh if reverse flag is set
    if reverse and hs_time < sh_time:
        hs_time, sh_time = sh_time, hs_time
    elif not reverse and hs_time > sh_time:
        hs_time, sh_time = sh_time, hs_time

    with open(os.path.join(directory, 'hs_time'), 'w') as f:
        f.write(f"{hs_time}\n")
    with open(os.path.join(directory, 'sh_time'), 'w') as f:
        f.write(f"{sh_time}\n")

# Walk the output directory and generate dummy outputs for all directories and subdirectories
for root, dirs, files in os.walk(args.output_dir):
    write_dummy_times(root, args.min, args.max, args.reverse)

print('Dummy benchmark outputs have been generated.')
