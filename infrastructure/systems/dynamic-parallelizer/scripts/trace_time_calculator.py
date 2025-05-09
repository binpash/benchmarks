import re
from datetime import datetime
from collections import defaultdict

# Function to parse log timestamp
def parse_timestamp(timestamp):
    return datetime.strptime(timestamp, "%Y-%m-%d %H:%M:%S,%f")

# Initialize dictionaries to hold start times, total times, and parse counts for each node ID
start_times = {}
total_times = defaultdict(int)
parse_counts = defaultdict(int)

# Regex pattern to match relevant log lines and capture the action, node ID, timestamp, and trace file
pattern = re.compile(
    r"INFO\|(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},\d{3})\|\[STATE_LOG\] "
    r"Read_trace_(start|end) - node (.+?) - trace (.+)"
)

input_file = "/srv/hs/report/output/benchmarks/max_temp/hs_log"

with open(input_file, "r") as log_file:
    for line in log_file:
        match = pattern.match(line)
        if match:
            timestamp, action, node_id, trace_file = match.groups()
            datetime_timestamp = parse_timestamp(timestamp)
            
            # Create a unique key for each node and trace combination
            key = f"{node_id}|{trace_file}"

            if action == 'start':
                start_times[key] = datetime_timestamp
            elif action == 'end' and key in start_times:
                # Calculate duration in seconds
                duration = (datetime_timestamp - start_times[key]).total_seconds()
                total_times[node_id] += duration
                parse_counts[node_id] += 1
                del start_times[key]  # Clean up

# Combine total times and counts into a single dict for sorting and printing
combined_stats = {node_id: {'time': total_times[node_id], 'count': parse_counts[node_id]}
                  for node_id in total_times}

# Sort the nodes by total time spent
sorted_combined_stats = sorted(combined_stats.items(), key=lambda x: x[1]['time'], reverse=True)

# Output the sorted results
for node_id, stats in sorted_combined_stats:
    print(f"Node ID: {node_id}, Times Parsed: {stats['count']}, Total Time Spent: {stats['time']:.3f} seconds")
