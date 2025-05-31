#!/usr/bin/env python3

import re
import argparse
import hashlib
import os
import sys

input_file = "vps-audit-report.txt"
input_file_2 = "vps-audit-negate-report.txt"
error_file = "vps-audit.err"
output_file = "vps-audit-processed.txt"
output_file_2 = "vps-audit-negate-processed.txt"
hash_folder = "hash"

ansi_escape = re.compile(r'\x1b\[[0-9;]*m')

audit_markers = [
    "Current uptime",
    "System up since",
    "System Uptime",
    "System Restart",
    "SSH Root Login",
    "SSH Password Auth",
    "SSH Port",
    "Firewall Status (UFW)",
    "Firewall Status (firewalld)",
    "Firewall Status (iptables)",
    "Firewall Status (nftables)",
    "Firewall Status",
    "Unattended Upgrades",
    "Fail2ban",
    "Failed Logins",
    "System Updates",
    "Running Services",
    "Port Scanning",
    "Port Security",
    "Disk Usage",
    "Memory Usage",
    "CPU Usage",
    "Sudo Logging",
    "Password Policy",
    "SUID Files"
]

sys_info_markers = [
    "Hostname",
    "Operating System",
    "Kernel Version",
    "Kernel",
    "OS",
    "Uptime",
    "CPU Model",
    "CPU Cores",
    "Total Memory",
    "Total Disk Space",
    "Public IP",
    "Load Average"
]

def clean_line(line):
    if line.startswith('-e'):
        line = line[2:].lstrip()  

    line = ansi_escape.sub('', line)

    line = re.sub(r"^Starting audit at [^\x1b]+", "Starting audit at", line)

    line = re.sub(r'vps-audit-report-\S+', '', line)
    line = re.sub(r'vps-audit-negate-report-\S+', '', line)

    if ':' in line and any(marker in line for marker in sys_info_markers):
        line = line.split(':', 1)[0] + ':'

    if any(marker in line for marker in audit_markers):
        line = next((marker for marker in audit_markers if marker in line), line)

    return line.strip() if line.strip() else None

def process_file(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            cleaned_line = clean_line(line)
            if cleaned_line:
                outfile.write(cleaned_line + '\n')

def generate_hash(file_path, hash_folder):
    filename = os.path.splitext(os.path.basename(file_path))[0]
    
    with open(file_path, 'rb') as f:
        hash_value = hashlib.sha256(f.read()).hexdigest()
    
    hash_file_path = os.path.join(hash_folder, f"{filename}.hash")
    with open(hash_file_path, 'w') as hash_file:
        hash_file.write(hash_value)
    
    print(f"{hash_file_path} {hash_value}")

def compare_hashes(file_path, hash_folder):
    filename = os.path.splitext(os.path.basename(file_path))[0]
    
    with open(file_path, 'rb') as f:
        current_hash = hashlib.sha256(f.read()).hexdigest()
    
    hash_file_path = os.path.join(hash_folder, f"{filename}.hash")
    if not os.path.exists(hash_file_path):
        print(f"No hash file found for {file_path}. Please generate the hash first.")
        sys.exit(1)

    with open(hash_file_path, 'r') as hash_file:
        # Read all lines and strip whitespace
        hashes = [line.strip() for line in hash_file.readlines() if line.strip()]
    
    # Check if the current hash matches any of the hashes in the file
    if current_hash not in hashes:
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process and optionally generate hashes for files.")
    parser.add_argument(
        "--generate",
        action="store_true",
        help="If specified, generate hashes for the processed files instead of comparing them."
    )
    args = parser.parse_args()

    os.makedirs(hash_folder, exist_ok=True)

    process_file(input_file, output_file)
    process_file(input_file_2, output_file_2)

    if args.generate:
        generate_hash(output_file, hash_folder)
        generate_hash(output_file_2, hash_folder)
    else:
        compare_hashes(output_file, hash_folder)
        compare_hashes(output_file_2, hash_folder)
