#!/usr/bin/env python3
import sys
import random
from pathlib import Path
import ipaddress

def generate_unique_ips(count, exclude_set=None):
    exclude_set = exclude_set or set()
    generated = set()
    while len(generated) < count:
        ip = str(ipaddress.IPv4Address(random.randint(0x0B000000, 0xDF000000)))  # 11.0.0.0 to 223.0.0.0
        if ip not in exclude_set:
            generated.add(ip)
    return generated

def main():
    if len(sys.argv) != 3:
        print("Usage: gen_comm_input.py <N> <output_dir>")
        sys.exit(1)

    random.seed(42)
    N = int(sys.argv[1])
    output_dir = Path(sys.argv[2])
    output_dir.mkdir(parents=True, exist_ok=True)

    # Common and unique line counts
    common_ratio = 0.5
    n_common = int(N * common_ratio)
    n_unique_each = N - n_common

    # Generate IPs
    common_ips = generate_unique_ips(n_common)
    unique_ips1 = generate_unique_ips(n_unique_each, exclude_set=common_ips)
    unique_ips2 = generate_unique_ips(n_unique_each, exclude_set=common_ips | unique_ips1)

    file1_lines = sorted(common_ips | unique_ips1)
    file2_lines = sorted(common_ips | unique_ips2)

    # Write files
    (output_dir / "file1").write_text('\n'.join(file1_lines) + '\n')
    (output_dir / "file2").write_text('\n'.join(file2_lines) + '\n')

    print(f"[+] Generated IP lists in {output_dir}")
    print(f"    - Total lines per file: {N}")
    print(f"    - Common IPs: {n_common}")
    print(f"    - Unique IPs per file: {n_unique_each}")
    print(f"    - Files: file1, file2")

if __name__ == "__main__":
    main()
