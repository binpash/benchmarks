#!/usr/bin/env python3
import sys
import random
import json

def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <number of data points>", file=sys.stderr)
        sys.exit(1)

    random.seed(42)
    n = int(sys.argv[1])
    for _ in range(n):
        last_octet = random.randint(0, 255)
        print(json.dumps({"ip": f"192.0.2.{last_octet}"}))

if __name__ == "__main__":
    main()
