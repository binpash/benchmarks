#!/usr/bin/env python3
import sys
import random

def generate_ip():
    if random.random() < 0.5:
        return f"128.151.150.{random.randint(0, 255)}"
    else:
        return ".".join(str(random.randint(0, 255)) for _ in range(4))

def generate_data():
    ip = generate_ip()
    num = random.randint(1, 200)
    line = f"{ip} {num}"
    return line

def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <number of data points>", file=sys.stderr)
        sys.exit(1)

    random.seed(42)
    n = sys.argv[1]
    n = int(n)
    for _ in range(n):
        print(generate_data())

if __name__ == "__main__":
    main()
