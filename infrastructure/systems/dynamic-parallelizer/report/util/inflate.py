import sys
import os
import shutil

def convert_to_bytes(size_str):
    """Convert a size string (like '5M') to bytes."""
    size_unit = size_str[-1].upper()
    size_value = int(size_str[:-1])

    if size_unit == 'K':
        return size_value * 1024
    elif size_unit == 'M':
        return size_value * 1024 * 1024
    elif size_unit == 'G':
        return size_value * 1024 * 1024 * 1024
    else:
        return int(size_str)  # If no unit, assume it's in bytes

def inflate_file(input_file, target_size_str):
    """Inflate the input file to the target size in memory and write it to disk."""
    target_size_bytes = convert_to_bytes(target_size_str)

    # Check if input file exists
    if not os.path.isfile(input_file):
        print(f"Error: Input file '{input_file}' not found.")
        return

    # Read the content of the input file
    with open(input_file, 'rb') as f:
        content = f.read()
    input_size = len(content)

    # Check if target size is smaller than input file size
    if target_size_bytes <= input_size:
        output_file = f"{target_size_str}-{os.path.basename(input_file)}"
        shutil.copy(input_file, output_file)
        print(f"Warning: Target size {target_size_str} is smaller than or equal to input file size. "
              f"Copied the file to {output_file} instead.")
        return

    # Duplicate content in memory to reach the target size
    while len(content) * 2 <= target_size_bytes:
        content += content  # Exponentially double the content

    # Final adjustment to reach the exact target size
    remaining_bytes = target_size_bytes - len(content)
    if remaining_bytes > 0:
        content += content[:remaining_bytes]

    # Write the inflated content to a new file
    output_file = f"{target_size_str}-{os.path.basename(input_file)}"
    with open(output_file, 'wb') as f:
        f.write(content)

    print(f"{input_file} inflated successfully to {target_size_str} as {output_file}.")

if __name__ == "__main__":
    # Check if enough arguments are provided
    if len(sys.argv) < 3:
        print("Usage: python inflate_file.py <input_file> <target_size_1> [<target_size_2> ...]")
        sys.exit(1)

    input_file = sys.argv[1]
    target_sizes = sys.argv[2:]

    for target_size in target_sizes:
        inflate_file(input_file, target_size)
