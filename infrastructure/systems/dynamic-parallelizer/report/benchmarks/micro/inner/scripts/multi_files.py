import os
import sys

def create_files(directory, file_count, file_size_in_bytes=64):
    """Create many tiny files filled with zeros."""
    zeros_data = b'\x00' * file_size_in_bytes
    if not os.path.exists(directory):
        os.makedirs(directory)
    for i in range(file_count):
        file_path = os.path.join(directory, f'tiny_file_{i}.dat')
        with open(file_path, 'wb') as f:
            f.write(zeros_data)

def read_files_in_directory(directory):
    """Read all files in the given directory."""
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        with open(file_path, 'rb') as f:
            data = f.read()
            # print(f"Read {len(data)} bytes from {filename}")

# Example usage
directory = sys.argv[1]
create_files(directory, 50000)  # Create 1000 tiny files
read_files_in_directory(directory)  # Read all the files in the directory
