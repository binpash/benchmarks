import sys

def create_large_file_with_zeros(file_path, size_in_mb):
    """Create a large file filled with zeros."""
    chunk_size = 1024 * 1024  # 1 MB
    zeros_data = b'\x00' * chunk_size
    with open(file_path, 'wb') as f:
        for _ in range(size_in_mb):
            f.write(zeros_data)

def read_large_file_in_chunks(file_path, chunk_size=1024 * 1024):
    """Read a large file in chunks."""
    with open(file_path, 'rb') as f:
        while True:
            chunk = f.read(chunk_size)
            if not chunk:
                break
            # Process the chunk (here, we'll just print its size)
            # print(f"Read chunk of size: {len(chunk)} bytes")

# Example usage
file_path = sys.argv[1]
size = int(sys.argv[2])
create_large_file_with_zeros(file_path, size)  # Create a 100 MB file filled with zeros
read_large_file_in_chunks(file_path)  # Read the file in 1 MB chunks
