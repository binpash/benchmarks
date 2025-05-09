import sys
import os

def inflate_file(input_file, times, inplace=False):
    """Inflate the input file by duplicating its content a specified number of times."""
    # Validate times parameter
    try:
        times = int(times)
        if times <= 0:
            raise ValueError("The times parameter must be a positive integer.")
    except ValueError as e:
        print(f"Error: {e}")
        return

    # Check if input file exists
    if not os.path.isfile(input_file):
        print(f"Error: Input file '{input_file}' not found.")
        return

    # Read the content of the input file
    with open(input_file, 'rb') as f:
        content = f.read()
    input_size = len(content)

    if input_size == 0:
        print(f"Error: Input file '{input_file}' is empty.")
        return

    # Duplicate content in memory
    inflated_content = content * times

    # Determine output file name
    if inplace:
        output_file = input_file
    else:
        output_file = f"{times}x-{os.path.basename(input_file)}"

    # Write the inflated content to the file
    with open(output_file, 'wb') as f:
        f.write(inflated_content)

    if inplace:
        print(f"{input_file} was inflated in-place {times} times.")
    else:
        print(f"{input_file} inflated {times} times as {output_file}.")

if __name__ == "__main__":
    # Check if enough arguments are provided
    if len(sys.argv) < 3:
        print("Usage: python inflate_file.py <input_file> <times> [-i]")
        sys.exit(1)

    input_file = sys.argv[1]
    times = sys.argv[2]
    inplace_flag = '-i' in sys.argv

    inflate_file(input_file, times, inplace=inplace_flag)
