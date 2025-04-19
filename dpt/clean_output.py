import re
import sys

def clean_output_file(input_path, output_path):
    png_line_pattern = re.compile(r'^"\S+\.png"(?:\s.+)?$')

    ansi_escape = re.compile(r'\x1B\[[0-?]*[ -/]*[@-~]')
    backspace = re.compile(r'\x08')

    with open(input_path, 'r', encoding='utf-8', errors='ignore') as infile, \
         open(output_path, 'w', encoding='utf-8') as outfile:

        for line in infile:
            line_clean = ansi_escape.sub('', line)
            line_clean = backspace.sub('', line_clean).strip()

            if png_line_pattern.match(line_clean):
                outfile.write(line_clean + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python clean_output.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    clean_output_file(input_file, output_file)
