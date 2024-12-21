#!/bin/bash

cd "$(realpath $(dirname "$0"))"
mkdir -p inputs

COMPRESSIONS=("gzip" "bzip2" "xz" "lzo" "lz4" "zstd")

for comp in "${COMPRESSIONS[@]}"; do
  input_dir="./inputs/$comp"
  mkdir -p "$input_dir"
  echo "Sample data for $comp compression" > "$input_dir/example.txt"
  echo "#!/bin/bash" > "$input_dir/startup.sh"
  echo "echo 'Startup script executed for $comp'" >> "$input_dir/startup.sh"
  chmod +x "$input_dir/startup.sh"
done

echo "Input directories prepared."
