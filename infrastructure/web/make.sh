#!/usr/bin/env bash
set -euo pipefail
set -x

src="index.pd.html"
dest="index.html"

header_content=$(<"header.html")
footer_content=$(<"footer.html")

pandoc "$src"                         \
  --from html                         \
  --template "$src"                   \
  --to html5                          \
  --variable title="benchmarks.sh"    \
  --output "$dest"
