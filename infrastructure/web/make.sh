#!/usr/bin/env bash
set -euo pipefail
set -x

src="html/index.html"
dest="index.html"
repo="binpash/benchmarks"

header_content=$(<"header.html")
footer_content=$(<"footer.html")
updated="$(LANG=en_us_88591; date +'%R'; date +'%m/%d/%Y')"
version=$(curl -s "https://api.github.com/repos/$repo/tags" | jq -r '.[0].name' || echo "0.1"); [ "$version" = "null" ] && version="0.1"
revision=$(curl -s "https://api.github.com/repos/$repo/commits" | jq -r '.[0].sha' | cut -c1-7 || echo "unknown")
commitmsg=$(curl -s "https://api.github.com/repos/$repo/commits/$revision" | jq -r '.commit.message' || echo "unknown")

pandoc "$src"                         \
  --standalone                        \
  --from html                         \
  --template "$src"                   \
  --to html5                          \
  --variable title="benchmarks.sh"    \
  --variable version="$version"       \
  --variable revision="$revision"     \
  --variable updated="$updated"       \
  --variable commitmsg="$commitmsg"   \
  --output "$dest"
