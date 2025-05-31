#!/bin/bash
# This is the main entry point of the search engine.
cd "$(dirname "$0")" || exit 1

while read -r url; do

  if [[ "$url" == "stop" ]]; then
    # stop the engine if it sees the string "stop" 
    exit;
  fi

  echo "[engine] crawling $url">/dev/stderr
  ./crawl.sh "$url" >${OUT}/content.txt
  echo "[engine] indexing $url">/dev/stderr
  ./index.sh ${OUT}/content.txt "$url"

  if  [[ "$(cat ${OUT}/visited.txt | wc -l)" -ge "$(cat ${OUT}/urls.txt | wc -l)" ]]; then
      # stop the engine if it has seen all available URLs
      break;
  fi

done < <(tail -f ${OUT}/urls.txt)
