#°/bin/bash
cat items.grouped.json | jq -r '.[] | .[] | select(.selected==1) | "\(.id)\t\(.title)"' | sed 's/\"//g'
