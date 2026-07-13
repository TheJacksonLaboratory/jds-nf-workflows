#!/usr/bin/env bash
set -euo pipefail

bad_reads="$1"

cat "$bad_reads" \
| awk '{print $5}' \
| sed -r 's/\,//g' \
| sort -n \
| uniq -c \
| awk '{print $2}' \
> ReadName_unique
