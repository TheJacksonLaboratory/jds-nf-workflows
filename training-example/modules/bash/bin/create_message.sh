#!/usr/bin/env bash

set -euo pipefail

message="$1"
output_file="$2"

printf '%s\n' "$message" > "$output_file"