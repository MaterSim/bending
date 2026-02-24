#!/bin/bash
set -euo pipefail

BASE_DIR="$(cd "$(dirname "$0")" && pwd)"


# Parameter grids
layers=(10 20 30)
k1s=(100)
k2s=(1 10 100)
k3s=(5)
v=(2e-7)

for layer in "${layers[@]}"; do
  for k1 in "${k1s[@]}"; do
    for k2 in "${k2s[@]}"; do
      for k3 in "${k3s[@]}"; do
        for v in "${v[@]}"; do
            job_name="L${layer}-s${k1}-b${k3}-l${k2}-v${v}"
            job_dir="${BASE_DIR}/Benchmark/${job_name}"
            dump_file="bend-L${layer}-s${k1}-b${k3}-l${k2}.dump"

            cd "${job_dir}"
            python "${BASE_DIR}/plot_combined.py" -d "${dump_file}"
            mv combined_plot.png "../bend-${job_name}.png"
        done
      done
    done
  done
done
