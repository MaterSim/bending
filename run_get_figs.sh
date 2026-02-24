#!/bin/bash
set -euo pipefail

BASE_DIR="$(cd "$(dirname "$0")" && pwd)"
LOG_DIR="${BASE_DIR}/logs"
FIG_ROOT="${BASE_DIR}/Figs"

mkdir -p "$LOG_DIR" "$FIG_ROOT"

# Activate Python environment (conda htocsp) for cron
CONDA_SH="/users/qzhu8/miniconda3/etc/profile.d/conda.sh"
if [ -f "$CONDA_SH" ]; then
  # shellcheck disable=SC1091
  source "$CONDA_SH"
  conda activate htocsp || true
elif command -v conda >/dev/null 2>&1; then
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate htocsp || true
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
OUT_DIR="${FIG_ROOT}/${STAMP}"
mkdir -p "$OUT_DIR"

# Generate figures
bash "${BASE_DIR}/get_fig.sh" | tee "${LOG_DIR}/get_fig_${STAMP}.log"

# Move generated figures into timestamped folder
find "${BASE_DIR}/Benchmark" -maxdepth 1 -type f -name "bend-*.png" -print0 \
  | xargs -0 -I{} mv "{}" "${OUT_DIR}/"

# Generate HTML gallery and index
python "${BASE_DIR}/generate_html.py" --outdir "${OUT_DIR}" --index

# Commit and push to GitHub (for GitHub Pages auto-deploy)
if command -v git >/dev/null 2>&1; then
  cd "${BASE_DIR}"
  git add Figs/
  git commit -m "Update figures: ${STAMP}" || true
  git push || true
fi

echo "Done: ${OUT_DIR}"
