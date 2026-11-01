#!/usr/bin/env bash
# End-to-end quickstart on the bundled synthetic fixtures: all six steps plus a regression diff.
# Needs: clinval installed, bcftools + samtools on PATH, Docker with the images
# clinval-mason:2.4.0 and clinval-happy:0.3.15 (README, "Quickstart"). Takes about a minute.
#   examples/quickstart/run_quickstart.sh [OUT_DIR]      (default: runs/quickstart)
set -euo pipefail
repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
target="${1:-$repo/runs/quickstart}"
if [[ -e "$target" ]]; then
  echo "$target exists; remove it or pass another directory" >&2
  exit 1
fi
mkdir -p "$target"
out="$(cd "$target" && pwd)"
suite="$repo/tests/fixtures/suite_erepo.yaml"
runtime="$repo/configs/runtime/docker-local-build.yaml"
acceptance="$repo/configs/acceptance/article-default.yaml"

# The systems under test: a toy caller + a toy ACMG classifier (baseline), and the same
# pipeline with the classifier's PP5/BP6 shortcut switched off (candidate).
pipeline() {  # NAME FLAGS
  sed -e "s|@NAME@|$1|g" -e "s|@REPO@|$repo|g" -e "s|@RUN@|$out|g" \
      -e "s|@PYTHON@|$(command -v python3)|g" -e "s|@FLAGS@|$2|g" \
      "$repo/examples/quickstart/toy_pipeline.yaml.in" > "$out/$1.yaml"
}
pipeline toy ""
pipeline toy-noPP5BP6 ', "--no-pp5-bp6"'
base="$out/toy.yaml"
cand="$out/toy-noPP5BP6.yaml"

set -x
clinval validate-suite "$suite"
clinval build-truth --suite "$suite" --out "$out"                        # steps 1-2
clinval simulate    --suite "$suite" --out "$out" --runtime "$runtime"   # step 3
for p in "$base" "$cand"; do
  clinval run-sut --suite "$suite" --out "$out" --pipeline "$p"                       # steps 4-5
  clinval compare --suite "$suite" --out "$out" --pipeline "$p" --runtime "$runtime"  # step 6
done
clinval metrics --pipeline "$base" --run "$out" --out "$out/metrics" --acceptance "$acceptance"
clinval diff --baseline "$base" --candidate "$cand" --run "$out" --out "$out/diff"
clinval report --metrics "$out/metrics" --diff "$out/diff" --out "$out/report.md"
set +x
echo
echo "Done. Read $out/report.md, then $out/sut/toy/compare/outcomes.tsv"
