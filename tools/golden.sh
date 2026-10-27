#!/usr/bin/env bash
# writes the golden stage dumps of the C++ metriko at the commit in rust/PORTED_FROM:
#   testdata/golden/<commit>/<stem>@<scale>/   (tools/dump/FORMAT.md)
# the library is taken from that commit through a temporary `git worktree`; tools/dump itself is taken from the
# working tree (it is checked against the library: it fails if it no longer mirrors that commit's lib.h).
#
#   tools/golden.sh [case-filter]        e.g. tools/golden.sh spot
# env: METRIKO_MODELS (models dir, default from testdata/cases.json)
#      METRIKO_DEPS   (dir holding libigl-src v2.6.4 and eigen-src 5.0.1, default metriko_web's native build deps;
#                      unset it to an empty string to let CMake fetch them)
#      METRIKO_GOLDEN_BUILD (build dir of the dump tool, default ${TMPDIR:-/tmp}/metriko-golden-<commit>)
set -euo pipefail

root=$(git -C "$(dirname "$0")/.." rev-parse --show-toplevel)
commit=$(tr -d '[:space:]' < "$root/rust/PORTED_FROM")
short=${commit:0:7}
filter=${1:-}
cases="$root/testdata/cases.json"
models=${METRIKO_MODELS:-$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["models_default"])' "$cases")}
deps=${METRIKO_DEPS-/Users/saki/dev/metriko_web/wasm/build-native/_deps}
work=${METRIKO_GOLDEN_BUILD:-${TMPDIR:-/tmp}/metriko-golden-$short}
out="$root/testdata/golden/$short"
bin="$work/build/metriko_dump"

# the library at the commit, through a temporary worktree; the build dir is kept, so rebuilding is incremental
src="$work/src"
mkdir -p "$work"
git -C "$root" worktree remove --force "$src" >/dev/null 2>&1 || rm -rf "$src"
git -C "$root" worktree add --detach "$src" "$commit" >/dev/null
trap 'git -C "$root" worktree remove --force "$src" >/dev/null 2>&1 || true' EXIT
fetch=()
if [[ -n "$deps" ]]; then
    fetch=(-DFETCHCONTENT_SOURCE_DIR_LIBIGL="$deps/libigl-src" -DFETCHCONTENT_SOURCE_DIR_EIGEN="$deps/eigen-src" -DFETCHCONTENT_FULLY_DISCONNECTED=ON)
fi
echo "building tools/dump against $short ..."
log="$work/build.log"
{ cmake -S "$root/tools/dump" -B "$work/build" -DMETRIKO_DIR="$src" "${fetch[@]}" && cmake --build "$work/build" -j8; } > "$log" 2>&1 \
    || { tail -30 "$log"; echo "build failed, see $log"; exit 1; }

mkdir -p "$out"
python3 -c '
import json, sys
for c in json.load(open(sys.argv[1]))["cases"]:
    print(c["mesh"], c["scale"])
' "$cases" | while read -r mesh scale; do
    stem=${mesh%.*}
    [[ -n "$filter" && "$stem" != *"$filter"* ]] && continue
    dir="$out/$stem@$scale"
    rm -rf "$dir"
    mkdir -p "$dir"
    t0=$(date +%s)
    if "$bin" "$models/$mesh" "$scale" "$dir" > "$dir/dump.log" 2>&1; then
        printf '%-24s %4ss  %6s  %s\n' "$stem@$scale" "$(( $(date +%s) - t0 ))" "$(du -sh "$dir" | cut -f1)" "$(grep -m1 '^hash' "$dir/dump.log")"
    else
        printf '%-24s FAILED (see %s)\n' "$stem@$scale" "$dir/dump.log"
        tail -3 "$dir/dump.log"
    fi
done
