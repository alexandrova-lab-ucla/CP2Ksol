#!/bin/bash
# Copy the model sources into a CP2K 2026.2 source tree.
#   usage: tools/install.sh <cp2k-2026.2 root>
# Every file that is replaced is kept as <file>.orig. Refuses anything that is not CP2K 2026.2.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ROOT="${1:?usage: tools/install.sh <cp2k-2026.2 root>}"
SRC="$ROOT/src"
[ -f "$SRC/cp2k_info.F" ] || { echo "no $SRC/cp2k_info.F: not a CP2K source tree" >&2; exit 1; }
grep -q 'cp2k_version = "CP2K version 2026.2"' "$SRC/cp2k_info.F" || {
  echo "$SRC is not CP2K 2026.2; this model is written against 2026.2 only" >&2; exit 1; }
n_new=0; n_mod=0
# files keep their path below src/ (src/motion/md_energies.F, src/start/cp2k.F, ...)
while IFS= read -r -d '' f; do
  rel="${f#"$HERE"/src/}"
  mkdir -p "$(dirname "$SRC/$rel")"
  if [ -e "$SRC/$rel" ]; then
    [ -e "$SRC/$rel.orig" ] || cp -p "$SRC/$rel" "$SRC/$rel.orig"
    n_mod=$((n_mod + 1))
  else
    n_new=$((n_new + 1))
  fi
  cp "$f" "$SRC/$rel"
done < <(find "$HERE/src" -type f -print0)
echo "installed into $SRC: $n_mod files replaced (originals kept as *.orig), $n_new files added"
