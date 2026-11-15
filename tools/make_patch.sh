#!/bin/bash
# Regenerate patches/cp2k-2026.2.diff from src/ against a pristine CP2K 2026.2 tree.
#   usage: tools/make_patch.sh <pristine cp2k-2026.2 root>
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PRISTINE="${1:?usage: tools/make_patch.sh <pristine cp2k-2026.2 root>}"
grep -q 'cp2k_version = "CP2K version 2026.2"' "$PRISTINE/src/cp2k_info.F" || {
  echo "$PRISTINE is not CP2K 2026.2" >&2; exit 1; }
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
mkdir -p "$TMP/a/src" "$TMP/b/src"
while IFS= read -r -d '' f; do
  rel="${f#"$HERE"/src/}"
  mkdir -p "$(dirname "$TMP/a/src/$rel")" "$(dirname "$TMP/b/src/$rel")"
  [ -e "$PRISTINE/src/$rel" ] && cp "$PRISTINE/src/$rel" "$TMP/a/src/$rel"
  cp "$f" "$TMP/b/src/$rel"
done < <(find "$HERE/src" -type f -print0)
( cd "$TMP" && diff -ruN a/src b/src ) > "$HERE/patches/cp2k-2026.2.diff" || true
echo "wrote patches/cp2k-2026.2.diff ($(wc -l < "$HERE/patches/cp2k-2026.2.diff") lines)"
