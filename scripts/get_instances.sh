#!/usr/bin/env bash
# Downloads the Anjos-large and Sko-large SRFLP instances (Rubio-Sanchez et al., 2016)
# from the OptSicom benchmark page and extracts the 40 instances used in the paper
# (n = 200, 300, 400, 500) into instances/Anjos-large and instances/Sko-large.
# The files are not modified; their SHA-256 checksums are verified against
# instances/large_instances.sha256 (the files used in the paper).
set -euo pipefail

URL="https://grafo.etsii.urjc.es/optsicom/srflp/srflp_instances.zip"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

echo "Downloading $URL ..."
curl -fsSL -o "$TMP/srflp_instances.zip" "$URL"

for n in 200 300 400 500; do
  unzip -q -o -j "$TMP/srflp_instances.zip" "Anjos-large/SRFLP_Rubio_AnKeVa_${n}_*.txt" -d "$ROOT/instances/Anjos-large"
  unzip -q -o -j "$TMP/srflp_instances.zip" "sko-large/SRFLP_Rubio_Sko_${n}_*.txt" -d "$ROOT/instances/Sko-large"
done

cd "$ROOT"
if command -v sha256sum >/dev/null 2>&1; then
  sha256sum -c --quiet instances/large_instances.sha256
else
  shasum -a 256 -c --quiet instances/large_instances.sha256
fi
echo "Done: 40 instances verified (instances/Anjos-large, instances/Sko-large)."
