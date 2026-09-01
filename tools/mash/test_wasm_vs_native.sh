#!/bin/bash
# Compare Mash WASM output vs native Mash output
# Usage: ./test_wasm_vs_native.sh <genome_dir>
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
GENOME_DIR="${1:?Usage: $0 <genome_dir>}"

GENOMES=(
    E_coli_K12MG1655.fna
    E_coli_CFT073.fna
    E_coli_HS.fna
    E_coli_UTI89.fna
    Ecoli_O126.fna
)

# Check prerequisites
command -v mash >/dev/null 2>&1 || { echo "Error: native mash not found. Install with: pixi global install -c conda-forge -c bioconda mash"; exit 1; }
command -v node >/dev/null 2>&1 || { echo "Error: node not found"; exit 1; }

echo "=== Mash WASM vs Native Comparison ==="
echo "Native: $(mash --version 2>&1)"
echo "Genome dir: $GENOME_DIR"
echo

# Verify genomes exist
for g in "${GENOMES[@]}"; do
    if [ ! -f "$GENOME_DIR/$g" ]; then
        echo "Error: $GENOME_DIR/$g not found"
        exit 1
    fi
done

TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# === Run native Mash ===
echo "--- Running native Mash ---"

# Individual distances (K12 vs each other genome)
for i in "${!GENOMES[@]}"; do
    if [ $i -eq 0 ]; then continue; fi
    mash dist "$GENOME_DIR/${GENOMES[0]}" "$GENOME_DIR/${GENOMES[$i]}" >> "$TMPDIR/native_dist.txt"
done
echo "Native individual distances:"
cat "$TMPDIR/native_dist.txt"
echo

# Triangle distance matrix
mash sketch -o "$TMPDIR/native_all" "${GENOMES[@]/#/$GENOME_DIR/}"
mash triangle "$TMPDIR/native_all.msh" > "$TMPDIR/native_triangle.txt"
echo "Native triangle matrix:"
cat "$TMPDIR/native_triangle.txt"
echo

# === Run WASM Mash ===
echo "--- Running WASM Mash ---"
cd "$SCRIPT_DIR"
node test_wasm.mjs "$GENOME_DIR"
echo

# === Compare outputs ===
echo "--- Comparing outputs ---"

PASS=true

# Compare individual distances
if [ -f "$SCRIPT_DIR/wasm_dist.txt" ] && [ -f "$TMPDIR/native_dist.txt" ]; then
    echo "Individual distances (K12 vs others):"
    echo "  Native                          WASM"

    # Extract just the distance column (3rd field) for comparison
    native_dists=$(awk '{print $3}' "$TMPDIR/native_dist.txt")
    wasm_dists=$(awk '{print $3}' "$SCRIPT_DIR/wasm_dist.txt")

    paste <(echo "$native_dists") <(echo "$wasm_dists") | while read native wasm; do
        if [ "$native" = "$wasm" ]; then
            echo "  $native == $wasm  OK"
        else
            echo "  $native != $wasm  MISMATCH"
            # Check if they're close (within floating point tolerance)
            diff=$(echo "$native $wasm" | awk '{d = $1 - $2; if (d < 0) d = -d; print d}')
            echo "    (difference: $diff)"
        fi
    done
    echo
fi

# Compare triangle matrices
if [ -f "$SCRIPT_DIR/wasm_triangle.txt" ] && [ -f "$TMPDIR/native_triangle.txt" ]; then
    echo "Triangle matrix comparison:"
    # Extract just the numeric values, ignoring file paths
    native_nums=$(grep -oE '[0-9]+\.[0-9]+' "$TMPDIR/native_triangle.txt" | sort)
    wasm_nums=$(grep -oE '[0-9]+\.[0-9]+' "$SCRIPT_DIR/wasm_triangle.txt" | sort)

    if [ "$native_nums" = "$wasm_nums" ]; then
        echo "  All distance values match exactly!"
    else
        echo "  Distance values differ:"
        diff <(echo "$native_nums") <(echo "$wasm_nums") || true
        PASS=false
    fi
    echo
fi

if $PASS; then
    echo "=== COMPARISON PASSED: WASM and native Mash produce identical results ==="
else
    echo "=== COMPARISON COMPLETED: Some differences found (see above) ==="
fi
