#!/bin/bash
set -euo pipefail

# ==============================================================================
# Post-build workarounds (runs after "npm run build")
# ==============================================================================

WORKER=.cloudflare/worker.mjs

# Adapter versions name the default export differently, so find it instead of guessing
ENTRY=$(grep -oE '^var [A-Za-z0-9_]+_default = \{' "$WORKER" | head -1 | sed 's/^var //; s/ = {$//')
if [ -z "$ENTRY" ]; then
	echo "ERROR: no default export in $WORKER; the Cloudflare adapter's output changed." >&2
	exit 1
fi
echo "const entry_default = $ENTRY;" >> "$WORKER"

# Append Durable Object Class "CDNStats" to Cloudflare Worker script
cat src/lib/stats.js >> "$WORKER"

# Append Cron logic to Cloudflare Worker script
cat src/lib/cron.js >> "$WORKER"

# Cron needs access to biowasm.json info
echo "const BIOWASM_CONFIG = "$(jq -rc . ../biowasm.json) >> "$WORKER"
