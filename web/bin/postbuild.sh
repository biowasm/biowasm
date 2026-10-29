#!/bin/bash
set -euo pipefail

# ==============================================================================
# Post-build workarounds (runs after "npm run build")
# ==============================================================================

WORKER=.cloudflare/worker.mjs

# The Cloudflare adapter names its default export differently across versions
# (`entry_default` on adapter-cloudflare-workers, `worker_default` on
# adapter-cloudflare), so look it up instead of hardcoding it. Without this the
# append below references an undefined variable and the Worker throws on startup.
ENTRY=$(grep -oE '^var [A-Za-z0-9_]+_default = \{' "$WORKER" | head -1 | sed 's/^var //; s/ = {$//')
if [ -z "$ENTRY" ]; then
	echo "ERROR: could not find the Worker's default export in $WORKER." >&2
	echo "The Cloudflare adapter's output changed; update bin/postbuild.sh." >&2
	exit 1
fi
echo "const entry_default = $ENTRY;" >> "$WORKER"

# Append Durable Object Class "CDNStats" to Cloudflare Worker script
cat src/lib/stats.js >> "$WORKER"

# Append Cron logic to Cloudflare Worker script
cat src/lib/cron.js >> "$WORKER"

# Cron needs access to biowasm.json info
echo "const BIOWASM_CONFIG = "$(jq -rc . ../biowasm.json) >> "$WORKER"
