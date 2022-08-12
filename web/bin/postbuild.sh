#!/bin/bash

# ==============================================================================
# Post-build workarounds (runs after "npm run build")
# ==============================================================================

# Append Durable Object Class "CDNStats" to Cloudflare Worker script
cat src/lib/stats.js >> .cloudflare/worker.mjs

# Append Cron logic to Cloudflare Worker script
cat src/lib/cron.js >> .cloudflare/worker.mjs

# Cron needs access to biowasm.json info
echo "const BIOWASM_CONFIG = "$(jq -rc . ../biowasm.json) >> .cloudflare/worker.mjs
