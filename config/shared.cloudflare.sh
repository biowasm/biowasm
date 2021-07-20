# Commonly used Emscripten compilation settings for WebAssembly modules that run on Cloudflare Workers
# Flags from https://github.com/cloudflare/worker-emscripten-template

. ./config/shared.default.sh

EM_FLAGS=$(cat <<EOF
    ${EM_FLAGS}
    -s DYNAMIC_EXECUTION=0
    -s TEXTDECODER=0
    -s MODULARIZE=1
    -s ENVIRONMENT="web"
    -s EXPORT_NAME="emscripten"
    --pre-js ../../../config/shared.cloudflare.js
EOF
)

# Remove extraneous whitespace
EM_FLAGS=$(echo $EM_FLAGS)

export EM_FLAGS;
