# Commonly used emcc compilation variables for Cloudflare Workers

. ./config/shared.default.sh

# Flags from https://github.com/cloudflare/worker-emscripten-template
EM_FLAGS=$(cat <<EOF
    ${EM_FLAGS}
    -s DYNAMIC_EXECUTION=0
    -s TEXTDECODER=0
    -s MODULARIZE=1
    -s ENVIRONMENT="web"
    -s EXPORT_NAME="emscripten"
    --pre-js ../../../config/pre.cloudflare.js
EOF
)

export EM_FLAGS;
export EM_FLAGS_THREADS;
