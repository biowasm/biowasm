# Commonly used Emscripten compilation settings for WebAssembly modules that run on
# Cloudflare Workers. Flags from https://github.com/cloudflare/worker-emscripten-template
# and https://github.com/robertaboukhalil/cf-workers-emscripten.

. ./config/shared.default.sh

# To customize the export name, use `-s EXPORT_NAME="Module"`
EM_FLAGS=$(cat <<EOF
    ${EM_FLAGS}
    -s TEXTDECODER=0
    -s ENVIRONMENT="web"
EOF
)

# Remove extraneous whitespace
EM_FLAGS=$(echo $EM_FLAGS)

export EM_FLAGS;
