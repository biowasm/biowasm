# Commonly used Emscripten compilation settings for WebAssembly modules that run in the browser

EM_FLAGS_BASE=$(cat <<EOF
    -s USE_ZLIB=1
    -s INVOKE_RUN=0
    -s FORCE_FILESYSTEM=1
    -s EXPORTED_RUNTIME_METHODS=["callMain"]
    -lworkerfs.js
EOF
)

EM_FLAGS=$(cat <<EOF
    ${EM_FLAGS_BASE}
    -s ALLOW_MEMORY_GROWTH=1
EOF
)

# Remove extraneous whitespace
EM_FLAGS=$(echo $EM_FLAGS)

export EM_FLAGS;
