# Commonly used emcc compilation variables

EM_FLAGS=$(cat <<EOF
    -s USE_ZLIB=1
    -s INVOKE_RUN=0
    -s ALLOW_MEMORY_GROWTH=1
    -s FORCE_FILESYSTEM=1
    -s EXTRA_EXPORTED_RUNTIME_METHODS=["callMain"]
    -lworkerfs.js
EOF
)

export EM_FLAGS;
