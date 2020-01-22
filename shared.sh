# Commonly used emcc compilation variables

EM_FLAGS_BASE=$(cat <<EOF
    -s USE_ZLIB=1
    -s INVOKE_RUN=0
    -s FORCE_FILESYSTEM=1
    -s EXTRA_EXPORTED_RUNTIME_METHODS=["callMain"]
    -lworkerfs.js
EOF
)

EM_FLAGS=$(cat <<EOF
    ${EM_FLAGS_BASE}
    -s ALLOW_MEMORY_GROWTH=1
EOF
)

# Note: avoid using memory growth with threads (see https://emscripten.org/docs/porting/pthreads.html)
EM_FLAGS_THREADS=$(cat <<EOF
    ${EM_FLAGS_BASE}
    -s USE_PTHREADS=1
    -s PTHREAD_POOL_SIZE=10
    -s TOTAL_MEMORY=100MB
    --threadprofiler
EOF
)

export EM_FLAGS;
export EM_FLAGS_THREADS;
