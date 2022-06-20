# Commonly used Emscripten compilation settings for WebAssembly modules that run in the browser
EM_FLAGS=$(cat <<EOF
    -s USE_ZLIB=1
    -s INVOKE_RUN=0
    -s FORCE_FILESYSTEM=1
    -s EXPORTED_RUNTIME_METHODS=["callMain","FS","PROXYFS","WORKERFS"]
    -s MODULARIZE=1
    -s ENVIRONMENT="web,worker"
    -s ALLOW_MEMORY_GROWTH=1
    -lworkerfs.js -lproxyfs.js
EOF
)
# Remove whitespace
EM_FLAGS=$(echo $EM_FLAGS)

# Nanosleep not supported in Emscripten
function EM_GNU_NANOSLEEP() {
    sed -i 's|if ${gl_cv_func_sleep_works+:} false|if true|g' configure
    sed -i 's|if ${ac_cv_search_nanosleep+:} false|if true|g' configure
    sed -i 's|if ${gl_cv_func_nanosleep+:} false|if true|g' configure
}

# Export
export EM_FLAGS;
export EM_GNU_NANOSLEEP;
