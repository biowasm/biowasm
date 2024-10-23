# Commonly used Emscripten compilation settings for WebAssembly modules that run in the browser.
# To load settings for development: . bin/shared.sh
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
EM_FLAGS=$(echo $EM_FLAGS)  # Remove whitespace

# Function to remove Nanosleep requirements for GNU tools (Nanosleep not supported in Emscripten)
function EM_GNU_NANOSLEEP() {
	echo "Running EM_GNU_NANOSLEEP..."
	sed -i 's|if ${gl_cv_func_sleep_works+:} false|if true|g' configure
	sed -i 's|if ${ac_cv_search_nanosleep+:} false|if true|g' configure
	sed -i 's|if ${gl_cv_func_nanosleep+:} false|if true|g' configure

	sed -i 's|if test ${gl_cv_func_sleep_works+y}|if true|g' configure
	sed -i 's|if test ${ac_cv_search_nanosleep+y}|if true|g' configure
	sed -i 's|if test ${gl_cv_func_nanosleep+y}|if true|g' configure
}

# Prevent ./configure from stalling at "checking whether strcasestr works in linear time..."
function EM_GNU_STRCASESTR_LINEAR() {
	echo "Running EM_GNU_STRCASESTR_LINEAR..."
	sed -i 's|if ${gl_cv_func_strcasestr_linear+:} false|if true|g' configure

	sed -i 's|if test ${gl_cv_func_strcasestr_linear+y}|if true|g' configure
}

# Shared LZMA flags
LZMA_VERSION="5.6.3"
DIR_LZMA=../../htslib/src/xz-${LZMA_VERSION}/src/liblzma
CFLAGS_LZMA="-I${DIR_LZMA}/api -I${DIR_LZMA}/api/lzma"
LDFLAGS_LZMA="-L${DIR_LZMA}/.libs"

# Export vars and functions
export EM_FLAGS CFLAGS_LZMA LDFLAGS_LZMA;
export -f EM_GNU_NANOSLEEP EM_GNU_STRCASESTR_LINEAR;
