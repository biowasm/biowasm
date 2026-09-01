#!/bin/bash
# Mash v2.3 -> WebAssembly compile script for biowasm
# Dependencies: Cap'n Proto (runtime), Boost (headers only)
set -e

# =============================================================================
# Step 1: Build Cap'n Proto runtime library with Emscripten
# =============================================================================

CAPNP_VERSION="0.7.0"
CAPNP_DIR="$(pwd)/capnp-wasm"
CAPNP_SRC="capnproto-c++-${CAPNP_VERSION}"

if [ ! -f "$CAPNP_DIR/lib/libcapnp.a" ]; then
    echo "=== Building Cap'n Proto ${CAPNP_VERSION} runtime for WASM ==="

    # Download source
    if [ ! -d "$CAPNP_SRC" ]; then
        curl -L -o capnp.tar.gz "https://capnproto.org/capnproto-c++-${CAPNP_VERSION}.tar.gz"
        tar xzf capnp.tar.gz
        rm capnp.tar.gz
    fi

    mkdir -p "$CAPNP_DIR/lib" "$CAPNP_DIR/include"

    cd "$CAPNP_SRC/src"

    # Fix arenaSpace size for wasm32 (ReaderArena is larger than arenaSpace on 32-bit)
    sed -i 's/void\* arenaSpace\[18 + sizeof(kj::MutexGuarded<void\*>) \/ sizeof(void\*)\]/void* arenaSpace[20 + sizeof(kj::MutexGuarded<void*>) \/ sizeof(void*)]/' capnp/message.h

    # Compile kj library (excluding tests and compiler-specific code)
    echo "--- Compiling kj library ---"
    KJ_SRCS=$(find kj -name '*.c++' \
        ! -name '*-test*' ! -name '*_test*' ! -name 'test*' \
        ! -path '*/compat/*' \
        ! -name 'main.c++' ! -name 'parse/*' \
        | sort)

    mkdir -p kj_objs
    for src in $KJ_SRCS; do
        obj="kj_objs/$(basename ${src%.c++}).o"
        echo "  Compiling $src -> $obj"
        em++ -O2 -std=c++14 -fexceptions -I. -c "$src" -o "$obj" 2>/dev/null || true
    done
    emar rcs "$CAPNP_DIR/lib/libkj.a" kj_objs/*.o

    # Compile capnp library (excluding tests, compiler, and RPC code)
    echo "--- Compiling capnp library ---"
    CAPNP_SRCS=$(find capnp -maxdepth 1 -name '*.c++' \
        ! -name '*-test*' ! -name '*_test*' ! -name 'test*' \
        ! -name 'compiler*' \
        | sort)

    mkdir -p capnp_objs
    for src in $CAPNP_SRCS; do
        obj="capnp_objs/$(basename ${src%.c++}).o"
        echo "  Compiling $src -> $obj"
        em++ -O2 -std=c++14 -fexceptions -I. -c "$src" -o "$obj" 2>/dev/null || true
    done
    emar rcs "$CAPNP_DIR/lib/libcapnp.a" capnp_objs/*.o

    # Copy headers
    cp -r kj "$CAPNP_DIR/include/"
    cp -r capnp "$CAPNP_DIR/include/"

    cd ../..
    echo "=== Cap'n Proto WASM build complete ==="
else
    echo "=== Cap'n Proto WASM libraries already built, skipping ==="
fi

# =============================================================================
# Step 2: Download Boost headers (for -DUSE_BOOST)
# =============================================================================

BOOST_VERSION="1_73_0"
BOOST_DIR="$(pwd)/boost-headers"

if [ ! -d "$BOOST_DIR/boost" ]; then
    echo "=== Downloading Boost headers ==="
    curl -L -o boost.tar.gz \
        "https://sourceforge.net/projects/boost/files/boost/1.73.0/boost_${BOOST_VERSION}.tar.gz/download"
    tar xzf boost.tar.gz "boost_${BOOST_VERSION}/boost/"
    mv "boost_${BOOST_VERSION}" "$BOOST_DIR"
    rm boost.tar.gz
    echo "=== Boost headers downloaded ==="
else
    echo "=== Boost headers already present, skipping ==="
fi

# =============================================================================
# Step 3: Generate Cap'n Proto schema C++ files (native capnp compiler)
# =============================================================================

echo "=== Generating Cap'n Proto schema ==="
cd src/mash/capnp
capnp compile -I "$CAPNP_DIR/include" -oc++ MinHash.capnp
cd ../../..

# =============================================================================
# Step 4: Compile all Mash source files
# =============================================================================

echo "=== Compiling Mash source files ==="

CXXFLAGS="-O2 -std=c++14 -DUSE_BOOST -s USE_ZLIB=1 -fexceptions"
INCLUDES="-Isrc -I${CAPNP_DIR}/include -I${BOOST_DIR}"

MASH_SRCS="
    src/mash/mash.cpp
    src/mash/Command.cpp
    src/mash/CommandBounds.cpp
    src/mash/CommandContain.cpp
    src/mash/CommandDistance.cpp
    src/mash/CommandFind.cpp
    src/mash/CommandInfo.cpp
    src/mash/CommandList.cpp
    src/mash/CommandPaste.cpp
    src/mash/CommandScreen.cpp
    src/mash/CommandSketch.cpp
    src/mash/CommandTaxScreen.cpp
    src/mash/CommandTriangle.cpp
    src/mash/hash.cpp
    src/mash/HashList.cpp
    src/mash/HashPriorityQueue.cpp
    src/mash/HashSet.cpp
    src/mash/MinHashHeap.cpp
    src/mash/MurmurHash3.cpp
    src/mash/Sketch.cpp
    src/mash/sketchParameterSetup.cpp
"

OBJS=""
for src in $MASH_SRCS; do
    obj="${src%.cpp}.o"
    echo "  Compiling $src"
    em++ $CXXFLAGS $INCLUDES -c "$src" -o "$obj"
    OBJS="$OBJS $obj"
done

# Compile capnp generated schema file
echo "  Compiling src/mash/capnp/MinHash.capnp.c++"
em++ $CXXFLAGS $INCLUDES -c src/mash/capnp/MinHash.capnp.c++ -o src/mash/capnp/MinHash.capnp.o
OBJS="$OBJS src/mash/capnp/MinHash.capnp.o"

# Compile C++ ABI exception stubs (needed for kj library in Emscripten 2.0.25)
cat > cxa_stubs.cpp << 'STUBEOF'
namespace __cxxabiv1 {
struct __cxa_eh_globals { void* caughtExceptions; unsigned int uncaughtExceptions; };
static __cxa_eh_globals eh_globals = {0, 0};
}
extern "C" {
__cxxabiv1::__cxa_eh_globals* __cxa_get_globals() { return &__cxxabiv1::eh_globals; }
__cxxabiv1::__cxa_eh_globals* __cxa_get_globals_fast() { return &__cxxabiv1::eh_globals; }
void* __cxa_current_exception_type() { return 0; }
}
STUBEOF
echo "  Compiling cxa_stubs.cpp"
em++ $CXXFLAGS -c cxa_stubs.cpp -o cxa_stubs.o
OBJS="$OBJS cxa_stubs.o"

# =============================================================================
# Step 5: Link into mash.js + mash.wasm
# =============================================================================

mkdir -p ../build
echo "=== Linking mash.js + mash.wasm ==="
em++ -O2 -fexceptions $OBJS \
    -L"${CAPNP_DIR}/lib" -lcapnp -lkj \
    $EM_FLAGS \
    -s DISABLE_EXCEPTION_CATCHING=0 \
    -o ../build/mash.js

echo "=== Build complete ==="
ls -lh ../build/mash.*
