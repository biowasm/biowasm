#!/bin/bash

# Remove extra space in front of Emscripten "-s" variables to avoid errors
emcmake cmake -DCMAKE_EXE_LINKER_FLAGS="-sTOTAL_STACK=2097152 -02 -sASSERTIONS=1 -sMODULARIZE=1 -sALLOW_MEMORY_GROWTH -sFORCE_FILESYSTEM=1 -sEXIT_RUNTIME=0 -s EXPORTED_RUNTIME_METHODS=["callMain","FS","PROXYFS","WORKERFS","UTF8ToString","getValue","AsciiToString"] -lworkerfs.js -lproxyfs.js -s INVOKE_RUN=0 -s ENVIRONMENT="web,worker" ${EM_FLAGS//-s /-s} -fwasm-exceptions --preload-file res@/hyphy --preload-file tests/hbltests@/tests"
emmake make hyphy
emmake make install
mv hyphy.{js,wasm,data} ../build/
