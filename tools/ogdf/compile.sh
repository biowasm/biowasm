# Dependencies
apt-get install doxygen -y

# Compile
emcmake cmake .
emmake make

# CMakeFiles/ex-energybased.dir/doc/examples/layout/energybased.cpp.o
emmake make doc/examples/layout/energybased.o
em++ -o ../examples/energybased.js \
    -L./ CMakeFiles/ex-energybased.dir/doc/examples/layout/energybased.cpp.o \
    -lOGDF -lCOIN $EM_FLAGS

    # --preload-file doc/examples/layout/sierpinski_04.gml@/ogdf/examples/sierpinski_04.gml
    # --preload-file doc/examples/layout/ERDiagram.gml@/ogdf/examples/ERDiagram.gml \
    # --preload-file doc/examples/layout/uk_Pack_Bary_EC_FRENC.gml@/ogdf/examples/uk_Pack_Bary_EC_FRENC.gml \
    # --preload-file doc/examples/layout/unix-history-time.gml@/ogdf/examples/unix-history-time.gml \
    # --preload-file doc/examples/layout/unix-history.gml@/ogdf/examples/unix-history.gml
