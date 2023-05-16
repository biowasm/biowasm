# Dependencies
apt-get install doxygen -y

# Compile
emcmake cmake .
emmake make

# CMakeFiles/ex-energybased.dir/doc/examples/layout/energybased.cpp.o
emmake make doc/examples/layout/energybased.o
em++ -o deleteme.js -L./ CMakeFiles/ex-energybased.dir/doc/examples/layout/energybased.cpp.o -lOGDF -lCOIN $EM_FLAGS
