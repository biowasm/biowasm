
# Support the old 3.8 version of Muscle where the code is not on GitHub.
# New versions point to a tag corresponding to a version; old version points to main
BRANCH_OR_TAG=$(git symbolic-ref -q --short HEAD || git describe --tags --exact-match)

# New versions
if [[ "$BRANCH_OR_TAG" != "main" ]]; then
    emmake make CC=emcc CXX=em++
# Muscle v3
else
    cd ..  # go to parent of src/ to make a new folder there
    wget "https://drive5.com/muscle/muscle_src_3.8.1551.tar.gz"
    mkdir src-v3
    tar -xf "muscle_src_3.8.1551.tar.gz" -C src-v3
    rm "muscle_src_3.8.1551.tar.gz"
    cd src-v3
    # Avoid error "sed: couldn't open temporary file ./sed.....: Permission denied"
    sed -r -e 's/strip muscle/#strip muscle/' -e 's/ -o muscle / -o ..\/build\/muscle.js \$\(EM_FLAGS\) /' Makefile > Makefile.tmp
    mv Makefile.tmp Makefile
    emmake make GPP=em++
    cd ../src
fi
