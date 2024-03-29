# Based on https://github.com/nika-begiashvili/libarchivejs/blob/master/lib/tools/docker/Dockerfile
# TODO: Run autoconf to generate configure files needed (but need Autoconf > 2.71)
curl -LO https://github.com/libarchive/libarchive/releases/download/v3.7.2/libarchive-3.7.2.zip
unzip libarchive-3.7.2.zip
cd libarchive-3.7.2

emconfigure ./configure --enable-static --disable-shared \
    --enable-bsdunzip=static --enable-largefile \
    --enable-bsdtar=static --enable-bsdcat=static --enable-bsdcpio=static \
    --enable-posix-regex-lib=libc \
    --disable-xattr --disable-acl \
    --without-xml2 --without-nettle --without-lzo2 --without-cng  --without-lz4 --without-expat \
    CFLAGS="-s USE_ZLIB=1 -s USE_BZIP2=1"

emmake make bsdunzip.js \
    EXEEXT=".js" \
    CFLAGS="-O2 -s USE_ZLIB=1 -s USE_BZIP2=1 -ffunction-sections -fdata-sections -fvisibility=hidden -D__LIBARCHIVE_ENABLE_VISIBILITY" \
    LDFLAGS="-O2 $EM_FLAGS"

mv bsdunzip.* ../../build/
