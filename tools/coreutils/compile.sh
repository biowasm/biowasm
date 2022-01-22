cd src/

# Install dependencies
sudo apt-get install -y autopoint gperf help2man gettext texinfo bison
./bootstrap

# Nanosleep not supported in Emscripten
sed -i 's|if ${gl_cv_func_sleep_works+:} false|if true|g' configure
sed -i 's|if ${ac_cv_search_nanosleep+:} false|if true|g' configure
sed -i 's|if ${gl_cv_func_nanosleep+:} false|if true|g' configure

# Disable assembly code
PREPROC="#if (defined __i386__ || defined __x86_64__) && defined __GNUC__"
PREPROC_ADD=" \&\& defined __NO_WAY_THIS_VAR_IS_DEFINED__"

# Configure
emconfigure ./configure \
  CC=emcc \
  FORCE_UNSAFE_CONFIGURE=1 \
  TIME_T_32_BIT_OK=yes \
  --host=wasm32

# This program needs gcc and should not be compiled to WebAssembly
sed -i 's|$(MAKE) src/make-prime-list$(EXEEXT)|gcc src/make-prime-list.c -o src/make-prime-list$(EXEEXT) -Ilib/|' Makefile

# Make all commands and skip "man" errors
# When update this list, need to update tools.json
emmake make all CC=emcc -k WERROR_CFLAGS=""

emmake make src/{hostname,basename,cat,chmod,comm,cp,cut,date,echo,env,fold,head,join,ls,md5sum,mkdir,mktemp,mv,nproc,paste,pwd,rm,rmdir,seq,shuf,sort,tail,tr,uniq,wc}.js \
  CC=emcc EXEEXT=.js \
  CFLAGS="-O2" \
  -k WERROR_CFLAGS=""

mv src/{hostname,basename,cat,chmod,comm,cp,cut,date,echo,env,fold,head,join,ls,md5sum,mkdir,mktemp,mv,nproc,paste,pwd,rm,rmdir,seq,shuf,sort,tail,tr,uniq,wc}.{js,wasm} ../build/
