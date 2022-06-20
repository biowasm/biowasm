# Contributing

See [CONTRIBUTING.md](https://github.com/biowasm/biowasm/blob/main/CONTRIBUTING.md).

## Setup

Tools listed in biowasm were compiled to WebAssembly using `Emscripten 2.0.25`.

```bash
# Fetch Emscripten docker image
docker pull emscripten/emsdk:2.0.25

# Create the container and mount ~/wasm to /src in the container
docker run \
    -it -d \
    -p 80:80 \
    --name wasm \
    --volume ~/wasm:/src \
    emscripten/emsdk:2.0.25

# Go into the container
docker exec -u root -it wasm bash
# While inside the container, install dependencies
apt-get update
apt-get install -y autoconf liblzma-dev less vim
# Create small web server for testing
cat << EOF > server.py
import http.server
import socketserver

handler = http.server.SimpleHTTPRequestHandler
handler.extensions_map['.wasm'] = 'application/wasm'
httpd = socketserver.TCPServer(('', 80), handler)
httpd.serve_forever()
EOF
chmod +x server.py
# Launch the web server
python3.7 /src/server.py &
```


## Compile an existing biowasm tool

```bash
# Go into your container
docker exec -it wasm bash

# Set up biowasm (only need to do this once)
cd biowasm/
make init

# Compile seqtk
VERSION=1.2 BRANCH=v1.2 make seqtk

# This will create tools/<tool name>/build with .js/.wasm files
ls tools/seqtk/build
```


## Add a new tool to compile

First, add the tool as a git module:

```bash
TOOL=minimap2
BRANCH=v2.22
REPO=https://github.com/lh3/minimap2.git

# Fetch codebase
mkdir -p tools/$TOOL
git submodule add $REPO tools/$TOOL/src

# Get specific version of the tool
cd tools/$TOOL/src
git checkout $BRANCH
git submodule update --init --recursive
cd -

# Stage changes for git
git add tools/$TOOL/src .gitmodules

# Other files needed
echo "TODO" > tools/$TOOL/README.md
echo "# TODO" > tools/$TOOL/compile.sh
chmod +x tools/$TOOL/compile.sh
```

You should also create the following files:

```bash
tools/<tool>/
    README.md        Details about the tool and dependencies
    compile.sh       Script that will run to compile the tool to WebAssembly (can use `$EM_FLAGS` for common flags)
    patches/    
        <tag>.patch  Patch applied to the code to compile it to WebAssembly; branch- or tag-specific (optional)
    configs/
        <tag>.json   Configuration file with info about which WebAssembly features are needed (see ssw for an example); branch- or tag-specific (optional)
```

Finally, you can edit:

* `config/tools.json` to make sure the new tool gets deployed
* `cloudflare/cdn/public/index.html` to list the new tool on the CDN's [packages page](https://cdn.biowasm.com/v2/)


## Future candidates

- [Sniffles](https://github.com/fritzsedlazeck/Sniffles)
- [FreeBayes](https://github.com/freebayes/freebayes)
- [HISAT2](https://github.com/DaehwanKimLab/hisat2)
