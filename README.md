# biowasm

![cdn-stg.biowasm.com](https://github.com/biowasm/biowasm/workflows/Deploy%20biowasm-stg/badge.svg) ![cdn.biowasm.com](https://github.com/biowasm/biowasm/workflows/Deploy%20biowasm-prd/badge.svg)

A repository of genomics tools, compiled from C/C++ to WebAssembly so they can run in a web browser:

* [samtools/htslib v1.10](tools/samtools/README.md)
* [bedtools v2.29](tools/bedtools/README.md)
* [fastp v0.20.1](tools/fastp/README.md)
* [seqtk v1.3](tools/seqtk/README.md)
* [wgsim](tools/wgsim/README.md)
* [bhtsne](tools/bhtsne/README.md)
* [seq-align](tools/seq-align/README.md)


## Get Started

biowasm modules are hosted on [cdn.biowasm.com](https://cdn.biowasm.com/index).

### Simple usage


```html
<script src="https://cdn.biowasm.com/aioli/latest/aioli.js"></script>
<script>

</script>
```


### Using Aioli


```html
<input id="myfile" type="file" multiple>
<script src="https://cdn.sandbox.bio/aioli/latest/aioli.js"></script>

<script>
let samtools = new Aioli("samtools/1.10");

// Initialize samtools and output the version
samtools
    .init()
    .then(() => samtools.exec("--version"))
    .then(d => console.log(d.stdout));

// When a user selects a .sam file from their computer,
// run `samtools view -q20` on the file
function loadFile(event)
{
    Aioli
        // First mount the file
        .mount(event.target.files[0])
        // Once it's mounted, run samtools view
        .then(file => samtools.exec(`view -q20 ${file.path}`))
        // Capture output
        .then(d => console.log(d.stdout));
}
document.getElementById("myfile").addEventListener("change", loadFile, false);
</script>
```


For a simple starter example, see [Aioli](https://github.com/biowasm/aioli#getting-started).



## Development

### Setup

Tools listed in biowasm were compiled to WebAssembly using `Emscripten 2.0.0`.

```bash
# Fetch Emscripten docker image
docker pull emscripten/emsdk:2.0.0

# Create the container and mount ~/wasm to /src in the container
docker run \
    -it -d \
    -p 80:80 \
    --name wasm \
    --volume ~/wasm:/src \
    emscripten/emsdk:2.0.0

docker exec -u root -it bash
apt-get install -y autoconf liblzma-dev less vim
cat << EOF > server.py
import http.server
import socketserver

handler = http.server.SimpleHTTPRequestHandler
handler.extensions_map['.wasm'] = 'application/wasm'
httpd = socketserver.TCPServer(('', 80), handler)
httpd.serve_forever()
EOF
chmod +x server.py
python3.7 /src/server.py &
```


### Compile a tool

```bash
# Go into your container
docker exec -it wasm bash

# Compile seqtk
cd biowasm/
make init  # only need to do this once
make seqtk

# This will create tools/<tool name>/build with .js/.wasm files
ls tools/seqtk/build
```


### Add a new tool

First, add the tool as a git module:

```bash
# Fetch codebase
mkdir -p tools/seqtk
git submodule add https://github.com/lh3/seqtk.git tools/seqtk/src

# Get specific version of the tool
cd tools/seqtk/src
git checkout v1.3
cd -

# Stage changes for git
git add tools/seqtk/src .gitmodules
```

You should also create the following files:

```bash
tools/<tool>/
    README.md   Details about the tool and dependencies
    compile.sh  Script that will run to compile the tool to WebAssembly (can use `$EM_FLAGS` for common flags)
    patch       Patches that need to be applied to the code to compile it to WebAssembly (optional)
```

## Todo

## To do

- Run each tool's tests?
- Support for Rust bioinformatics tools such as [sourmash](https://github.com/dib-lab/sourmash/tree/v3.2.2/src/core) and [rust-bio](https://github.com/rust-bio/rust-bio)
