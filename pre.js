// Fetch the corresponding .wasm and .data files relative to where the .js file is stored,
// whether it's locally or on a CDN. Works in the main thread and in WebWorkers.
Module['locateFile'] = (path, dir) => {
        var dirJS = "", dirRoot = "";
        if(typeof WorkerGlobalScope !== 'undefined' && self instanceof WorkerGlobalScope)
                dirJS = self.location.href;
        else if(document.currentScript)
                dirJS = document.currentScript.src;

        dirRoot = dirJS.match(/.*\//)[0];

        console.log(`path=${path}, dirRoot=${dirRoot}, dir=${dir}`)
        return dirRoot + path
}

