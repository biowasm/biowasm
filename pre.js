// Fetch the corresponding .wasm and .data files relative to where the .js file is stored,
// whether it's locally or on a CDN. Works in the main thread and in WebWorkers.
Module['locateFile'] = function(path, dir)
{
	var dirJS = "", dirRoot = "";
	// Inside a WebWorker
	if(typeof WorkerGlobalScope !== 'undefined' && self instanceof WorkerGlobalScope)
		dirJS = self.location.href;
	// Inside the main thread
	else if(document.currentScript)
		dirJS = document.currentScript.src;

	// Get folder where main .js file is located
	dirRoot = dirJS.match(/.*\//)[0];
	return dirRoot + path;
}
