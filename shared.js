// Fetch the corresponding .wasm and .data files relative to where the .js file is stored,
// whether it's locally or on a CDN. Works in the main thread and in WebWorkers.
var Module = typeof Module !== 'undefined' ? Module : {};

Module['locateFile'] = function(path, dir)
{
	var dirRoot = "";

	// Use hardcoded path for all files
	if(typeof BIOWASM_URL !== 'undefined')
		dirRoot = BIOWASM_URL;
	
	// Or infer it from the path to the JS file
	else
	{
		var dirJS = "";
		// Inside a WebWorker
		if(typeof WorkerGlobalScope !== 'undefined' && self instanceof WorkerGlobalScope)
			dirJS = self.location.href;
		// In the main thread
		else if(document.currentScript)
			dirJS = document.currentScript.src;
	
		dirRoot = dirJS.substring(0, dirJS.lastIndexOf("/") + 1);
	}

	return dirRoot + path;
}
