// Cloudflare Worker Script.
// Must specify `wasm_modules = { WASM_MODULE = "./wgsim.wasm" }` in wrangler.toml.
// See <https://github.com/robertaboukhalil/cf-workers-emscripten/> for an example.

import Module from "wgsim.js";

addEventListener("fetch", event => {
	event.respondWith(handleRequest(event.request));
});

async function handleRequest(request) {
	// Initialize WebAssembly module
	let output = "";
	const m = await Module({
		// By default, stdout/stderr is output to console.log/warn
		print: text => output += `${text}\n`,
		printErr: text => output += `${text}\n`,

		// Instead of downloading the .wasm file, fetch it from a global var
		instantiateWasm: (imports, callback) => {
			// Note that "PI_WASM" is defined in wrangler.toml
			const instance = new WebAssembly.Instance(PI_WASM, imports);
			callback(instance);
			return instance.exports;
		}
	});

	// Call main function
	m.callMain([ "--help" ]);

	// Return output
	return new Response(output);
}
