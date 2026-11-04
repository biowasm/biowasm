// biowasm SIMD shim v1 (bin/simd.js)
(() => {
	const ModuleNoSIMD = Module;
	const simd = WebAssembly.validate(new Uint8Array([0,97,115,109,1,0,0,0,1,5,1,96,0,1,123,3,2,1,0,10,10,1,8,0,65,0,253,15,253,98,11]));
	Module = (args = {}) => {
		const locateFile = args.locateFile;
		if(simd && locateFile && args.simd !== false) args = { ...args, locateFile: (path, prefix) => locateFile(path.replace(/\.wasm$/, "-simd.wasm"), prefix) };
		return ModuleNoSIMD(args);
	};
})();
