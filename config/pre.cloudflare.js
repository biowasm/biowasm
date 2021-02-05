// Emscripten pre.js file for Cloudflare Workers environment
const document = this

// A way to retrieve stdout/stderr
Module.OUT = { stdout: "", stderr: "" };
Module.print = (function(d) { Module.OUT.stdout += d + "\n"; });
Module.printErr = (function(d) { Module.OUT.stderr += + "\n"; });
Module.getOutput = (function(d) { return Module.OUT; });
