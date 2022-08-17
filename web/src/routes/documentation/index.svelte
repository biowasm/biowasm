<script>
import { browser } from "$app/env";
import { Alert, ListGroup, ListGroupItem } from "sveltestrap";
import CONFIG from "@/biowasm.json";
import CodePen from "$components/CodePen.svelte";
import CodeBlock from "$components/CodeBlock.svelte";

// Sample code
import codeSimple from "$examples/simple.html?raw";
import codeUserFiles from "$examples/files-local.html?raw";
import codeRemoteFiles from "$examples/files-remote.html?raw";
import codeMultipleTools from "$examples/multiple-tools.html?raw";
</script>

<svelte:head>
	<title>Documentation</title>
</svelte:head>

<h4>Documentation</h4>

<div class="col-xl-8">
	<!-- ================================================================================ -->
	<!-- Table of Contents -->
	<!-- ================================================================================ -->
	<ListGroup flush class="small">
		<ListGroupItem></ListGroupItem>
		<ListGroupItem tag="a" href="#introduction">Introduction</ListGroupItem>
		<ListGroupItem tag="a" href="#examples">Code Examples</ListGroupItem>
		<ListGroupItem tag="a" href="#aioli">API documentation</ListGroupItem>
		<ListGroupItem tag="a" href="#no-cdn">Biowasm without a CDN</ListGroupItem>
		<ListGroupItem></ListGroupItem>
	</ListGroup>

	<!-- ================================================================================ -->
	<!-- Introduction -->
	<!-- ================================================================================ -->
	<h4 id="introduction">Introduction</h4>
	<p>
		The biowasm project compiles popular C/C++ genomics tools to <a href="https://developer.mozilla.org/en-US/docs/WebAssembly" target="_blank">WebAssembly</a> so they can run in a web browser. Biowasm includes a JavaScript library called <strong>Aioli</strong> that helps you run these tools with very little setup, as we'll see below.
	</p>

	<br />

	<!-- ================================================================================ -->
	<!-- Code examples -->
	<!-- ================================================================================ -->
	<h4 id="examples">Code Examples</h4>

	<!-- Simple example -->
	<h5 id="simple-example">A simple example</h5>
	<p>
		Say you're building a web app that analyzes sequencing data by parsing user-provided <a href="https://en.wikipedia.org/wiki/Binary_Alignment_Map" target="_blank">BAM files</a>. On the command-line, you would use the program <a href="http://www.htslib.org/" target="_blank">samtools</a> to parse these files; with biowasm, you can write a few lines of code and run <code>samtools</code> directly in the browser:
	</p>
	{#if browser}
		<CodePen code={codeSimple} /><br />
	{/if}
	<p>
		Note that we ran <code>samtools view</code> on a sample file <code>/samtools/examples/toy.sam</code> that comes preloaded with biowasm for testing purposes.
	</p>

	<!-- Process user files -->
	<h5 id="local-files">Process local user files</h5>
	<p>
		To run <code>samtools</code> on a user-provided file (download one <a href="https://raw.githubusercontent.com/samtools/samtools/develop/examples/toy.sam" target="_blank">here</a>), we need to mount the user's file onto a virtual file system before we can use it within Aioli:
	</p>
	<Alert color="info">
		<strong>Note:</strong> Mounting a file does <strong>not</strong> load its contents in memory! For example, selecting a 500GB BAM file below will work just fine since <code>samtools view -H</code> only parses the header of the file.<br /><br />One common use case for biowasm is to efficiently process small subsets of large files in the browser, e.g. parse the very top/bottom of a file, or randomly subsample a large file.
	</Alert>
	{#if browser}
		<CodePen code={codeUserFiles} /><br />
	{/if}

	<!-- Process remote files -->
	<h5 id="remote-files">Process remote user files</h5>
	<p>
		You can even mount URLs as long as they are <a href="https://developer.mozilla.org/en-US/docs/Web/HTTP/CORS" target="_blank">CORS-enabled</a>:
	</p>
	<Alert color="info">
		<strong>Note:</strong> Similar to local files, mounting URLs does not load their entire contents into memory. Below, we mount a 300MB BAM file and ask <code>samtools</code> for a small subset of chromosome 20. As you can see in your browser's <code>Network</code> tab, only a few megabytes of relevant data are downloaded for this step.
	</Alert>
	{#if browser}
		<CodePen code={codeRemoteFiles} /><br />
	{/if}

	<!-- Multiple tools -->
	<h5 id="multiple-tools">Run multiple tools</h5>
	<p>
		Biowasm also supports running multiple tools at once. Below, we first run <code>seqtk</code> on the output of a <code>samtools</code> command:
	</p>
	{#if browser}
		<CodePen code={codeMultipleTools} /><br />
	{/if}

	<br />

	<!-- ================================================================================ -->
	<!-- Aioli API -->
	<!-- ================================================================================ -->
	<h4 id="aioli">Aioli API</h4>

	<!-- Initialization -->
	<h5 id="aioli-init">Initialization</h5>
	<h6 id="aioli-init-simple">Simple Initialization</h6>
	<p>
		Most tools can be initialized using a name and a version (see the <a href={CONFIG.url}>Packages</a> page for supported versions):
	</p>

	<CodeBlock lang="javascript" code={`
		new Aioli(["samtools/1.10", "seqtk/1.3"]);
	`} />
	<p>
		For tools that include multiple sub-tools, you need to specify which ones you want to load:
	</p>
	<CodeBlock lang="javascript" code={`
		new Aioli(["coreutils/head/8.32", "coreutils/tail/8.32"]);
	`} />

	<h6 id="aioli-init-advanced">Advanced Initialization</h6>
	<p>
		If you need to customize the URL where the WebAssembly assets are stored, or the lazy-loading behavior of a module, you can use the following initialization format:
	</p>
	<CodeBlock lang="javascript" code={`
		new Aioli([{
		    tool: "coreutils",
		    version: "8.32",
		    program: "head",         // Optional: sub-tool name; not needed for most tools (default: same as tool name)
		    urlPrefix: "./assets/",  // Optional: custom path to .wasm assets (default: biowasm CDN)
		    loading: "lazy",         // Optional: if set to "lazy", only downloads WebAssembly modules when needed, not at initialization (default: eager)
		    reinit: false,           // Optional: if set to true, will reinitialize module after each invocation; not needed for most tools
		}], {
		    printInterleaved: true,  // Optional: whether to return interleaved stdout/stderr; if false, returns object with stdout/stderr keys (default: true)
		    debug: false,            // Optional: set to true to see console log messages for debugging (default: false)
		});
	`} />

	<h5 id="aioli-api">Utility functions</h5>
	<h6 id="aioli-api-exec">Run a command</h6>
	<CodeBlock lang="javascript" code={`
		// Simple command
		await CLI.exec("samtools view");

		// To escape arguments, use an array instead of a string
		await CLI.exec("jq", [
		    ".some.data",
		    "data.json"
		]);
	`} />

	<h6 id="aioli-api-mount">Mount local and remote files</h6>
	<p>See examples above for how to mount <a href="#local-files">local files</a> and <a href="#remote-files">remote URLs</a>.</p>
	<p>The function <code>CLI.mount()</code> can mount one (or an array of) <code>File</code> objects, <code>Blob</code> objects, a <code>FileList</code>, strings, or URLs:</p>
	<CodeBlock lang="javascript" code={`
		// Mount a list of user-provided File objects
		document.getElementById("myfiles").addEventListener("change", async event => {
		    await CLI.mount(event.target.files);
		});

		// Mount URLs
		await CLI.mount([
		    { name: "filename.txt", url: "https://url..." },
		    { name: "filename2.txt", url: "https://url2..." },
		]);

		// Mount a string to path filename.txt
		await CLI.mount([{
		    name: "filename.txt",
		    data: "This is a data file"
		}]);

		// Mount a Blob to path filename.txt
		await CLI.mount([{
		    name: "filename.txt",
		    data: <Blob>
		}]);
	`} />

	<h6 id="aioli-api-fs">File system utilities</h6>
	<CodeBlock lang="javascript" code={`
		// Returns a blob URL so the user can download a file out of the virtual file system
		const url = await CLI.download("/path/to/a/file");

		// Basic ls, cd, mkdir, pwd utilities
		await CLI.mkdir("/some/path");
		await CLI.ls("/some/path");
		await CLI.cd("/some/path");
		await CLI.pwd();
	`} />

	<Alert color="info">
		<strong>Note:</strong> File system utilities operate on a virtual file system; they cannot operate directly on user files.
	</Alert>

	<h6 id="aioli-api-stdin">Standard IO</h6>
	<p>If a tool needs to interact with <code>stdin</code>, set <code>CLI.stdin</code> before making the call to the tool:</p>
	<CodeBlock lang="javascript" code={`
		CLI = await new Aioli(["coreutils/cat/8.32"]);
		CLI.stdin = "Hello World";
		await CLI.exec("cat");  // Will return Hello World
	`} />

	<br />

	<!-- No CDN -->
	<h4 id="no-cdn">Biowasm without the CDN</h4>
	<p>
		The examples above use the biowasm CDN, which is where the WebAssembly assets for the genomics tools are hosted for free. Most tools that use biowasm rely on the CDN, but if you would like to store your assets alongside your app, here's how to do so:
	</p>

	<h6>Step 1: Install the Aioli package</h6>
	<p>
		Instead of importing Aioli using
	</p>
	<CodeBlock lang="xml" code={`
		<script src="https://biowasm.com${CONFIG.url}/aioli.js"></script>
	`} />
	<p>
		you can install Aioli with npm:
	</p>
	<CodeBlock code={`
		npm install --save "@biowasm/aioli"
	`} />
	<p>
		Then you can import Aioli as follows:
	</p>
	<CodeBlock lang="javascript" code={`
		import Aioli from "@biowasm/aioli";
	`} />

	<h6>Step 2: Download Biowasm assets</h6>
	<p>
		Note that even if you import Aioli locally with <code>npm</code>, the WebAssembly modules will still be downloaded from the biowasm CDN unless you download those assets locally as well.
	</p>
	<p>
		To do so, navigate to the <a href={CONFIG.url}>Packages page</a> and download the files for each tool of interest.
	</p>

	<h6>Step 3: Set up Aioli</h6>
	<p>
		Once your assets are downloaded, you need to let Aioli know the relative URL path where it can find those assets from your server. Use the config parameter <code>urlPrefix</code> as shown in the <a href="#aioli-init-advanced:~:text=urlPrefix%3A%20%22./assets/%22%2C">Initialization</a> section above.
	</p>

	<br /><br /><br />
</div>

<style>
h4:not(:first-child) {
	margin-top: 40px;
}

h5, h6 {
	margin-top: 20px;
}

p {
	line-height: 1.7em;
}
</style>
