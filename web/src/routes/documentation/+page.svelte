<script>
import { browser } from "$app/environment";
import { page } from "$app/stores";
import { Alert, Icon } from "@sveltestrap/sveltestrap";
import CONFIG from "@/biowasm.json";
import Playground from "$components/Playground.svelte";
import CodeBlock from "$components/CodeBlock.svelte";

// Sample code
import codeSimple from "$examples/simple.html?raw";
import codeUserFiles from "$examples/files-local.html?raw";
import codeRemoteFiles from "$examples/files-remote.html?raw";
import codeMultipleTools from "$examples/multiple-tools.html?raw";

// Sidebar navigation: these are anchors into the sections below
const NAV = [
	{ href: "#introduction", label: "Introduction", icon: "book" },
	{ href: "#examples", label: "Code Examples", icon: "code-slash", children: [
		{ href: "#simple-example", label: "A simple example" },
		{ href: "#local-files", label: "Process local user files" },
		{ href: "#remote-files", label: "Process remote user files" },
		{ href: "#multiple-tools", label: "Run multiple tools" }
	] },
	{ href: "#aioli", label: "API documentation", icon: "braces", children: [
		{ href: "#aioli-init", label: "Initialization" },
		{ href: "#aioli-api", label: "Utility functions" }
	] },
	{ href: "#no-cdn", label: "Biowasm without a CDN", icon: "cloud-slash", children: [
		{ href: "#no-cdn-install", label: "Install the Aioli package" },
		{ href: "#no-cdn-assets", label: "Download Biowasm assets" },
		{ href: "#no-cdn-setup", label: "Set up Aioli" }
	] }
];

// Highlight the section we're in; the top of the page counts as the first one
$: activeHash = $page.url.hash || NAV[0].href;
</script>

<svelte:head>
	<title>Documentation | Biowasm</title>
</svelte:head>

<div class="docs-layout mb-5">
	<aside class="docs-nav">
		<div class="docs-nav-inner">
			<div class="docs-nav-title text-uppercase text-muted">Documentation</div>
			<nav>
				{#each NAV as item (item.href)}
					<a href={item.href} class="docs-nav-link" class:active={activeHash === item.href}>
						<Icon name={item.icon} />
						<span>{item.label}</span>
					</a>

					{#each item.children || [] as child (child.href)}
						<a href={child.href} class="docs-nav-link docs-nav-sublink" class:active={activeHash === child.href}>
							<span>{child.label}</span>
						</a>
					{/each}
				{/each}
			</nav>
		</div>
	</aside>

	<main class="docs-main">
	<!-- ================================================================================ -->
	<!-- Introduction -->
	<!-- ================================================================================ -->
	<h4 id="introduction">Introduction</h4>
	<p>
		The biowasm project compiles popular C/C++ genomics tools to <a href="https://developer.mozilla.org/en-US/docs/WebAssembly" target="_blank">WebAssembly</a> so they can run in a web browser. Biowasm includes a JavaScript library called <strong>Aioli</strong> that helps you run these tools with very little setup, as we'll see below.
	</p>


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
		<Playground code={codeSimple} /><br />
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
		<Playground code={codeUserFiles} /><br />
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
		<Playground code={codeRemoteFiles} /><br />
	{/if}

	<!-- Multiple tools -->
	<h5 id="multiple-tools">Run multiple tools</h5>
	<p>
		Biowasm also supports running multiple tools at once. Below, we first run <code>seqtk</code> on the output of a <code>samtools</code> command:
	</p>
	{#if browser}
		<Playground code={codeMultipleTools} /><br />
	{/if}


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

		// Basic file system utilities
		await CLI.mkdir("/some/path");
		await CLI.ls("/some/path");
		await CLI.cd("/some/path");
		await CLI.pwd();

		// Call Emscripten's virtual file system utilities (docs: https://emscripten.org/docs/api_reference/Filesystem-API.html)
		await CLI.fs.unlink("/some/path/here"); // Delete a file
		await CLI.fs.stat("/some/path/here"); // Get file size/timestamps
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


	<!-- No CDN -->
	<h4 id="no-cdn">Biowasm without the CDN</h4>
	<p>
		The examples above use the biowasm CDN, which is where the WebAssembly assets for the genomics tools are hosted for free. Most tools that use biowasm rely on the CDN, but if you would like to store your assets alongside your app, here's how to do so:
	</p>

	<h6 id="no-cdn-install">Step 1: Install the Aioli package</h6>
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

	<h6 id="no-cdn-assets">Step 2: Download Biowasm assets</h6>
	<p>
		Note that even if you import Aioli locally with <code>npm</code>, the WebAssembly modules will still be downloaded from the biowasm CDN unless you download those assets locally as well.
	</p>
	<p>
		To do so, navigate to the <a href={CONFIG.url}>Packages page</a> and download the files for each tool of interest.
	</p>

	<h6 id="no-cdn-setup">Step 3: Set up Aioli</h6>
	<p>
		Once your assets are downloaded, you need to let Aioli know the relative URL path where it can find those assets from your server. Use the config parameter <code>urlPrefix</code> as shown in the <a href="#aioli-init-advanced:~:text=urlPrefix%3A%20%22./assets/%22%2C">Initialization</a> section above.
	</p>

	</main>
</div>

<style>
.docs-layout {
	display: flex;
	gap: 2.5rem;
	max-width: 1120px;
	align-items: flex-start;
}

.docs-nav {
	flex: 0 0 260px;
	width: 260px;
	/* Span the full row, otherwise the sticky child has no room to travel */
	align-self: stretch;
}

.docs-nav-inner {
	position: sticky;
	top: 1.5rem;
	max-height: calc(100vh - 3rem);
	overflow-y: auto;
}

.docs-nav-title {
	font-size: 0.72rem;
	font-weight: 700;
	letter-spacing: 0.06em;
	margin-bottom: 0.75rem;
	padding-left: 0.75rem;
}

.docs-nav-link {
	display: flex;
	align-items: center;
	gap: 0.6rem;
	padding: 0.45rem 0.75rem;
	border-radius: 6px;
	color: var(--bs-body-color);
	text-decoration: none;
	font-size: 0.95rem;
	font-weight: 600;
	line-height: 1.3;
}

.docs-nav-sublink {
	/* 0.75rem padding + 1em icon + 0.6rem gap, so text lines up with the headings */
	padding-left: 2.3rem;
	padding-top: 0.3rem;
	padding-bottom: 0.3rem;
	font-size: 0.875rem;
	font-weight: 400;
	color: var(--bs-secondary-color);
}

.docs-nav-link:hover {
	background: var(--bs-tertiary-bg);
	color: var(--bs-primary-text-emphasis);
}

.docs-nav-link.active {
	color: var(--bs-primary-text-emphasis);
	background: var(--bs-primary-bg-subtle);
	font-weight: 600;
}

.docs-main {
	flex: 1 1 auto;
	min-width: 0;
}

/* Stack the sidebar on top on narrow screens */
@media (max-width: 767.98px) {
	.docs-layout {
		flex-direction: column;
		gap: 1rem;
	}

	.docs-nav {
		flex-basis: auto;
		width: 100%;
	}

	.docs-nav-inner {
		position: static;
	}

	.docs-nav nav {
		display: flex;
		flex-wrap: wrap;
		gap: 0.25rem;
	}

	/* Top-level sections only, so the stacked nav stays short */
	.docs-nav-sublink {
		display: none;
	}
}

/* The navbar scrolls away, so the sticky sidebar shifts up by its height the moment it
   pins. Animating the jump to an anchor makes that glide instead of snapping. */
:global(html) {
	scroll-behavior: smooth;
}

@media (prefers-reduced-motion: reduce) {
	:global(html) {
		scroll-behavior: auto;
	}
}

/* Stop anchored sections flush against the viewport edge, matching the sidebar's offset */
h4, h5, h6 {
	scroll-margin-top: 1.5rem;
}

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
