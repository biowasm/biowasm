<script>
import { browser } from "$app/env";
import { Alert, ListGroup, ListGroupItem } from "sveltestrap";
import CodePen from "$components/CodePen.svelte";
import codeSimple from "$examples/simple.html?raw";
import codeUserFiles from "$examples/files-local.html?raw";
import codeRemoteFiles from "$examples/files-remote.html?raw";
import codeMultipleTools from "$examples/multiple-tools.html?raw";
</script>

<svelte:head>
	<title>Get Started</title>
</svelte:head>

<h4>Get Started</h4>

<div class="col-xl-8">
	<!-- Table of Contents -->
	<ListGroup flush class="small">
		<ListGroupItem></ListGroupItem>
		<ListGroupItem tag="a" href="#introduction">Introduction</ListGroupItem>
		<ListGroupItem tag="a" href="#simple-example">A simple example</ListGroupItem>
		<ListGroupItem tag="a" href="#local-files">Process local user files</ListGroupItem>
		<ListGroupItem tag="a" href="#remote-files">Process remote user files</ListGroupItem>
		<ListGroupItem tag="a" href="#multiple-tools">Run multiple tools</ListGroupItem>
		<ListGroupItem></ListGroupItem>
	</ListGroup>

	<!-- Introduction -->
	<h5 id="introduction">Introduction</h5>
	<p>
		The biowasm project compiles popular C/C++ genomics tools to <a href="https://developer.mozilla.org/en-US/docs/WebAssembly" target="_blank">WebAssembly</a> so they can run in a web browser. Biowasm includes a JavaScript library called <strong>Aioli</strong> that helps you run these tools with very little setup, as we'll see below.
	</p>

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
</div>

<style>
h5 {
	margin-top: 40px;
}

p {
	line-height: 1.7em;
}
</style>
