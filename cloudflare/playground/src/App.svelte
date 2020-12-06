<script>
// Imports
import { onMount } from "svelte";
import popper from "popper.js";
import Bootstrap from "bootstrap";
import "bootstrap/dist/css/bootstrap.min.css";
import introJs from "intro.js";
import "intro.js/introjs.css";

import CommandLine from "./CommandLine.svelte";
import { TOOLS } from "./config.js";

// -----------------------------------------------------------------------------
// Globals
// -----------------------------------------------------------------------------

// State
let output = { stdout: "", stderr: "Loading..." };
let tool = new URL(window.location).searchParams.get("tool") || "samtools";

// Function from <CommandLine /> for running a command
let launch;

// -----------------------------------------------------------------------------
// On load
// -----------------------------------------------------------------------------

// Run default command on page load
onMount(() => loadTool());

// -----------------------------------------------------------------------------
// Utilities
// -----------------------------------------------------------------------------

async function loadTool(newTool=null) {
	// Launch new tool
	newTool = newTool || tool;
	tool = newTool;
	launch(`${tool} --help`);

	// Update URL
	window.history.pushState(tool, `biowasm playground - ${tool}`, `?tool=${tool}`);
}
</script>

<div class="jumbotron mt-5 mt-md-2 pb-1">
	<div class="container mt-2 mb-2">
		<div class="row">
			<div class="col-8">
				<h2>Playground</h2>
			</div>
			<div class="col-4 text-right" data-step="1" data-intro="These are the tools supported in this playground">
				{#each Object.keys(TOOLS).sort() as t}
					<button class="btn-tool mr-2 btn {t == tool ? "btn-secondary" : "btn-outline-secondary"}" on:click={() => loadTool(t)}>{t}</button>
				{/each}
			</div>
			<div class="col-12">
				<p class="text-secondary">
					Uses WebAssembly to run bioinformatics tools directly in your browser
					<button class="btn btn-sm btn-outline-secondary" on:click={() => introJs().start()}>Take a quick tour</button>
				</p>
			</div>
		</div>
	</div>
</div>

<div class="container">
	<CommandLine
		bind:launch={launch}
		on:output={msg => output = msg.detail} />

	<pre class="border rounded border-primary p-3" style="height:55vh"><span id="stdout">{output.stdout}</span><span id="stderr" style="color:#bbb">{output.stderr}</span></pre>
</div>
