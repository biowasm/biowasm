<script>
// Imports
import { onMount } from "svelte";
import popper from "popper.js";
import Bootstrap from "bootstrap";
import "bootstrap/dist/css/bootstrap.min.css";

import CommandLine from "./CommandLine.svelte";
import { TOOLS } from "./config.js";

// -----------------------------------------------------------------------------
// Globals
// -----------------------------------------------------------------------------

// State
let output = "";
let tool = new URL(window.location).searchParams.get("tool") || "samtools";

// Function from <CommandLine /> for running current command
let run = () => {};

// -----------------------------------------------------------------------------
// On load
// -----------------------------------------------------------------------------

// Run default command on page load
onMount(() => run());

// -----------------------------------------------------------------------------
// Utilities
// -----------------------------------------------------------------------------

async function loadTool(newTool) {
	// Launch new tool
	tool = newTool;
	run();

	// Update URL
	window.history.pushState(tool, `biowasm playground - ${tool}`, `?tool=${tool}`);
}
</script>

<div class="jumbotron mt-5 mt-md-2 pb-1">
	<div class="mt-4 mb-4 row">
		<div class="col-6">
			<h2>Playground</h2>
		</div>
		<div class="col-6 text-right">
			{#each Object.keys(TOOLS) as t}
				<button class="mr-2 btn {t == tool ? "btn-secondary" : "btn-outline-secondary"}" on:click={() => loadTool(t)}>{t}</button>
			{/each}
		</div>
	</div>
</div>

<div class="container">
	<CommandLine
		bind:run={run}
		command={`${tool} --version`}
		on:output={msg => output = msg.detail.stdout.trim() + msg.detail.stderr.trim()} />

	<pre class="border rounded border-primary p-3" style="height:55vh">{output}</pre>
</div>
