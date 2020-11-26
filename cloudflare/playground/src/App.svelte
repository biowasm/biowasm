<script>
// Imports
import { onMount } from "svelte";
import popper from "popper.js";
import Bootstrap from "bootstrap";
import "bootstrap/dist/css/bootstrap.min.css";
import CommandLine from "./CommandLine.svelte";

// State
let output = "";
let tool = new URL(window.location).searchParams.get("tool") || "samtools";
let run = null;

// -----------------------------------------------------------------------------
// On load
// -----------------------------------------------------------------------------

onMount(() => run());

// -----------------------------------------------------------------------------
// Utilities
// -----------------------------------------------------------------------------

function loadTool(newTool) {
	tool = newTool;
	console.log(run);
	run(tool);
}
</script>

<div class="jumbotron mt-5 mt-md-2 pb-1">
	<div class="container mt-4 mb-4">
		<h2>Playground</h2>
		<p class="lead">
			Launch: 
			<button class="btn btn-sm btn-outline-secondary" on:click={() => loadTool("samtools")}>samtools</button>
			<button class="btn btn-sm btn-outline-secondary" on:click={() => loadTool("bedtools2")}>bedtools2</button>
		</p>
	</div>
</div>

<div class="container">
	<CommandLine
		bind:run={run}
		command={`${tool} --version`}
		on:output={msg => output = msg.detail.stdout.trim() + msg.detail.stderr.trim()} />

	<pre class="border rounded border-primary p-3" style="height:55vh">{output}</pre>
</div>
