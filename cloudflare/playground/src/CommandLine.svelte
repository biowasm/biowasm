<script>
// TODO:
// - Have Aioli running commands in here and just pushing result to parent component?


// Exports
export let command = "";        // Command to execute (e.g. samtools --version)
export let disabled = false;    // Whether to disable the input or not

// Imports
import { afterUpdate, createEventDispatcher } from "svelte";
import { Aioli } from "@biowasm/aioli";


// -----------------------------------------------------------------------------
// Globals
// -----------------------------------------------------------------------------

// Constants
const TOOLS = {
    "samtools": { module: "samtools", version: "1.10" },
    "bedtools2": { module: "bedtools2", version: "2.29.2" }
};

// 
let dispatch = createEventDispatcher();  // used to communicate with parent component

// User input
let program = null;                      // e.g. bedtools
let args = null;                         // e.g. intersect -a a.bed -b b.bed

// UI
let msgInfo = "Supported: " + Object.values(TOOLS).map(d => `<code>${d.module}</code>`).join(", ");
let msgError = "";
let elTextbox = null;                    // DOM element for the CLI input


// -----------------------------------------------------------------------------
// Reactive Statements
// -----------------------------------------------------------------------------

// Split program name from args
$: program = command.split(" ").shift();
$: args = command.replace(`${program} `, "").trim();


// -----------------------------------------------------------------------------
// Utility functions
// -----------------------------------------------------------------------------

// Execute a command
function run()
{
	// Send a message to parent component about what to execute
	disabled = true;

	// dispatch("execute", {
	// 	program: program,
	// 	args: args,
	// 	done: () => disabled = false
	// });
}

// Focus on command line once the DOM settles
afterUpdate(() => elTextbox.focus());


// -----------------------------------------------------------------------------
// HTML
// -----------------------------------------------------------------------------
</script>

<style>
input {
	font-family: monospace;
}
</style>

<!-- Info message -->
<div class="row">
	<div class="col-12">
		<span class="text-muted">
			<span class="text-info">{@html msgInfo}&nbsp;</span>
		</span>
	</div>
</div>

<!-- CLI -->
<div class="row mt-2">
	<div class="col-12">
		<div class="input-group">
			<div class="input-group-prepend">
				<span class="input-group-text">$</span>
			</div>
			<input
				type="text"
				class="form-control form-control-lg"
				disabled={disabled}
				bind:this={elTextbox}
				bind:value={command}
				on:keydown={event => event.key == "Enter" ? run() : null}
			/>
			<div class="input-group-append">
				<button
					class="btn btn-lg {disabled ? 'btn-secondary' : 'btn-info'}"
					disabled={disabled}
					on:click={run}
				>
					Run
				</button>
			</div>
		</div>
	</div>
</div>

<!-- Error message -->
<div class="row">
	<div class="col-12">
		<small class="text-muted">
			<span class="text-danger">{@html msgError}&nbsp;</span>
		</small>
	</div>
</div>
