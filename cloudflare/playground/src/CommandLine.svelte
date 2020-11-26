<script>
// Exports
export let command = "";        // Command to execute (e.g. samtools --version)
export let execute = false;     // Set this to true when changing the command to execute it on load
export let disabled = false;    // Whether to disable the input or not

// Imports
import jQuery from "jquery";
import { onMount, afterUpdate, createEventDispatcher } from "svelte";
import { Aioli } from "@biowasm/aioli";
import { TOOLS } from "./config.js";


// -----------------------------------------------------------------------------
// Globals
// -----------------------------------------------------------------------------

// State
let aiolis = {};         // e.g. { samtools: Aioli("samtools/1.10"), ... }
let dispatch = createEventDispatcher();

// User input
let program = null;      // e.g. "bedtools"
let args = null;         // e.g. "intersect -a a.bed -b b.bed"

// UI
let msgInfo = "";
let msgError = "";
let elTextbox = null;    // DOM element for the CLI input


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
async function run()
{
	msgError = "";

	// Is this a valid program?
	if(program == "") {
		msgError = `Please enter a command`;
		throw msgError;
	}
	if(!(program in TOOLS)) {
		msgError = `Program <code>${program}</code> is not supported`;
		throw msgError;
	}

	// Initialize Aioli object if not already initialized
	disabled = true;
	let aioli = null;
	if(program in aiolis) {
		aioli = aiolis[program];
	} else {
		msgInfo = `Initializing ${program}...`;
		aioli = new Aioli(TOOLS[program].aioli);
		aiolis[program] = aioli;
		await aioli.init();

		// Mount sample files
		for(let file of TOOLS[program].files)
			await Aioli.mount(file.url, name=file.name);
		console.log(await aioli.ls("/urls"));

		msgInfo = "";
	}

	// Run command and send output to parent component
	let output = await aioli.exec(args);
	dispatch("output", output);

	disabled = false;
}


// On component mount
onMount(() => {
	// Run the command provided now?
	if(execute)
		run();
});

// Once the DOM settles
afterUpdate(() => {
	// Focus on command line
	elTextbox.focus();

	// Enable jQuery tooltips. Needs to be here because tooltips are dynamically generated based on `program`
	jQuery("[data-toggle='tooltip']").tooltip();
});


// -----------------------------------------------------------------------------
// HTML
// -----------------------------------------------------------------------------
</script>

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
				style="font-size:100%; font-family:'Courier New',Courier,monospace"
			/>
			<div class="input-group-append">
				<button
					class="btn btn-md btn-outline-secondary dropdown-toggle"
					data-toggle="dropdown"
					aria-expanded="false"
				>
					Examples
				</button>

				<div class="dropdown-menu">
					{#each TOOLS[program]?.queries || [] as queries}
						<h6 class="dropdown-header">{queries.header}</h6>

						{#each queries.items as item}
							<a
								class="dropdown-item"
								data-toggle="tooltip"
								data-placement="left"
								href="#"
								title="{item.tooltip}"
								on:click={() => {
									command = item.command;
									setTimeout(run, 200);
									msgInfo = item.description;
								}}
							>
								&nbsp;&nbsp;<code>{item.label}</code>
							</a>
						{/each}
					{/each}
				</div>
			</div>

			<div class="input-group-append">
				<button
					class="btn btn-md {disabled ? 'btn-secondary' : 'btn-primary'}"
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
