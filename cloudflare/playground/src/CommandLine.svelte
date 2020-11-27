<script>
// Exports
export let command = "";        // Command to execute (e.g. samtools --version)

// Imports
import { afterUpdate, createEventDispatcher } from "svelte";
import jQuery from "jquery";
import { Aioli } from "@biowasm/aioli";

import { TOOLS, UTILITIES, BIOWASM_URL } from "./config.js";

// -----------------------------------------------------------------------------
// Globals
// -----------------------------------------------------------------------------

let dispatch = createEventDispatcher();

// State
let State = {
	aiolis: {},       // e.g. { samtools: Aioli("samtools/1.10"), ... }
	program: null,    // e.g. "samtools"
}

// UI
let UI = {
	msgInfo: "",      // Info text above the input box
	msgError: "",     // Error text below the input box
	textbox: null,    // DOM element for the input box
	disabled: false,  // Disable CommandLine
}

// -----------------------------------------------------------------------------
// On load
// -----------------------------------------------------------------------------

// Once the DOM settles
afterUpdate(() => {
	// Focus on command line
	UI.textbox.focus();

	// Enable jQuery tooltips. Needs to be here because tooltips are dynamically generated based on `program`
	jQuery("[data-toggle='tooltip']").tooltip();
});

// -----------------------------------------------------------------------------
// Utility functions
// -----------------------------------------------------------------------------

// Execute a command and update UI
export async function launch(cmd=null)
{
	// By default, execute what's in the input box
	cmd = cmd || command;

	// Prep UI
	UI.msgInfo = `Running...`;
	UI.disabled = true;

	// Make UI match command being run
	if(cmd != command)
		command = cmd;

	// Run command and send output to parent component
	let output = {};
	try {
		output = await run(cmd);		
	} catch (error) {
		console.warn(error);
		output = { stdout: "", stderr: String(error) };
	}
	dispatch("output", output);

	// Revert UI
	UI.msgInfo = "Ready.";
	UI.disabled = false;
}

// Execute a command and return value
export async function run(cmd)
{
	let program = cmd.split(" ").shift();
	let args = cmd.replace(`${program} `, "").trim();

	// Validate user input
	let error = "";
	if(program == "")
		error = `Please enter a command`;
	if(!(program in TOOLS) && !(program in UTILITIES))
		error = `Program <code>${program}</code> is not supported`;
	if(args.includes("|"))
		error = "Piping is not currently supported";
	if(error != "")
		throw error;

	// Process non-WebAssembly utilities
	if(program in UTILITIES)
	{
		// Retrieve first Aioli object we see; file system utilities need existing object to work
		let aioli = Object.values(State.aiolis).shift();
		// Run utility and send output to parent component
		let output = { stdout: "", stderr: "" };
		try {
			if(cmd == program)
				args = null;
			output.stdout = await UTILITIES[program](aioli, args);
		} catch (error) {
			console.error(error);
			output.stderr = "Error: invalid command";
		}
		return output;
	}

	// Process WebAssembly utilities. Initialize Aioli object if needed.
	let aioli = null;
	if(program in State.aiolis) {
		aioli = State.aiolis[program];
	} else {
		// Create Aioli object for program
		let config = TOOLS[program].aioli;
		config.urlModule = `${BIOWASM_URL}/${config.module}/${config.version}`;
		config.urlAioli = `${BIOWASM_URL}/aioli/latest/aioli.worker.js`;

		aioli = new Aioli(config);
		State.aiolis[program] = aioli;
		await aioli.init();

		// Mount sample files
		if(TOOLS[program].files != null)
			for(let file of TOOLS[program].files)
				await Aioli.mount(file.url, name=file.name);
		
		// Set working directory
		aioli.fs("chdir", "/urls");
	}
	State.program = program;
	return await aioli.exec(args);
}

// -----------------------------------------------------------------------------
// HTML
// -----------------------------------------------------------------------------
</script>

<!-- Info message -->
<div class="row">
	<div class="col-12">
		<span class="text-muted">
			<span class="text-info">{@html UI.msgInfo}&nbsp;</span>
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
				disabled={UI.disabled}
				bind:this={UI.textbox}
				bind:value={command}
				on:keydown={event => event.key == "Enter" ? launch(command) : null}
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
					{#each TOOLS[State.program]?.queries || [] as queries}
						<h6 class="dropdown-header">{queries.header}</h6>

						{#each queries.items as item}
							<a
								class="dropdown-item"
								data-toggle="tooltip"
								data-placement="left"
								href="#"
								title="{item.tooltip}"
								on:click={async () => {
									await launch(item.command);
									UI.msgInfo = item.description;
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
					class="btn btn-md {UI.disabled ? 'btn-secondary' : 'btn-primary'}"
					disabled={UI.disabled}
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
			<span class="text-danger">{@html UI.msgError}&nbsp;</span>
		</small>
	</div>
</div>
