<script>
// Exports
export let command = "";

// Imports
import { afterUpdate, createEventDispatcher } from "svelte";
import jQuery from "jquery";
import { Aioli } from "@biowasm/aioli";

import { TOOLS, UTILITIES, PIPING, BIOWASM_URL } from "./config.js";

// -----------------------------------------------------------------------------
// Globals
// -----------------------------------------------------------------------------

let dispatch = createEventDispatcher();

// State
let State = {
	aiolis: {},       // e.g. { samtools: Aioli("samtools/1.10"), ... }
	program: null,    // e.g. "samtools"
	timeout: null,    // Timeout object for making the "Stop" button visible
}

// UI
let UI = {
	msgInfo: "",      // Info text above the input box
	msgError: "",     // Error text below the input box
	textbox: null,    // DOM element for the input box
	fileInput: null,  // DOM element for the hidden file input box
	disabled: false,  // Disable CommandLine
	stopBtnShow: false,      // If true, show "Stop" instead of "Run" btn
	stopBtnDisabled: false,  // If true, "Stop" button is disabled
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

	// Enable .terminal buttons to launch a command
	jQuery(".terminal").off("click").click(function() {
		launch(this.value);
	});
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
	UI.msgError = "";
	UI.disabled = true;
	State.timeout = setTimeout(() => UI.stopBtnShow = true, 2000);

	// Make UI match command being run
	if(cmd != command)
		command = cmd;

	// Support very basic output redirection.
	// Assumes ">" is only used for redirection (and not say part of a string argument)
	let redirection = null;
	let redirections = cmd.split(">").map(d => d.trim());
	if(redirections.length > 2) {
		UI.msgError = "Unsupported command.";
		UI.disabled = false;
		return;
	} else if(redirections.length == 2) {
		cmd = redirections[0];
		redirection = redirections[1];
	}

	// Support a few convenient piping functions.
	// Assumes "|" is only used for piping (and not say part of a filename).
	let piping = null;
	let pipes = cmd.split("|").map(d => d.trim());
	if(pipes.length > 2 || (pipes.length == 2 && !Object.keys(PIPING).includes(pipes[1]))) {
		UI.msgError = "Piping not supported.";
		UI.disabled = false;
		return;
	} else if(pipes.length == 2) {
		cmd = pipes[0];
		piping = pipes[1];
	}

	// Run command
	let output = {};
	try {
		output = await run(cmd);
		if(piping != null)
			output.stdout = PIPING[piping](output.stdout);
	} catch (error) {
		console.warn(error);
		output = { stdout: "", stderr: String(error) };
	}

	// Send result to parent component unless using redirection
	if(redirection != null) {
		let blob = new Blob([output.stdout], { type: "text/plain" });
		await Aioli.mount(blob, redirection);
		output.stdout = "ok";
	}
	dispatch("output", output);

	// Revert UI
	UI.msgInfo = "Ready.";
	UI.disabled = false;
	UI.stopBtnShow = UI.stopBtnDisabled = false;
	clearTimeout(State.timeout);
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
	let aioli = State.aiolis[program];
	if(aioli != null) {
		aioli = State.aiolis[program];
	} else {
		// Create Aioli object for program
		let config = TOOLS[program].aioli;
		config.urlModule = `${BIOWASM_URL}/${config.module}/${config.version}`;
		config.urlAioli = `${BIOWASM_URL}/aioli/latest/aioli.worker.js`;

		aioli = new Aioli(config);
		State.aiolis[program] = aioli;
		await aioli.init();

		// Mount URLs with sample data
		if(Array.isArray(TOOLS[program].files))
			for(let file of TOOLS[program].files)
				await Aioli.mount(file.url, file.name);
		// Mount user files that were previously mounted in other workers
		for(let file of Aioli.files)
			if(file.source == "file")
				await Aioli.mount(file.file, file.name);

		// Set working directory
		aioli.fs("chdir", "/urls");
	}
	State.program = program;
	return await aioli.exec(args);
}

// Stop what's currently running
export async function stop()
{
	console.group("Stop");

	let aioli = State.aiolis[State.program];
	if(aioli != null)
	{
		// Delete worker from Aioli.workers
		for(let i in Aioli.workers)
			if(Aioli.workers[i].worker == aioli.worker) {
				Aioli.workers.splice(i, 1);
				console.log("Deleted from Aioli.workers");
				break;
			}

		// Terminate the WebWorker
		await aioli.worker.terminate();
		State.aiolis[State.program] = null;
		console.log("Aioli Worker terminated");

		// Reload the worker
		await run(`${State.program} --help`);
		console.log("Aioli Worker reloaded");
	}
	console.groupEnd("Stop");

	UI.stopBtnShow = UI.stopBtnDisabled = false;
	UI.msgInfo = "";
	UI.msgError = "Stopped.";
	UI.disabled = false;
}

// Mount selected files
async function mountFiles()
{
	let files = Array.from(UI.fileInput.files);
	for(let file of files) {
		let fileMounted = await Aioli.mount(file);
		console.log(`Mounted ${file.name} to ${fileMounted.path}`);
	}
	launch("ls /data");
}


// -----------------------------------------------------------------------------
// HTML
// -----------------------------------------------------------------------------
</script>

<div class="row">
	<!-- Info message -->
	<div class="col-10">
		<span class="text-muted">
			<span class="text-info">{@html UI.msgInfo}</span>
		</span>
	</div>
	<!-- Choose a file locally -->
	<div class="col-2 text-right">
		<span class="text-muted">
			<button
				on:click={() => UI.fileInput.click()}
				type="button" class="btn btn-link p-0 text-info"
				data-step="4"
				data-position="left"
				data-intro="Load files from your computer and run bioinformatics commands on them. Your files are mounted (virtually) to <code>/data</code> (this app can only read, not modify your files)"
				style="vertical-align: baseline;">
				Load a local file
			</button>
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
				data-step="2"
				data-intro="This is where you enter commands. Pressing 'Enter' executes the command and displays the result below."
				disabled={UI.disabled}
				bind:this={UI.textbox}
				bind:value={command}
				on:keydown={event => event.key == "Enter" ? launch(command) : null}
				style="font-size:100%; font-family:'Courier New',Courier,monospace" />
			<div class="input-group-append">
				<button
					class="btn btn-md btn-outline-secondary dropdown-toggle"
					data-toggle="dropdown"
					aria-expanded="false"
					data-step="3"
					data-position="left"
					data-intro="Sample queries to get you started">
					Examples
				</button>

				<div class="dropdown-menu">
					{#each TOOLS[State.program]?.queries || [] as queries}
						<h6 class="dropdown-header">{queries.header}</h6>

						{#each queries.items as item}
							<a
								href="#"
								class="dropdown-item"
								data-toggle="tooltip"
								data-placement="left"
								data-html="true"
								data-original-title="{item.tooltip}"
								title="{item.tooltip}"
								on:click={async () => {
									await launch(item.command);
									UI.msgInfo = item.description;
								}}>
								&nbsp;&nbsp;<code>{item.label}</code>
							</a>
						{/each}
					{/each}
				</div>
			</div>

			{#if !UI.stopBtnShow}
			<div class="input-group-append">
				<button
					class="btn btn-md {UI.disabled ? 'btn-secondary' : 'btn-primary'}"
					disabled={UI.disabled}
					on:click={launch(command)}>
					Run
				</button>
			</div>

			{:else}
			<div class="input-group-append">
				<button
					class="btn btn-md {UI.stopBtnDisabled ? 'btn-secondary' : 'btn-danger'}"
					disabled={UI.stopBtnDisabled}
					on:click={() => {
						UI.stopBtnDisabled = true;
						stop();
					}}>
					{ UI.stopBtnDisabled ? "Stopping..." : "Stop" }
				</button>
			</div>
			{/if}
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

<!-- Hidden file input -->
<input
	bind:this={UI.fileInput}
	on:change={mountFiles}
	type="file"
	hidden multiple>
