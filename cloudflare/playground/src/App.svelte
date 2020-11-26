<script>
// Bootstrap
import popper from "popper.js";
import Bootstrap from "bootstrap";
import "bootstrap/dist/css/bootstrap.min.css";

// Imports
import CommandLine from "./CommandLine.svelte";

// State
let output = "";
</script>

<CommandLine
	command="samtools --version"
	execute={false}
	on:output={msg => output = msg.detail.stdout.trim() + msg.detail.stderr.trim()} />

<pre class="border rounded border-primary p-3" style="height:55vh">{output}</pre>
