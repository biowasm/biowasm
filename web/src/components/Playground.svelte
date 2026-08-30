<script>
import { Icon } from "sveltestrap";
import hljs from "highlight.js/lib/core";
import xml from "highlight.js/lib/languages/xml";
import javascript from "highlight.js/lib/languages/javascript";
import "highlight.js/styles/github.css";

hljs.registerLanguage("xml", xml);
hljs.registerLanguage("javascript", javascript);

export let code = "<strong>Hello</strong> World.";
export let lang = "xml";
export let autorun = true;
export let height = 250;

// Trim extra whitespace before highlighting the code (same as CodeBlock.svelte)
const codeTrimmed = code.split("\n").map(d => d.replace("\t\t", "")).join("\n").trim();

// Drop the top margin the browser puts on the first element (e.g. an <h4>) so the output
// starts at the top of the pane instead of ~20px down
const STYLE = "<style>body > :first-child { margin-top: 0 }</style>";

let source = codeTrimmed;                          // what's in the editor
let srcdoc = autorun ? STYLE + codeTrimmed : "";   // what the iframe is running
let editor, highlighted;

$: codeHighlighted = hljs.highlight(source, { language: lang }).value;

// Blank srcdoc first so re-running resets the iframe's state
function run() {
	srcdoc = "";
	setTimeout(() => (srcdoc = STYLE + source), 0);
}

// Keep the highlighted layer aligned with the textarea while scrolling
function syncScroll() {
	highlighted.scrollTop = editor.scrollTop;
	highlighted.scrollLeft = editor.scrollLeft;
}
</script>

<div class="border rounded overflow-hidden mb-3">
	<div class="d-flex align-items-center justify-content-between bg-light border-bottom px-2 py-1">
		<span class="small text-muted">Edit the code and click Run</span>
		<button class="btn btn-sm btn-outline-primary py-0" on:click={run}>Run <Icon name="play-fill" /></button>
	</div>

	<div class="d-flex flex-column flex-md-row">
		<div class="editor position-relative w-100 border-bottom" style="height: {height}px">
			<pre class="position-absolute top-0 start-0 w-100 h-100 m-0 p-2" bind:this={highlighted} aria-hidden="true">{@html codeHighlighted}</pre>
			<textarea
				class="position-relative w-100 h-100 m-0 p-2 border-0"
				bind:this={editor}
				bind:value={source}
				on:scroll={syncScroll}
				spellcheck="false"
				aria-label="Code editor"
			></textarea>
		</div>

		<div class="w-100" style="height: {height}px">
			{#if srcdoc}
				<iframe class="w-100 h-100 border-0" {srcdoc} title="Result" sandbox="allow-scripts allow-downloads allow-modals"></iframe>
			{:else}
				<div class="d-flex align-items-center justify-content-center h-100 bg-light">
					<button class="btn btn-sm btn-outline-primary" on:click={run}>Run <Icon name="play-fill" /></button>
				</div>
			{/if}
		</div>
	</div>
</div>

<style>
/* The highlighted <pre> and the <textarea> are stacked exactly on top of each other, so
   these must render identically for the text to line up. */
.editor pre,
.editor textarea {
	font-family: var(--bs-font-monospace);
	font-size: 12px;
	line-height: 1.5;
	tab-size: 4;
	white-space: pre;
	overflow: auto;
}

.editor pre {
	pointer-events: none;
}

/* Transparent text (with a visible caret) so the highlighted code shows through. The
   selection colour must be translucent too, else it hides the code underneath it. */
.editor textarea {
	background: transparent;
	color: transparent;
	caret-color: #000;
	resize: none;
	outline: none;
}

.editor textarea::selection {
	background: rgba(103, 169, 245, 0.35);
}

@media (min-width: 768px) {
	.editor {
		border-bottom: 0 !important;
		border-right: var(--bs-border-width) solid var(--bs-border-color);
	}
}
</style>
