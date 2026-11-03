<script>
// Live minimap2 demo aligning editable reads to the human mitochondrial genome bundled in its .data file
import { onMount } from "svelte";
import { Icon } from "@sveltestrap/sveltestrap";
import CodeBlock from "$components/CodeBlock.svelte";

const TOOL = "minimap2/2.22";
const REFERENCE = "/minimap2/MT-human.fa";
// Smaller minimizers than the sr preset's -k21 -w11, which can't seed 60bp reads that carry an edit
const READ_LENGTH = 60;
const FLAGS = "-ax sr -k11 -w5";
const BASES = "ACGT";
// Same check Aioli makes (via wasm-feature-detect) to decide whether to load minimap2-simd
const SIMD_TEST = new Uint8Array([0,97,115,109,1,0,0,0,1,5,1,96,0,1,123,3,2,1,0,10,10,1,8,0,65,0,253,15,253,98,11]);

// Slices of MT-human.fa with one kind of edit each, checked against a native minimap2
const READS = `
>exact
AGCCTGTTCTGTAATCGATAAACCCCGATCAACCTCACCACCTCTTGCTCAGCCTATATA
>snp
TAAAAATAAAATGACAGTTTGAACATACAAGACCCACCCCATTCCTCCCCACACTCATCG
>insertion
GGACATCCCGATGGTGCAGCCGCTATTAAAGATGGTTCGTTTGTTCAACGATTAAAGTCC
>deletion
CGCCTAACCGCTAACATTACTGCAGGCCACCTCATGCACCTAATTGGAAGCGCCACCCTA
>reverse
ACATAAATAGTATGGCTTTGAAGAAGGCGTGGGTACAGATGTGCAGGAATGCTAGGTGTG
`.trim();

const CODE = `
<script src="https://biowasm.com/cdn/v3/aioli.js"><\/script>
<script type="module">
const CLI = await new Aioli(["${TOOL}"]);
const [reads] = await CLI.mount([{ name: "reads.fa", data: ">read1\\nACGT..." }]);
const sam = await CLI.exec(\`minimap2 ${FLAGS} ${REFERENCE} \${reads}\`);
<\/script>
`;

let reads = READS;
let status = "loading";  // loading, ready, running, error
let error = "";
let sam = "";
let samLines = [];  // [{ key, fields: [{ text, changed }] }]
let nbAligned = 0;
let elapsed = null;
let showCode = false;
let copied = false;
let CLI;
let simd = false;
let slow = false;  // only show a busy state for slow runs, else it flickers
let genome = "";
let nbRuns = 0;

onMount(async () => {
	simd = WebAssembly.validate(SIMD_TEST);
	try {
		await loadScript("https://biowasm.com/cdn/v3/aioli.js");
		CLI = await new Aioli([TOOL], { printInterleaved: false });
		const fasta = await CLI.fs.readFile(REFERENCE, { encoding: "utf8" });
		genome = fasta.split("\n").filter(line => !line.startsWith(">")).join("").toUpperCase();
		status = "ready";
	} catch (e) {
		status = "error";
		error = e.message || String(e);
	}
});

async function run() {
	if(status === "running")
		return;
	status = "running";
	error = "";
	const timer = setTimeout(() => (slow = true), 300);
	try {
		// Remounting the same name doesn't overwrite it
		const [path] = await CLI.mount([{ name: `reads-${++nbRuns}.fa`, data: reads }]);
		const start = performance.now();
		const output = await CLI.exec(`minimap2 ${FLAGS} ${REFERENCE} ${path}`);
		elapsed = performance.now() - start;
		sam = output.stdout.split("\n").filter(line => !line.startsWith("@PG")).join("\n").trim();
		if(!sam.split("\n").some(line => !line.startsWith("@")))
			throw new Error("No reads found. Paste them in FASTA format: a >name line, then the sequence.");
		samLines = diffSAM(samLines, sam, ++nbAligned);
		status = "ready";
	} catch (e) {
		status = "error";
		error = e.message || String(e);
	} finally {
		clearTimeout(timer);
		slow = false;
	}
}

// Cut new reads from random spots in the genome, keeping edits mid-read so minimap2 can anchor them
function shuffle() {
	const randInt = (min, max) => Math.round(min) + Math.floor(Math.random() * (Math.round(max) - Math.round(min) + 1));
	const randBases = n => Array.from({ length: n }, () => BASES[randInt(0, 3)]).join("");
	const slice = n => {
		const start = randInt(0, genome.length - n);
		return genome.slice(start, start + n);
	};
	const complement = { A: "T", C: "G", G: "C", T: "A" };

	const snp = slice(READ_LENGTH).split("");
	const i = randInt(READ_LENGTH * 0.2, READ_LENGTH * 0.8);
	snp[i] = BASES.replace(snp[i], "")[randInt(0, 2)];

	const insertionLength = randInt(1, 4);
	const insertion = slice(READ_LENGTH - insertionLength);
	const j = randInt(READ_LENGTH * 0.3, READ_LENGTH * 0.6);

	const deletionLength = randInt(1, 4);
	const deletion = slice(READ_LENGTH + deletionLength);
	const k = randInt(READ_LENGTH * 0.3, READ_LENGTH * 0.6);

	reads = Object.entries({
		exact: slice(READ_LENGTH),
		snp: snp.join(""),
		insertion: insertion.slice(0, j) + randBases(insertionLength) + insertion.slice(j),
		deletion: deletion.slice(0, k) + deletion.slice(k + deletionLength),
		reverse: slice(READ_LENGTH).split("").reverse().map(b => complement[b] || "N").join("")
	}).map(([name, seq]) => `>${name}\n${seq}`).join("\n");

	run();
}

async function copyCode() {
	await navigator.clipboard.writeText(CODE.trim());
	copied = true;
	setTimeout(() => (copied = false), 1500);
}

// Flag each field that differs from the same read's line in the previous run
function diffSAM(previous, sam, run) {
	const before = Object.fromEntries(previous.map(line => [line.fields[0].text, line.fields.map(f => f.text)]));
	return sam.split("\n").map((line, i) => {
		const fields = line.split("\t");
		const old = before[fields[0]];
		return {
			key: `${run}-${i}`,
			fields: fields.map((text, j) => ({ text, changed: old !== undefined && old[j] !== text }))
		};
	});
}

function loadScript(src) {
	if(window.Aioli)
		return Promise.resolve();
	return new Promise((resolve, reject) => {
		const script = document.createElement("script");
		script.src = src;
		script.onload = resolve;
		script.onerror = () => reject(new Error(`Could not load ${src}`));
		document.head.appendChild(script);
	});
}

function handleKeydown(evt) {
	if((evt.metaKey || evt.ctrlKey) && evt.key === "Enter" && status !== "loading" && status !== "running") {
		evt.preventDefault();
		run();
	}
}

</script>

<div class="card shadow-sm">
	<div class="card-header d-flex flex-wrap align-items-center justify-content-between gap-2 bg-white">
		<span class="small fw-semibold text-nowrap">
			<Icon name="lightning-charge-fill" class="text-warning" />
			<span class="d-none d-sm-inline">Live demo:</span><span class="d-sm-none">Demo:</span> <code>minimap2</code> <span class="d-none d-sm-inline">in your browser</span>
			{#if simd}
				<span class="badge text-bg-success fw-normal ms-1" title="Your browser supports WebAssembly SIMD, so this runs the vectorized build of minimap2">SIMD</span>
			{/if}
		</span>
		<div class="d-flex gap-1 flex-shrink-0">
			<button class="btn btn-sm btn-outline-secondary py-0 text-nowrap" on:click={shuffle} disabled={!genome || slow} title="New random reads">
				<span class="d-none d-sm-inline me-1">Shuffle</span><Icon name="shuffle" />
			</button>
			<button class="btn btn-sm btn-primary py-0 text-nowrap" on:click={run} disabled={status === "loading" || slow}>
				Align <Icon name="play-fill" />
			</button>
		</div>
	</div>

	<div class="card-body p-0">
		<div class="d-flex flex-column flex-sm-row">
			<div class="pane d-flex flex-column">
				<label class="small text-muted px-3 pt-2 mb-1" for="demo-reads">Reads (edit them!)</label>
				<textarea
					id="demo-reads"
					class="reads form-control border-0 rounded-0 px-3 flex-grow-1"
					spellcheck="false"
					bind:value={reads}
					on:keydown={handleKeydown}
				></textarea>
			</div>
			<div class="pane output d-flex flex-column">
				<span class="small text-muted px-3 pt-2 mb-1">SAM output (vs. human mtDNA)</span>
				{#if sam}
					<pre class="sam flex-grow-1 m-0 px-3 pb-2">{#each samLines as line (line.key)}{#each line.fields as field, j}{j ? "\t" : ""}<span class:changed={field.changed}>{field.text}</span>{/each}{"\n"}{/each}</pre>
				{:else}
					<div class="flex-grow-1 small text-muted fst-italic px-3">Click Align to run minimap2</div>
				{/if}
			</div>
		</div>

		{#if status === "error"}
			<div class="text-danger small border-top px-3 py-2">{error}</div>
		{/if}
	</div>

	<div class="card-footer bg-white d-flex align-items-center justify-content-between small">
		<span class="text-muted">
			{#if elapsed !== null}
				Aligned in {Math.round(elapsed)} ms. Your data never left your browser.
			{/if}
		</span>
		<button class="btn btn-link btn-sm p-0 ms-2 text-nowrap" on:click={() => (showCode = !showCode)}>
			<Icon name="code-slash" /> {showCode ? "Hide" : "Show"} code
		</button>
	</div>
	{#if showCode}
		<div class="code position-relative px-3 pb-1">
			<CodeBlock code={CODE} />
			<button class="copy btn btn-sm btn-light border py-0 px-1 position-absolute" on:click={copyCode} title="Copy to clipboard" aria-label="Copy to clipboard">
				<Icon name={copied ? "check-lg" : "clipboard"} class={copied ? "text-success" : ""} />
			</button>
		</div>
	{/if}
</div>

<style>
.reads {
	font-family: var(--bs-font-monospace);
	font-size: 12px;
	resize: none;
	/* Wrap rather than scroll on narrow screens */
	white-space: pre-wrap;
	word-break: break-all;
}

/* iOS Safari zooms the page when focusing an input with text under 16px */
@media (max-width: 767.98px) {
	.reads {
		font-size: 16px;
	}
}

.pane {
	flex: 1 1 0;
	min-width: 0;
	height: 160px;
}

.output {
	border-top: var(--bs-border-width) solid var(--bs-border-color);
}

@media (min-width: 576px) {
	.output {
		border-top: 0;
		border-left: var(--bs-border-width) solid var(--bs-border-color);
	}
}

/* Long SAM lines scroll sideways instead of wrapping, so columns stay readable */
/* Fades out, and the keyed #each restarts it on every run */
.changed {
	border-radius: 2px;
	animation: flash 2.5s ease-out;
}

@keyframes flash {
	from { background: rgba(255, 193, 7, 0.6); }
	to { background: transparent; }
}

.sam {
	font-size: 12px;
	overflow: auto;
	white-space: pre;
}

.reads:focus {
	box-shadow: none;
}

.copy {
	top: 0.35rem;
	right: calc(1rem + 0.35rem);
}

.code :global(pre) {
	font-size: 11px;
	/* CodeBlock's p-3 is a utility class, hence the !important */
	padding: 0.5rem 0.75rem !important;
	margin-bottom: 0.5rem;
}
</style>
