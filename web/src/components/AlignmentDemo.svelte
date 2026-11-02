<script>
// Live minimap2 demo aligning editable reads to the human mitochondrial genome bundled in its .data file
import { onMount } from "svelte";
import { Icon } from "@sveltestrap/sveltestrap";
import CodeBlock from "$components/CodeBlock.svelte";

const TOOL = "minimap2/2.22";
const REFERENCE = "/minimap2/MT-human.fa";
const REFERENCE_LENGTH = 16569;
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
let alignments = [];
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
		await run();
		const fasta = await CLI.fs.readFile(REFERENCE, { encoding: "utf8" });
		genome = fasta.split("\n").filter(line => !line.startsWith(">")).join("").toUpperCase();
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
		alignments = parseSAM(output.stdout);
		if(alignments.length === 0)
			throw new Error("No reads found. Paste them in FASTA format: a >name line, then the sequence.");
		status = "ready";
	} catch (e) {
		status = "error";
		error = e.message || String(e);
	} finally {
		clearTimeout(timer);
		slow = false;
	}
}

function parseSAM(sam) {
	return sam.split("\n").filter(line => line && !line.startsWith("@")).map(line => {
		const [name, flag, , pos, , cigar, , , , seq] = line.split("\t");
		return {
			name,
			flag: +flag,
			mapped: !(flag & 4),
			strand: flag & 16 ? "-" : "+",
			pos: +pos,
			cigar,
			length: seq.length
		};
	}).filter(a => !(a.flag & 256) && !(a.flag & 2048));
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

$: mapped = alignments.filter(a => a.mapped);
</script>

<div class="card shadow-sm">
	<div class="card-header d-flex align-items-center justify-content-between bg-white">
		<span class="small fw-semibold text-nowrap">
			<Icon name="lightning-charge-fill" class="text-warning" />
			<span class="d-none d-sm-inline">Live demo:</span><span class="d-sm-none">Demo:</span> <code>minimap2</code> <span class="d-none d-sm-inline">in your browser</span>
			{#if simd}
				<span class="badge text-bg-success fw-normal ms-1" title="Your browser supports WebAssembly SIMD, so this runs the vectorized build of minimap2">SIMD</span>
			{/if}
		</span>
		<div class="d-flex gap-1 flex-shrink-0 ms-2">
			<button class="btn btn-sm btn-outline-secondary py-0" on:click={shuffle} disabled={!genome || slow} title="New random reads" aria-label="New random reads">
				<Icon name="shuffle" />
			</button>
			<button class="btn btn-sm btn-primary py-0 text-nowrap" on:click={run} disabled={status === "loading" || slow}>
				{#if status === "loading" || slow}
					<span class="spinner-border spinner-border-sm" role="status" aria-hidden="true"></span>
				{:else}
					Align <Icon name="play-fill" />
				{/if}
			</button>
		</div>
	</div>

	<div class="card-body p-0">
		<label class="small text-muted px-3 pt-2 mb-1 d-block" for="demo-reads">
			Reads to align to the human mitochondrial genome (edit them!)
		</label>
		<textarea
			id="demo-reads"
			class="reads form-control border-0 rounded-0 px-3"
			rows="5"
			spellcheck="false"
			bind:value={reads}
			on:keydown={handleKeydown}
		></textarea>

		<div class="border-top px-3 py-2">
			<svg class="w-100 d-block" viewBox="0 0 1000 28" preserveAspectRatio="none" aria-label="Read positions along the genome">
				<rect x="0" y="12" width="1000" height="4" rx="2" class="genome" />
				{#each mapped as a}
					<rect
						x={(a.pos / REFERENCE_LENGTH) * 1000}
						y="6"
						width={Math.max(6, (a.length / REFERENCE_LENGTH) * 1000)}
						height="16"
						rx="2"
						class="read"
					>
						<title>{a.name}: {a.pos.toLocaleString()} ({a.strand}) {a.cigar}</title>
					</rect>
				{/each}
			</svg>
			<div class="d-flex justify-content-between small text-muted"><span>1</span><span>{REFERENCE_LENGTH.toLocaleString()} bp</span></div>
		</div>

		{#if status === "error"}
			<div class="text-danger small px-3 pb-2">{error}</div>
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
	resize: vertical;
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

.reads:focus {
	box-shadow: none;
}

.genome {
	fill: var(--bs-gray-300);
}

.read {
	fill: var(--bs-primary);
	opacity: 0.75;
	transition: x 0.3s;
}

.read:hover {
	opacity: 1;
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
