<script>
import { Icon } from "@sveltestrap/sveltestrap";
import LinkPackage from "$components/LinkPackage.svelte";
import AlignmentDemo from "$components/AlignmentDemo.svelte";
import CONFIG from "@/biowasm.json";

const NB_PACKAGES = CONFIG.tools.filter(t => t.listed !== false).length;

const WHY = [
	{ icon: "piggy-bank", title: "Reduce cloud costs", text: "Run computations in each user's browser; no cloud infra to manage, and costs don't grow with usage." },
	{ icon: "shield-lock", title: "Data privacy", text: "Files are analyzed on the user's device and never uploaded, which is essential for clinical usage." },
	{ icon: "lightning-charge", title: "Better UX", text: "For small interactive analyses, results come back faster, with no uploads or job queues to wait for." },
	{ icon: "mortarboard", title: "Teaching & offline use", text: "No environment to set up, and tools keep working even with intermittent Internet connection." }
];

// Apps that use biowasm; logos live in static/users/. `label` is shown as text next to (or instead of) the logo,
// `icon` is a Bootstrap icon shown before the label, `dim` darkens pale logos so they stay legible
// in grayscale, `size` scales down logos that look too heavy, and `text` is shown as a tooltip
const USERS = [
	{ name: "Nanopore", url: "https://bonito.epi2me.io", logo: "/users/nanopore.svg", label: null, text: "In-browser basecaller, uses <code>minimap2</code> and <code>samtools</code> after basecalling" },
	{ name: "sandbox.bio", url: "https://sandbox.bio", logo: null, icon: "terminal-fill", label: "sandbox.bio", text: "Runs command-line tools in the browser to power interactive tutorials" },
	{ name: "SeqtoID", url: "https://seqtoid.org/", logo: "/users/seqtoid.png", dim: true, label: null, text: "Runs <code>htsfile</code> and <code>seqtk</code> to identify data issues before file upload" },
	{ name: "ViralWasm", url: "https://niema-lab.github.io/ViralWasm", logo: null, label: "🧪 ViralWasm", text: "Runs <code>minimap2</code> and <code>ViralConsensus</code> for viral molecular epidemiology analysis" },
	{ name: "Ribbon", url: "https://genomeribbon.com", logo: null, icon: "bezier", label: "Ribbon", text: "Runs <code>samtools</code> in the browser to parse, estimate coverage and subsample .bam files" },
	{ name: "42basepairs", url: "https://42basepairs.com", logo: "/users/42basepairs.svg", label: "42basepairs", text: "Runs <code>samtools</code>, <code>bedtools</code>, <code>bcftools</code> and other tools to preview genomic files" },
	{ name: "Wasm4Learn", url: "https://ifb-elixirfr.github.io/Wasm4Learn/", logo: "/users/wasm4learn.svg", dim: true, label: null, text: "Runs <code>gawk</code>, <code>grep</code> and <code>sed</code> in the browser for bioinformatics teaching" },
	{ name: "Datamonkey", url: "https://v3.datamonkey.org", logo: "/users/datamonkey.svg", dim: true, size: 0.85, label: null, text: "Runs <code>hyphy</code> in the browser for molecular evolution analysis" },
	{ name: "Datagrok", url: "https://datagrok.ai", logo: "/users/datagrok.svg", label: "datagrok", text: "Runs <code>kalign</code> in the browser for multiple-sequence alignment analysis" },
	{ name: "Diana", url: "https://diana.cpg.org.au", logo: "/users/diana.webp", label: "Diana", text: "Runs <code>samtools</code>, <code>bgzip</code> and <code>tabix</code> in the browser to index genomes and annotations for TraDIS analysis" },
	{ name: "GBIF", url: "https://www.gbif.org/occurrence/search?view=dashboard&nucleotideSequence.sequence=TTATATTTTATATTCGGAATTTGGGCAGGAATAATTGGGACATCATTAAGATTATTAATTCGAGCTGAATTAGGAAACCCAGGATCTCTAATTGGGGATGATCAAATTTATAACTCAATTGTAACAGCTCATGCTTTTATTATAATTTTCTTTATAGTAATACCTATTATAATTGGAG&nucleotideSequence.similarity=%5B%2298%22%2C%2299%22%5D", logo: "/users/gbif.svg", label: null, text: "Runs <code>kalign</code> to align sequences for its sequence phylogeny dashboard charts" },
	{ name: "Sequence Manipulation Suite", url: "https://paulstothard.github.io/sequence-manipulation-suite/", logo: null, icon: "tools", label: "SMS3", text: "The Sequence Manipulation Suite runs <code>minimap2</code>, <code>muscle</code>, <code>samtools</code> and other tools in the browser alongside its classic sequence analysis utilities" },
	{ name: "json4u", url: "https://json4u.com", logo: "/users/json4u.svg", label: "json4u", text: "Runs <code>jq</code> in the browser to power its JSON query feature" },
	{ name: "CRISPRLungo", url: "https://pinellolab.github.io/CRISPRLungo/", logo: "/users/crisprlungo.png", size: 0.75, label: null, text: "Runs <code>minimap2</code> as part of a long-read CRISPR analysis pipeline" },
	{ name: "bedqc", url: "https://quinlan-lab.github.io/bedqc", logo: null, icon: "bar-chart", label: "bedqc", text: "Runs <code>bedtools</code> in the browser to validate .bed files" }
];
</script>

<svelte:head>
	<title>Biowasm | Bioinformatics in the browser</title>
</svelte:head>

<!-- Jumbotron -->
<div class="px-3 px-md-5 py-5 mb-4 bg-light rounded-3">
	<div class="row align-items-center gx-0 gx-lg-5 gy-4">
		<div class="col-lg-5">
			<h1 class="fs-2 fw-bold mb-3">Run bioinformatics tools in your browser</h1>
			<p class="fs-5 mb-4">Use command-line tools in your web apps with just a few lines of code. Powered by WebAssembly.</p>
			<div class="d-flex flex-wrap gap-2">
				<a href="/documentation" class="btn btn-primary" type="button">
					Get Started <Icon name="arrow-right" />
				</a>
				<a href={CONFIG.url} class="btn btn-outline-primary" type="button">
					Browse packages ({NB_PACKAGES}) <Icon name="arrow-right" />
				</a>
			</div>
		</div>
		<div class="col-lg-7">
			<AlignmentDemo />
		</div>
	</div>
</div>

<!-- Who uses biowasm -->
<div class="px-3 px-md-5 pt-3 pb-4 mb-4">
	<h2 class="users-title small fw-normal text-uppercase text-body-secondary text-center mb-3">Who uses biowasm?</h2>
	<div class="d-flex flex-wrap align-items-center justify-content-center gap-2 gap-md-3">
		{#each USERS as app}
			<a class="app position-relative px-3 py-2 rounded-3 text-reset text-decoration-none" href={app.url}>
				<span class="app-logo d-flex align-items-center gap-1">
					{#if app.logo}
						<img src={app.logo} alt={app.label ? "" : app.name} class:dim={app.dim} class:icon-img={app.label} style:height={app.size && `${2 * app.size}rem`} />
					{/if}
					{#if app.icon}
						<Icon name={app.icon} class="fs-5" />
					{/if}
					{#if app.label}
						<span class="fs-6 fw-bold">{app.label}</span>
					{/if}
				</span>
				<span class="app-tip">{@html app.text}</span>
			</a>
		{/each}
	</div>
</div>

<div class="row g-4 align-items-md-stretch mb-5">
	<!-- Why the browser -->
	<div class="col-lg-8">
		<div class="h-100 p-5 bg-light rounded-3">
			<h2 class="fs-4 fw-bold mb-4">Why bioinformatics in the browser?</h2>
			<div class="row g-4">
				{#each WHY as reason}
					<div class="col-md-6">
						<div class="d-flex align-items-center mb-2">
							<span class="why-icon d-inline-flex align-items-center justify-content-center rounded-3 me-2">
								<Icon name={reason.icon} />
							</span>
							<h3 class="fs-6 fw-semibold mb-0">{reason.title}</h3>
						</div>
						<p class="text-muted mb-0">{reason.text}</p>
					</div>
				{/each}
			</div>
		</div>
	</div>

	<!-- Supported tools -->
	<div class="col-lg-4">
		<div class="h-100 p-5 rounded-3">
			<h2 class="fs-4 fw-bold">Supported Tools</h2>
			<p>
				Biowasm supports many C/C++ bioinformatics tools (<LinkPackage tool="samtools" />, <LinkPackage tool="bedtools" />, <LinkPackage tool="bcftools" />, <LinkPackage tool="bowtie2" />, <LinkPackage tool="minimap2" />, <LinkPackage tool="fastp" />, <LinkPackage tool="wgsim" />, etc.)
				and general data wrangling utilities (<LinkPackage tool="jq" />, <LinkPackage tool="gawk" />, <LinkPackage tool="grep" />, <LinkPackage tool="sed" />, etc.)
			</p>
			<a href={CONFIG.url} class="btn btn-outline-primary" type="button">
				Browse all packages <Icon name="arrow-right" />
			</a>
		</div>
	</div>
</div>

<style>
.why-icon {
	width: 2rem;
	height: 2rem;
	color: var(--bs-primary);
	background: var(--bs-primary-bg-subtle);
}
.users-title {
	letter-spacing: 0.05em;
}
.app:hover {
	background: var(--bs-tertiary-bg);
}
.app-logo {
	min-height: 2rem;
	filter: grayscale(1);
	opacity: 0.75;
}
.app:hover .app-logo {
	filter: none;
	opacity: 1;
}
.app-logo img {
	height: 2rem;
	max-width: 12rem;
	object-fit: contain;
}
/* Icons paired with a label match the size of the Bootstrap icons */
.app-logo img.icon-img {
	height: 1.5rem;
}
.app:not(:hover) .app-logo img.dim {
	filter: brightness(0.45);
}
/* Custom tooltip, since the native `title` one only shows up after a delay */
.app-tip {
	display: none;
	position: absolute;
	top: 100%;
	left: 50%;
	transform: translateX(-50%);
	z-index: 10;
	width: max-content;
	max-width: 16rem;
	margin-top: 0.25rem;
	padding: 0.375rem 0.625rem;
	font-size: 0.8125rem;
	text-align: center;
	color: var(--bs-white);
	background: var(--bs-dark);
	border-radius: var(--bs-border-radius);
	pointer-events: none;
}
.app:hover .app-tip,
.app:focus-visible .app-tip {
	display: block;
}
.app-tip :global(code) {
	color: inherit;
	font-weight: 600;
}
</style>
