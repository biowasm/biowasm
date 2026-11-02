<script>
import { Icon, Table } from "@sveltestrap/sveltestrap";
import LinkPackage from "$components/LinkPackage.svelte";
import AlignmentDemo from "$components/AlignmentDemo.svelte";
import CONFIG from "@/biowasm.json";

const WHY = [
	{ icon: "piggy-bank", title: "Reduce cloud costs", text: "Run computations in each user's browser; no cloud infra to manage, and costs don't grow with usage." },
	{ icon: "shield-lock", title: "Data privacy", text: "Files are analyzed on the user's device and never uploaded, which is essential for clinical usage." },
	{ icon: "lightning-charge", title: "Better UX", text: "For small interactive analyses, results come back faster, with no uploads or job queues to wait for." },
	{ icon: "mortarboard", title: "Teaching & offline use", text: "No environment to set up, and tools keep working even with intermittent Internet connection." }
];
</script>

<svelte:head>
	<title>Biowasm | Bioinformatics in the browser</title>
</svelte:head>

<!-- Jumbotron -->
<div class="px-3 px-md-5 px-lg-4 py-4 mb-4 bg-light rounded-3">
	<div class="row align-items-center gx-0 gx-lg-4 gy-4">
		<div class="col-lg-5">
			<h1 class="fs-2 fw-bold">Run bioinformatics tools in your browser</h1>
			<p class="fs-5">Use command-line tools in your web apps with just a few lines of code. Powered by WebAssembly.</p>
			<a href="/documentation" class="btn btn-primary" type="button">
				Get Started <Icon name="arrow-right" />
			</a>
		</div>
		<div class="col-lg-7">
			<AlignmentDemo />
		</div>
	</div>
</div>

<!-- Why the browser -->
<div class="px-3 px-md-5 pb-4 mb-4">
	<h2 class="fs-4 fw-bold mb-4">Why run bioinformatics tools in the browser?</h2>
	<div class="row g-4">
		{#each WHY as reason}
			<div class="col-md-6 col-lg-3">
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

<!-- More info -->
<div class="row align-items-md-stretch">
	<div class="col-md-4">
		<div class="h-100 p-5 text-bg-light rounded-3">
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
	<div class="col-md-8">
		<div class="h-100 p-5 rounded-3">
			<h2 class="fs-4 fw-bold">Who Uses Biowasm?</h2>
			<Table>
				<thead>
					<tr>
						<th>App</th>
						<th>Why Biowasm?</th>
					</tr>
				</thead>
				<tbody>
					<tr>
						<td><a class="fw-bold" href="https://sandbox.bio">sandbox.bio</a></td>
						<td>Runs command-line tools in the browser to power interactive tutorials</td>
					</tr>
					<tr>
						<td><a class="fw-bold" href="https://42basepairs.com">42basepairs</a></td>
						<td>Runs <code>samtools</code>, <code>bedtools</code>, <code>bcftools</code> and other tools to preview genomic files</td>
					</tr>
					<tr>
						<td><a class="fw-bold" href="https://czid.org">CZ ID</a></td>
						<td>Runs <code>htsfile</code> and <code>seqtk</code> to identify data issues before file upload</td>
					</tr>
					<tr>
						<td><a class="fw-bold" href="https://bonito.epi2me.io">Nanopore</a></td>
						<td>In-browser basecaller, uses <code>minimap2</code> and <code>samtools</code> after basecalling</td>
					</tr>
					<tr>
						<td><a class="fw-bold" href="https://niema-lab.github.io/ViralWasm">ViralWasm</a></td>
						<td>Runs <code>minimap2</code> and <code>ViralConsensus</code> for viral molecular epidemiology analysis</td>
					</tr>
					<tr>
						<td><a class="fw-bold" href="https://pinellolab.github.io/CRISPRLungo/">CRISPRLungo</a></td>
						<td>Runs `minimap2` as part of a long-read CRISPR analysis pipeline</td>
					</tr>
					<tr>
						<td><a class="fw-bold" href="https://datagrok.ai">Datagrok</a></td>
						<td>Runs <code>kalign</code> in the browser for multiple-sequence alignment analysis</td>
					</tr>
					<tr>
						<td><a class="fw-bold" href="https://quinlan-lab.github.io/bedqc">bedqc</a></td>
						<td>Runs <code>bedtools</code> in the browser to validate .bed files</td>
					</tr>
					<tr>
						<td><a class="fw-bold" href="https://genomeribbon.com">Ribbon</a></td>
						<td>Runs <code>samtools</code> in the browser to parse, estimate coverage and subsample .bam files</td>
					</tr>
					<tr>
						<td><a class="fw-bold" href="https://fastq.bio">fastq.bio</a></td>
						<td>Runs <code>fastp</code> in the browser to evaluate sequencing data quality</td>
					</tr>
				</tbody>
			</Table>
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
</style>
