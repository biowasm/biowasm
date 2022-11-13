<script>
import { Icon, Table } from "sveltestrap";
import LinkPackage from "$components/LinkPackage.svelte";
import CONFIG from "@/biowasm.json";
</script>

<svelte:head>
	<title>Biowasm</title>
</svelte:head>

<!-- Jumbotron -->
<div class="p-5 mb-4 bg-light rounded-3">
	<div class="container-fluid py-5">
		<h1 class="display-6 fw-bold col-md-8">Run C/C++ genomics tools in the browser</h1>
		<p class="fs-4 col-md-8">Use bioinformatics command-line tools in your web apps with just a few lines of code. Powered by WebAssembly.</p>
		<a href="/documentation" class="btn btn-primary btn-lg" type="button">
			Get Started <Icon name="arrow-right" />
		</a>
	</div>
</div>

<!-- More info -->
<div class="row align-items-md-stretch">
	<div class="col-md-4">
		<div class="h-100 p-5 text-bg-dark rounded-3">
			<h2 class="fw-bold">Supported Tools</h2>
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
			<h2 class="fw-bold">Who Uses Biowasm?</h2>
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
				</tbody>
			</Table>
		</div>
	</div>
</div>
